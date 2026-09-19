import pytest
from fastapi.testclient import TestClient

from app.config import settings
from app.enrichment import enriched_distribution
from app.envelope import convolve, merge_peaks, power_polynomial, top_peaks
from app.isotopes import ISOTOPES
from app.parser import parse_formula


def synthetic_peaks(formula, symbol, mass_number, fraction, top=10):
    """Observed-style peaks with ``symbol``'s isotope ``mass_number``
    enriched to ``fraction`` (all other elements natural)."""
    parsed = parse_formula(formula)
    dist = enriched_distribution(symbol, mass_number, fraction)
    env = [(0.0, 1.0)]
    for sym, count in parsed.composition.items():
        base = list(dist) if sym == symbol else list(ISOTOPES[sym])
        env = convolve(env, power_polynomial(base, count, 0.0001, 1e-12), 0.0001, 1e-12)
    return top_peaks(merge_peaks(env, 0.0001), top)


def synthetic_peaks_dual(formula, enrichments, top=200):
    """Observed-style peaks with two elements enriched.

    ``enrichments`` maps element symbol to (mass_number, fraction).
    """
    parsed = parse_formula(formula)
    distributions = {
        sym: enriched_distribution(sym, mass_number, fraction)
        for sym, (mass_number, fraction) in enrichments.items()
    }
    env = [(0.0, 1.0)]
    for sym, count in parsed.composition.items():
        base = list(distributions.get(sym, ISOTOPES[sym]))
        env = convolve(env, power_polynomial(base, count, 0.0001, 1e-12), 0.0001, 1e-12)
    return top_peaks(merge_peaks(env, 0.0001), top)


@pytest.fixture()
def client(tmp_path, monkeypatch):
    # Isolate each test run with its own SQLite file.
    monkeypatch.setattr(settings, "db_path", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def test_health(client):
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_envelope_glucose(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6H12O6", "top_n": 10})
    assert resp.status_code == 200
    body = resp.json()
    assert body["monoisotopic_mass"] == pytest.approx(180.063388, abs=1e-5)
    assert body["peak_count"] == 10
    assert len(body["peaks"]) == 10
    assert max(p["intensity"] for p in body["peaks"]) == pytest.approx(100.0)
    assert body["cached"] is False
    masses = [p["mass"] for p in body["peaks"]]
    assert masses == sorted(masses)


def test_envelope_with_label(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6H6[13C]2", "top_n": 5})
    assert resp.status_code == 200
    body = resp.json()
    assert body["labels"] == {"[13C]": 2}
    assert body["monoisotopic_mass"] == pytest.approx(
        6 * 12.0 + 6 * 1.00782503223 + 2 * 13.00335483507, abs=1e-5
    )


def test_envelope_cached_on_second_call(client):
    payload = {"formula": "C2H6O", "top_n": 5}
    first = client.post("/api/v1/envelope", json=payload).json()
    second = client.post("/api/v1/envelope", json=payload).json()
    assert first["cached"] is False
    assert second["cached"] is True
    assert first["peaks"] == second["peaks"]


def test_invalid_formula_reports_position(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6Xx2"})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["position"] == 2
    assert "Xx" in detail["message"]


def test_invalid_label_mass_reports_position(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6[20C]"})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["position"] == 2
    assert "20" in detail["message"]


def test_top_n_bounds(client):
    resp = client.post("/api/v1/envelope", json={"formula": "H2O", "top_n": 0})
    assert resp.status_code == 422
    resp = client.post(
        "/api/v1/envelope", json={"formula": "H2O", "top_n": settings.top_n_max + 1}
    )
    assert resp.status_code == 422


def test_missing_formula(client):
    resp = client.post("/api/v1/envelope", json={})
    assert resp.status_code == 422


def test_elements_endpoint(client):
    resp = client.get("/api/v1/elements")
    assert resp.status_code == 200
    body = resp.json()
    assert "C" in body and "O" in body
    assert body["C"][0]["mass"] == pytest.approx(12.0)


def _fit_payload(**overrides):
    payload = {
        "formula": "C6",
        "label": "13C",
        "peaks": [
            {"mass": 72.0, "intensity": 100.0},
            {"mass": 73.00335, "intensity": 6.5},
        ],
        "step_percent": 1.0,
        "mass_tolerance": 0.05,
    }
    payload.update(overrides)
    return payload


def test_enrichment_fit_recovers_abundance(client):
    # Natural-abundance glucose envelope -> best 13C abundance near the
    # natural 1.07%, i.e. the 1.0% grid point.
    env = client.post("/api/v1/envelope", json={"formula": "C6H12O6", "top_n": 10})
    peaks = env.json()["peaks"]
    resp = client.post(
        "/api/v1/enrichment-fit",
        json=_fit_payload(formula="C6H12O6", peaks=peaks, step_percent=0.5),
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["best_abundance_percent"] == pytest.approx(1.0)
    assert body["search_points"] == 201
    assert body["mse"] >= 0.0
    assert body["pair_count"] == len(body["pairs"])
    # every observed peak appears exactly once in the pairing
    seen = [p["observed_mass"] for p in body["pairs"] if p["observed_mass"] is not None]
    assert len(seen) == len(peaks)


def test_enrichment_fit_perfect_match_zero_mse(client):
    # C6 at 50% 13C: binomial envelope, 7 isotopologue peaks.
    peaks = [
        {"mass": m, "intensity": i}
        for m, i in synthetic_peaks("C6", "C", 13, 0.5, top=7)
    ]
    resp = client.post(
        "/api/v1/enrichment-fit",
        json=_fit_payload(peaks=peaks, step_percent=0.5),
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["best_abundance_percent"] == pytest.approx(50.0)
    assert body["mse"] == pytest.approx(0.0, abs=1e-9)


def test_enrichment_fit_search_points_with_coarse_step(client):
    resp = client.post(
        "/api/v1/enrichment-fit", json=_fit_payload(step_percent=10.0)
    )
    assert resp.status_code == 200
    assert resp.json()["search_points"] == 11


def test_enrichment_fit_invalid_label(client):
    for label, position in [("C13", 0), ("13Zz", 2), ("20C", 0), ("[13C]", 0), ("", 0)]:
        resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(label=label))
        assert resp.status_code == 422, label
        detail = resp.json()["detail"]
        assert detail["error"] == "invalid_label"
        assert detail["position"] == position


def test_enrichment_fit_invalid_formula(client):
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(formula="C6Xx2"))
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_formula"


def test_enrichment_fit_empty_peaks(client):
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(peaks=[]))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_peaks"


def test_enrichment_fit_non_increasing_masses(client):
    peaks = [
        {"mass": 72.0, "intensity": 100.0},
        {"mass": 71.0, "intensity": 5.0},
    ]
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(peaks=peaks))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_peaks"
    assert detail["index"] == 1


def test_enrichment_fit_duplicate_masses(client):
    peaks = [
        {"mass": 72.0, "intensity": 100.0},
        {"mass": 72.0, "intensity": 5.0},
    ]
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(peaks=peaks))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_peaks"
    assert detail["index"] == 1
    assert "Duplicate" in detail["message"]


def test_enrichment_fit_parameter_bounds(client):
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(step_percent=0.4))
    assert resp.status_code == 422
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(step_percent=10.5))
    assert resp.status_code == 422
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(mass_tolerance=0.0))
    assert resp.status_code == 422
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(mass_tolerance=-0.1))
    assert resp.status_code == 422
    peaks = [{"mass": float(m), "intensity": 1.0} for m in range(201)]
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(peaks=peaks))
    assert resp.status_code == 422
    peaks = [{"mass": 72.0, "intensity": -1.0}]
    resp = client.post("/api/v1/enrichment-fit", json=_fit_payload(peaks=peaks))
    assert resp.status_code == 422


def _dual_payload(**overrides):
    payload = {
        "formula": "C4O2",
        "labels": [
            {"label": "13C", "step_percent": 10.0},
            {"label": "18O", "step_percent": 10.0},
        ],
        "peaks": [
            {"mass": m, "intensity": i}
            for m, i in synthetic_peaks_dual(
                "C4O2", {"C": (13, 0.5), "O": (18, 0.4)}
            )
        ],
        "mass_tolerance": 0.005,
    }
    payload.update(overrides)
    return payload


def test_dual_enrichment_fit_recovers_both_abundances(client):
    resp = client.post("/api/v1/dual-enrichment-fit", json=_dual_payload())
    assert resp.status_code == 200
    body = resp.json()
    assert body["labels"] == ["13C", "18O"]
    assert body["best_abundances_percent"] == pytest.approx([50.0, 40.0])
    assert body["mse"] == pytest.approx(0.0, abs=1e-9)
    assert body["search_points"] == 121
    assert body["pair_count"] == len(body["pairs"])
    # every observed peak appears exactly once in the pairing
    seen = [p["observed_mass"] for p in body["pairs"] if p["observed_mass"] is not None]
    assert len(seen) == len(body["pairs"])
    # pairs are returned in ascending mass order
    masses = [
        p["observed_mass"] if p["observed_mass"] is not None else p["predicted_mass"]
        for p in body["pairs"]
    ]
    assert masses == sorted(masses)


def test_dual_enrichment_fit_independent_step_sizes(client):
    payload = _dual_payload(
        labels=[
            {"label": "13C", "step_percent": 5.0},
            {"label": "18O", "step_percent": 10.0},
        ]
    )
    resp = client.post("/api/v1/dual-enrichment-fit", json=payload)
    assert resp.status_code == 200
    body = resp.json()
    assert body["search_points"] == 21 * 11
    assert body["best_abundances_percent"] == pytest.approx([50.0, 40.0])


def test_dual_enrichment_fit_same_element_labels_are_locatable_422(client):
    payload = _dual_payload(
        labels=[{"label": "13C"}, {"label": "12C"}],
    )
    resp = client.post("/api/v1/dual-enrichment-fit", json=payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "duplicate_label_element"
    # index locates the second (offending) label in the labels list
    assert detail["index"] == 1
    assert "'C'" in detail["message"]


def test_dual_enrichment_fit_malformed_label_still_reports_position(client):
    payload = _dual_payload(
        labels=[{"label": "13C"}, {"label": "20O"}],
    )
    resp = client.post("/api/v1/dual-enrichment-fit", json=payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_label"
    assert detail["position"] == 0


def test_dual_enrichment_fit_requires_exactly_two_labels(client):
    payload = _dual_payload()
    one = {**payload, "labels": [{"label": "13C"}]}
    resp = client.post("/api/v1/dual-enrichment-fit", json=one)
    assert resp.status_code == 422
    three = {
        **payload,
        "labels": [{"label": "13C"}, {"label": "18O"}, {"label": "15N"}],
    }
    resp = client.post("/api/v1/dual-enrichment-fit", json=three)
    assert resp.status_code == 422


def test_dual_enrichment_fit_invalid_formula(client):
    resp = client.post(
        "/api/v1/dual-enrichment-fit", json=_dual_payload(formula="C4Xx2")
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_formula"


def test_dual_enrichment_fit_invalid_peaks(client):
    resp = client.post(
        "/api/v1/dual-enrichment-fit", json=_dual_payload(peaks=[])
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_peaks"


def test_dual_enrichment_fit_parameter_bounds(client):
    bad_steps = [
        [{"label": "13C", "step_percent": 0.4}, {"label": "18O"}],
        [{"label": "13C"}, {"label": "18O", "step_percent": 10.5}],
    ]
    for labels in bad_steps:
        resp = client.post(
            "/api/v1/dual-enrichment-fit", json=_dual_payload(labels=labels)
        )
        assert resp.status_code == 422
    resp = client.post(
        "/api/v1/dual-enrichment-fit", json=_dual_payload(mass_tolerance=0.0)
    )
    assert resp.status_code == 422
    peaks = [{"mass": float(m), "intensity": 1.0} for m in range(201)]
    resp = client.post(
        "/api/v1/dual-enrichment-fit", json=_dual_payload(peaks=peaks)
    )
    assert resp.status_code == 422


def test_dual_enrichment_fit_tie_picks_lexicographic_minimum(client):
    # Neither C nor N appears in H2O: every grid point scores the same.
    peaks = [
        {"mass": m, "intensity": i}
        for m, i in synthetic_peaks("H2O", "C", 13, 0.5, top=5)
    ]
    payload = {
        "formula": "H2O",
        "labels": [{"label": "13C", "step_percent": 10.0}, {"label": "15N", "step_percent": 5.0}],
        "peaks": peaks,
        "mass_tolerance": 0.05,
    }
    resp = client.post("/api/v1/dual-enrichment-fit", json=payload)
    assert resp.status_code == 200
    body = resp.json()
    assert body["best_abundances_percent"] == [0.0, 0.0]
    assert body["search_points"] == 11 * 21
