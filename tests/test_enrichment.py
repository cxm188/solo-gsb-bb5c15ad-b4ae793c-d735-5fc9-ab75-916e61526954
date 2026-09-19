import pytest

from app.enrichment import (
    LabelPairError,
    PeakListError,
    abundance_grid,
    enriched_distribution,
    fit_dual_enrichment,
    fit_enrichment,
    normalise_base_peak,
    pair_and_score,
    validate_observed_peaks,
)
from app.envelope import convolve, merge_peaks, power_polynomial, top_peaks
from app.isotopes import ISOTOPES
from app.parser import parse_formula

TOL = 0.0001
PRUNE = 1e-12


def synthetic_peaks(formula, symbol, mass_number, fraction, top=10):
    """Observed-style peaks for ``formula`` with ``symbol``'s isotope
    ``mass_number`` enriched to ``fraction`` (all other elements natural)."""
    parsed = parse_formula(formula)
    dist = enriched_distribution(symbol, mass_number, fraction)
    env = [(0.0, 1.0)]
    for sym, count in parsed.composition.items():
        base = list(dist) if sym == symbol else list(ISOTOPES[sym])
        env = convolve(env, power_polynomial(base, count, TOL, PRUNE), TOL, PRUNE)
    return top_peaks(merge_peaks(env, TOL), top)


def synthetic_peaks_dual(formula, enrichments, top=10):
    """Observed-style peaks for ``formula`` with two elements enriched.

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
        env = convolve(env, power_polynomial(base, count, TOL, PRUNE), TOL, PRUNE)
    return top_peaks(merge_peaks(env, TOL), top)


def test_enriched_distribution_sums_to_one():
    for fraction in (0.0, 0.0107, 0.25, 0.5, 1.0):
        dist = enriched_distribution("C", 13, fraction)
        assert sum(a for _, a in dist) == pytest.approx(1.0)
        assert all(a > 0.0 for _, a in dist)
    # endpoints collapse to a single isotope
    assert enriched_distribution("C", 13, 0.0) == ((12.0, 1.0),)
    assert enriched_distribution("C", 13, 1.0) == ((13.00335483507, 1.0),)


def test_enriched_distribution_at_natural_abundance_matches_natural():
    dist = enriched_distribution("C", 13, 0.0107)
    assert [m for m, _ in dist] == [m for m, _ in ISOTOPES["C"]]
    assert [a for _, a in dist] == pytest.approx([a for _, a in ISOTOPES["C"]])


def test_enriched_distribution_keeps_other_isotope_ratios():
    # 18O enriched to 50%: 16O and 17O keep their natural ratio.
    dist = enriched_distribution("O", 18, 0.5)
    by_mass = {round(m): a for m, a in dist}
    assert by_mass[18] == pytest.approx(0.5)
    natural_16, natural_17 = 0.99757, 0.00038
    assert by_mass[16] / by_mass[17] == pytest.approx(natural_16 / natural_17)


def test_enriched_distribution_single_isotope_element():
    # Fluorine has only 19F: enrichment cannot change anything.
    assert enriched_distribution("F", 19, 0.3) == ((18.99840316273, 1.0),)


def test_abundance_grid_endpoints_and_counts():
    grid = abundance_grid(0.5)
    assert len(grid) == 201 and grid[0] == 0.0 and grid[-1] == 100.0
    assert abundance_grid(10.0) == [i * 10.0 for i in range(11)]
    grid = abundance_grid(3.0)
    assert len(grid) == 35 and grid[-2] == 99.0 and grid[-1] == 100.0
    grid = abundance_grid(7.0)
    assert grid[-2] == 98.0 and grid[-1] == 100.0


def test_validate_observed_peaks_empty():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed_peaks([])
    assert excinfo.value.index is None


def test_validate_observed_peaks_duplicate_mass():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed_peaks([(100.0, 1.0), (100.0, 2.0)])
    assert excinfo.value.index == 1
    assert "Duplicate" in excinfo.value.message


def test_validate_observed_peaks_decreasing_mass():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed_peaks([(100.0, 1.0), (101.0, 2.0), (100.5, 3.0)])
    assert excinfo.value.index == 2
    assert "increasing" in excinfo.value.message


def test_validate_observed_peaks_all_zero_intensity():
    with pytest.raises(PeakListError):
        validate_observed_peaks([(100.0, 0.0), (101.0, 0.0)])


def test_normalise_base_peak():
    assert normalise_base_peak([(10.0, 50.0), (11.0, 100.0)]) == [
        (10.0, 50.0),
        (11.0, 100.0),
    ]


def test_pair_and_score_exact_match():
    obs = [(100.0, 100.0), (101.0, 50.0)]
    pairs, mse = pair_and_score(obs, list(obs), 0.01)
    assert mse == pytest.approx(0.0)
    assert pairs == [(0, 0), (1, 1)]


def test_pair_and_score_unmatched_both_sides():
    obs = [(100.0, 80.0)]
    pred = [(100.0, 50.0), (200.0, 30.0)]
    pairs, mse = pair_and_score(obs, pred, 0.01)
    assert pairs == [(0, 0), (None, 1)]
    # (80-50)^2 + (0-30)^2 over 2 pairs
    assert mse == pytest.approx((30.0**2 + 30.0**2) / 2)


def test_pair_and_score_respects_tolerance_and_nearest():
    obs = [(100.0, 100.0)]
    pred = [(100.04, 1.0), (100.05, 2.0), (101.0, 3.0)]
    pairs, _ = pair_and_score(obs, pred, 0.06)
    # nearest within tolerance (100.04) wins; the rest are unmatched
    assert pairs[0] == (0, 0)
    assert (None, 1) in pairs and (None, 2) in pairs
    # nothing within tolerance of 100.0 here:
    pairs, _ = pair_and_score(obs, [(105.0, 1.0)], 0.06)
    assert pairs == [(0, None), (None, 0)]


def test_pair_and_score_one_to_one_matching():
    # Two observed peaks compete for one predicted peak: nearest wins.
    obs = [(100.00, 10.0), (100.03, 20.0)]
    pred = [(100.01, 15.0)]
    pairs, _ = pair_and_score(obs, pred, 0.05)
    assert pairs == [(0, 0), (1, None)]


def test_fit_enrichment_recovers_abundance():
    observed = synthetic_peaks("C6", "C", 13, 0.5, top=7)
    fit = fit_enrichment(
        parse_formula("C6"), "C", 13, observed,
        step_percent=0.5, mass_tolerance=0.05,
        merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundance_percent == pytest.approx(50.0)
    assert fit.mse == pytest.approx(0.0, abs=1e-9)
    assert fit.search_points == 201


def test_fit_enrichment_coarse_grid_picks_nearest():
    observed = synthetic_peaks("C6", "C", 13, 0.42, top=7)
    fit = fit_enrichment(
        parse_formula("C6"), "C", 13, observed,
        step_percent=10.0, mass_tolerance=0.05,
        merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundance_percent == pytest.approx(40.0)
    assert fit.search_points == 11


def test_fit_enrichment_tie_picks_lower_abundance():
    # H2O has no carbon: every candidate gives the same envelope.
    observed = [(18.01056468, 100.0), (19.01394, 0.05)]
    fit = fit_enrichment(
        parse_formula("H2O"), "C", 13, observed,
        step_percent=5.0, mass_tolerance=0.05,
        merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundance_percent == 0.0


def test_fit_dual_enrichment_recovers_both_abundances():
    observed = synthetic_peaks_dual(
        "C4O2", {"C": (13, 0.5), "O": (18, 0.4)}, top=200
    )
    fit = fit_dual_enrichment(
        parse_formula("C4O2"), ("C", 13), ("O", 18), observed,
        step_a_percent=10.0, step_b_percent=10.0,
        mass_tolerance=0.005, merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundances_percent == pytest.approx((50.0, 40.0))
    assert fit.mse == pytest.approx(0.0, abs=1e-9)
    assert fit.search_points == 11 * 11


def test_fit_dual_enrichment_coarse_grids_pick_nearest():
    observed = synthetic_peaks_dual(
        "C4O2", {"C": (13, 0.42), "O": (18, 0.62)}, top=200
    )
    fit = fit_dual_enrichment(
        parse_formula("C4O2"), ("C", 13), ("O", 18), observed,
        step_a_percent=10.0, step_b_percent=10.0,
        mass_tolerance=0.005, merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundances_percent == pytest.approx((40.0, 60.0))
    assert fit.search_points == 121


def test_fit_dual_enrichment_tie_picks_lexicographic_minimum():
    # Formula contains neither enriched element: every grid point scores
    # identically, so (0, 0) must win.
    observed = synthetic_peaks("H2O", "C", 13, 0.5, top=5)
    fit = fit_dual_enrichment(
        parse_formula("H2O"), ("C", 13), ("N", 15), observed,
        step_a_percent=5.0, step_b_percent=10.0,
        mass_tolerance=0.05, merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundances_percent == (0.0, 0.0)
    assert fit.search_points == 21 * 11


def test_fit_dual_enrichment_one_element_absent_from_formula():
    # 13C on C6 fitted jointly with 15N on a nitrogen-free molecule:
    # the 15N dimension is degenerate, so its lexicographically smallest
    # value (0) is taken while 13C is recovered.
    observed = synthetic_peaks("C6", "C", 13, 0.5, top=7)
    fit = fit_dual_enrichment(
        parse_formula("C6"), ("C", 13), ("N", 15), observed,
        step_a_percent=0.5, step_b_percent=5.0,
        mass_tolerance=0.005, merge_tolerance=TOL, prune=PRUNE,
    )
    assert fit.best_abundances_percent[0] == pytest.approx(50.0)
    assert fit.best_abundances_percent[1] == 0.0
    assert fit.search_points == 201 * 21


def test_fit_dual_enrichment_rejects_same_element_labels():
    observed = synthetic_peaks("C6", "C", 13, 0.5, top=7)
    with pytest.raises(LabelPairError) as excinfo:
        fit_dual_enrichment(
            parse_formula("C6"), ("C", 12), ("C", 13), observed,
            step_a_percent=5.0, step_b_percent=5.0,
            mass_tolerance=0.05, merge_tolerance=TOL, prune=PRUNE,
        )
    assert excinfo.value.index == 1
    assert "'C'" in excinfo.value.message


def test_fit_dual_enrichment_matches_single_fit_when_second_is_natural():
    # With the second element pinned at its natural isotope (15N at 0% on a
    # nitrogen-free molecule) the best 13C value matches the single fit.
    observed = synthetic_peaks("C6", "C", 13, 0.5, top=7)
    single = fit_enrichment(
        parse_formula("C6"), "C", 13, observed,
        step_percent=5.0, mass_tolerance=0.005,
        merge_tolerance=TOL, prune=PRUNE,
    )
    dual = fit_dual_enrichment(
        parse_formula("C6"), ("C", 13), ("N", 15), observed,
        step_a_percent=5.0, step_b_percent=50.0,
        mass_tolerance=0.005, merge_tolerance=TOL, prune=PRUNE,
    )
    assert dual.best_abundances_percent[0] == single.best_abundance_percent
