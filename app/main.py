"""FastAPI application: mass-spectrometry isotope-envelope precompute API."""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from . import __version__, db
from .config import settings
from .enrichment import (
    PeakListError,
    fit_enrichment,
    normalise_base_peak,
    validate_observed_peaks,
)
from .envelope import isotope_envelope, monoisotopic_mass, top_peaks
from .isotopes import ISOTOPES
from .parser import FormulaError, LabelError, parse_formula, parse_label


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Dependencies are installed at image build time; the SQLite schema is
    # created here so a fresh container is ready to serve immediately.
    db.init_db(settings.db_path)
    yield
    db.close_db()


app = FastAPI(
    title="MS Precompute API",
    version=__version__,
    description=(
        "Computes natural-abundance isotope envelopes (mass/intensity "
        "polynomial convolution) and monoisotopic masses for molecular "
        "formulas, with optional pure-isotope labels such as [13C]; also "
        "estimates single-isotope enrichment by fitting candidate "
        "abundances against observed peaks."
    ),
    lifespan=lifespan,
)


class EnvelopeRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    top_n: int = Field(
        settings.top_n_default,
        ge=1,
        le=settings.top_n_max,
        description="Number of most intense peaks to return",
    )


class PeakOut(BaseModel):
    mass: float
    intensity: float


class EnvelopeResponse(BaseModel):
    formula: str
    composition: dict[str, int]
    labels: dict[str, int]
    monoisotopic_mass: float
    base_peak_intensity: float
    peak_count: int
    peaks: list[PeakOut]
    cached: bool


class ObservedPeak(BaseModel):
    mass: float = Field(..., description="Observed peak mass (Da)")
    intensity: float = Field(..., ge=0, description="Observed peak intensity")


class EnrichmentFitRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    label: str = Field(
        ...,
        description="Single isotope label to estimate, e.g. '13C'",
        examples=["13C"],
    )
    peaks: list[ObservedPeak] = Field(
        ...,
        max_length=200,
        description="Observed peaks, sorted by strictly increasing mass",
    )
    step_percent: float = Field(
        1.0,
        ge=0.5,
        le=10.0,
        description="Abundance grid step in percentage points (0.5-10)",
    )
    mass_tolerance: float = Field(
        ...,
        gt=0,
        description="Mass tolerance (Da) for pairing observed/predicted peaks",
    )


class PeakPairOut(BaseModel):
    observed_mass: float | None
    observed_intensity: float | None
    predicted_mass: float | None
    predicted_intensity: float | None


class EnrichmentFitResponse(BaseModel):
    formula: str
    label: str
    best_abundance_percent: float
    mse: float
    search_points: int
    pair_count: int
    pairs: list[PeakPairOut]


@app.exception_handler(FormulaError)
async def formula_error_handler(request: Request, exc: FormulaError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_formula",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.exception_handler(LabelError)
async def label_error_handler(request: Request, exc: LabelError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_label",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.exception_handler(PeakListError)
async def peak_list_error_handler(request: Request, exc: PeakListError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_peaks",
                "message": exc.message,
                "index": exc.index,
            }
        },
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@app.get("/api/v1/elements")
def elements() -> dict:
    """List supported element symbols and their natural isotopes."""
    return {
        symbol: [{"mass": m, "abundance": a} for m, a in isotopes]
        for symbol, isotopes in sorted(ISOTOPES.items())
    }


@app.post("/api/v1/envelope", response_model=EnvelopeResponse)
def envelope(req: EnvelopeRequest) -> dict:
    """Compute the isotope envelope and monoisotopic mass of a formula."""
    try:
        parsed = parse_formula(req.formula)
    except FormulaError as exc:
        db.log_request(req.formula, req.top_n, ok=False, message=exc.message)
        raise

    cached = db.get_cached(req.formula, req.top_n)
    if cached is not None:
        db.log_request(req.formula, req.top_n, ok=True, message="cache hit")
        return {**cached, "cached": True}

    env = isotope_envelope(parsed, settings.merge_tolerance, settings.prune_threshold)
    peaks = top_peaks(env, req.top_n)
    result = {
        "formula": req.formula,
        "composition": parsed.composition,
        "labels": {f"[{mass_no}{symbol}]": count for (symbol, mass_no), count in parsed.labels.items()},
        "monoisotopic_mass": round(monoisotopic_mass(parsed), 6),
        "base_peak_intensity": 100.0,
        "peak_count": len(peaks),
        "peaks": [
            {"mass": round(m, 6), "intensity": round(i, 6)} for m, i in peaks
        ],
        "cached": False,
    }
    db.put_cache(req.formula, req.top_n, {k: v for k, v in result.items() if k != "cached"})
    db.log_request(req.formula, req.top_n, ok=True, message="computed")
    return result


@app.post("/api/v1/enrichment-fit", response_model=EnrichmentFitResponse)
def enrichment_fit(req: EnrichmentFitRequest) -> dict:
    """Estimate the abundance of a single isotope label by fitting the
    computed envelope against observed peaks over a 0-100% grid."""
    parsed = parse_formula(req.formula)
    symbol, mass_number = parse_label(req.label)

    observed = [(p.mass, p.intensity) for p in req.peaks]
    validate_observed_peaks(observed)

    fit = fit_enrichment(
        parsed,
        symbol,
        mass_number,
        observed,
        step_percent=req.step_percent,
        mass_tolerance=req.mass_tolerance,
        merge_tolerance=settings.merge_tolerance,
        prune=settings.prune_threshold,
    )

    obs = normalise_base_peak(observed)
    pairs = [
        {
            "observed_mass": round(obs[i][0], 6) if i is not None else None,
            "observed_intensity": round(obs[i][1], 6) if i is not None else None,
            "predicted_mass": round(fit.predicted[j][0], 6) if j is not None else None,
            "predicted_intensity": round(fit.predicted[j][1], 6) if j is not None else None,
        }
        for i, j in fit.pairs
    ]
    pairs.sort(
        key=lambda p: p["observed_mass"]
        if p["observed_mass"] is not None
        else p["predicted_mass"]
    )
    return {
        "formula": req.formula,
        "label": req.label,
        "best_abundance_percent": round(fit.best_abundance_percent, 6),
        "mse": fit.mse,
        "search_points": fit.search_points,
        "pair_count": len(pairs),
        "pairs": pairs,
    }
