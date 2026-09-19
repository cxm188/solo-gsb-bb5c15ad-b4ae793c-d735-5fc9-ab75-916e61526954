"""Isotope-enrichment fitting by grid search over candidate abundances.

For each candidate abundance of a single isotope label (e.g. 13C), the
labelled element's natural isotope distribution is replaced by one in which
the labelled isotope carries the candidate abundance while the remaining
isotopes keep their natural ratios. The molecular envelope is rebuilt with
the existing convolution machinery, normalised to base peak = 100 (as are
the observed peaks), and paired with the observed peaks by nearest mass
within a tolerance. The abundance minimising the intensity mean squared
error wins; ties go to the lower abundance.

A dual-label variant (:func:`fit_dual_enrichment`) grid-searches two labels
on *different* elements jointly: the abundance-independent envelope (every
other element) and each labelled element's per-grid-point polynomial are
each computed once and then combined on the 2-D grid. MSE ties are broken
by lexicographic (first label, second label) abundance order.
"""

from __future__ import annotations

from dataclasses import dataclass

from .envelope import Peak, convolve, label_mass_shift, merge_peaks, power_polynomial
from .isotopes import ISOTOPES, isotope_by_mass_number
from .parser import ParsedFormula


class PeakListError(ValueError):
    """Validation error for the observed peak list, carrying the 0-based
    index of the offending peak (``None`` when not index-specific)."""

    def __init__(self, message: str, index: int | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.index = index


class LabelPairError(ValueError):
    """Validation error for a pair of dual-fit labels, carrying the 0-based
    index of the offending label in the request's ``labels`` list."""

    def __init__(self, message: str, index: int) -> None:
        super().__init__(message)
        self.message = message
        self.index = index


def validate_observed_peaks(peaks: list[Peak]) -> None:
    """Require a non-empty list of strictly mass-increasing observed peaks."""
    if not peaks:
        raise PeakListError("Observed peak list must contain at least one peak")
    for i in range(1, len(peaks)):
        prev, cur = peaks[i - 1][0], peaks[i][0]
        if cur == prev:
            raise PeakListError(
                f"Duplicate peak mass {cur}; masses must be strictly increasing", i
            )
        if cur < prev:
            raise PeakListError(
                f"Peak masses must be sorted in increasing order; "
                f"mass {cur} at index {i} follows {prev}",
                i,
            )
    if max(i for _, i in peaks) <= 0:
        raise PeakListError(
            "At least one observed peak must have a positive intensity"
        )


def enriched_distribution(
    symbol: str, mass_number: int, fraction: float
) -> tuple[tuple[float, float], ...]:
    """Isotope distribution of ``symbol`` with the ``mass_number`` isotope
    enriched to ``fraction``; the other natural isotopes keep their natural
    abundance ratios and share the remaining ``1 - fraction``. Zero-abundance
    entries are dropped so the result keeps the positive-abundance invariant
    of the convolution machinery."""
    natural = ISOTOPES[symbol]
    labelled_mass = isotope_by_mass_number(symbol, mass_number)[0]
    others = [(m, a) for m, a in natural if round(m) != mass_number]
    if not others:
        # Single-isotope element: enrichment cannot change the distribution.
        return ((labelled_mass, 1.0),)
    scale = (1.0 - fraction) / sum(a for _, a in others)
    dist = [(m, a * scale) for m, a in others if a * scale > 0.0]
    if fraction > 0.0:
        dist.append((labelled_mass, fraction))
    dist.sort()
    return tuple(dist)


def abundance_grid(step_percent: float) -> list[float]:
    """Candidate abundances in percent: 0, step, 2*step, ..., 100.

    The 100% endpoint is always included (appended when the step does not
    land on it exactly).
    """
    n = int(100.0 / step_percent + 1e-9)
    grid = [i * step_percent for i in range(n + 1)]
    if grid[-1] < 100.0 - 1e-9:
        grid.append(100.0)
    return grid


def normalise_base_peak(peaks: list[Peak]) -> list[Peak]:
    """Rescale intensities so the base (most intense) peak is 100."""
    base = max(i for _, i in peaks)
    return [(m, 100.0 * i / base) for m, i in peaks]


def pair_and_score(
    observed: list[Peak], predicted: list[Peak], tolerance: float
) -> tuple[list[tuple[int | None, int | None]], float]:
    """Pair observed and predicted peaks by nearest mass within ``tolerance``.

    Greedy global matching: all within-tolerance (observed, predicted)
    candidate pairs are considered, smallest mass difference first; each
    peak takes part in at most one pair. Unmatched peaks are paired with a
    zero-intensity placeholder. Returns the pairs as
    ``(observed_index | None, predicted_index | None)`` plus the intensity
    mean squared error over all pairs.
    """
    candidates: list[tuple[float, int, int]] = []
    j0 = 0
    n_pred = len(predicted)
    for i, (m_obs, _) in enumerate(observed):
        lo, hi = m_obs - tolerance, m_obs + tolerance
        while j0 < n_pred and predicted[j0][0] < lo:
            j0 += 1
        j = j0
        while j < n_pred and predicted[j][0] <= hi:
            candidates.append((abs(predicted[j][0] - m_obs), i, j))
            j += 1
    candidates.sort()

    matched_obs: dict[int, int] = {}
    used_pred: set[int] = set()
    for _, i, j in candidates:
        if i not in matched_obs and j not in used_pred:
            matched_obs[i] = j
            used_pred.add(j)

    pairs: list[tuple[int | None, int | None]] = []
    squared_error = 0.0
    for i, (_, i_obs) in enumerate(observed):
        j = matched_obs.get(i)
        if j is None:
            pairs.append((i, None))
            squared_error += i_obs * i_obs
        else:
            pairs.append((i, j))
            diff = i_obs - predicted[j][1]
            squared_error += diff * diff
    for j, (_, i_pred) in enumerate(predicted):
        if j not in used_pred:
            pairs.append((None, j))
            squared_error += i_pred * i_pred
    return pairs, squared_error / len(pairs)


@dataclass
class FitResult:
    best_abundance_percent: float
    mse: float
    search_points: int
    # (observed_index | None, predicted_index | None) pairs at the best fit
    pairs: list[tuple[int | None, int | None]]
    predicted: list[Peak]  # normalised envelope at the best abundance


def _envelope_excluding(parsed: ParsedFormula, symbol: str, tolerance: float, prune: float) -> list[Peak]:
    """Natural-abundance envelope of every element except ``symbol``
    (computed once and reused across all candidate abundances)."""
    envelope: list[Peak] = [(0.0, 1.0)]
    for sym, count in parsed.composition.items():
        if sym == symbol:
            continue
        element_poly = power_polynomial(list(ISOTOPES[sym]), count, tolerance, prune)
        envelope = convolve(envelope, element_poly, tolerance, prune)
    return envelope


def _candidate_envelope(
    parsed: ParsedFormula,
    symbol: str,
    distribution: tuple[tuple[float, float], ...],
    others_envelope: list[Peak],
    tolerance: float,
    prune: float,
) -> list[Peak]:
    """Full envelope with ``symbol``'s distribution replaced by ``distribution``."""
    envelope = others_envelope
    count = parsed.composition.get(symbol, 0)
    if count:
        element_poly = power_polynomial(list(distribution), count, tolerance, prune)
        # Convolution is commutative; iterating the smaller operand first
        # gives the merge sort fewer, longer runs to fuse.
        if len(element_poly) <= len(envelope):
            envelope = convolve(element_poly, envelope, tolerance, prune)
        else:
            envelope = convolve(envelope, element_poly, tolerance, prune)
    shift = label_mass_shift(parsed)
    if shift:
        envelope = [(m + shift, i) for m, i in envelope]
    return merge_peaks(envelope, tolerance)


def fit_enrichment(
    parsed: ParsedFormula,
    symbol: str,
    mass_number: int,
    observed: list[Peak],
    step_percent: float,
    mass_tolerance: float,
    merge_tolerance: float,
    prune: float,
) -> FitResult:
    """Grid-search the label abundance (0–100%) best matching the observed
    peaks. ``observed`` must be sorted by strictly increasing mass with
    non-negative intensities, at least one of them positive."""
    obs = normalise_base_peak(observed)
    others = _envelope_excluding(parsed, symbol, merge_tolerance, prune)
    grid = abundance_grid(step_percent)

    best: FitResult | None = None
    for pct in grid:
        dist = enriched_distribution(symbol, mass_number, pct / 100.0)
        env = _candidate_envelope(parsed, symbol, dist, others, merge_tolerance, prune)
        pred = normalise_base_peak(env)
        pairs, mse = pair_and_score(obs, pred, mass_tolerance)
        # Strict improvement only: on a tie the lower abundance already
        # visited is kept.
        if best is None or mse < best.mse:
            best = FitResult(pct, mse, len(grid), pairs, pred)
    # The grid always spans 0–100%, so a best candidate always exists.
    assert best is not None
    return best


@dataclass
class DualFitResult:
    best_abundances_percent: tuple[float, float]
    mse: float
    search_points: int
    # (observed_index | None, predicted_index | None) pairs at the best fit
    pairs: list[tuple[int | None, int | None]]
    predicted: list[Peak]  # normalised envelope at the best abundances


def _envelope_excluding_two(
    parsed: ParsedFormula, symbol_a: str, symbol_b: str, tolerance: float, prune: float
) -> list[Peak]:
    """Natural-abundance envelope of every element except the two labelled
    elements (computed once and reused over the whole 2-D grid)."""
    envelope: list[Peak] = [(0.0, 1.0)]
    for sym, count in parsed.composition.items():
        if sym == symbol_a or sym == symbol_b:
            continue
        element_poly = power_polynomial(list(ISOTOPES[sym]), count, tolerance, prune)
        envelope = convolve(envelope, element_poly, tolerance, prune)
    return envelope


def _element_grid_polynomials(
    parsed: ParsedFormula,
    symbol: str,
    mass_number: int,
    grid: list[float],
    tolerance: float,
    prune: float,
) -> list[list[Peak]]:
    """Raised isotope polynomial of ``symbol`` for every grid abundance
    (computed once per grid point and reused across the other dimension)."""
    count = parsed.composition.get(symbol, 0)
    if not count:
        # Element absent from the formula: a unit factor for every point.
        return [[(0.0, 1.0)] for _ in grid]
    return [
        power_polynomial(
            list(enriched_distribution(symbol, mass_number, pct / 100.0)),
            count,
            tolerance,
            prune,
        )
        for pct in grid
    ]


def fit_dual_enrichment(
    parsed: ParsedFormula,
    label_a: tuple[str, int],
    label_b: tuple[str, int],
    observed: list[Peak],
    step_a_percent: float,
    step_b_percent: float,
    mass_tolerance: float,
    merge_tolerance: float,
    prune: float,
) -> DualFitResult:
    """Grid-search two label abundances (each 0–100%) jointly. The two
    labels must belong to **different** elements; same-element pairs raise
    :class:`LabelPairError`.

    ``observed`` must be sorted by strictly increasing mass with
    non-negative intensities, at least one of them positive. On an MSE tie
    the lexicographically smallest ``(abundance_a, abundance_b)`` wins;
    iterating ``a`` in the outer loop and ``b`` in the inner loop visits
    points in exactly that order.
    """
    symbol_a, mass_number_a = label_a
    symbol_b, mass_number_b = label_b
    if symbol_a == symbol_b:
        raise LabelPairError(
            f"The two labels must target different elements; both target "
            f"element '{symbol_a}'",
            1,
        )

    obs = normalise_base_peak(observed)
    rest = _envelope_excluding_two(
        parsed, symbol_a, symbol_b, merge_tolerance, prune
    )
    grid_a = abundance_grid(step_a_percent)
    grid_b = abundance_grid(step_b_percent)
    polys_a = _element_grid_polynomials(
        parsed, symbol_a, mass_number_a, grid_a, merge_tolerance, prune
    )
    polys_b = _element_grid_polynomials(
        parsed, symbol_b, mass_number_b, grid_b, merge_tolerance, prune
    )
    shift = label_mass_shift(parsed)

    best: DualFitResult | None = None
    for ia, pct_a in enumerate(grid_a):
        for ib, pct_b in enumerate(grid_b):
            # rest is independent of both abundances; each element
            # polynomial is reused across the other grid dimension.
            envelope = convolve(rest, polys_a[ia], merge_tolerance, prune)
            envelope = convolve(envelope, polys_b[ib], merge_tolerance, prune)
            if shift:
                envelope = [(m + shift, i) for m, i in envelope]
            envelope = merge_peaks(envelope, merge_tolerance)
            pred = normalise_base_peak(envelope)
            pairs, mse = pair_and_score(obs, pred, mass_tolerance)
            # Strict improvement only: ties keep the first point visited,
            # i.e. the lexicographically smallest abundance pair.
            if best is None or mse < best.mse:
                best = DualFitResult((pct_a, pct_b), mse, len(grid_a) * len(grid_b), pairs, pred)
    assert best is not None
    return best
