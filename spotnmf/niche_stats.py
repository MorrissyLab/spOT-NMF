"""Statistically principled niche inference.

This module replaces the two fixed cutoffs that used to define spatial niches
(a hard-coded 90th-percentile usage threshold and a 0.199 co-occurrence
threshold) with data-driven presence calls and an explicit null model.

The three problems with the fixed-threshold approach were:

1. **The two thresholds were coupled.** Zeroing every program below its own
   90th percentile forces *every* program to be present in ~10% of spots. Under
   independence the conditional probability ``P(B | A)`` is then ~0.10, so an
   edge cutoff of 0.199 silently means "keep pairs with ~2x the co-occurrence
   expected by chance". Move the quantile to 0.95 and the same 0.199 becomes a
   4x cutoff. The intent (2-fold enrichment) was reasonable; encoding it as a
   raw probability cutoff was not. :func:`cooccurrence_stats` reports
   ``log2_oe`` -- log2 observed/expected -- which is invariant to the presence
   calling rule, so the enrichment target is stated directly.

2. **Fixed-rate binarization flattens real differences in program extent.** A
   focal program truly occupying 2% of a section is padded out to 10% with
   background spots; a broad stromal program covering 40% is truncated to its
   hottest quarter. :func:`call_presence` calls presence per program from the
   shape of its own usage distribution instead.

3. **Spots are not independent observations.** Spatial autocorrelation means the
   effective sample size is far below the spot count, so a Fisher/chi-square test
   on raw spots is badly anticonservative. :func:`add_permutation_pvalues`
   compares the observed statistic against a null that preserves each program's
   own spatial structure and destroys only the alignment *between* programs.

Typical use::

    B, presence_info = call_presence(usage, method="otsu")
    coords = coords_from_index(usage.index)
    stats = cooccurrence_stats(B, sample="CRC")
    stats = add_permutation_pvalues(stats, B, coords=coords, null="torus")
    edges = select_edges(stats, fdr=0.05, min_log2_oe=1.0)
"""

from __future__ import annotations

import re
import warnings
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

__all__ = [
    "otsu_threshold",
    "call_presence",
    "coords_from_index",
    "cooccurrence_stats",
    "add_permutation_pvalues",
    "select_edges",
    "stability_selection",
    "parameter_sweep",
    "clr_transform",
    "compositional_check",
]

# Presence calls outside this range are almost certainly a failed threshold
# search rather than biology, so they are clamped back to the boundary quantile.
_MIN_PREVALENCE = 0.005
_MAX_PREVALENCE = 0.60


# --------------------------------------------------------------------------
# 1. Presence calling
# --------------------------------------------------------------------------

def otsu_threshold(x: np.ndarray, n_candidates: int = 256) -> float:
    """Find the usage threshold that maximizes between-class variance (Otsu).

    Otsu's method splits a distribution into two classes without assuming a
    parametric form, which suits NMF usages: they are zero-inflated and
    right-skewed, so a background mode and a "program is really here" mode are
    typically visible but not Gaussian. Candidate thresholds are drawn at evenly
    spaced quantiles rather than evenly spaced values so the search is invariant
    to the scale of the usage column.

    Args:
        x (numpy.ndarray): 1-D array of usage values for one program.
        n_candidates (int): Number of candidate thresholds to evaluate.
            Defaults to 256.

    Returns:
        float: The threshold maximizing between-class variance. Values strictly
        greater than the threshold are considered present. Returns ``inf`` when
        the input is constant (no split is possible).
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or np.allclose(x.max(), x.min()):
        return float("inf")

    xs = np.sort(x)
    n = xs.size

    # Candidate cut positions, evenly spaced in rank, strictly interior so both
    # classes are non-empty.
    idx = np.unique(
        np.linspace(1, n - 1, num=min(n_candidates, max(n - 1, 1)), dtype=int)
    )
    if idx.size == 0:
        return float("inf")

    csum = np.cumsum(xs)
    total = csum[-1]

    w0 = idx / n                      # weight of the "absent" class
    w1 = 1.0 - w0
    mu0 = csum[idx - 1] / idx         # mean of the "absent" class
    mu1 = (total - csum[idx - 1]) / (n - idx)

    between = w0 * w1 * (mu0 - mu1) ** 2
    best = idx[int(np.argmax(between))]
    # Threshold sits between the two flanking order statistics.
    return float(0.5 * (xs[best - 1] + xs[best]))


def _mixture_threshold(x: np.ndarray, seed: int = 0, max_weight: float = 0.40) -> float:
    """Threshold from a two-component Gaussian mixture on log usage.

    A background/present mixture is only meaningful when the two components are
    actually separable. Consensus spOT-NMF usages are continuous and unimodal on
    the log scale -- a right-skewed bulk with a long tail and, on the CRC test
    matrix, no exact zeros at all -- so a two-component fit tends to bisect the
    bulk rather than isolate the tail, assigning half the section to "present".

    Rather than return that silently, the fit is rejected when the high
    component claims more than ``max_weight`` of the mass, and the call falls
    back to :func:`otsu_threshold`.

    Args:
        x (numpy.ndarray): 1-D array of usage values for one program.
        seed (int): Random seed for the mixture fit. Defaults to 0.
        max_weight (float): Reject the fit if the high component carries more
            than this fraction of the mass. Defaults to 0.40.

    Returns:
        float: Usage threshold above which a spot counts as present.
    """
    from sklearn.mixture import GaussianMixture

    x = np.asarray(x, dtype=float)
    nz = x[x > 0]
    if nz.size < 20 or np.allclose(nz.max(), nz.min()):
        return otsu_threshold(x)

    lx = np.log(nz).reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=seed, n_init=3).fit(lx)
    hi = int(np.argmax(gm.means_.ravel()))

    # Separation checks: the "present" component must be a genuine minority mode,
    # and the two means must be at least one pooled SD apart.
    sds = np.sqrt(gm.covariances_.ravel())
    means = gm.means_.ravel()
    separation = abs(means[1] - means[0]) / max(np.sqrt((sds ** 2).mean()), 1e-12)
    if gm.weights_[hi] > max_weight or separation < 1.0:
        return otsu_threshold(x)

    order = np.argsort(lx.ravel())
    post = gm.predict_proba(lx[order])[:, hi]
    crossing = np.flatnonzero(post > 0.5)
    if crossing.size == 0:
        return otsu_threshold(x)
    return float(np.exp(lx.ravel()[order][crossing[0]]))


def call_presence(
    usage: pd.DataFrame,
    method: str = "otsu",
    quantile: float = 0.90,
    min_prevalence: float = _MIN_PREVALENCE,
    max_prevalence: float = _MAX_PREVALENCE,
    seed: int = 0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Binarize a usage matrix into per-spot program presence calls.

    Each program is thresholded independently using its own usage distribution,
    so programs are free to differ in spatial extent. The legacy ``"quantile"``
    method is retained -- unchanged in behaviour -- so the fixed-threshold
    results remain reproducible and can be included in a robustness sweep.

    Args:
        usage (pandas.DataFrame): Spots (rows) by programs (columns) usage
            matrix.
        method (str): Presence calling rule. ``"otsu"`` maximizes between-class
            variance (default, no distributional assumption); ``"mixture"`` fits
            a two-component Gaussian mixture to log non-zero usage;
            ``"quantile"`` reproduces the legacy fixed-quantile cutoff.
        quantile (float): Quantile used when ``method="quantile"``. Defaults to
            0.90, the value that was previously hard-coded.
        min_prevalence (float): Presence calls yielding a smaller fraction of
            spots are clamped back to this fraction. Defaults to 0.005.
        max_prevalence (float): Presence calls yielding a larger fraction of
            spots are clamped back to this fraction. Defaults to 0.60.
        seed (int): Seed for the mixture fit. Defaults to 0.

    Returns:
        Tuple[pandas.DataFrame, pandas.DataFrame]: ``(B, info)`` where ``B`` is a
        boolean spots-by-programs presence matrix and ``info`` records, per
        program, the chosen ``threshold``, the realised ``prevalence``, and
        whether the call was ``clamped``.

    Raises:
        ValueError: If ``method`` is not one of the supported rules.
    """
    if method not in {"otsu", "mixture", "quantile"}:
        raise ValueError(
            f"Unknown presence method {method!r}; expected 'otsu', 'mixture' or 'quantile'."
        )

    thresholds, prevalences, clamped = {}, {}, {}
    presence = np.empty((len(usage), usage.shape[1]), dtype=bool)

    for j, col in enumerate(usage.columns):
        x = usage[col].to_numpy(dtype=float)

        if method == "quantile":
            thr = float(np.quantile(x, quantile))
            # The legacy rule zeroed values *below* the quantile and then tested
            # "> 0", so a spot exactly at the threshold counts as present.
            present = (x >= thr) & (x > 0)
        else:
            thr = otsu_threshold(x) if method == "otsu" else _mixture_threshold(x, seed=seed)
            present = x > thr

        prev = float(present.mean())
        was_clamped = False

        # A degenerate threshold search (all-present or all-absent) is not
        # informative; fall back to the nearest admissible quantile.
        if prev < min_prevalence or prev > max_prevalence:
            target = min_prevalence if prev < min_prevalence else max_prevalence
            thr = float(np.quantile(x, 1.0 - target))
            present = (x > thr) & (x > 0)
            prev = float(present.mean())
            was_clamped = True

        thresholds[col] = thr
        prevalences[col] = prev
        clamped[col] = was_clamped
        presence[:, j] = present

    B = pd.DataFrame(presence, index=usage.index, columns=usage.columns)
    info = pd.DataFrame(
        {
            "threshold": pd.Series(thresholds),
            "prevalence": pd.Series(prevalences),
            "clamped": pd.Series(clamped),
        }
    ).loc[usage.columns]
    info.index.name = "program"

    n_clamped = int(info["clamped"].sum())
    if n_clamped:
        warnings.warn(
            f"{n_clamped} of {len(info)} programs had their presence call clamped into "
            f"[{min_prevalence}, {max_prevalence}]; inspect the returned info table.",
            RuntimeWarning,
            stacklevel=2,
        )
    return B, info


# --------------------------------------------------------------------------
# 2. Spatial coordinates
# --------------------------------------------------------------------------

_HD_BARCODE = re.compile(r"^s_\d+um_(\d+)_(\d+)-\d+$")


def coords_from_index(index: Sequence[str]) -> Optional[np.ndarray]:
    """Recover integer grid coordinates from Visium HD spot barcodes.

    Visium HD barcodes encode the bin's grid position, e.g.
    ``s_024um_00131_00210-1`` is row 131, column 210. This lets the spatial null
    models run directly from ``topics_per_spot_*.csv`` without threading an
    AnnData object through the network step. For other platforms, pass
    coordinates explicitly instead.

    Args:
        index (Sequence[str]): Spot barcodes, typically ``usage.index``.

    Returns:
        Optional[numpy.ndarray]: An ``(n_spots, 2)`` integer array of
        ``(row, col)`` positions, or None if the barcodes are not Visium HD
        grid barcodes.
    """
    rows, cols = [], []
    for bc in index:
        m = _HD_BARCODE.match(str(bc))
        if m is None:
            return None
        rows.append(int(m.group(1)))
        cols.append(int(m.group(2)))
    return np.column_stack([np.asarray(rows), np.asarray(cols)])


def _as_grid(coords: np.ndarray) -> Tuple[np.ndarray, int, int]:
    """Map arbitrary spatial coordinates onto a dense integer lattice.

    Real-valued coordinates (e.g. ``adata.obsm["spatial"]``) are ranked onto a
    regular lattice so torus shifts are well defined. Coordinates that are
    already integer grid indices pass through unchanged apart from an offset.

    Args:
        coords (numpy.ndarray): ``(n_spots, 2)`` array of spatial coordinates.

    Returns:
        Tuple[numpy.ndarray, int, int]: ``(grid, n_rows, n_cols)`` where ``grid``
        holds integer lattice positions.
    """
    coords = np.asarray(coords, dtype=float)
    grid = np.empty(coords.shape, dtype=np.int64)
    for axis in (0, 1):
        vals = coords[:, axis]
        uniq = np.unique(vals)
        # Dense rank: ties (spots sharing a row/column) keep the same index.
        grid[:, axis] = np.searchsorted(uniq, vals)
    return grid, int(grid[:, 0].max()) + 1, int(grid[:, 1].max()) + 1


# --------------------------------------------------------------------------
# 3. Co-occurrence statistics
# --------------------------------------------------------------------------

def _contingency(Bv: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Compute the full pairwise 2x2 co-occurrence counts in one pass.

    Args:
        Bv (numpy.ndarray): ``(n_spots, n_programs)`` binary presence matrix.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray, int]: ``(n11, marginal, n_spots)``
        where ``n11[i, j]`` counts spots with both programs present and
        ``marginal[i]`` counts spots with program ``i`` present.
    """
    X = np.asarray(Bv, dtype=np.float32)
    n11 = X.T @ X
    marginal = X.sum(axis=0)
    return n11, marginal, X.shape[0]


def _log2_oe(n11: np.ndarray, marginal: np.ndarray, n_spots: int) -> np.ndarray:
    """Log2 observed-over-expected co-occurrence, with a Laplace correction.

    Expected counts assume the two programs are spatially independent. Because
    the ratio is taken against each program's *own* realised marginal, the
    statistic does not depend on how presence was called -- which is what
    decouples the enrichment target from the binarization rule.

    Args:
        n11 (numpy.ndarray): Pairwise counts of joint presence.
        marginal (numpy.ndarray): Per-program presence counts.
        n_spots (int): Total number of spots.

    Returns:
        numpy.ndarray: ``(n_programs, n_programs)`` matrix of log2(O/E).
    """
    expected = np.outer(marginal, marginal) / max(n_spots, 1)
    return np.log2((n11 + 0.5) / (expected + 0.5))


def cooccurrence_stats(B: pd.DataFrame, sample: str, include_self: bool = False) -> pd.DataFrame:
    """Compute pairwise co-occurrence enrichment for every program pair.

    Replaces the ``joblib``-parallel O(P^2) Python loop in
    ``niche_networks.compute_pairwise_stats`` with a single matrix product, and
    reports enrichment statistics that carry a meaningful null value rather than
    a bare conditional probability.

    The returned frame keeps the legacy column names (``program_one``,
    ``program_two``, ``n``, ``<sample>_P2pos``, ``<sample>_P2neg``) so existing
    downstream code continues to work, and adds:

    * ``log2_oe`` -- log2 observed/expected joint presence. 0 under
      independence, 1 at two-fold enrichment. Symmetric.
    * ``log_or`` -- log odds ratio with a Haldane-Anscombe correction.
      Symmetric.
    * ``jaccard`` -- intersection over union of the two presence sets.
    * ``cond_prob`` -- the legacy ``(P2pos + 1) / (n + 2)`` statistic, retained
      for comparison against previously published networks.

    Args:
        B (pandas.DataFrame): Boolean spots-by-programs presence matrix from
            :func:`call_presence`.
        sample (str): Sample identifier used to name the legacy columns.
        include_self (bool): Whether to keep self-pairs. Defaults to False.

    Returns:
        pandas.DataFrame: One row per ordered program pair.
    """
    programs = list(B.columns)
    n11, marginal, n_spots = _contingency(B.to_numpy())

    n10 = marginal[:, None] - n11          # program_one present, program_two absent
    log2_oe = _log2_oe(n11, marginal, n_spots)

    n01 = marginal[None, :] - n11
    n00 = n_spots - n11 - n10 - n01
    # Haldane-Anscombe correction keeps the odds ratio finite for empty cells.
    log_or = np.log(
        ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
    )
    union = marginal[:, None] + marginal[None, :] - n11
    jaccard = np.divide(n11, union, out=np.zeros_like(n11), where=union > 0)

    i_idx, j_idx = np.meshgrid(
        np.arange(len(programs)), np.arange(len(programs)), indexing="ij"
    )
    keep = np.ones(i_idx.shape, dtype=bool) if include_self else (i_idx != j_idx)

    out = pd.DataFrame(
        {
            "group1": sample,
            "program_one": np.asarray(programs)[i_idx[keep]],
            "program_two": np.asarray(programs)[j_idx[keep]],
            "n": marginal[i_idx[keep]].astype(int),
            f"{sample}_P2pos": n11[keep].astype(int),
            f"{sample}_P2neg": n10[keep].astype(int),
            "n_both": n11[keep].astype(int),
            "n_one": marginal[i_idx[keep]].astype(int),
            "n_two": marginal[j_idx[keep]].astype(int),
            "n_spots": n_spots,
            "log2_oe": log2_oe[keep],
            "log_or": log_or[keep],
            "jaccard": jaccard[keep],
        }
    )
    out["cond_prob"] = (out[f"{sample}_P2pos"] + 1) / (out["n"] + 2)
    out["prevalence_one"] = out["n_one"] / n_spots
    out["prevalence_two"] = out["n_two"] / n_spots
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------
# 4. Null models
# --------------------------------------------------------------------------

def _torus_shift_matrix(Bv: np.ndarray, grid: np.ndarray, n_rows: int, n_cols: int,
                        rng: np.random.Generator) -> np.ndarray:
    """Rigidly translate each program's spatial map on a torus.

    The torus-translation (rigid-shift) test comes from spatial ecology, where
    it is the standard way to test species association without destroying each
    species' own aggregation. Applied here, it preserves each program's spatial
    autocorrelation *exactly* -- the map is moved, not reshuffled -- and
    destroys only the alignment between programs. Each program receives an
    independent shift, so every pair sees a random relative offset.

    Unlike block permutation it introduces no block-size parameter, which
    matters when the whole point is to stop justifying arbitrary constants: a
    block smaller than the autocorrelation range leaves the null anticonservative,
    and the correct block size is not knowable in advance.

    Args:
        Bv (numpy.ndarray): ``(n_spots, n_programs)`` binary presence matrix.
        grid (numpy.ndarray): ``(n_spots, 2)`` integer lattice positions.
        n_rows (int): Number of lattice rows.
        n_cols (int): Number of lattice columns.
        rng (numpy.random.Generator): Random source.

    Returns:
        numpy.ndarray: Presence matrix with each program independently shifted.
    """
    n_spots, n_prog = Bv.shape
    # Scatter onto a dense lattice so a shift is a simple index roll, then
    # gather back at the observed (in-tissue) positions only.
    dense = np.zeros((n_rows, n_cols, n_prog), dtype=Bv.dtype)
    dense[grid[:, 0], grid[:, 1], :] = Bv

    dr = rng.integers(0, n_rows, size=n_prog)
    dc = rng.integers(0, n_cols, size=n_prog)

    out = np.empty((n_spots, n_prog), dtype=Bv.dtype)
    src_r = (grid[:, 0][:, None] - dr[None, :]) % n_rows
    src_c = (grid[:, 1][:, None] - dc[None, :]) % n_cols
    for j in range(n_prog):
        out[:, j] = dense[src_r[:, j], src_c[:, j], j]
    return out


def _block_permute_matrix(Bv: np.ndarray, block_id: np.ndarray,
                          rng: np.random.Generator) -> np.ndarray:
    """Permute contiguous spatial blocks, preserving within-block structure.

    Retained as a secondary null for platforms without a regular lattice. It is
    weaker than the torus shift: correlation on scales larger than one block is
    destroyed, so if the block is smaller than the true autocorrelation range
    the null remains anticonservative.

    Args:
        Bv (numpy.ndarray): ``(n_spots, n_programs)`` binary presence matrix.
        block_id (numpy.ndarray): Block assignment per spot.
        rng (numpy.random.Generator): Random source.

    Returns:
        numpy.ndarray: Presence matrix with blocks independently permuted per
        program.
    """
    out = np.empty_like(Bv)
    blocks = [np.flatnonzero(block_id == b) for b in np.unique(block_id)]
    for j in range(Bv.shape[1]):
        order = rng.permutation(len(blocks))
        # Reassign block contents to shuffled block positions, matching by rank
        # within block so unequal block sizes are handled by recycling.
        for dest, src in enumerate(order):
            d_idx, s_idx = blocks[dest], blocks[src]
            take = np.resize(Bv[s_idx, j], d_idx.size)
            out[d_idx, j] = take
    return out


def _assign_blocks(grid: np.ndarray, n_blocks: int) -> np.ndarray:
    """Assign spots to a roughly square grid of contiguous spatial blocks.

    Args:
        grid (numpy.ndarray): ``(n_spots, 2)`` integer lattice positions.
        n_blocks (int): Target number of blocks.

    Returns:
        numpy.ndarray: Integer block id per spot.
    """
    side = max(int(np.sqrt(max(n_blocks, 1))), 1)
    r = grid[:, 0].astype(float)
    c = grid[:, 1].astype(float)
    # np.ptp(arr) rather than arr.ptp(): the method was removed in NumPy 2.0.
    rb = np.minimum((r - r.min()) / max(np.ptp(r) + 1e-9, 1e-9) * side, side - 1).astype(int)
    cb = np.minimum((c - c.min()) / max(np.ptp(c) + 1e-9, 1e-9) * side, side - 1).astype(int)
    return rb * side + cb


def add_permutation_pvalues(
    stats: pd.DataFrame,
    B: pd.DataFrame,
    coords: Optional[np.ndarray] = None,
    null: str = "torus",
    n_perm: int = 1000,
    n_blocks: int = 64,
    seed: int = 0,
    statistic: str = "log2_oe",
    tail_approx: bool = True,
) -> pd.DataFrame:
    """Attach spatial-permutation p-values and BH-adjusted q-values.

    This is the step that turns the module from a thresholding heuristic into a
    test with a controlled error rate. Spots in a tissue section are strongly
    autocorrelated, so a Fisher or chi-square test that treats them as
    independent draws is severely anticonservative -- on tens of thousands of
    Visium HD bins it will call almost every pair significant. Comparing the
    observed statistic against maps that retain each program's own spatial
    structure removes that inflation.

    Args:
        stats (pandas.DataFrame): Output of :func:`cooccurrence_stats`.
        B (pandas.DataFrame): The presence matrix those statistics came from.
        coords (Optional[numpy.ndarray]): ``(n_spots, 2)`` spatial coordinates
            aligned to ``B.index``. Required for the spatial nulls; when None
            the function falls back to ``"label"`` with a warning.
        null (str): ``"torus"`` (default, recommended), ``"block"``, or
            ``"label"``. ``"label"`` ignores spatial structure and should be
            used only as a sanity floor -- it is the anticonservative null this
            module exists to avoid.
        n_perm (int): Number of permutations. Defaults to 1000, which with
            the tail approximation resolves p-values well past what
            Benjamini-Hochberg needs across P(P-1) pairs.
        n_blocks (int): Target block count for ``null="block"``. Defaults to 64.
        seed (int): Random seed, so results are reproducible. Defaults to 0.
        statistic (str): Which column to test. Defaults to ``"log2_oe"``.
        tail_approx (bool): Extend p-values below the ``1 / (n_perm + 1)``
            sampling floor by fitting a generalized Pareto distribution to the
            null tail. Without this, the strongest edges tie at the floor and
            Benjamini-Hochberg cannot resolve them across thousands of pairs.
            Defaults to True.

    Returns:
        pandas.DataFrame: ``stats`` with added columns ``null_mean``,
        ``null_sd``, ``z`` (standardized effect against the null),
        ``p_empirical`` (add-one corrected, floored at ``1/(n_perm+1)``),
        ``p_perm`` (tail-approximated when ``tail_approx``, else identical to
        ``p_empirical``) and ``q_perm`` (Benjamini-Hochberg).

    Raises:
        ValueError: If ``null`` is not a supported null model.
    """
    if null not in {"torus", "block", "label"}:
        raise ValueError(f"Unknown null {null!r}; expected 'torus', 'block' or 'label'.")

    if coords is None and null != "label":
        warnings.warn(
            f"No spatial coordinates supplied; falling back from null={null!r} to 'label'. "
            "Label permutation ignores spatial autocorrelation and is anticonservative -- "
            "treat the resulting p-values as a floor, not a result.",
            RuntimeWarning,
            stacklevel=2,
        )
        null = "label"

    rng = np.random.default_rng(seed)
    Bv = B.to_numpy().astype(np.float32)
    programs = list(B.columns)
    pos = {p: i for i, p in enumerate(programs)}
    i_idx = stats["program_one"].map(pos).to_numpy()
    j_idx = stats["program_two"].map(pos).to_numpy()

    if null in {"torus", "block"}:
        grid, n_rows, n_cols = _as_grid(np.asarray(coords))
        block_id = _assign_blocks(grid, n_blocks) if null == "block" else None

    null_draws = np.empty((n_perm, len(stats)), dtype=np.float64)
    for k in range(n_perm):
        if null == "torus":
            Bp = _torus_shift_matrix(Bv, grid, n_rows, n_cols, rng)
        elif null == "block":
            Bp = _block_permute_matrix(Bv, block_id, rng)
        else:
            Bp = Bv[rng.permutation(Bv.shape[0]), :]
            # Permuting all columns together would preserve every pairing, so
            # each program must be shuffled independently.
            for j in range(Bp.shape[1]):
                Bp[:, j] = Bv[rng.permutation(Bv.shape[0]), j]

        n11_p, marg_p, n_spots_p = _contingency(Bp)
        if statistic == "log2_oe":
            mat = _log2_oe(n11_p, marg_p, n_spots_p)
        else:
            n10_p = marg_p[:, None] - n11_p
            n01_p = marg_p[None, :] - n11_p
            n00_p = n_spots_p - n11_p - n10_p - n01_p
            mat = np.log(((n11_p + 0.5) * (n00_p + 0.5)) / ((n10_p + 0.5) * (n01_p + 0.5)))
        null_draws[k] = mat[i_idx, j_idx]

    obs = stats[statistic].to_numpy()
    null_mean = null_draws.mean(axis=0)
    null_sd = null_draws.std(axis=0, ddof=1)
    # Add-one correction: an empirical p-value can never be exactly zero, which
    # keeps the BH step well behaved at small n_perm.
    exceed = (null_draws >= obs[None, :]).sum(axis=0)
    p_empirical = (exceed + 1.0) / (n_perm + 1.0)

    if tail_approx:
        p_use = _gpd_tail_pvalues(null_draws, obs)
    else:
        p_use = p_empirical

    out = stats.copy()
    out["null_mean"] = null_mean
    out["null_sd"] = null_sd
    out["z"] = np.divide(
        obs - null_mean, null_sd, out=np.zeros_like(obs), where=null_sd > 0
    )
    out["p_empirical"] = p_empirical
    out["p_perm"] = p_use
    out["q_perm"] = _bh(p_use)
    out.attrs["null"] = null
    out.attrs["n_perm"] = n_perm
    out.attrs["tail_approx"] = tail_approx
    return out


def _gpd_tail_pvalues(
    null_draws: np.ndarray, obs: np.ndarray, n_exceed: int = 100
) -> np.ndarray:
    """Approximate permutation p-values beyond the resolution of the sampler.

    An empirical permutation p-value cannot fall below ``1 / (n_perm + 1)``. With
    P=50 programs there are 2450 ordered pairs, so Benjamini-Hochberg needs
    p-values down to ~2e-5 at the top rank -- unreachable with a few hundred
    permutations, which leaves the strongest edges tied at the floor and
    unresolvable. On the CRC matrix this alone drove the selected edge count to
    zero despite 22 pairs sitting above z = 5.

    Following Knijnenburg et al. (2009), the upper tail of the permutation null
    is modelled with a generalized Pareto distribution fitted to the largest
    ``n_exceed`` draws, which extends the usable p-value range by orders of
    magnitude without running orders of magnitude more permutations.

    Args:
        null_draws (numpy.ndarray): ``(n_perm, n_pairs)`` null statistic draws.
        obs (numpy.ndarray): ``(n_pairs,)`` observed statistics.
        n_exceed (int): Number of upper-tail draws used for the fit. Defaults
            to 100.

    Returns:
        numpy.ndarray: Tail-approximated one-sided p-values, floored three
        orders of magnitude below the sampling floor so the fit is never
        extrapolated further than the data can support. Falls back to the
        empirical value wherever the fit fails or the observation lies inside
        the body of the null.
    """
    from scipy.stats import genpareto

    n_perm, n_pairs = null_draws.shape
    empirical = ((null_draws >= obs[None, :]).sum(axis=0) + 1.0) / (n_perm + 1.0)
    n_exceed = int(min(n_exceed, max(n_perm // 4, 10)))
    # Extrapolation is trustworthy for roughly three orders of magnitude beyond
    # the sampling floor; beyond that the fit is unconstrained by data, so
    # p-values are reported at the floor rather than at a fabricated precision.
    p_floor = (1.0 / (n_perm + 1.0)) * 1e-3
    out = empirical.copy()

    # Only pairs whose empirical p-value is at (or near) the sampling floor need
    # extrapolation; elsewhere the empirical estimate is already well resolved.
    needs_tail = np.flatnonzero((null_draws >= obs[None, :]).sum(axis=0) <= 10)
    for j in needs_tail:
        draws = np.sort(null_draws[:, j])
        tail = draws[-n_exceed:]
        thresh = tail[0]
        if obs[j] <= thresh:
            continue
        exceedances = tail - thresh
        try:
            shape, loc, scale = genpareto.fit(exceedances, floc=0)
            if not np.isfinite(scale) or scale <= 0:
                continue
            sf = genpareto.sf(obs[j] - thresh, shape, loc=loc, scale=scale)
        except Exception:
            continue
        if np.isfinite(sf):
            # Scale the conditional tail probability by the exceedance rate.
            out[j] = min(empirical[j], max((n_exceed / n_perm) * float(sf), p_floor))
    return out


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values.

    Implemented directly rather than via ``statsmodels`` so the network path
    stays importable in minimal environments.

    Args:
        p (numpy.ndarray): Raw p-values.

    Returns:
        numpy.ndarray: BH-adjusted q-values, in the input order.
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity from the largest p-value down.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0, 1)
    return out


# --------------------------------------------------------------------------
# 5. Edge selection
# --------------------------------------------------------------------------

def select_edges(
    stats: pd.DataFrame,
    fdr: float = 0.05,
    min_log2_oe: float = 1.0,
    min_prevalence_frac: float = 0.01,
    keep_depletion: bool = False,
) -> pd.DataFrame:
    """Select niche edges by error rate and effect size instead of a magic number.

    Two criteria replace the single 0.199 cutoff:

    * **Statistical** -- the pair must survive multiple-testing correction
      against the spatial null (``q_perm < fdr``). With P programs there are
      P(P-1) ordered pairs (2450 at P=50), and the previous implementation
      applied no correction at all.
    * **Biological** -- the pair must show at least ``min_log2_oe`` enrichment.
      The default of 1.0 is exactly two-fold, which is what ``edge_threshold =
      0.199`` was implicitly encoding at ``quantile = 0.90``. Stated this way it
      no longer changes meaning when the presence rule changes.

    Args:
        stats (pandas.DataFrame): Output of :func:`add_permutation_pvalues`.
        fdr (float): Benjamini-Hochberg false discovery rate. Defaults to 0.05.
        min_log2_oe (float): Minimum log2 observed/expected enrichment.
            Defaults to 1.0 (two-fold).
        min_prevalence_frac (float): Both programs must be present in at least
            this fraction of spots. Replaces the old absolute ``n_bins=1000``,
            which required >=10,000 spots to ever pass and so silently produced
            empty graphs on standard Visium. Defaults to 0.01.
        keep_depletion (bool): If True, also keep significantly *depleted*
            pairs, tagged ``sign = -1``. They are excluded from community
            detection -- mutual exclusion is evidence against shared niche
            membership, not for it. Defaults to False.

    Returns:
        pandas.DataFrame: The surviving edges, with a ``sign`` column
        (+1 enriched, -1 depleted).

    Raises:
        KeyError: If ``stats`` has no ``q_perm`` column, i.e.
            :func:`add_permutation_pvalues` was not run.
    """
    if "q_perm" not in stats.columns:
        raise KeyError(
            "stats has no 'q_perm' column; run add_permutation_pvalues() before select_edges()."
        )

    prevalent = (
        (stats["prevalence_one"] >= min_prevalence_frac)
        & (stats["prevalence_two"] >= min_prevalence_frac)
    )
    significant = stats["q_perm"] < fdr

    enriched = significant & prevalent & (stats["log2_oe"] >= min_log2_oe)
    if keep_depletion:
        depleted = significant & prevalent & (stats["log2_oe"] <= -min_log2_oe)
        keep = enriched | depleted
    else:
        keep = enriched

    out = stats.loc[keep].copy()
    out["sign"] = np.where(out["log2_oe"] >= 0, 1, -1)
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------
# 6. Stability selection
# --------------------------------------------------------------------------

def stability_selection(
    usage: pd.DataFrame,
    coords: Optional[np.ndarray] = None,
    n_subsamples: int = 50,
    subsample_frac: float = 0.5,
    n_blocks: int = 64,
    presence_method: str = "otsu",
    presence_quantile: float = 0.90,
    fdr: float = 0.05,
    min_log2_oe: float = 1.0,
    n_perm: int = 100,
    null: str = "torus",
    seed: int = 0,
    sample: str = "sample",
) -> Tuple[pd.DataFrame, float]:
    """Rank edges by how often they survive on spatial subsamples.

    Implements Meinshausen-Buhlmann stability selection. Subsampling is done in
    contiguous spatial *blocks* rather than individual spots: resampling spots
    independently would break the autocorrelation the null models are built to
    respect, and would make every subsample look like the full section.

    The returned bound is the Meinshausen-Buhlmann inequality
    ``E[V] <= q^2 / ((2*pi_thr - 1) * p)`` on the expected number of falsely
    selected edges, valid for ``pi_thr > 0.5``. This gives the methods section
    an actual error guarantee to quote in place of a chosen threshold.

    Args:
        usage (pandas.DataFrame): Spots-by-programs usage matrix.
        coords (Optional[numpy.ndarray]): Spatial coordinates aligned to
            ``usage.index``.
        n_subsamples (int): Number of spatial subsamples. Defaults to 50.
        subsample_frac (float): Fraction of blocks retained per subsample.
            Defaults to 0.5, as assumed by the Meinshausen-Buhlmann bound.
        n_blocks (int): Number of contiguous blocks to subsample over.
            Defaults to 64.
        presence_method (str): Passed to :func:`call_presence`.
        presence_quantile (float): Passed to :func:`call_presence`.
        fdr (float): Passed to :func:`select_edges`.
        min_log2_oe (float): Passed to :func:`select_edges`.
        n_perm (int): Permutations per subsample. Defaults to 100.
        null (str): Null model per subsample. Defaults to ``"torus"``.
        seed (int): Random seed. Defaults to 0.
        sample (str): Sample label. Defaults to ``"sample"``.

    Returns:
        Tuple[pandas.DataFrame, float]: ``(freq, bound_per_pi)`` where ``freq``
        holds each program pair's selection frequency and mean enrichment, and
        ``bound_per_pi`` is ``q^2 / p`` -- divide by ``(2*pi_thr - 1)`` to get
        the expected number of false selections at a chosen threshold.
    """
    rng = np.random.default_rng(seed)
    if coords is None:
        coords = coords_from_index(usage.index)

    if coords is None:
        block_id = rng.integers(0, n_blocks, size=len(usage))
        warnings.warn(
            "No coordinates available; subsampling spots at random rather than in "
            "spatial blocks. Selection frequencies will be optimistic.",
            RuntimeWarning,
            stacklevel=2,
        )
    else:
        grid, _, _ = _as_grid(np.asarray(coords))
        block_id = _assign_blocks(grid, n_blocks)

    blocks = np.unique(block_id)
    counts: dict = {}
    effects: dict = {}
    n_selected_per_run = []

    for it in range(n_subsamples):
        chosen = rng.choice(blocks, size=max(int(len(blocks) * subsample_frac), 1), replace=False)
        mask = np.isin(block_id, chosen)
        sub = usage.loc[mask]
        sub_coords = None if coords is None else np.asarray(coords)[mask]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            Bs, _ = call_presence(sub, method=presence_method, quantile=presence_quantile, seed=seed + it)
            st = cooccurrence_stats(Bs, sample=sample)
            st = add_permutation_pvalues(
                st, Bs, coords=sub_coords, null=null, n_perm=n_perm, seed=seed + it
            )
            sel = select_edges(st, fdr=fdr, min_log2_oe=min_log2_oe)

        n_selected_per_run.append(len(sel))
        for row in sel.itertuples(index=False):
            key = (row.program_one, row.program_two)
            counts[key] = counts.get(key, 0) + 1
            effects.setdefault(key, []).append(row.log2_oe)

    if not counts:
        return pd.DataFrame(
            columns=["program_one", "program_two", "selection_freq", "mean_log2_oe"]
        ), 0.0

    freq = pd.DataFrame(
        [
            {
                "program_one": a,
                "program_two": b,
                "selection_freq": c / n_subsamples,
                "mean_log2_oe": float(np.mean(effects[(a, b)])),
            }
            for (a, b), c in counts.items()
        ]
    ).sort_values("selection_freq", ascending=False, ignore_index=True)

    n_programs = usage.shape[1]
    n_candidate_pairs = n_programs * (n_programs - 1)
    q_avg = float(np.mean(n_selected_per_run))
    bound_per_pi = (q_avg ** 2) / max(n_candidate_pairs, 1)
    freq.attrs["q_avg"] = q_avg
    freq.attrs["n_candidate_pairs"] = n_candidate_pairs
    return freq, bound_per_pi


# --------------------------------------------------------------------------
# 7. Robustness sweep
# --------------------------------------------------------------------------

def _partition_from_edges(edges: pd.DataFrame, programs: Sequence[str], seed: int = 0) -> dict:
    """Cluster programs into niches from a selected edge list.

    Only positively enriched edges contribute, and Infomap is run against a
    seeded random source so the partition is reproducible -- the previous
    implementation left it unseeded, so niche labels could change between runs
    of identical input.

    Args:
        edges (pandas.DataFrame): Selected edges with ``program_one`` /
            ``program_two`` / ``log2_oe``.
        programs (Sequence[str]): Full program list, so isolated programs still
            receive a label.
        seed (int): Random seed for Infomap. Defaults to 0.

    Returns:
        dict: Mapping of program name to niche id (-1 when unassigned).
    """
    import random

    import igraph as ig

    labels = {p: -1 for p in programs}
    pos_edges = edges[edges.get("sign", 1) > 0] if "sign" in edges.columns else edges
    if pos_edges.empty:
        return labels

    tuples = list(
        pos_edges[["program_one", "program_two", "log2_oe"]].itertuples(index=False, name=None)
    )
    g = ig.Graph.TupleList(tuples, directed=True, vertex_name_attr="name", edge_attrs=["weight"])

    state = random.Random(seed)
    try:
        ig.set_random_number_generator(state)
        comms = g.community_infomap(trials=10, edge_weights="weight")
    finally:
        ig.set_random_number_generator(random)

    for name, member in zip(g.vs["name"], comms.membership):
        labels[name] = int(member)
    return labels


def parameter_sweep(
    usage: pd.DataFrame,
    coords: Optional[np.ndarray] = None,
    presence_methods: Sequence[str] = ("otsu", "mixture", "quantile"),
    quantiles: Sequence[float] = (0.80, 0.85, 0.90, 0.95),
    fdrs: Sequence[float] = (0.01, 0.05, 0.10),
    min_log2_oes: Sequence[float] = (0.5, 1.0, 1.5),
    n_perm: int = 100,
    null: str = "torus",
    seed: int = 0,
    sample: str = "sample",
) -> pd.DataFrame:
    """Evaluate how the inferred niches respond to every free parameter.

    This is the evidence the reviewer asked for: rather than defending a single
    setting, it reports what the network does across the whole parameter grid.
    For each combination it records the edge count, the niche count, and the
    agreement of both the edge set (Jaccard) and the partition (adjusted Rand
    index) with the reference setting.

    **Read ``ari`` and ``jaccard`` jointly with ``n_edges``.** Both agreement
    measures saturate trivially in the sparse corner -- a setting that keeps two
    edges will reproduce those two edges perfectly. The result to report is the
    widest contiguous region where agreement is high *and* the network is
    non-trivial, not the maximum of either column. The ``degenerate`` flag marks
    rows with too few edges to interpret.

    Args:
        usage (pandas.DataFrame): Spots-by-programs usage matrix.
        coords (Optional[numpy.ndarray]): Spatial coordinates aligned to
            ``usage.index``. Parsed from the index when None.
        presence_methods (Sequence[str]): Presence rules to sweep.
        quantiles (Sequence[float]): Quantiles to sweep for the ``"quantile"``
            rule. Ignored by the other rules, which take no quantile.
        fdrs (Sequence[float]): FDR levels to sweep.
        min_log2_oes (Sequence[float]): Enrichment cutoffs to sweep.
        n_perm (int): Permutations per grid cell. Defaults to 100.
        null (str): Null model. Defaults to ``"torus"``.
        seed (int): Random seed. Defaults to 0.
        sample (str): Sample label. Defaults to ``"sample"``.

    Returns:
        pandas.DataFrame: One row per parameter combination, with columns
        ``presence_method``, ``quantile``, ``fdr``, ``min_log2_oe``,
        ``n_edges``, ``n_niches``, ``median_prevalence``, ``jaccard`` and
        ``ari`` (both against the reference setting), ``degenerate``, and a
        ``key`` identifying each row.
    """
    from sklearn.metrics import adjusted_rand_score

    if coords is None:
        coords = coords_from_index(usage.index)
    programs = list(usage.columns)

    # Presence calls only depend on the rule, so compute each once.
    presence_cache = {}
    for method in presence_methods:
        qs = quantiles if method == "quantile" else (None,)
        for q in qs:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                presence_cache[(method, q)] = call_presence(
                    usage, method=method, quantile=(q if q is not None else 0.90), seed=seed
                )

    # Statistics likewise depend only on the presence call, not on fdr/effect size.
    stats_cache = {}
    for key, (B, _info) in presence_cache.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            st = cooccurrence_stats(B, sample=sample)
            stats_cache[key] = add_permutation_pvalues(
                st, B, coords=coords, null=null, n_perm=n_perm, seed=seed
            )

    rows = []
    edge_sets, partitions = {}, {}
    # Rows are keyed by an integer id rather than the parameter tuple: a None
    # quantile (the methods that take none) becomes NaN once it round-trips
    # through the frame, and NaN != NaN would break the lookup.
    for (method, q), st in stats_cache.items():
        _B, info = presence_cache[(method, q)]
        for fdr in fdrs:
            for mo in min_log2_oes:
                sel = select_edges(st, fdr=fdr, min_log2_oe=mo)
                key = len(rows)
                edge_sets[key] = set(zip(sel["program_one"], sel["program_two"]))
                partitions[key] = _partition_from_edges(sel, programs, seed=seed)
                rows.append(
                    {
                        "key": key,
                        "presence_method": method,
                        "quantile": q,
                        "fdr": fdr,
                        "min_log2_oe": mo,
                        "n_edges": len(sel),
                        "n_niches": len({v for v in partitions[key].values() if v >= 0}),
                        "median_prevalence": float(info["prevalence"].median()),
                        "_is_reference": (method == "otsu" and fdr == 0.05 and mo == 1.0),
                    }
                )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    # Reference: the recommended default, or the densest cell if it is absent.
    ref_rows = out.index[out["_is_reference"]]
    ref_key = int(out.loc[ref_rows[0], "key"]) if len(ref_rows) else max(
        edge_sets, key=lambda k: len(edge_sets[k])
    )

    ref_edges = edge_sets[ref_key]
    ref_part = partitions[ref_key]
    ref_labels = [ref_part[p] for p in programs]

    jac, ari = [], []
    for key in out["key"]:
        es = edge_sets[key]
        union = len(ref_edges | es)
        jac.append(len(ref_edges & es) / union if union else np.nan)
        ari.append(adjusted_rand_score(ref_labels, [partitions[key][p] for p in programs]))

    out["jaccard"] = jac
    out["ari"] = ari
    # Agreement metrics saturate trivially on near-empty graphs; flag rather
    # than drop, so the sparse corner stays visible instead of looking optimal.
    out["degenerate"] = out["n_edges"] < max(len(programs) // 2, 3)
    out.attrs["reference"] = ref_key
    out = out.drop(columns=["_is_reference"])
    return out.sort_values(
        ["presence_method", "quantile", "fdr", "min_log2_oe"], ignore_index=True
    )


# --------------------------------------------------------------------------
# 8. Compositional cross-check
# --------------------------------------------------------------------------

def clr_transform(usage: pd.DataFrame, pseudocount: Optional[float] = None) -> pd.DataFrame:
    """Centred log-ratio transform of a compositional usage matrix.

    ``topics_per_spot`` rows sum to 1, so the usages are compositional: a
    program can appear absent purely because others are abundant. The induced
    negative dependence is of order ``-1/(K-1)`` (about -0.02 at K=50), so it is
    a second-order effect here rather than a reason to abandon co-occurrence --
    and it biases *against* detecting positive association, which means edges
    that survive are not artefacts of the constraint. This transform supports
    reporting that check rather than asserting it.

    Args:
        usage (pandas.DataFrame): Row-compositional usage matrix.
        pseudocount (Optional[float]): Added before the log. Defaults to half
            the smallest non-zero usage.

    Returns:
        pandas.DataFrame: CLR-transformed usages.
    """
    X = usage.to_numpy(dtype=float)
    if pseudocount is None:
        nz = X[X > 0]
        pseudocount = float(nz.min() / 2) if nz.size else 1e-9
    L = np.log(X + pseudocount)
    return pd.DataFrame(L - L.mean(axis=1, keepdims=True), index=usage.index, columns=usage.columns)


def compositional_check(usage: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    """Compare selected edges against a CLR-based association measure.

    Reports, for every selected edge, the Spearman correlation of the two
    programs' CLR-transformed usages. Edges driven by genuine co-localization
    should show positive CLR association; edges that appear only in the binary
    analysis warrant a closer look.

    Args:
        usage (pandas.DataFrame): The original compositional usage matrix.
        edges (pandas.DataFrame): Selected edges from :func:`select_edges`.

    Returns:
        pandas.DataFrame: ``edges`` with an added ``clr_spearman`` column.
    """
    from scipy.stats import rankdata

    if edges.empty:
        out = edges.copy()
        out["clr_spearman"] = pd.Series(dtype=float)
        return out

    clr = clr_transform(usage)
    R = np.apply_along_axis(rankdata, 0, clr.to_numpy())
    R = (R - R.mean(axis=0)) / R.std(axis=0)
    corr = (R.T @ R) / R.shape[0]
    pos = {p: i for i, p in enumerate(usage.columns)}

    out = edges.copy()
    out["clr_spearman"] = [
        corr[pos[a], pos[b]] for a, b in zip(out["program_one"], out["program_two"])
    ]
    return out
