"""Tests for the statistically principled niche inference module.

These cover the behaviours the fixed-threshold implementation got wrong:
presence calls that were identical for every program, an edge statistic with no
null, a null that ignored spatial autocorrelation, and thresholds that changed
meaning when another parameter moved.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spotnmf import niche_stats as ns


def _toy_usage(n_side: int = 24, seed: int = 0):
    """Build a small spatial usage matrix with two planted co-occurring pairs.

    Returns:
        Tuple[pandas.DataFrame, numpy.ndarray]: usage matrix and coordinates.
    """
    rng = np.random.default_rng(seed)
    rr, cc = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    coords = np.column_stack([rr.ravel(), cc.ravel()])

    # Two disjoint spatial quadrants carry the planted signal.
    quad_a = (rr < n_side // 2) & (cc < n_side // 2)
    quad_b = (rr >= n_side // 2) & (cc >= n_side // 2)

    cols = {}
    cols["A1"] = quad_a.ravel() * 1.0 + rng.random(n_side ** 2) * 0.05
    cols["A2"] = quad_a.ravel() * 1.0 + rng.random(n_side ** 2) * 0.05
    cols["B1"] = quad_b.ravel() * 1.0 + rng.random(n_side ** 2) * 0.05
    cols["B2"] = quad_b.ravel() * 1.0 + rng.random(n_side ** 2) * 0.05
    for k in range(4):
        cols[f"N{k}"] = rng.random(n_side ** 2)

    usage = pd.DataFrame(cols)
    usage.index = [f"s_024um_{r:05d}_{c:05d}-1" for r, c in coords]
    return usage, coords


def test_quantile_presence_forces_identical_prevalence():
    """The legacy rule makes every program present in the same share of spots."""
    usage, _ = _toy_usage()
    _B, info = ns.call_presence(usage, method="quantile", quantile=0.90)
    # Every program lands on exactly the same prevalence -- the flattening that
    # motivated replacing this rule.
    assert info["prevalence"].nunique() == 1
    assert info["prevalence"].iloc[0] == pytest.approx(0.10, abs=0.01)


def test_otsu_presence_recovers_differing_extent():
    """Data-driven presence lets programs differ in spatial extent."""
    usage, _ = _toy_usage()
    _B, info = ns.call_presence(usage, method="otsu")
    # The planted quadrant programs cover ~25% of the section; the noise
    # programs should not be forced onto that same value.
    assert info["prevalence"].nunique() > 1
    for prog in ("A1", "A2", "B1", "B2"):
        assert info.loc[prog, "prevalence"] == pytest.approx(0.25, abs=0.05)


def test_log2_oe_is_invariant_to_presence_rule():
    """The enrichment statistic does not shift when the binarization changes.

    This is the property the legacy 0.199 cutoff lacked: its meaning was tied to
    the 90th-percentile step, so changing the quantile silently changed the
    effective enrichment demanded.
    """
    usage, _ = _toy_usage()
    scores = {}
    for q in (0.80, 0.90):
        B, _info = ns.call_presence(usage, method="quantile", quantile=q)
        st = ns.cooccurrence_stats(B, sample="S")
        pair = st[(st.program_one == "A1") & (st.program_two == "A2")]
        scores[q] = float(pair["log2_oe"].iloc[0])
    # Both binarizations see the same planted pair as strongly enriched.
    assert scores[0.80] > 1.0 and scores[0.90] > 1.0
    # The legacy conditional probability, by contrast, tracks the prevalence.
    conds = {}
    for q in (0.80, 0.90):
        B, _info = ns.call_presence(usage, method="quantile", quantile=q)
        st = ns.cooccurrence_stats(B, sample="S")
        pair = st[(st.program_one == "A1") & (st.program_two == "A2")]
        conds[q] = float(pair["cond_prob"].iloc[0])
    assert conds[0.80] != pytest.approx(conds[0.90], abs=1e-6)


def test_cooccurrence_stats_matches_direct_counts():
    """Vectorized contingency counts agree with a direct computation."""
    usage, _ = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S").set_index(["program_one", "program_two"])
    both = int((B["A1"] & B["A2"]).sum())
    assert st.loc[("A1", "A2"), "n_both"] == both
    assert st.loc[("A1", "A2"), "n"] == int(B["A1"].sum())
    # log2(O/E) is symmetric even though the frame stores ordered pairs.
    assert st.loc[("A1", "A2"), "log2_oe"] == pytest.approx(st.loc[("A2", "A1"), "log2_oe"])


def test_spatial_null_is_wider_than_label_permutation():
    """Ignoring autocorrelation understates the null spread."""
    usage, coords = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    lab = ns.add_permutation_pvalues(st, B, coords=coords, null="label", n_perm=60, seed=0)
    tor = ns.add_permutation_pvalues(st, B, coords=coords, null="torus", n_perm=60, seed=0)
    assert tor["null_sd"].median() > lab["null_sd"].median()


def test_select_edges_recovers_planted_pairs():
    """The planted co-occurring pairs survive selection; noise pairs do not."""
    usage, coords = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    st = ns.add_permutation_pvalues(st, B, coords=coords, null="torus", n_perm=100, seed=0)
    sel = ns.select_edges(st, fdr=0.10, min_log2_oe=1.0)
    found = {tuple(sorted(p)) for p in zip(sel.program_one, sel.program_two)}
    assert ("A1", "A2") in found
    assert ("B1", "B2") in found
    # A planted pair should never link across the two disjoint quadrants.
    assert ("A1", "B1") not in found


def test_select_edges_requires_pvalues():
    """Selection refuses to run without a null model having been fitted."""
    usage, _ = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    with pytest.raises(KeyError, match="add_permutation_pvalues"):
        ns.select_edges(st, fdr=0.05)


def test_coords_from_index_parses_and_rejects():
    """Visium HD barcodes yield coordinates; other indices return None."""
    usage, coords = _toy_usage()
    parsed = ns.coords_from_index(usage.index)
    assert parsed is not None
    assert np.array_equal(parsed, coords)
    assert ns.coords_from_index(["spot1", "spot2"]) is None


def test_missing_coords_falls_back_to_label_with_warning():
    """A missing coordinate set must degrade loudly, not silently."""
    usage, _ = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    with pytest.warns(RuntimeWarning, match="anticonservative"):
        out = ns.add_permutation_pvalues(st, B, coords=None, null="torus", n_perm=20, seed=0)
    assert out.attrs["null"] == "label"


def test_bh_matches_statsmodels():
    """The in-tree BH implementation agrees with the reference one."""
    statsmodels = pytest.importorskip("statsmodels.stats.multitest")
    rng = np.random.default_rng(0)
    p = rng.random(200) ** 3
    _rej, expected, _a, _b = statsmodels.multipletests(p, method="fdr_bh")
    assert np.allclose(ns._bh(p), expected)


def test_permutation_pvalues_are_reproducible():
    """A fixed seed gives identical p-values."""
    usage, coords = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    a = ns.add_permutation_pvalues(st, B, coords=coords, null="torus", n_perm=40, seed=3)
    b = ns.add_permutation_pvalues(st, B, coords=coords, null="torus", n_perm=40, seed=3)
    assert np.allclose(a["p_perm"], b["p_perm"])


def test_negative_control_yields_no_edges():
    """Destroying real alignment while keeping autocorrelation should find nothing."""
    usage, coords = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    grid, n_rows, n_cols = ns._as_grid(coords)
    rng = np.random.default_rng(11)
    shifted = ns._torus_shift_matrix(B.to_numpy().astype(np.float32), grid, n_rows, n_cols, rng)
    Bs = pd.DataFrame(shifted.astype(bool), index=B.index, columns=B.columns)

    st = ns.cooccurrence_stats(Bs, sample="S")
    st = ns.add_permutation_pvalues(st, Bs, coords=coords, null="torus", n_perm=100, seed=0)
    sel = ns.select_edges(st, fdr=0.05, min_log2_oe=1.0)
    assert len(sel) == 0


def test_clr_transform_is_centred():
    """CLR rows sum to zero, as the transform requires."""
    usage, _ = _toy_usage()
    comp = usage.div(usage.sum(axis=1), axis=0)
    clr = ns.clr_transform(comp)
    assert np.allclose(clr.to_numpy().sum(axis=1), 0.0, atol=1e-8)


def test_call_presence_rejects_unknown_method():
    """An unsupported presence rule fails loudly."""
    usage, _ = _toy_usage()
    with pytest.raises(ValueError, match="Unknown presence method"):
        ns.call_presence(usage, method="nope")


def test_add_permutation_pvalues_rejects_unknown_null():
    """An unsupported null model fails loudly."""
    usage, coords = _toy_usage()
    B, _ = ns.call_presence(usage, method="otsu")
    st = ns.cooccurrence_stats(B, sample="S")
    with pytest.raises(ValueError, match="Unknown null"):
        ns.add_permutation_pvalues(st, B, coords=coords, null="nope")
