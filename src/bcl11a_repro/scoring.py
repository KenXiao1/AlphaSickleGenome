"""Independent quantitative scoring and numerical cross-validation against AlphaGenome 2."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from alphagenome.data import genome
from alphagenome.interpretation import ism as ism_lib
from alphagenome.models.v2 import variant_scorers
import anndata
import numpy as np
import pandas as pd

from .constants import SC_ATAC_CENTER_MASK_WIDTH, SC_ATAC_PSEUDOCOUNT


def compute_center_mask_log2_sum_score(
    ref_values: np.ndarray,
    alt_values: np.ndarray,
    interval_start: int,
    variant_zero_based_start: int,
    width: int = SC_ATAC_CENTER_MASK_WIDTH,
    pseudocount: float = SC_ATAC_PSEUDOCOUNT,
) -> dict[str, Any]:
    """Independently computes the AlphaGenome 2 CenterMaskScorer (DIFF_LOG2_SUM) from track arrays.

    Paper Methods (p. 35):
        Delta_ATAC = log2(sum_{i in W_1201} y_ALT,i + 1) - log2(sum_{i in W_1201} y_REF,i + 1)

    Window alignment:
        For an odd width W = 2*r + 1 centered on a 1-bp SNV at 0-based `variant_zero_based_start`,
        `variant.reference_interval.resize(W)` spans `[variant_start - r, variant_start + r + 1)`.
    """
    if ref_values.shape != alt_values.shape:
        raise ValueError(
            f"REF shape {ref_values.shape} != ALT shape {alt_values.shape}."
        )
    if ref_values.ndim != 2:
        raise ValueError(
            f"Expected 2D track array (positions, tracks), got shape {ref_values.shape}."
        )
    if width <= 0 or width % 2 != 1:
        raise ValueError(f"Center-mask width must be positive odd integer, got {width}.")

    half_width = width // 2
    rel_center = variant_zero_based_start - interval_start
    win_start = rel_center - half_width
    win_end = rel_center + half_width + 1

    if win_start < 0 or win_end > ref_values.shape[0]:
        raise ValueError(
            f"Scoring window [{win_start}, {win_end}) out of bounds for length {ref_values.shape[0]}."
        )

    ref_win = ref_values[win_start:win_end, :].astype(np.float64)
    alt_win = alt_values[win_start:win_end, :].astype(np.float64)

    ref_sum = ref_win.sum(axis=0)
    alt_sum = alt_win.sum(axis=0)

    if np.any(ref_sum + pseudocount <= 0) or np.any(alt_sum + pseudocount <= 0):
        raise ValueError("Summed accessibility plus pseudocount must be strictly positive.")

    delta_log2 = np.log2(alt_sum + pseudocount) - np.log2(ref_sum + pseudocount)
    fold_change_p1 = (alt_sum + pseudocount) / (ref_sum + pseudocount)
    pct_change_sum = np.divide(
        100.0 * (alt_sum - ref_sum),
        ref_sum,
        out=np.zeros_like(ref_sum, dtype=np.float64),
        where=ref_sum > 0,
    )

    return {
        "window_0based": (interval_start + win_start, interval_start + win_end),
        "relative_slice": (win_start, win_end),
        "ref_sum": ref_sum,
        "alt_sum": alt_sum,
        "delta_log2": delta_log2,
        "fold_change_p1": fold_change_p1,
        "pct_change_sum": pct_change_sum,
    }


def cross_validate_scores(
    api_tidy_df: pd.DataFrame,
    ref_values: np.ndarray,
    alt_values: np.ndarray,
    interval_start: int,
    variant_zero_based_start: int,
    width: int = SC_ATAC_CENTER_MASK_WIDTH,
    pseudocount: float = SC_ATAC_PSEUDOCOUNT,
    atol: float = 1e-4,
) -> pd.DataFrame:
    """Cross-validates official `score_variant()` output against independent array recomputation."""
    if len(api_tidy_df) != ref_values.shape[1]:
        raise ValueError(
            f"Track count mismatch: api_tidy_df has {len(api_tidy_df)} rows, "
            f"track array has {ref_values.shape[1]} columns."
        )

    manual = compute_center_mask_log2_sum_score(
        ref_values=ref_values,
        alt_values=alt_values,
        interval_start=interval_start,
        variant_zero_based_start=variant_zero_based_start,
        width=width,
        pseudocount=pseudocount,
    )

    out = api_tidy_df.copy().reset_index(drop=True)
    api_scores = out["raw_score"].to_numpy(dtype=np.float64)
    manual_scores = manual["delta_log2"]
    abs_diff = np.abs(api_scores - manual_scores)

    out["ref_sum_1201bp"] = manual["ref_sum"]
    out["alt_sum_1201bp"] = manual["alt_sum"]
    out["api_raw_score"] = api_scores
    out["manual_log2_score"] = manual_scores
    out["abs_diff"] = abs_diff
    out["fold_change_p1"] = manual["fold_change_p1"]
    out["pct_change_sum"] = manual["pct_change_sum"]

    max_diff = float(np.max(abs_diff))
    if max_diff > atol:
        raise ValueError(
            f"Numerical cross-validation FAILED: max |api - manual| = {max_diff:.6g} > atol={atol}. "
            "Check REF/ALT ordering, window centering, log base, or pseudocount."
        )
    return out


def ism_anndata_list_to_tidy_df(
    ism_scores: Sequence[Sequence[anndata.AnnData]],
    scorers: Sequence[variant_scorers.VariantScorerTypes],
) -> pd.DataFrame:
    """Attaches `variant_scorer` to ISM AnnData outputs (omitted by SDK v2.0.0) and tidies."""
    for per_variant_list in ism_scores:
        if len(per_variant_list) != len(scorers):
            raise ValueError(
                f"ISM scorer count mismatch: {len(per_variant_list)} vs {len(scorers)}."
            )
        for ix, adata in enumerate(per_variant_list):
            if "variant_scorer" not in adata.uns:
                adata.uns["variant_scorer"] = scorers[ix]

    tidy_df = variant_scorers.tidy_scores(ism_scores)
    if tidy_df is None or tidy_df.empty:
        raise ValueError("ISM tidy_scores produced an empty DataFrame.")
    return tidy_df


def build_ism_contribution_matrix(
    ism_tidy_df: pd.DataFrame,
    ism_interval: genome.Interval,
    sequence: str,
    track_name: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Builds the (width, 4) mean-centered ISM sequence attribution matrix for a specific track."""
    sub = ism_tidy_df[ism_tidy_df["track_name"] == track_name].copy()
    if sub.empty:
        raise ValueError(f"No ISM rows found for track_name='{track_name}'.")

    variants = [genome.Variant.from_str(str(v)) for v in sub["variant_id"]]
    scores = sub["raw_score"].to_numpy(dtype=np.float64).tolist()

    mat = ism_lib.ism_matrix(
        variant_scores=scores,
        variants=variants,
        interval=ism_interval,
        multiply_by_sequence=True,
        vocabulary="ACGT",
        require_fully_filled=True,
    )

    rows = []
    for idx in range(ism_interval.width):
        pos_0 = ism_interval.start + idx
        pos_1 = pos_0 + 1
        base = sequence[idx]
        rows.append(
            {
                "position_0based": pos_0,
                "position_1based": pos_1,
                "sequence_base": base,
                "attr_A": float(mat[idx, 0]),
                "attr_C": float(mat[idx, 1]),
                "attr_G": float(mat[idx, 2]),
                "attr_T": float(mat[idx, 3]),
                "active_base_attribution": float(mat[idx].sum()),
            }
        )
    return mat, pd.DataFrame(rows)
