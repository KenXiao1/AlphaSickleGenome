"""End-to-end cached execution & offline regeneration pipeline for BCL11A AG2 reproduction."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone
import hashlib
import json
import os
import pathlib
import platform
import subprocess
from typing import Any

from alphagenome import __version__ as ag_version
from alphagenome.models.v2 import dna_client
from alphagenome.models.v2 import variant_scorers
import numpy as np
import pandas as pd

from .constants import (
    ALPHAGENOME_GENCODE_V46_FEATHER_URL,
    BCL11A_GENE_ID,
    CONDITION_BLOOD_NAIVE_B,
    CONDITION_BM_ERYTHROBLAST,
    CONDITION_BM_NAIVE_B_TUTORIAL,
    ENHANCER_55KB_DHS_TUTORIAL,
    ENHANCER_55KB_PAPER,
    ENHANCER_58KB_CORE,
    ENHANCER_58KB_DHS_TUTORIAL,
    ENHANCER_62KB_DHS_TUTORIAL,
    ENHANCER_62KB_PAPER,
    OBSERVED_SC_ATAC_TRACK_SUFFIX,
    REFERENCE_ASSEMBLY,
    SAVED_LOCUS_SLICE_0BASED,
    SC_ATAC_CENTER_MASK_WIDTH,
    SC_ATAC_PSEUDOCOUNT,
    SG1617_GRNA_TARGET_0BASED,
    UCSC_HG38_SEQUENCE_URL,
    VARIANT_STR,
)
from .plotting import generate_all_figures
from .scoring import (
    build_ism_contribution_matrix,
    cross_validate_scores,
    ism_anndata_list_to_tidy_df,
)
from .validation import (
    assert_no_secret_leak,
    build_bcl11a_variant_and_windows,
    load_api_key_safely,
    verify_enhancer_and_gene_annotation,
    verify_reference_sequence,
    verify_supported_conditions,
)


def _download_with_curl(url: str, dest: pathlib.Path, timeout_sec: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    curl_bin = "curl.exe" if os.name == "nt" else "curl"
    subprocess.check_call(
        [curl_bin, "-fsSL", "--retry", "2", "--max-time", str(timeout_sec), url, "-o", str(dest)]
    )


def _sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ensure_ucsc_and_gtf_cache(
    cache_dir: pathlib.Path,
    *,
    offline: bool,
    refresh_cache: bool,
) -> tuple[dict[str, Any], pd.DataFrame]:
    ucsc_path = cache_dir / "ucsc_hg38_chr2_60495000_60495600.json"
    bcl11a_gtf_parquet = cache_dir / "bcl11a_gencode_v46.parquet"
    full_gtf_feather = cache_dir / "gencode.v46.annotation.gtf.gz.feather"

    if refresh_cache or not ucsc_path.exists():
        if offline:
            raise FileNotFoundError(
                f"Offline mode enabled, but UCSC cache file is missing: {ucsc_path}"
            )
        _download_with_curl(UCSC_HG38_SEQUENCE_URL, ucsc_path, timeout_sec=60)

    if refresh_cache or not bcl11a_gtf_parquet.exists():
        if offline:
            raise FileNotFoundError(
                f"Offline mode enabled, but GENCODE v46 BCL11A cache is missing: {bcl11a_gtf_parquet}"
            )
        if refresh_cache or not full_gtf_feather.exists():
            _download_with_curl(ALPHAGENOME_GENCODE_V46_FEATHER_URL, full_gtf_feather, timeout_sec=180)
        gtf_all = pd.read_feather(full_gtf_feather)
        bcl11a_sub = gtf_all[
            gtf_all["gene_id"].fillna("").str.startswith(BCL11A_GENE_ID)
        ].copy()
        bcl11a_sub.to_parquet(bcl11a_gtf_parquet, index=False)

    ucsc_payload = json.loads(ucsc_path.read_text(encoding="utf-8"))
    bcl11a_gtf_df = pd.read_parquet(bcl11a_gtf_parquet)
    return ucsc_payload, bcl11a_gtf_df


def regenerate_figures_from_processed(
    processed_dir: pathlib.Path,
    figures_dir: pathlib.Path,
) -> dict[str, dict[str, pathlib.Path]]:
    """Regenerates Figures A–D strictly from saved files in `processed_dir` (no API or raw cache needed)."""
    meta = json.loads((processed_dir / "run_metadata.json").read_text(encoding="utf-8"))
    tracks_df = pd.read_parquet(processed_dir / "sc_atac_window_10kb_tracks.parquet")
    scores_df = pd.read_csv(processed_dir / "variant_effect_scores_comparison.csv")
    bcl11a_gtf_df = pd.read_parquet(processed_dir / "bcl11a_gencode_v46.parquet")
    ism_df = pd.read_csv(processed_dir / "ism_attribution_chr2_60495262_60495285.csv")

    ref_ism = ism_df[ism_df["allele_background"] == "REF"][
        ["attr_A", "attr_C", "attr_G", "attr_T"]
    ].to_numpy(dtype=np.float32)
    alt_ism = ism_df[ism_df["allele_background"] == "ALT"][
        ["attr_A", "attr_C", "attr_G", "attr_T"]
    ].to_numpy(dtype=np.float32)

    return generate_all_figures(
        tracks_df=tracks_df,
        scores_df=scores_df,
        bcl11a_gtf_df=bcl11a_gtf_df,
        ism_ref_matrix=ref_ism,
        ism_alt_matrix=alt_ism,
        metadata=meta,
        output_dir=figures_dir,
    )


def run_reproduction_pipeline(
    *,
    dotenv_path: pathlib.Path,
    cache_dir: pathlib.Path,
    processed_dir: pathlib.Path,
    figures_dir: pathlib.Path,
    raw_export_dir: pathlib.Path | None = None,
    offline: bool = False,
    refresh_cache: bool = False,
    model_version: dna_client.ModelVersion = dna_client.ModelVersion.ALL_FOLDS,
    data_provenance_label: str = "REAL_ALPHAGENOME_2_API",
    client_factory: Callable[[str], Any] | None = None,
) -> dict[str, Any]:
    """Runs validation, cached/live AlphaGenome 2 API inference, cross-validation, and plotting."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    if raw_export_dir is not None:
        raw_export_dir.mkdir(parents=True, exist_ok=True)

    # Step 2 (Part A): Build variant & intervals, and validate GRCh38 sequence + GENCODE v46 BEFORE any API call
    spec = build_bcl11a_variant_and_windows()
    ucsc_payload, bcl11a_gtf_df = _ensure_ucsc_and_gtf_cache(
        cache_dir, offline=offline, refresh_cache=refresh_cache
    )
    seq_verification = verify_reference_sequence(ucsc_payload, spec.variant)
    gene_verification = verify_enhancer_and_gene_annotation(bcl11a_gtf_df, spec.variant)

    # Lazy API client initialization: ONLY created if at least one cache artifact is missing or refresh_cache=True
    api_key_secret = ""
    if not offline:
        api_key_secret = load_api_key_safely(dotenv_path, environ=os.environ)

    _client_instance: Any = None
    api_calls_made = 0

    def _get_client() -> Any:
        nonlocal _client_instance
        if offline:
            raise RuntimeError("Cannot invoke AlphaGenome 2 API while `offline=True`.")
        if _client_instance is None:
            if client_factory is not None:
                _client_instance = client_factory(api_key_secret)
            else:
                _client_instance = dna_client.create(
                    api_key_secret,
                    model_version=model_version,
                    timeout=15.0,
                )
        return _client_instance

    # Step 1: Query or load `supported_conditions()` and verify exact condition CURIEs
    cond_cache_csv = cache_dir / "supported_conditions_human.csv"
    if refresh_cache or not cond_cache_csv.exists():
        cond_df = _get_client().supported_conditions()
        api_calls_made += 1
        cond_df.to_csv(cond_cache_csv, index=False)
    else:
        cond_df = pd.read_csv(cond_cache_csv)

    target_conditions = [CONDITION_BM_ERYTHROBLAST, CONDITION_BLOOD_NAIVE_B]
    verified_conditions = verify_supported_conditions(cond_df, target_conditions)

    # Also record tutorial B-cell condition if present in supported_df for documentation transparency
    tutorial_b_info = None
    if CONDITION_BM_NAIVE_B_TUTORIAL in set(cond_df["condition"]):
        tutorial_b_info = verify_supported_conditions(
            cond_df, [CONDITION_BM_NAIVE_B_TUTORIAL]
        )[CONDITION_BM_NAIVE_B_TUTORIAL]

    # Step 3 & Step 4: `predict_variant()` and `score_variant()` for Condition A & Condition B
    sc_atac_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_ATAC"]
    condition_arrays: dict[str, dict[str, Any]] = {}
    cross_val_tables: list[pd.DataFrame] = []
    raw_track_meta_tables: list[pd.DataFrame] = []

    for cond_id in target_conditions:
        slug = cond_id.replace(":", "_")
        npz_path = cache_dir / f"predict_variant_{slug}.npz"
        meta_csv_path = cache_dir / f"predict_variant_{slug}_metadata.csv"
        score_csv_path = cache_dir / f"score_variant_{slug}.csv"

        if refresh_cache or not (npz_path.exists() and meta_csv_path.exists()):
            v_out = _get_client().predict_variant(
                interval=spec.input_interval_1mb,
                variant=spec.variant,
                requested_outputs=[dna_client.OutputType.SC_ATAC],
                condition=cond_id,
            )
            api_calls_made += 1
            ref_sc = v_out.reference.sc_atac
            alt_sc = v_out.alternate.sc_atac
            if ref_sc is None or alt_sc is None:
                raise RuntimeError(f"SC_ATAC output was None for condition {cond_id}.")
            np.savez_compressed(
                npz_path,
                reference_sc_atac=ref_sc.values,
                alternate_sc_atac=alt_sc.values,
                interval_start=np.int64(spec.input_interval_1mb.start),
                interval_end=np.int64(spec.input_interval_1mb.end),
            )
            ref_sc.metadata.to_csv(meta_csv_path, index=False)

        if refresh_cache or not score_csv_path.exists():
            s_out = _get_client().score_variant(
                interval=spec.input_interval_1mb,
                variant=spec.variant,
                condition=cond_id,
                variant_scorers=[sc_atac_scorer],
            )
            api_calls_made += 1
            tidy_s = variant_scorers.tidy_scores([s_out])
            if tidy_s is None or tidy_s.empty:
                raise RuntimeError(f"Empty tidy_scores returned for {cond_id}.")
            tidy_s.to_csv(score_csv_path, index=False)

        with np.load(npz_path) as npz_data:
            ref_arr = npz_data["reference_sc_atac"]
            alt_arr = npz_data["alternate_sc_atac"]
        track_meta_df = pd.read_csv(meta_csv_path)
        score_tidy_df = pd.read_csv(score_csv_path)

        cv_df = cross_validate_scores(
            api_tidy_df=score_tidy_df,
            ref_values=ref_arr,
            alt_values=alt_arr,
            interval_start=spec.input_interval_1mb.start,
            variant_zero_based_start=spec.variant.start,
            width=SC_ATAC_CENTER_MASK_WIDTH,
            pseudocount=SC_ATAC_PSEUDOCOUNT,
            atol=1e-4,
        )
        cv_df.insert(0, "condition", cond_id)
        cross_val_tables.append(cv_df)
        raw_track_meta_tables.append(track_meta_df.assign(condition=cond_id))

        obs_mask = (track_meta_df["imputation_mode"] == "observed").to_numpy()
        if int(obs_mask.sum()) != 1:
            raise ValueError(
                f"Expected exactly 1 observed SC_ATAC track in {cond_id}, found {int(obs_mask.sum())}."
            )
        obs_idx = int(np.where(obs_mask)[0][0])
        condition_arrays[cond_id] = {
            "ref_all": ref_arr,
            "alt_all": alt_arr,
            "obs_idx": obs_idx,
            "obs_track_name": str(track_meta_df["name"].iloc[obs_idx]),
            "metadata": track_meta_df,
        }

    # Saturated ISM across Paper Fig. 4j window (chr2:60495262-60495285) for BM erythroblast REF & ALT
    ery_slug = CONDITION_BM_ERYTHROBLAST.replace(":", "_")
    ism_ref_csv = cache_dir / f"ism_ref_{ery_slug}.csv"
    ism_alt_csv = cache_dir / f"ism_alt_{ery_slug}.csv"

    for tag, interval_var, out_csv in [
        ("ref", None, ism_ref_csv),
        ("alt", spec.variant, ism_alt_csv),
    ]:
        if refresh_cache or not out_csv.exists():
            ism_scores = _get_client().score_ism_variants(
                interval=spec.input_interval_1mb,
                ism_interval=spec.ism_interval_23bp,
                condition=CONDITION_BM_ERYTHROBLAST,
                variant_scorers=[sc_atac_scorer],
                interval_variant=interval_var,
                progress_bar=False,
            )
            api_calls_made += 1
            ism_tidy = ism_anndata_list_to_tidy_df(ism_scores, [sc_atac_scorer])
            ism_tidy.to_csv(out_csv, index=False)

    ism_ref_raw_df = pd.read_csv(ism_ref_csv)
    ism_alt_raw_df = pd.read_csv(ism_alt_csv)
    ery_obs_track = condition_arrays[CONDITION_BM_ERYTHROBLAST]["obs_track_name"]

    ism_ref_mat, ism_ref_attr_df = build_ism_contribution_matrix(
        ism_tidy_df=ism_ref_raw_df,
        ism_interval=spec.ism_interval_23bp,
        sequence=seq_verification["ism_ref_sequence"],
        track_name=ery_obs_track,
    )
    ism_alt_mat, ism_alt_attr_df = build_ism_contribution_matrix(
        ism_tidy_df=ism_alt_raw_df,
        ism_interval=spec.ism_interval_23bp,
        sequence=seq_verification["ism_alt_sequence"],
        track_name=ery_obs_track,
    )

    # Export compact raw API records to `raw_export_dir` if configured
    if raw_export_dir is not None:
        (raw_export_dir / "ucsc_hg38_chr2_60495000_60495600.json").write_text(
            json.dumps(ucsc_payload, indent=2), encoding="utf-8"
        )
        cond_subset = cond_df[
            cond_df["condition"].isin(
                [CONDITION_BM_ERYTHROBLAST, CONDITION_BLOOD_NAIVE_B, CONDITION_BM_NAIVE_B_TUTORIAL]
            )
        ]
        cond_subset.to_csv(raw_export_dir / "supported_conditions_selected.csv", index=False)
        pd.concat(raw_track_meta_tables, ignore_index=True).to_csv(
            raw_export_dir / "predict_variant_sc_atac_track_metadata.csv", index=False
        )
        for cond_id in target_conditions:
            slug = cond_id.replace(":", "_")
            pd.read_csv(cache_dir / f"score_variant_{slug}.csv").to_csv(
                raw_export_dir / f"score_variant_{slug}.csv", index=False
            )
        ism_ref_raw_df.to_csv(raw_export_dir / f"ism_ref_{ery_slug}.csv", index=False)
        ism_alt_raw_df.to_csv(raw_export_dir / f"ism_alt_{ery_slug}.csv", index=False)

    # Build reusable 10-kb processed locus slice (`chr2:60,490,000-60,500,000`) & 1201-bp CSV
    slice_s, slice_e = SAVED_LOCUS_SLICE_0BASED
    rel_s = slice_s - spec.input_interval_1mb.start
    rel_e = slice_e - spec.input_interval_1mb.start
    pos_0 = np.arange(slice_s, slice_e, dtype=np.int64)

    ery_info = condition_arrays[CONDITION_BM_ERYTHROBLAST]
    b_info = condition_arrays[CONDITION_BLOOD_NAIVE_B]
    tracks_10kb_cols: dict[str, Any] = {
        "position_0based": pos_0,
        "position_1based": pos_0 + 1,
        "bm_erythroblast_ref": ery_info["ref_all"][rel_s:rel_e, ery_info["obs_idx"]],
        "bm_erythroblast_alt": ery_info["alt_all"][rel_s:rel_e, ery_info["obs_idx"]],
        "blood_naive_b_ref": b_info["ref_all"][rel_s:rel_e, b_info["obs_idx"]],
        "blood_naive_b_alt": b_info["alt_all"][rel_s:rel_e, b_info["obs_idx"]],
    }
    # Also retain all individual channels in the 10-kb parquet for complete reproducibility
    for prefix, info in [("bm_erythroblast", ery_info), ("blood_naive_b", b_info)]:
        for ch_idx, row in info["metadata"].iterrows():
            short_ch = str(row["name"]).removeprefix(f"{row['biosample_curie']}_")
            col_slug = short_ch.replace(" ", "_").replace("-", "_")
            tracks_10kb_cols[f"{prefix}_ch{ch_idx}_{col_slug}_ref"] = info["ref_all"][rel_s:rel_e, ch_idx]
            tracks_10kb_cols[f"{prefix}_ch{ch_idx}_{col_slug}_alt"] = info["alt_all"][rel_s:rel_e, ch_idx]

    tracks_df = pd.DataFrame(tracks_10kb_cols)
    tracks_parquet_path = processed_dir / "sc_atac_window_10kb_tracks.parquet"
    tracks_df.to_parquet(tracks_parquet_path, index=False)

    # Also save 1201-bp observed tracks as plain CSV for immediate inspection
    tracks_1201_df = tracks_df[
        (tracks_df["position_0based"] >= spec.scoring_window_1201bp.start)
        & (tracks_df["position_0based"] < spec.scoring_window_1201bp.end)
    ][
        [
            "position_0based",
            "position_1based",
            "bm_erythroblast_ref",
            "bm_erythroblast_alt",
            "blood_naive_b_ref",
            "blood_naive_b_alt",
        ]
    ]
    tracks_1201_df.to_csv(processed_dir / "sc_atac_1201bp_observed.csv", index=False)

    scores_comparison_df = pd.concat(cross_val_tables, ignore_index=True)
    scores_csv_path = processed_dir / "variant_effect_scores_comparison.csv"
    scores_comparison_df.to_csv(scores_csv_path, index=False)

    ism_combined_df = pd.concat(
        [
            ism_ref_attr_df.assign(allele_background="REF"),
            ism_alt_attr_df.assign(allele_background="ALT"),
        ],
        ignore_index=True,
    )
    ism_csv_path = processed_dir / "ism_attribution_chr2_60495262_60495285.csv"
    ism_combined_df.to_csv(ism_csv_path, index=False)

    bcl11a_gtf_df.to_parquet(processed_dir / "bcl11a_gencode_v46.parquet", index=False)

    metadata: dict[str, Any] = {
        "data_provenance": data_provenance_label,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "alphagenome_package_version": ag_version,
        "alphagenome_version": ag_version,
        "api_client_class": "alphagenome.models.v2.dna_client.DnaClient",
        "grpc_service_package": "google.gdm.gdmscience.alphagenome.v2main",
        "model_version": model_version.name,
        "api_key_configured": bool(api_key_secret) if not offline else "offline_mode",
        "api_calls_executed_in_this_run": api_calls_made,
        "reference_assembly": REFERENCE_ASSEMBLY,
        "variant": VARIANT_STR,
        "variant_coordinates": {
            "chromosome": spec.variant.chromosome,
            "position_1based": spec.variant.position,
            "interval_0based_half_open": [spec.variant.start, spec.variant.end],
            "reference_bases_plus_strand": spec.variant.reference_bases,
            "alternate_bases_plus_strand": spec.variant.alternate_bases,
        },
        "windows_0based_half_open": {
            "input_1mb": [spec.input_interval_1mb.start, spec.input_interval_1mb.end],
            "scoring_center_mask_1201bp": [
                spec.scoring_window_1201bp.start,
                spec.scoring_window_1201bp.end,
            ],
            "enhancer_62kb_paper": list(ENHANCER_62KB_PAPER),
            "enhancer_58kb_core_paper_fig4k": list(ENHANCER_58KB_CORE),
            "enhancer_55kb_paper": list(ENHANCER_55KB_PAPER),
            "enhancer_62kb_dhs_tutorial": list(ENHANCER_62KB_DHS_TUTORIAL),
            "enhancer_58kb_dhs_tutorial": list(ENHANCER_58KB_DHS_TUTORIAL),
            "enhancer_55kb_dhs_tutorial": list(ENHANCER_55KB_DHS_TUTORIAL),
            "sg1617_grna_target": list(SG1617_GRNA_TARGET_0BASED),
            "ism_fig4j_23bp": [spec.ism_interval_23bp.start, spec.ism_interval_23bp.end],
            "saved_locus_slice_10kb": list(SAVED_LOCUS_SLICE_0BASED),
        },
        "sequence_verification": seq_verification,
        "gene_and_enhancer_verification": gene_verification,
        "conditions": {
            "condition_a_main": verified_conditions[CONDITION_BM_ERYTHROBLAST],
            "condition_b_control": verified_conditions[CONDITION_BLOOD_NAIVE_B],
            "tutorial_b_cell_comparison_note": tutorial_b_info,
            "primary_comparable_channel_suffix": OBSERVED_SC_ATAC_TRACK_SUFFIX,
        },
        "scorer_parameters": {
            "scorer_key": "SC_ATAC",
            "scorer_repr": str(sc_atac_scorer),
            "requested_output": sc_atac_scorer.requested_output.name,
            "width_bp": sc_atac_scorer.width,
            "aggregation_type": sc_atac_scorer.aggregation_type.name,
            "formula": "log2(sum_{i in W_1201}(y_ALT,i) + 1) - log2(sum_{i in W_1201}(y_REF,i) + 1)",
            "pseudocount": SC_ATAC_PSEUDOCOUNT,
            "max_abs_cross_validation_diff": float(scores_comparison_df["abs_diff"].max()),
        },
        "checksums_sha256": {
            "sc_atac_window_10kb_tracks.parquet": _sha256_file(tracks_parquet_path),
            "variant_effect_scores_comparison.csv": _sha256_file(scores_csv_path),
            "ism_attribution_chr2_60495262_60495285.csv": _sha256_file(ism_csv_path),
        },
    }

    assert_no_secret_leak(metadata, api_key_secret)
    meta_path = processed_dir / "run_metadata.json"
    meta_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    figures = generate_all_figures(
        tracks_df=tracks_df,
        scores_df=scores_comparison_df,
        bcl11a_gtf_df=bcl11a_gtf_df,
        ism_ref_matrix=ism_ref_mat,
        ism_alt_matrix=ism_alt_mat,
        metadata=metadata,
        output_dir=figures_dir,
    )

    return {
        "data_provenance": data_provenance_label,
        "metadata_path": str(meta_path),
        "scores_csv_path": str(scores_csv_path),
        "tracks_parquet_path": str(tracks_parquet_path),
        "ism_csv_path": str(ism_csv_path),
        "figures": {k: {fmt: str(p) for fmt, p in v.items()} for k, v in figures.items()},
        "api_calls_made": api_calls_made,
    }
