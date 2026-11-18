"""Offline unit tests for AlphaGenome 2 BCL11A enhancer reproduction.

NOTE ON DATA PROVENANCE:
All arrays, sequences, GTF rows, and score tables constructed in this test module
are SYNTHETIC TEST FIXTURES (`SYNTHETIC_*`) used strictly to test coordinate
conventions, guardrails, scoring math, secret redaction, and plotting code
offline without an API key. They are never mixed with real AlphaGenome 2 outputs.
"""

from __future__ import annotations

import json
import pathlib
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bcl11a_repro.constants import (  # noqa: E402
    BCL11A_GENE_ID,
    BCL11A_MANE_TSS_1BASED_OR_0POINT,
    CONDITION_BLOOD_NAIVE_B,
    CONDITION_BM_ERYTHROBLAST,
    ENHANCER_58KB_CORE,
    ISM_WINDOW_0BASED,
    OBSERVED_SC_ATAC_TRACK_SUFFIX,
    VARIANT_ALT,
    VARIANT_CHROM,
    VARIANT_POS_1BASED,
    VARIANT_REF,
)
from bcl11a_repro.plotting import generate_all_figures  # noqa: E402
from bcl11a_repro.pipeline import run_reproduction_pipeline  # noqa: E402
from bcl11a_repro.scoring import (  # noqa: E402
    build_ism_contribution_matrix,
    compute_center_mask_log2_sum_score,
    cross_validate_scores,
    ism_anndata_list_to_tidy_df,
)
from bcl11a_repro.validation import (  # noqa: E402
    assert_no_secret_leak,
    build_bcl11a_variant_and_windows,
    load_api_key_safely,
    verify_enhancer_and_gene_annotation,
    verify_reference_sequence,
    verify_supported_conditions,
)


def _make_synthetic_ucsc_payload(
    *,
    genome_build: str = "hg38",
    ref_base_at_target: str = "G",
    motif_tail: str = "ATAA",
) -> dict[str, object]:
    """Builds a SYNTHETIC UCSC sequence response spanning chr2:60495000-60495600."""
    start = 60495000
    end = 60495600
    seq = ["C"] * (end - start)
    rel = VARIANT_POS_1BASED - 1 - start  # 265
    seq[rel] = ref_base_at_target
    for i, ch in enumerate(motif_tail, start=1):
        seq[rel + i] = ch
    return {
        "genome": genome_build,
        "chrom": "chr2",
        "start": start,
        "end": end,
        "dna": "".join(seq),
    }


def _make_synthetic_bcl11a_gtf(strand: str = "-") -> pd.DataFrame:
    """Builds a SYNTHETIC minimal GTF DataFrame matching GENCODE v46 schema."""
    return pd.DataFrame(
        [
            {
                "Chromosome": "chr2",
                "Start": 60450519,
                "End": 60554467,
                "Strand": strand,
                "Feature": "gene",
                "gene_id": "ENSG00000119866.22",
                "gene_name": "BCL11A",
                "transcript_id": None,
                "tag": None,
            },
            {
                "Chromosome": "chr2",
                "Start": 60457193,
                "End": 60553654,
                "Strand": strand,
                "Feature": "transcript",
                "gene_id": "ENSG00000119866.22",
                "gene_name": "BCL11A",
                "transcript_id": "ENST00000642384.2",
                "tag": "basic,Ensembl_canonical,MANE_Select",
            },
            {
                "Chromosome": "chr2",
                "Start": 60553152,
                "End": 60553654,
                "Strand": strand,
                "Feature": "exon",
                "gene_id": "ENSG00000119866.22",
                "gene_name": "BCL11A",
                "transcript_id": "ENST00000642384.2",
                "tag": "basic,Ensembl_canonical,MANE_Select",
            },
            {
                "Chromosome": "chr2",
                "Start": 60457193,
                "End": 60460600,
                "Strand": strand,
                "Feature": "exon",
                "gene_id": "ENSG00000119866.22",
                "gene_name": "BCL11A",
                "transcript_id": "ENST00000642384.2",
                "tag": "basic,Ensembl_canonical,MANE_Select",
            },
        ]
    )


class TestCoordinatesAndValidation(unittest.TestCase):
    def test_variant_and_windows_coordinates(self) -> None:
        spec = build_bcl11a_variant_and_windows()
        self.assertEqual(str(spec.variant), "chr2:60495266:G>A")
        self.assertEqual(spec.variant.chromosome, VARIANT_CHROM)
        self.assertEqual(spec.variant.position, VARIANT_POS_1BASED)
        self.assertEqual(spec.variant.reference_bases, VARIANT_REF)
        self.assertEqual(spec.variant.alternate_bases, VARIANT_ALT)
        # 0-based half-open single-base interval
        self.assertEqual(spec.variant.start, 60495265)
        self.assertEqual(spec.variant.end, 60495266)
        # 1-Mb input interval (2^20 = 1,048,576 bp)
        self.assertEqual(spec.input_interval_1mb.start, 59970978)
        self.assertEqual(spec.input_interval_1mb.end, 61019554)
        self.assertEqual(spec.input_interval_1mb.width, 1_048_576)
        # Symmetric 1201-bp center-mask window [pos_0based - 600, pos_0based + 601)
        self.assertEqual(spec.scoring_window_1201bp.start, 60494665)
        self.assertEqual(spec.scoring_window_1201bp.end, 60495866)
        self.assertEqual(spec.scoring_window_1201bp.width, 1201)

    def test_verify_reference_sequence_passes_on_valid_hg38_g(self) -> None:
        spec = build_bcl11a_variant_and_windows()
        payload = _make_synthetic_ucsc_payload(genome_build="hg38", ref_base_at_target="G")
        info = verify_reference_sequence(payload, spec.variant)
        self.assertEqual(info["verified_ref_base"], "G")
        self.assertEqual(info["plus_strand_motif_60495266_60495270"], "GATAA")

    def test_verify_reference_sequence_rejects_wrong_ref_base_or_assembly(self) -> None:
        spec = build_bcl11a_variant_and_windows()
        bad_base = _make_synthetic_ucsc_payload(genome_build="hg38", ref_base_at_target="C")
        with self.assertRaises(ValueError):
            verify_reference_sequence(bad_base, spec.variant)

        bad_genome = _make_synthetic_ucsc_payload(genome_build="hg19", ref_base_at_target="G")
        with self.assertRaises(ValueError):
            verify_reference_sequence(bad_genome, spec.variant)

    def test_verify_enhancer_and_gene_annotation_minus_strand(self) -> None:
        spec = build_bcl11a_variant_and_windows()
        gtf_df = _make_synthetic_bcl11a_gtf(strand="-")
        res = verify_enhancer_and_gene_annotation(gtf_df, spec.variant)
        self.assertEqual(res["gene_id_unversioned"], BCL11A_GENE_ID)
        self.assertEqual(res["gene_strand"], "-")
        self.assertEqual(res["mane_tss"], BCL11A_MANE_TSS_1BASED_OR_0POINT)
        # Downstream distance along minus strand: 60553654 - 60495266 = +58388 bp
        self.assertEqual(res["signed_distance_from_tss_bp"], 58388)
        self.assertTrue(res["in_58kb_core_enhancer"])
        self.assertTrue(res["in_58kb_dhs_enhancer"])

        # Plus-strand corruption must fail rather than flip alleles or coordinates
        bad_gtf = _make_synthetic_bcl11a_gtf(strand="+")
        with self.assertRaises(ValueError):
            verify_enhancer_and_gene_annotation(bad_gtf, spec.variant)

    def test_verify_supported_conditions_strict_lookup(self) -> None:
        cond_df = pd.DataFrame(
            [
                {
                    "condition": CONDITION_BM_ERYTHROBLAST,
                    "biosample_type": "primary_cell",
                    "biosample_name": "bone marrow erythroblast",
                    "sc_atac": "observed",
                    "sc_gene_expression": "observed",
                },
                {
                    "condition": CONDITION_BLOOD_NAIVE_B,
                    "biosample_type": "primary_cell",
                    "biosample_name": "blood naive B cell",
                    "sc_atac": "observed",
                    "sc_gene_expression": "observed",
                },
            ]
        )
        verified = verify_supported_conditions(
            cond_df, [CONDITION_BM_ERYTHROBLAST, CONDITION_BLOOD_NAIVE_B]
        )
        self.assertEqual(len(verified), 2)
        self.assertEqual(verified[CONDITION_BM_ERYTHROBLAST]["sc_atac"], "observed")
        self.assertEqual(verified[CONDITION_BLOOD_NAIVE_B]["biosample_name"], "blood naive B cell")

        with self.assertRaises(ValueError):
            verify_supported_conditions(cond_df, ["UBERON:9999999_CL:9999999"])


class TestScoringAndCrossValidation(unittest.TestCase):
    def test_center_mask_scorer_math_sign_order_and_pseudocount(self) -> None:
        # Synthetic 2000-bp interval [1000, 3000) with variant at 0-based 2000
        interval_start = 1000
        var_start = 2000  # relative index 1000; 1201-bp window is [400, 1601)
        ref = np.zeros((2000, 2), dtype=np.float32)
        alt = np.zeros((2000, 2), dtype=np.float32)

        # Track 0: REF has total 999 inside window, ALT has total 249 inside window
        # log2(249 + 1) - log2(999 + 1) = log2(250) - log2(1000) = -2.0
        ref[1000, 0] = 999.0
        alt[1000, 0] = 249.0

        # Put huge noise OUTSIDE the 1201-bp window [400, 1601) to prove window isolation
        ref[100, 0] = 50000.0
        alt[1800, 0] = 90000.0

        # Track 1: Zero signal in both REF and ALT -> log2(0+1) - log2(0+1) = 0.0
        res = compute_center_mask_log2_sum_score(
            ref_values=ref,
            alt_values=alt,
            interval_start=interval_start,
            variant_zero_based_start=var_start,
            width=1201,
            pseudocount=1.0,
        )
        np.testing.assert_allclose(res["ref_sum"], [999.0, 0.0], rtol=0, atol=1e-6)
        np.testing.assert_allclose(res["alt_sum"], [249.0, 0.0], rtol=0, atol=1e-6)
        np.testing.assert_allclose(res["delta_log2"], [-2.0, 0.0], rtol=0, atol=1e-6)

    def test_cross_validate_scores_detects_sign_flip(self) -> None:
        interval_start = 1000
        var_start = 2000
        ref = np.zeros((2000, 1), dtype=np.float32)
        alt = np.zeros((2000, 1), dtype=np.float32)
        ref[1000, 0] = 31.0
        alt[1000, 0] = 7.0  # log2(8) - log2(32) = 3 - 5 = -2.0

        good_api_df = pd.DataFrame(
            [{"track_name": "synthetic_track", "imputation_mode": "observed", "raw_score": -2.0}]
        )
        cv = cross_validate_scores(
            api_tidy_df=good_api_df,
            ref_values=ref,
            alt_values=alt,
            interval_start=interval_start,
            variant_zero_based_start=var_start,
            width=1201,
            atol=1e-5,
        )
        self.assertEqual(len(cv), 1)
        self.assertLess(cv["abs_diff"].iloc[0], 1e-6)

        flipped_api_df = pd.DataFrame(
            [{"track_name": "synthetic_track", "imputation_mode": "observed", "raw_score": +2.0}]
        )
        with self.assertRaises(ValueError):
            cross_validate_scores(
                api_tidy_df=flipped_api_df,
                ref_values=ref,
                alt_values=alt,
                interval_start=interval_start,
                variant_zero_based_start=var_start,
                width=1201,
                atol=1e-5,
            )

    def test_ism_contribution_matrix_mean_centering(self) -> None:
        from alphagenome.data import genome

        ism_interval = genome.Interval("chr2", 60495265, 60495268)  # 3 bp: G, A, T
        seq = "GAT"
        # For position 60495266 (ref G), alt alleles A, C, T all reduce accessibility by -3.0
        # Mean alt effect = -3.0 => ref G contribution = 0 - (-3.0) = +3.0!
        rows = [
            {"variant_id": "chr2:60495266:G>A", "track_name": "t_obs", "raw_score": -3.0},
            {"variant_id": "chr2:60495266:G>C", "track_name": "t_obs", "raw_score": -3.0},
            {"variant_id": "chr2:60495266:G>T", "track_name": "t_obs", "raw_score": -3.0},
            {"variant_id": "chr2:60495267:A>C", "track_name": "t_obs", "raw_score": 0.0},
            {"variant_id": "chr2:60495267:A>G", "track_name": "t_obs", "raw_score": 0.0},
            {"variant_id": "chr2:60495267:A>T", "track_name": "t_obs", "raw_score": 0.0},
            {"variant_id": "chr2:60495268:T>A", "track_name": "t_obs", "raw_score": -1.5},
            {"variant_id": "chr2:60495268:T>C", "track_name": "t_obs", "raw_score": -1.5},
            {"variant_id": "chr2:60495268:T>G", "track_name": "t_obs", "raw_score": -1.5},
        ]
        mat, df = build_ism_contribution_matrix(
            ism_tidy_df=pd.DataFrame(rows),
            ism_interval=ism_interval,
            sequence=seq,
            track_name="t_obs",
        )
        self.assertEqual(mat.shape, (3, 4))
        # Alphabet ACGT -> G is index 2, A is index 0, T is index 3
        self.assertAlmostEqual(float(mat[0, 2]), 3.0, places=5)
        self.assertAlmostEqual(float(mat[0, 0]), 0.0, places=5)
        self.assertAlmostEqual(float(mat[1, 0]), 0.0, places=5)
        self.assertAlmostEqual(float(mat[2, 3]), 1.5, places=5)
        self.assertEqual(len(df), 3)

    def test_ism_anndata_list_to_tidy_handles_sdk_missing_uns_scorer(self) -> None:
        """Regression test: SDK v2.0.0 score_ism_variants omits uns['variant_scorer']."""
        from alphagenome.data import genome
        from alphagenome.models.v2 import variant_scorers
        import anndata

        scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_ATAC"]
        var_obj = genome.Variant("chr2", 60495266, "G", "A")
        interval_obj = genome.Interval("chr2", 59970978, 61019554)
        adata = anndata.AnnData(
            X=np.array([[-2.193]], dtype=np.float32),
            obs=None,
            var=pd.DataFrame(
                [{"name": "t_obs", "strand": ".", "imputation_mode": "observed"}],
                index=["0"],
            ),
            uns={"interval": interval_obj, "variant": var_obj},  # No 'variant_scorer'!
        )
        tidy_df = ism_anndata_list_to_tidy_df([[adata]], [scorer])
        self.assertEqual(len(tidy_df), 1)
        self.assertEqual(str(tidy_df["variant_id"].iloc[0]), "chr2:60495266:G>A")
        self.assertAlmostEqual(float(tidy_df["raw_score"].iloc[0]), -2.193, places=5)


class TestSecurityAndPlottingOffline(unittest.TestCase):
    def test_api_key_loading_and_secret_redaction(self) -> None:
        fake_secret = "SYNTHETIC_FAKE_KEY_987654321"
        with tempfile.TemporaryDirectory() as tmp:
            env_path = pathlib.Path(tmp) / ".env"
            env_path.write_text(f"ALPHAGENOME_API_KEY={fake_secret}\n", encoding="utf-8")
            loaded = load_api_key_safely(env_path, environ={})
            self.assertEqual(loaded, fake_secret)

            safe_meta = {"model_version": "ALL_FOLDS", "api_key_configured": True}
            assert_no_secret_leak(safe_meta, fake_secret)

            leaky_meta = {"header": fake_secret}
            with self.assertRaises(ValueError):
                assert_no_secret_leak(leaky_meta, fake_secret)

    def test_offline_figure_generation_from_synthetic_processed_data(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = pathlib.Path(tmp) / "figures"
            # Build minimal synthetic processed inputs clearly tagged as synthetic
            positions_0based = np.arange(60490000, 60500000, dtype=np.int64)
            n = len(positions_0based)
            peak = np.exp(-0.5 * ((positions_0based - 60495268) / 25.0) ** 2)
            tracks_df = pd.DataFrame(
                {
                    "position_0based": positions_0based,
                    "position_1based": positions_0based + 1,
                    "bm_erythroblast_ref": (25.0 * peak + 0.2).astype(np.float32),
                    "bm_erythroblast_alt": (5.0 * peak + 0.2).astype(np.float32),
                    "blood_naive_b_ref": np.full(n, 0.05, dtype=np.float32),
                    "blood_naive_b_alt": np.full(n, 0.05, dtype=np.float32),
                }
            )
            scores_df = pd.DataFrame(
                [
                    {
                        "condition": CONDITION_BM_ERYTHROBLAST,
                        "biosample_name": "bone marrow erythroblast",
                        "track_name": f"{CONDITION_BM_ERYTHROBLAST}_{OBSERVED_SC_ATAC_TRACK_SUFFIX}",
                        "Assay title": "scATAC-seq 10X",
                        "data_pipeline": "cellranger",
                        "imputation_mode": "observed",
                        "ref_sum_1201bp": 1100.0,
                        "alt_sum_1201bp": 240.0,
                        "api_raw_score": -2.19,
                        "manual_log2_score": -2.19,
                        "abs_diff": 0.0,
                        "pct_change_sum": -78.18,
                    },
                    {
                        "condition": CONDITION_BLOOD_NAIVE_B,
                        "biosample_name": "blood naive B cell",
                        "track_name": f"{CONDITION_BLOOD_NAIVE_B}_{OBSERVED_SC_ATAC_TRACK_SUFFIX}",
                        "Assay title": "scATAC-seq 10X",
                        "data_pipeline": "cellranger",
                        "imputation_mode": "observed",
                        "ref_sum_1201bp": 54.0,
                        "alt_sum_1201bp": 53.9,
                        "api_raw_score": -0.002,
                        "manual_log2_score": -0.002,
                        "abs_diff": 0.0,
                        "pct_change_sum": -0.18,
                    },
                ]
            )
            ism_len = ISM_WINDOW_0BASED[1] - ISM_WINDOW_0BASED[0]
            ism_ref_mat = np.zeros((ism_len, 4), dtype=np.float32)
            ism_ref_mat[3, 2] = 1.8  # synthetic G at 60495266
            ism_alt_mat = np.zeros((ism_len, 4), dtype=np.float32)

            gtf_df = _make_synthetic_bcl11a_gtf(strand="-")
            meta = {
                "data_source_type": "SYNTHETIC_UNIT_TEST_FIXTURE",
                "alphagenome_version": "2.0.0",
                "model_version": "ALL_FOLDS",
                "variant": "chr2:60495266:G>A",
                "enhancer_58kb_core": list(ENHANCER_58KB_CORE),
            }

            saved = generate_all_figures(
                tracks_df=tracks_df,
                scores_df=scores_df,
                bcl11a_gtf_df=gtf_df,
                ism_ref_matrix=ism_ref_mat,
                ism_alt_matrix=ism_alt_mat,
                metadata=meta,
                output_dir=out_dir,
            )
            self.assertEqual(len(saved), 4)
            for fig_key, paths in saved.items():
                self.assertTrue(paths["png"].exists(), f"Missing PNG for {fig_key}")
                self.assertTrue(paths["pdf"].exists(), f"Missing PDF for {fig_key}")
                self.assertGreater(paths["png"].stat().st_size, 5000)
                self.assertGreater(paths["pdf"].stat().st_size, 2000)


class TestPipelineCachingAndGuardrails(unittest.TestCase):
    def _populate_synthetic_raw_cache(self, cache_dir: pathlib.Path, ref_base: str = "G") -> None:
        """Populates cache_dir with SYNTHETIC raw artifacts to test offline/cache-hit execution."""
        from alphagenome.data import genome
        from alphagenome.interpretation import ism as ism_lib

        cache_dir.mkdir(parents=True, exist_ok=True)
        ucsc_payload = _make_synthetic_ucsc_payload(genome_build="hg38", ref_base_at_target=ref_base)
        (cache_dir / "ucsc_hg38_chr2_60495000_60495600.json").write_text(
            json.dumps(ucsc_payload), encoding="utf-8"
        )
        _make_synthetic_bcl11a_gtf(strand="-").to_parquet(
            cache_dir / "bcl11a_gencode_v46.parquet", index=False
        )
        pd.DataFrame(
            [
                {
                    "condition": CONDITION_BM_ERYTHROBLAST,
                    "biosample_type": "primary_cell",
                    "biosample_name": "bone marrow erythroblast",
                    "sc_atac": "observed",
                    "sc_gene_expression": "observed",
                },
                {
                    "condition": CONDITION_BLOOD_NAIVE_B,
                    "biosample_type": "primary_cell",
                    "biosample_name": "blood naive B cell",
                    "sc_atac": "observed",
                    "sc_gene_expression": "observed",
                },
            ]
        ).to_csv(cache_dir / "supported_conditions_human.csv", index=False)

        spec = build_bcl11a_variant_and_windows()
        seq_len = spec.input_interval_1mb.width
        rel_var = spec.variant.start - spec.input_interval_1mb.start

        for cond_id, bname, ref_peak, alt_peak in [
            (CONDITION_BM_ERYTHROBLAST, "bone marrow erythroblast", 31.0, 7.0),
            (CONDITION_BLOOD_NAIVE_B, "blood naive B cell", 0.0, 0.0),
        ]:
            slug = cond_id.replace(":", "_")
            ref_arr = np.zeros((seq_len, 1), dtype=np.float32)
            alt_arr = np.zeros((seq_len, 1), dtype=np.float32)
            ref_arr[rel_var, 0] = ref_peak
            alt_arr[rel_var, 0] = alt_peak
            np.savez_compressed(
                cache_dir / f"predict_variant_{slug}.npz",
                reference_sc_atac=ref_arr,
                alternate_sc_atac=alt_arr,
                interval_start=np.int64(spec.input_interval_1mb.start),
                interval_end=np.int64(spec.input_interval_1mb.end),
            )
            tname = f"{cond_id}_{OBSERVED_SC_ATAC_TRACK_SUFFIX}"
            pd.DataFrame(
                [
                    {
                        "name": tname,
                        "strand": ".",
                        "biosample_curie": cond_id,
                        "imputation_mode": "observed",
                        "biosample_type": "primary_cell",
                        "biosample_name": bname,
                        "Assay title": "scATAC-seq 10X",
                        "suspension_type": "nuclear",
                        "data_pipeline": "cellranger",
                    }
                ]
            ).to_csv(cache_dir / f"predict_variant_{slug}_metadata.csv", index=False)

            expected_score = float(np.log2(alt_peak + 1.0) - np.log2(ref_peak + 1.0))
            pd.DataFrame(
                [
                    {
                        "variant_id": str(spec.variant),
                        "scored_interval": str(spec.input_interval_1mb),
                        "output_type": "SC_ATAC",
                        "variant_scorer": "CenterMaskScorer(SC_ATAC, 1201, DIFF_LOG2_SUM)",
                        "track_name": tname,
                        "track_strand": ".",
                        "imputation_mode": "observed",
                        "biosample_curie": cond_id,
                        "biosample_type": "primary_cell",
                        "biosample_name": bname,
                        "Assay title": "scATAC-seq 10X",
                        "suspension_type": "nuclear",
                        "data_pipeline": "cellranger",
                        "raw_score": expected_score,
                    }
                ]
            ).to_csv(cache_dir / f"score_variant_{slug}.csv", index=False)

        # Synthetic ISM cache for BM erythroblast
        ism_int = genome.Interval("chr2", ISM_WINDOW_0BASED[0], ISM_WINDOW_0BASED[1])
        dna_full = ucsc_payload["dna"]
        ism_ref_seq = str(dna_full)[ISM_WINDOW_0BASED[0] - 60495000 : ISM_WINDOW_0BASED[1] - 60495000]
        rel_in_ism = VARIANT_POS_1BASED - 1 - ISM_WINDOW_0BASED[0]
        ism_alt_seq = ism_ref_seq[:rel_in_ism] + "A" + ism_ref_seq[rel_in_ism + 1 :]
        tname_ery = f"{CONDITION_BM_ERYTHROBLAST}_{OBSERVED_SC_ATAC_TRACK_SUFFIX}"

        for tag, seq_str in [("ref", ism_ref_seq), ("alt", ism_alt_seq)]:
            snvs = ism_lib.ism_variants(ism_int, seq_str)
            pd.DataFrame(
                [
                    {
                        "variant_id": str(v),
                        "track_name": tname_ery,
                        "imputation_mode": "observed",
                        "raw_score": -1.5 if (tag == "ref" and v.position == VARIANT_POS_1BASED) else -0.05,
                    }
                    for v in snvs
                ]
            ).to_csv(
                cache_dir / f"ism_{tag}_{CONDITION_BM_ERYTHROBLAST.replace(':', '_')}.csv",
                index=False,
            )

    def test_pipeline_offline_cache_hit_never_calls_api_or_needs_key(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            cache_dir = base / "cache"
            proc_dir = base / "processed"
            fig_dir = base / "figures"
            self._populate_synthetic_raw_cache(cache_dir, ref_base="G")

            def _forbidden_client_factory(_api_key: str):
                raise AssertionError("API client should NEVER be called on offline cache hit!")

            summary = run_reproduction_pipeline(
                dotenv_path=base / "nonexistent.env",
                cache_dir=cache_dir,
                processed_dir=proc_dir,
                figures_dir=fig_dir,
                offline=True,
                refresh_cache=False,
                data_provenance_label="SYNTHETIC_UNIT_TEST_FIXTURE",
                client_factory=_forbidden_client_factory,
            )
            self.assertTrue((proc_dir / "run_metadata.json").exists())
            self.assertTrue((proc_dir / "variant_effect_scores_comparison.csv").exists())
            self.assertTrue((proc_dir / "sc_atac_window_10kb_tracks.parquet").exists())
            self.assertTrue((fig_dir / "fig_a_bcl11a_locus.png").exists())
            self.assertTrue((fig_dir / "fig_d_paper_comparison.pdf").exists())
            self.assertEqual(summary["data_provenance"], "SYNTHETIC_UNIT_TEST_FIXTURE")

    def test_pipeline_aborts_before_inference_if_ref_check_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            cache_dir = base / "cache"
            self._populate_synthetic_raw_cache(cache_dir, ref_base="T")
            env_file = base / ".env"
            env_file.write_text("ALPHAGENOME_API_KEY=SYNTHETIC_KEY_DO_NOT_CALL\n", encoding="utf-8")

            client_called = False

            def _spy_client_factory(_api_key: str):
                nonlocal client_called
                client_called = True
                raise AssertionError("Should abort before client creation!")

            with self.assertRaises(ValueError):
                run_reproduction_pipeline(
                    dotenv_path=env_file,
                    cache_dir=cache_dir,
                    processed_dir=base / "processed",
                    figures_dir=base / "figures",
                    offline=False,
                    refresh_cache=False,
                    client_factory=_spy_client_factory,
                )
            self.assertFalse(client_called)


if __name__ == "__main__":
    unittest.main()
