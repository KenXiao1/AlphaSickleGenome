"""Strict scientific & security validation checks prior to AlphaGenome 2 inference."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import dataclasses
import json
import pathlib
from typing import Any

from alphagenome.data import gene_annotation
from alphagenome.data import genome
from dotenv import dotenv_values
import pandas as pd

from .constants import (
    BCL11A_GENE_ID,
    BCL11A_MANE_TSS_1BASED_OR_0POINT,
    BCL11A_STRAND,
    ENHANCER_58KB_CORE,
    ENHANCER_58KB_DHS_TUTORIAL,
    GATAA_MOTIF_0BASED,
    INPUT_SEQUENCE_LENGTH_1MB,
    ISM_WINDOW_0BASED,
    SC_ATAC_CENTER_MASK_WIDTH,
    UCSC_GENOME_ID,
    VARIANT_ALT,
    VARIANT_CHROM,
    VARIANT_POS_1BASED,
    VARIANT_REF,
)


@dataclasses.dataclass(frozen=True)
class VariantWindowSpec:
    """Validated variant and genomic intervals in 0-based half-open API coordinates."""

    variant: genome.Variant
    reference_interval_0based: genome.Interval
    input_interval_1mb: genome.Interval
    scoring_window_1201bp: genome.Interval
    enhancer_58kb_core_interval: genome.Interval
    enhancer_58kb_dhs_interval: genome.Interval
    ism_interval_23bp: genome.Interval


def build_bcl11a_variant_and_windows() -> VariantWindowSpec:
    """Builds the canonical chr2:60495266:G>A variant and its official analysis windows."""
    variant = genome.Variant(
        chromosome=VARIANT_CHROM,
        position=VARIANT_POS_1BASED,
        reference_bases=VARIANT_REF,
        alternate_bases=VARIANT_ALT,
    )
    ref_interval = variant.reference_interval
    if (ref_interval.start, ref_interval.end) != (VARIANT_POS_1BASED - 1, VARIANT_POS_1BASED):
        raise ValueError(
            f"Unexpected 0-based interval for {variant}: {ref_interval}"
        )

    input_1mb = ref_interval.resize(INPUT_SEQUENCE_LENGTH_1MB)
    scoring_1201 = ref_interval.resize(SC_ATAC_CENTER_MASK_WIDTH)
    core_enhancer = genome.Interval(
        VARIANT_CHROM, ENHANCER_58KB_CORE[0], ENHANCER_58KB_CORE[1], name="+58 kb core"
    )
    dhs_enhancer = genome.Interval(
        VARIANT_CHROM,
        ENHANCER_58KB_DHS_TUTORIAL[0],
        ENHANCER_58KB_DHS_TUTORIAL[1],
        name="+58 kb DHS",
    )
    ism_interval = genome.Interval(
        VARIANT_CHROM, ISM_WINDOW_0BASED[0], ISM_WINDOW_0BASED[1], name="sg1617 ISM window"
    )

    if not input_1mb.contains(scoring_1201):
        raise ValueError("1-Mb input window does not contain 1201-bp scoring window.")
    if not scoring_1201.contains(core_enhancer):
        raise ValueError("1201-bp scoring window does not contain +58 kb core enhancer.")
    if not core_enhancer.contains(ref_interval):
        raise ValueError("+58 kb core enhancer does not contain variant position.")

    return VariantWindowSpec(
        variant=variant,
        reference_interval_0based=ref_interval,
        input_interval_1mb=input_1mb,
        scoring_window_1201bp=scoring_1201,
        enhancer_58kb_core_interval=core_enhancer,
        enhancer_58kb_dhs_interval=dhs_enhancer,
        ism_interval_23bp=ism_interval,
    )


def verify_reference_sequence(
    ucsc_payload: Mapping[str, Any],
    variant: genome.Variant,
) -> dict[str, Any]:
    """Verifies GRCh38/hg38 reference allele and positive-strand GATAA motif without auto-fixing."""
    genome_name = str(ucsc_payload.get("genome", ""))
    chrom = str(ucsc_payload.get("chrom", ""))
    start = int(ucsc_payload["start"])
    end = int(ucsc_payload["end"])
    dna = str(ucsc_payload.get("dna", "")).upper()

    if genome_name != UCSC_GENOME_ID:
        raise ValueError(
            f"Reference genome mismatch: expected '{UCSC_GENOME_ID}', got '{genome_name}'."
        )
    if chrom != variant.chromosome:
        raise ValueError(
            f"Chromosome mismatch: expected '{variant.chromosome}', got '{chrom}'."
        )
    if len(dna) != (end - start):
        raise ValueError(
            f"DNA length {len(dna)} does not match interval [{start}, {end}) width {end - start}."
        )
    if not (start <= variant.start and variant.end <= end):
        raise ValueError(
            f"UCSC sequence [{start}, {end}) does not cover variant [{variant.start}, {variant.end})."
        )

    observed_ref = dna[variant.start - start : variant.end - start]
    if observed_ref != variant.reference_bases:
        raise ValueError(
            f"GRCh38 reference base check FAILED at {variant.chromosome}:{variant.position}: "
            f"expected '{variant.reference_bases}', observed '{observed_ref}'. Aborting."
        )

    motif_start, motif_end = GATAA_MOTIF_0BASED
    motif_seq = dna[motif_start - start : motif_end - start]
    if motif_seq != "GATAA":
        raise ValueError(
            f"Expected positive-strand 'GATAA' motif at chr2:{motif_start + 1}-{motif_end}, "
            f"observed '{motif_seq}'."
        )

    ism_start, ism_end = ISM_WINDOW_0BASED
    ism_ref_seq = dna[ism_start - start : ism_end - start]
    rel_var_in_ism = variant.start - ism_start
    ism_alt_seq = (
        ism_ref_seq[:rel_var_in_ism]
        + variant.alternate_bases
        + ism_ref_seq[rel_var_in_ism + len(variant.reference_bases) :]
    )

    return {
        "verified_genome": genome_name,
        "verified_chrom": chrom,
        "verified_position_1based": variant.position,
        "verified_position_0based": variant.start,
        "verified_ref_base": observed_ref,
        "plus_strand_motif_60495266_60495270": motif_seq,
        "ism_window_0based": [ism_start, ism_end],
        "ism_ref_sequence": ism_ref_seq,
        "ism_alt_sequence": ism_alt_seq,
    }


def verify_enhancer_and_gene_annotation(
    bcl11a_gtf_df: pd.DataFrame,
    variant: genome.Variant,
) -> dict[str, Any]:
    """Verifies BCL11A strand, MANE Select TSS, and +58 kb intronic enhancer membership."""
    genes = bcl11a_gtf_df[bcl11a_gtf_df["Feature"] == "gene"]
    if len(genes) != 1:
        raise ValueError(f"Expected exactly 1 BCL11A gene row, found {len(genes)}.")

    gene_row = genes.iloc[0]
    gene_id_raw = str(gene_row["gene_id"])
    gene_id_unversioned = gene_id_raw.split(".", maxsplit=1)[0]
    gene_strand = str(gene_row["Strand"])
    gene_start = int(gene_row["Start"])
    gene_end = int(gene_row["End"])

    if gene_id_unversioned != BCL11A_GENE_ID:
        raise ValueError(f"Expected gene {BCL11A_GENE_ID}, got {gene_id_unversioned}.")
    if gene_strand != BCL11A_STRAND:
        raise ValueError(
            f"BCL11A must lie on negative strand ('-'), got '{gene_strand}'."
        )

    mane_df = gene_annotation.filter_to_mane_select_transcript(bcl11a_gtf_df)
    mane_tx = mane_df[mane_df["Feature"] == "transcript"]
    if len(mane_tx) != 1:
        raise ValueError(f"Expected 1 MANE_Select transcript for BCL11A, got {len(mane_tx)}.")

    tss_df = gene_annotation.extract_tss(mane_tx)
    mane_tss = int(tss_df["Start"].iloc[0])
    if mane_tss != BCL11A_MANE_TSS_1BASED_OR_0POINT:
        raise ValueError(
            f"Expected BCL11A MANE Select TSS {BCL11A_MANE_TSS_1BASED_OR_0POINT}, got {mane_tss}."
        )

    if not (gene_start <= variant.start < gene_end):
        raise ValueError(
            f"Variant {variant} lies outside BCL11A gene body [{gene_start}, {gene_end})."
        )

    # For a minus-strand gene, downstream of TSS corresponds to smaller genomic coordinates:
    signed_distance_from_tss = mane_tss - variant.position
    in_core = ENHANCER_58KB_CORE[0] <= variant.start < ENHANCER_58KB_CORE[1]
    in_dhs = ENHANCER_58KB_DHS_TUTORIAL[0] <= variant.start < ENHANCER_58KB_DHS_TUTORIAL[1]
    if not (in_core and in_dhs):
        raise ValueError(
            f"Variant {variant} does not fall inside the BCL11A +58 kb enhancer intervals."
        )

    return {
        "gene_id": gene_id_raw,
        "gene_id_unversioned": gene_id_unversioned,
        "gene_name": str(gene_row["gene_name"]),
        "gene_strand": gene_strand,
        "gene_body_0based": [gene_start, gene_end],
        "mane_transcript_id": str(mane_tx["transcript_id"].iloc[0]),
        "mane_transcript_0based": [int(mane_tx["Start"].iloc[0]), int(mane_tx["End"].iloc[0])],
        "mane_tss": mane_tss,
        "signed_distance_from_tss_bp": signed_distance_from_tss,
        "in_58kb_core_enhancer": in_core,
        "in_58kb_dhs_enhancer": in_dhs,
    }


def verify_supported_conditions(
    supported_df: pd.DataFrame,
    required_conditions: Sequence[str],
) -> dict[str, dict[str, Any]]:
    """Confirms exact condition identifiers and scATAC support from official supported_conditions()."""
    if "condition" not in supported_df.columns or "sc_atac" not in supported_df.columns:
        raise ValueError("supported_df must contain 'condition' and 'sc_atac' columns.")

    valid_modes = {"observed", "imputed_supported", "imputed_exploratory"}
    out: dict[str, dict[str, Any]] = {}
    for cond_id in required_conditions:
        matches = supported_df[supported_df["condition"] == cond_id]
        if matches.empty:
            raise ValueError(
                f"Condition '{cond_id}' was not found in official supported_conditions()."
            )
        if len(matches) > 1:
            raise ValueError(f"Duplicate rows found for condition '{cond_id}'.")
        row = matches.iloc[0]
        sc_atac_mode = str(row["sc_atac"])
        if sc_atac_mode not in valid_modes:
            raise ValueError(
                f"Condition '{cond_id}' has invalid sc_atac mode '{sc_atac_mode}'."
            )
        out[cond_id] = {
            "condition": cond_id,
            "biosample_type": str(row.get("biosample_type", "")),
            "biosample_name": str(row.get("biosample_name", "")),
            "sc_atac": sc_atac_mode,
            "sc_gene_expression": str(row.get("sc_gene_expression", "")),
        }
    return out


def load_api_key_safely(
    dotenv_path: pathlib.Path,
    environ: Mapping[str, str] | None = None,
) -> str:
    """Loads ALPHAGENOME_API_KEY from .env or environment without ever echoing the value."""
    file_vals = dotenv_values(dotenv_path) if dotenv_path.exists() else {}
    env_map = environ if environ is not None else {}
    key = (
        (file_vals.get("ALPHAGENOME_API_KEY") or "")
        or (env_map.get("ALPHAGENOME_API_KEY") or "")
    ).strip()
    if not key:
        raise RuntimeError(
            "ALPHAGENOME_API_KEY is not configured in .env or environment variables."
        )
    return key


def assert_no_secret_leak(payload: Any, secret: str) -> None:
    """Raises ValueError if the secret string appears anywhere in serialized payload."""
    if not secret:
        return
    serialized = payload if isinstance(payload, str) else json.dumps(payload, default=str)
    if secret in serialized:
        raise ValueError("SECURITY VIOLATION: API secret detected in serialized output.")
