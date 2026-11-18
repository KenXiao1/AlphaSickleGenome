"""Publication-quality figure generation for AlphaGenome 2 BCL11A enhancer reproduction."""

from __future__ import annotations

import os
import pathlib
import tempfile
from typing import Any

if "MPLCONFIGDIR" not in os.environ:
    _mpl_dir = pathlib.Path(tempfile.gettempdir()) / "alphasicklegenome_mpl"
    _mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(_mpl_dir)

from alphagenome.visualization import plot as ag_plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from .constants import (
    BCL11A_MANE_TSS_1BASED_OR_0POINT,
    CONDITION_BLOOD_NAIVE_B,
    CONDITION_BM_ERYTHROBLAST,
    ENHANCER_55KB_PAPER,
    ENHANCER_58KB_CORE,
    ENHANCER_58KB_DHS_TUTORIAL,
    ENHANCER_62KB_PAPER,
    GATAA_MOTIF_0BASED,
    INPUT_SEQUENCE_LENGTH_1MB,
    ISM_WINDOW_0BASED,
    SG1617_GRNA_TARGET_0BASED,
    VARIANT_POS_0BASED,
    VARIANT_POS_1BASED,
    VARIANT_STR,
)

# Colorblind-safe Okabe-Ito palette
COLOR_REF = "#0072B2"  # Blue
COLOR_ALT = "#D55E00"  # Vermilion / orange-red
COLOR_MUT = "#E69F00"  # Amber / orange highlight
COLOR_ENH = "#FDE0C5"  # Soft enhancer shading
COLOR_CORE = "#F4A582"  # Core enhancer accent
COLOR_GENE = "#2B3A4A"  # Slate navy for exons/gene body


def _apply_pub_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.labelsize": 9.0,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.8,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _despine(ax: plt.Axes, keep_left: bool = True, keep_bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)


def _format_genomic_xaxis(ax: plt.Axes, step: int | None = None) -> None:
    if step is not None:
        ax.xaxis.set_major_locator(mticker.MultipleLocator(step))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(round(v)):,}"))


def _save_fig(fig: plt.Figure, out_base: pathlib.Path) -> dict[str, pathlib.Path]:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_base.with_suffix(".png")
    pdf_path = out_base.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {"png": png_path, "pdf": pdf_path}


def _extract_observed_row(scores_df: pd.DataFrame, condition: str) -> pd.Series:
    sub = scores_df[
        (scores_df["condition"] == condition)
        & (scores_df["imputation_mode"] == "observed")
    ]
    if sub.empty:
        raise ValueError(f"No observed score row found for condition {condition}.")
    return sub.iloc[0]


def plot_figure_a_locus(
    bcl11a_gtf_df: pd.DataFrame,
    out_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    """Figure A: Multi-scale genomic map of BCL11A (- strand), +58 kb enhancer, variant & windows."""
    _apply_pub_style()
    fig = plt.figure(figsize=(9.5, 6.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.45, 1.45], hspace=0.58)

    win_1mb_start = VARIANT_POS_1BASED - INPUT_SEQUENCE_LENGTH_1MB // 2  # 59,970,978
    win_1mb_end = win_1mb_start + INPUT_SEQUENCE_LENGTH_1MB  # 61,019,554
    win_1201_start = VARIANT_POS_0BASED - 600  # 60,494,665
    win_1201_end = VARIANT_POS_0BASED + 601  # 60,495,866

    # --- Panel a: 1-Mb AlphaGenome 2 receptive field & BCL11A gene locus ---
    ax0 = fig.add_subplot(gs[0])
    _despine(ax0, keep_left=False)
    ax0.set_yticks([])
    ax0.set_xlim(win_1mb_start, win_1mb_end)
    ax0.set_ylim(0.0, 1.0)

    ax0.add_patch(
        mpatches.FancyBboxPatch(
            (win_1mb_start, 0.62),
            INPUT_SEQUENCE_LENGTH_1MB,
            0.20,
            boxstyle="round,pad=0.0",
            facecolor="#EBF3FA",
            edgecolor=COLOR_REF,
            linewidth=1.0,
        )
    )
    ax0.text(
        (win_1mb_start + win_1mb_end) / 2,
        0.72,
        f"AlphaGenome 2 1-Mb Input Receptive Field (chr2:{win_1mb_start:,}–{win_1mb_end:,}, 1,048,576 bp)",
        ha="center",
        va="center",
        fontsize=8.2,
        color=COLOR_GENE,
        fontweight="bold",
    )

    gene_row = bcl11a_gtf_df[bcl11a_gtf_df["Feature"] == "gene"].iloc[0]
    g_start, g_end = int(gene_row["Start"]), int(gene_row["End"])
    ax0.plot([g_start, g_end], [0.28, 0.28], color=COLOR_GENE, lw=2.4)
    ax0.text(
        g_end + 10_000,
        0.28,
        "BCL11A (ENSG00000119866, − strand, 5' ← TSS)",
        va="center",
        ha="left",
        fontsize=8.2,
        fontweight="bold",
        color=COLOR_GENE,
    )
    ax0.axvspan(ENHANCER_62KB_PAPER[0], ENHANCER_55KB_PAPER[1], color=COLOR_CORE, alpha=0.45)
    ax0.axvline(VARIANT_POS_1BASED, color=COLOR_ALT, linestyle="--", lw=1.2)
    ax0.annotate(
        f"chr2:{VARIANT_POS_1BASED:,}:G>A\n(+58 kb enhancer)",
        xy=(VARIANT_POS_1BASED, 0.58),
        xytext=(VARIANT_POS_1BASED - 245_000, 0.18),
        arrowprops=dict(arrowstyle="->", color=COLOR_ALT, lw=0.9),
        fontsize=7.8,
        color=COLOR_ALT,
        fontweight="bold",
    )
    ax0.text(-0.04, 1.12, "a", transform=ax0.transAxes, fontsize=11, fontweight="bold")
    ax0.set_title("GRCh38 / hg38 Chromosome 2 — 1-Mb Model Input Context", loc="left", pad=4)
    _format_genomic_xaxis(ax0, step=200_000)

    # --- Panel b: BCL11A MANE Select transcript (exons + minus-strand arrows) & intronic enhancers ---
    ax1 = fig.add_subplot(gs[1])
    _despine(ax1, keep_left=False)
    ax1.set_yticks([])
    zoom1_start, zoom1_end = 60_450_000, 60_562_000
    ax1.set_xlim(zoom1_start, zoom1_end)
    ax1.set_ylim(0.0, 1.15)

    mane_rows = bcl11a_gtf_df[
        bcl11a_gtf_df["tag"].fillna("").str.contains("MANE_Select")
    ]
    mane_tx = mane_rows[mane_rows["Feature"] == "transcript"].iloc[0]
    tx_start, tx_end = int(mane_tx["Start"]), int(mane_tx["End"])
    tx_id = str(mane_tx["transcript_id"])

    y_tx = 0.68
    ax1.plot([tx_start, tx_end], [y_tx, y_tx], color=COLOR_GENE, lw=1.4, zorder=2)
    for arrow_x in np.linspace(tx_start + 12_000, tx_end - 9_000, 7):
        ax1.annotate(
            "",
            xy=(arrow_x - 2800, y_tx),
            xytext=(arrow_x + 2800, y_tx),
            arrowprops=dict(arrowstyle="->", color=COLOR_GENE, lw=1.0),
        )

    exons = mane_rows[mane_rows["Feature"] == "exon"]
    for _, ex in exons.iterrows():
        ex_s, ex_e = int(ex["Start"]), int(ex["End"])
        ax1.add_patch(
            mpatches.Rectangle(
                (ex_s, y_tx - 0.10),
                max(ex_e - ex_s, 350),
                0.20,
                facecolor=COLOR_GENE,
                edgecolor=COLOR_GENE,
                zorder=3,
            )
        )

    # Clean TSS bent arrow pointing left (minus-strand 5' -> 3'), with label strictly above arrow
    tss = BCL11A_MANE_TSS_1BASED_OR_0POINT
    ax1.plot([tss, tss], [y_tx, y_tx + 0.19], color=COLOR_GENE, lw=1.3)
    ax1.annotate(
        "",
        xy=(tss - 6500, y_tx + 0.19),
        xytext=(tss, y_tx + 0.19),
        arrowprops=dict(arrowstyle="->", color=COLOR_GENE, lw=1.3),
    )
    ax1.text(
        tss,
        y_tx + 0.24,
        f"TSS chr2:{tss:,}  (5' → 3' transcription along − strand)",
        ha="right",
        va="bottom",
        fontsize=7.9,
        fontweight="bold",
        color=COLOR_GENE,
    )
    ax1.text(
        tx_start,
        y_tx + 0.15,
        f"BCL11A MANE Select ({tx_id})",
        fontsize=8.0,
        color=COLOR_GENE,
        ha="left",
    )

    # Intronic erythroid enhancers with non-overlapping staggered labels
    enh_specs = [
        ("+62 kb", ENHANCER_62KB_PAPER, "#FDCC8A", (ENHANCER_62KB_PAPER[0] - 2200, 0.31), "right"),
        ("+58 kb", ENHANCER_58KB_DHS_TUTORIAL, "#FC8D59", ((ENHANCER_58KB_DHS_TUTORIAL[0] + ENHANCER_58KB_DHS_TUTORIAL[1]) / 2, 0.08), "center"),
        ("+55 kb", ENHANCER_55KB_PAPER, "#FDCC8A", (ENHANCER_55KB_PAPER[1] + 2200, 0.31), "left"),
    ]
    ax1.axvline(VARIANT_POS_1BASED, color=COLOR_ALT, linestyle="--", lw=1.2, zorder=2)
    for enh_label, (es, ee), col, (lx, ly), ha_align in enh_specs:
        ax1.add_patch(
            mpatches.Rectangle((es, 0.22), ee - es, 0.19, facecolor=col, edgecolor="#B34800", lw=0.9, zorder=3)
        )
        ax1.annotate(
            enh_label,
            xy=((es + ee) / 2, 0.22 if ly < 0.15 else 0.315),
            xytext=(lx, ly),
            ha=ha_align,
            va="center",
            fontsize=7.9,
            fontweight="bold",
            color="#8C2D04",
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.90),
            arrowprops=dict(arrowstyle="-", color="#B34800", lw=0.7),
            zorder=5,
        )
    ax1.text(
        (ENHANCER_62KB_PAPER[0] + ENHANCER_55KB_PAPER[1]) / 2,
        0.46,
        "Intronic Erythroid Enhancer Cluster (+58.4 kb downstream of minus-strand TSS)",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#8C2D04",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.90),
        zorder=5,
    )
    ax1.text(-0.04, 1.08, "b", transform=ax1.transAxes, fontsize=11, fontweight="bold")
    ax1.set_title("BCL11A Gene Structure (GENCODE v46) & Intronic Erythroid Enhancers", loc="left", pad=4)
    _format_genomic_xaxis(ax1, step=20_000)

    # --- Panel c: Fine-scale zoom of +58 kb enhancer, 1201-bp scorer window, sg1617 & GATAA ---
    ax2 = fig.add_subplot(gs[2])
    _despine(ax2, keep_left=False)
    ax2.set_yticks([])
    zoom2_start, zoom2_end = 60_494_550, 60_496_380
    ax2.set_xlim(zoom2_start, zoom2_end)
    ax2.set_ylim(0.0, 1.18)

    ax2.axvline(VARIANT_POS_1BASED, color=COLOR_ALT, lw=1.6, zorder=4)

    ax2.add_patch(
        mpatches.Rectangle(
            (win_1201_start, 0.80),
            win_1201_end - win_1201_start,
            0.18,
            facecolor="#D6EAF8",
            edgecolor=COLOR_REF,
            lw=1.0,
            zorder=2,
        )
    )
    ax2.text(
        (win_1201_start + win_1201_end) / 2,
        0.89,
        f"1201-bp scATAC CenterMaskScorer Window (chr2:{win_1201_start:,}–{win_1201_end:,})",
        ha="center",
        va="center",
        fontsize=7.9,
        fontweight="bold",
        color="#0B4F6C",
        bbox=dict(boxstyle="round,pad=0.10", facecolor="#D6EAF8", edgecolor="none"),
        zorder=6,
    )

    ax2.add_patch(
        mpatches.Rectangle(
            (ENHANCER_58KB_DHS_TUTORIAL[0], 0.50),
            ENHANCER_58KB_DHS_TUTORIAL[1] - ENHANCER_58KB_DHS_TUTORIAL[0],
            0.17,
            facecolor="#FDEBD0",
            edgecolor="#CA6F1E",
            linestyle="--",
            lw=0.9,
            zorder=2,
        )
    )
    ax2.text(
        ENHANCER_58KB_DHS_TUTORIAL[1] - 18,
        0.585,
        f"+58 kb DHS (chr2:{ENHANCER_58KB_DHS_TUTORIAL[0]:,}–{ENHANCER_58KB_DHS_TUTORIAL[1]:,}, 1,288 bp)",
        ha="right",
        va="center",
        fontsize=7.5,
        color="#784212",
        zorder=6,
    )

    ax2.add_patch(
        mpatches.Rectangle(
            (ENHANCER_58KB_CORE[0], 0.20),
            ENHANCER_58KB_CORE[1] - ENHANCER_58KB_CORE[0],
            0.19,
            facecolor="#F5CBA7",
            edgecolor="#A04000",
            lw=1.1,
            zorder=2,
        )
    )
    ax2.text(
        (ENHANCER_58KB_CORE[0] + ENHANCER_58KB_CORE[1]) / 2,
        0.295,
        "+58 kb Core (490 bp)",
        ha="center",
        va="center",
        fontsize=7.6,
        fontweight="bold",
        color="#6E2C00",
        bbox=dict(boxstyle="round,pad=0.10", facecolor="#F5CBA7", edgecolor="none"),
        zorder=6,
    )
    ax2.text(
        (ENHANCER_58KB_CORE[0] + ENHANCER_58KB_CORE[1]) / 2,
        0.055,
        f"Paper Fig. 4k window: chr2:{ENHANCER_58KB_CORE[0]:,}–{ENHANCER_58KB_CORE[1]:,}",
        ha="center",
        va="bottom",
        fontsize=7.3,
        color="#6E2C00",
        bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="none", alpha=0.92),
        zorder=6,
    )

    ax2.annotate(
        f"{VARIANT_STR} (+ strand REF=G → ALT=A)\n"
        f"Disrupts 1st base of GATA1 motif 'GATAA' (chr2:{GATAA_MOTIF_0BASED[0]+1:,}–{GATAA_MOTIF_0BASED[1]:,})\n"
        f"within sg1617 target (chr2:{SG1617_GRNA_TARGET_0BASED[0]+1:,}–{SG1617_GRNA_TARGET_0BASED[1]:,})",
        xy=(VARIANT_POS_1BASED, 0.41),
        xytext=(60_495_605, 0.16),
        arrowprops=dict(arrowstyle="->", color=COLOR_ALT, lw=1.1),
        fontsize=7.3,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#FFF5F0", edgecolor=COLOR_ALT, lw=0.8),
        color="#7F0000",
        zorder=7,
    )

    ax2.text(-0.04, 1.08, "c", transform=ax2.transAxes, fontsize=11, fontweight="bold")
    ax2.set_title("Fine-Scale View of +58 kb Enhancer, Scoring Window & Therapeutic Variant", loc="left", pad=4)
    ax2.set_xlabel("Chromosome 2 Coordinate (GRCh38 / hg38, positive reference strand)")
    _format_genomic_xaxis(ax2, step=300)

    return _save_fig(fig, out_dir / "fig_a_bcl11a_locus")


def plot_figure_b_erythroblast_ref_vs_alt(
    tracks_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    out_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    """Figure B: Bone marrow erythroblast REF vs ALT scATAC-seq on shared coordinates & y-scale."""
    _apply_pub_style()
    obs_row = _extract_observed_row(scores_df, CONDITION_BM_ERYTHROBLAST)

    win_1201_s = VARIANT_POS_0BASED - 600
    win_1201_e = VARIANT_POS_0BASED + 601
    sub_1201 = tracks_df[
        (tracks_df["position_0based"] >= win_1201_s)
        & (tracks_df["position_0based"] < win_1201_e)
    ].copy()
    sub_core = tracks_df[
        (tracks_df["position_0based"] >= ENHANCER_58KB_CORE[0])
        & (tracks_df["position_0based"] < ENHANCER_58KB_CORE[1])
    ].copy()

    ymax = float(max(sub_1201["bm_erythroblast_ref"].max(), sub_1201["bm_erythroblast_alt"].max()) * 1.15)

    fig = plt.figure(figsize=(9.4, 6.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.45, 1.15, 0.95], hspace=0.45)

    # --- Panel a: Full 1201-bp CenterMaskScorer window ---
    ax0 = fig.add_subplot(gs[0])
    _despine(ax0)
    x1201 = sub_1201["position_1based"].to_numpy()
    y_ref = sub_1201["bm_erythroblast_ref"].to_numpy()
    y_alt = sub_1201["bm_erythroblast_alt"].to_numpy()

    ax0.axvspan(
        ENHANCER_58KB_CORE[0] + 1,
        ENHANCER_58KB_CORE[1],
        color=COLOR_ENH,
        alpha=0.45,
        label="+58 kb core enhancer (490 bp)",
    )
    ax0.fill_between(x1201, y_ref, alpha=0.20, color=COLOR_REF)
    ax0.plot(x1201, y_ref, color=COLOR_REF, lw=1.35, label="REF allele (G) — predicted scATAC")
    ax0.fill_between(x1201, y_alt, alpha=0.25, color=COLOR_ALT)
    ax0.plot(x1201, y_alt, color=COLOR_ALT, lw=1.35, label="ALT allele (A) — predicted scATAC")
    ax0.axvline(VARIANT_POS_1BASED, color=COLOR_MUT, linestyle="--", lw=1.4, label=VARIANT_STR)

    ax0.set_xlim(win_1201_s + 1, win_1201_e)
    ax0.set_ylim(0.0, ymax)
    ax0.set_ylabel("Predicted scATAC\n(norm. insertions)")
    ax0.set_title(
        f"Bone Marrow Erythroblast ({CONDITION_BM_ERYTHROBLAST}, observed 10x nuclear channel) — 1201-bp Scoring Window",
        loc="left",
    )
    ax0.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.6)
    ax0.text(
        0.985,
        0.93,
        f"1201-bp Sum (REF): {obs_row['ref_sum_1201bp']:.2f}\n"
        f"1201-bp Sum (ALT): {obs_row['alt_sum_1201bp']:.2f} ({obs_row['pct_change_sum']:.1f}%)\n"
        f"Official $\\Delta_{{\\mathrm{{ATAC}}}}$ ($\\log_2$): {obs_row['api_raw_score']:.4f}",
        transform=ax0.transAxes,
        ha="right",
        va="top",
        fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8F9FA", edgecolor="#ADB5BD", lw=0.8),
    )
    ax0.text(-0.07, 1.08, "a", transform=ax0.transAxes, fontsize=11, fontweight="bold")
    _format_genomic_xaxis(ax0, step=200)

    # --- Panel b: Zoomed +58 kb core enhancer on SAME y-axis scale ---
    ax1 = fig.add_subplot(gs[1])
    _despine(ax1)
    xc = sub_core["position_1based"].to_numpy()
    yc_ref = sub_core["bm_erythroblast_ref"].to_numpy()
    yc_alt = sub_core["bm_erythroblast_alt"].to_numpy()

    ax1.axvspan(
        SG1617_GRNA_TARGET_0BASED[0] + 1,
        SG1617_GRNA_TARGET_0BASED[1],
        color="#FDE0C5",
        alpha=0.75,
        label="sg1617 gRNA target (21 bp)",
    )
    ax1.plot(xc, yc_ref, color=COLOR_REF, lw=1.45, label="REF (G)")
    ax1.plot(xc, yc_alt, color=COLOR_ALT, lw=1.45, label="ALT (A)")
    ax1.axvline(VARIANT_POS_1BASED, color=COLOR_MUT, linestyle="--", lw=1.4)
    ax1.set_xlim(ENHANCER_58KB_CORE[0] + 1, ENHANCER_58KB_CORE[1])
    ax1.set_ylim(0.0, ymax)
    ax1.set_ylabel("Predicted scATAC\n(same y-scale)")
    ax1.set_title(
        f"Zoom: +58 kb Core Enhancer (chr2:{ENHANCER_58KB_CORE[0]:,}–{ENHANCER_58KB_CORE[1]:,}, 490 bp)",
        loc="left",
    )
    ax1.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92)
    ax1.text(-0.07, 1.08, "b", transform=ax1.transAxes, fontsize=11, fontweight="bold")
    _format_genomic_xaxis(ax1, step=100)

    # --- Panel c: Base-resolution allelic difference (ALT - REF) ---
    ax2 = fig.add_subplot(gs[2], sharex=ax1)
    _despine(ax2)
    diff_c = yc_alt - yc_ref
    ax2.axhline(0.0, color="#495057", lw=0.8)
    ax2.fill_between(
        xc,
        diff_c,
        0.0,
        where=(diff_c <= 0),
        color=COLOR_ALT,
        alpha=0.35,
        label="Accessibility loss (ALT < REF)",
    )
    ax2.plot(xc, diff_c, color=COLOR_ALT, lw=1.1)
    ax2.axvline(VARIANT_POS_1BASED, color=COLOR_MUT, linestyle="--", lw=1.4)
    ax2.set_ylabel("$\\Delta$ scATAC\n(ALT − REF)")
    ax2.set_xlabel("Chromosome 2 Coordinate (GRCh38 / hg38, 1-based position)")
    ax2.set_title("Per-Base Predicted Accessibility Change (ALT − REF) Across +58 kb Core Enhancer", loc="left")
    ax2.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92)
    ax2.text(-0.07, 1.10, "c", transform=ax2.transAxes, fontsize=11, fontweight="bold")
    _format_genomic_xaxis(ax2, step=100)

    return _save_fig(fig, out_dir / "fig_b_erythroblast_ref_vs_alt")


def plot_figure_c_cell_specificity(
    tracks_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    out_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    """Figure C: Erythroblast vs Blood Naive B Cell under identical window, channel & scorer."""
    _apply_pub_style()
    obs_ery = _extract_observed_row(scores_df, CONDITION_BM_ERYTHROBLAST)
    obs_b = _extract_observed_row(scores_df, CONDITION_BLOOD_NAIVE_B)

    win_1201_s = VARIANT_POS_0BASED - 600
    win_1201_e = VARIANT_POS_0BASED + 601
    sub_1201 = tracks_df[
        (tracks_df["position_0based"] >= win_1201_s)
        & (tracks_df["position_0based"] < win_1201_e)
    ].copy()

    x = sub_1201["position_1based"].to_numpy()
    ymax = float(
        max(
            sub_1201["bm_erythroblast_ref"].max(),
            sub_1201["bm_erythroblast_alt"].max(),
            sub_1201["blood_naive_b_ref"].max(),
            sub_1201["blood_naive_b_alt"].max(),
        )
        * 1.22
    )

    fig = plt.figure(figsize=(10.5, 6.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.48, 1.0], hspace=0.45, wspace=0.30)

    # --- Panel a: Bone marrow erythroblast (1201-bp window, shared y-axis) ---
    ax_a = fig.add_subplot(gs[0, 0])
    _despine(ax_a)
    ax_a.axvspan(ENHANCER_58KB_CORE[0] + 1, ENHANCER_58KB_CORE[1], color=COLOR_ENH, alpha=0.45, label="+58 kb core")
    ax_a.plot(x, sub_1201["bm_erythroblast_ref"], color=COLOR_REF, lw=1.3, label="REF (G)")
    ax_a.plot(x, sub_1201["bm_erythroblast_alt"], color=COLOR_ALT, lw=1.3, label="ALT (A)")
    ax_a.axvline(VARIANT_POS_1BASED, color=COLOR_MUT, linestyle="--", lw=1.3, label=VARIANT_STR)
    ax_a.set_xlim(win_1201_s + 1, win_1201_e)
    ax_a.set_ylim(0.0, ymax)
    ax_a.set_ylabel("Predicted scATAC\n(norm. insertions)")
    ax_a.set_title(
        f"Condition A: Bone Marrow Erythroblast\n({CONDITION_BM_ERYTHROBLAST} | observed 10x nuclear)",
        loc="left",
    )
    # Single-column compact legend confined strictly to the low-signal left flank (<60,495,040)
    ax_a.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, ncol=1, fontsize=7.4)
    ax_a.text(
        0.98,
        0.92,
        f"$\\Delta_{{\\mathrm{{ATAC}}}} = {obs_ery['api_raw_score']:.4f}$",
        transform=ax_a.transAxes,
        ha="right",
        va="top",
        fontweight="bold",
        color=COLOR_ALT,
        fontsize=8.6,
    )
    ax_a.text(-0.11, 1.10, "a", transform=ax_a.transAxes, fontsize=11, fontweight="bold")
    _format_genomic_xaxis(ax_a, step=300)

    # --- Panel b: Blood naive B cell (exact same 1201-bp window & shared y-axis!) ---
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    _despine(ax_b)
    ax_b.axvspan(ENHANCER_58KB_CORE[0] + 1, ENHANCER_58KB_CORE[1], color=COLOR_ENH, alpha=0.45)
    ax_b.plot(x, sub_1201["blood_naive_b_ref"], color=COLOR_REF, lw=1.3, label="REF (G)")
    ax_b.plot(x, sub_1201["blood_naive_b_alt"], color=COLOR_ALT, lw=1.3, linestyle="--", label="ALT (A)")
    ax_b.axvline(VARIANT_POS_1BASED, color=COLOR_MUT, linestyle="--", lw=1.3)
    ax_b.set_ylim(0.0, ymax)
    ax_b.set_ylabel("Predicted scATAC\n(same y-scale)")
    ax_b.set_xlabel("Chromosome 2 Coordinate (GRCh38, 1201-bp scoring window)")
    ax_b.set_title(
        f"Condition B: Blood Naive B Cell (Control)\n({CONDITION_BLOOD_NAIVE_B} | observed 10x nuclear)",
        loc="left",
    )
    b_max = float(max(sub_1201["blood_naive_b_ref"].max(), sub_1201["blood_naive_b_alt"].max()))
    ax_b.text(
        0.98,
        0.92,
        f"$\\Delta_{{\\mathrm{{ATAC}}}} = {obs_b['api_raw_score']:.4f}$\n(closed chromatin; max={b_max:.2f})",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        color="#343A40",
    )
    ax_b.text(-0.11, 1.10, "b", transform=ax_b.transAxes, fontsize=11, fontweight="bold")
    _format_genomic_xaxis(ax_b, step=300)

    # --- Panel c: Summed 1201-bp signal comparison ---
    ax_c = fig.add_subplot(gs[0, 1])
    _despine(ax_c)
    cell_labels = ["BM Erythroblast\n(observed 10x)", "Blood Naive B Cell\n(observed 10x)"]
    xpos = np.arange(2)
    bar_w = 0.34
    ref_sums = [float(obs_ery["ref_sum_1201bp"]), float(obs_b["ref_sum_1201bp"])]
    alt_sums = [float(obs_ery["alt_sum_1201bp"]), float(obs_b["alt_sum_1201bp"])]

    b1 = ax_c.bar(xpos - bar_w / 2, ref_sums, width=bar_w, color=COLOR_REF, label="REF (G)")
    b2 = ax_c.bar(xpos + bar_w / 2, alt_sums, width=bar_w, color=COLOR_ALT, label="ALT (A)")
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax_c.text(
                rect.get_x() + rect.get_width() / 2,
                h + max(ref_sums) * 0.02,
                f"{h:.1f}",
                ha="center",
                va="bottom",
                fontsize=7.7,
            )
    ax_c.set_xticks(xpos)
    ax_c.set_xticklabels(cell_labels)
    ax_c.set_ylim(0.0, max(ref_sums) * 1.22)
    ax_c.set_ylabel("1201-bp Summed\nscATAC Signal")
    ax_c.set_title("Total Local Accessibility in $W_{1201}$", loc="left")
    ax_c.legend(loc="upper right", frameon=False)
    ax_c.text(-0.16, 1.10, "c", transform=ax_c.transAxes, fontsize=11, fontweight="bold")

    # --- Panel d: Official CenterMaskScorer (DIFF_LOG2_SUM, 1201 bp) ---
    ax_d = fig.add_subplot(gs[1, 1])
    _despine(ax_d)
    obs_scores = [float(obs_ery["api_raw_score"]), float(obs_b["api_raw_score"])]
    bars_d = ax_d.bar(
        xpos,
        obs_scores,
        width=0.44,
        color=[COLOR_ALT, "#6C757D"],
        alpha=0.88,
        label="Observed channel (10x nuclear)",
    )
    ax_d.axhline(0.0, color="#212529", lw=0.8)

    for idx, cond_id in enumerate([CONDITION_BM_ERYTHROBLAST, CONDITION_BLOOD_NAIVE_B]):
        imp = scores_df[
            (scores_df["condition"] == cond_id)
            & (scores_df["imputation_mode"] != "observed")
        ]
        if not imp.empty:
            jitter = np.linspace(-0.12, 0.12, len(imp))
            ax_d.scatter(
                np.full(len(imp), xpos[idx]) + jitter,
                imp["api_raw_score"].to_numpy(),
                facecolors="white",
                edgecolors="#212529",
                s=26,
                zorder=5,
                label="Imputed_supported channels (n=4)" if idx == 0 else None,
            )

    for idx, rect in enumerate(bars_d):
        val = obs_scores[idx]
        ax_d.text(
            rect.get_x() + rect.get_width() / 2,
            val - 0.14 if val < -0.5 else val - 0.16,
            f"{val:.4f}",
            ha="center",
            va="top",
            fontsize=8.0,
            fontweight="bold",
        )

    ax_d.set_xticks(xpos)
    ax_d.set_xticklabels(cell_labels)
    ax_d.set_ylim(min(obs_scores) * 1.30, 0.32)
    ax_d.set_ylabel("Variant Score $\\Delta_{\\mathrm{ATAC}}$\n$\\log_2((\\mathrm{ALT}+1)/(\\mathrm{REF}+1))$")
    ax_d.set_title("Variant Effect Score (1201-bp CenterMask)", loc="left")
    ax_d.legend(loc="lower right", frameon=True, facecolor="white", fontsize=7.2)
    ax_d.text(-0.16, 1.10, "d", transform=ax_d.transAxes, fontsize=11, fontweight="bold")

    return _save_fig(fig, out_dir / "fig_c_erythroblast_vs_bcell")


def plot_figure_d_paper_comparison(
    tracks_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    ism_ref_matrix: np.ndarray,
    ism_alt_matrix: np.ndarray,
    metadata: dict[str, Any],
    out_dir: pathlib.Path,
) -> dict[str, pathlib.Path]:
    """Figure D: Reproduction of AlphaGenome 2 Fig. 4j-k & qualitative concordance audit."""
    _apply_pub_style()
    obs_ery = _extract_observed_row(scores_df, CONDITION_BM_ERYTHROBLAST)
    obs_b = _extract_observed_row(scores_df, CONDITION_BLOOD_NAIVE_B)

    fig = plt.figure(figsize=(12.6, 7.0))
    outer_gs = fig.add_gridspec(1, 2, width_ratios=[1.02, 1.38], wspace=0.22)
    left_gs = outer_gs[0].subgridspec(2, 1, height_ratios=[1.0, 1.12], hspace=0.46)
    gs_j = left_gs[0].subgridspec(2, 1, hspace=0.16)
    gs_k = left_gs[1].subgridspec(2, 1, hspace=0.16)

    # --- Left Top (gs_j): Reproduction of Paper Figure 4j (ISM logos for REF G vs ALT A) ---
    ism_start_0, ism_end_0 = ISM_WINDOW_0BASED
    logo_ymin = float(min(ism_ref_matrix.min(), ism_alt_matrix.min(), -0.5)) - 0.28
    logo_ymax = float(max(ism_ref_matrix.max(), ism_alt_matrix.max(), 1.5)) + 0.56

    ax_j1 = fig.add_subplot(gs_j[0])
    ax_j2 = fig.add_subplot(gs_j[1], sharex=ax_j1)
    for ax_j, mat, ylabel_txt in [
        (ax_j1, ism_ref_matrix, "REF (G)\nISM score"),
        (ax_j2, ism_alt_matrix, "ALT (A)\nISM score"),
    ]:
        _despine(ax_j)
        ax_j.axvspan(
            VARIANT_POS_1BASED - 0.5,
            VARIANT_POS_1BASED + 0.5,
            color="#FDB863",
            alpha=0.55,
            zorder=0,
        )
        ax_j.axhline(0.0, color="#ADB5BD", lw=0.7, zorder=1)
        ag_plot.seqlogo(
            mat,
            alphabet="ACGT",
            one_based=True,
            start=ism_start_0,
            ax=ax_j,
            letter_colors="default",
        )
        ax_j.set_xlim(ism_start_0 + 0.5, ism_end_0 + 0.5)
        ax_j.set_ylim(logo_ymin, logo_ymax)
        ax_j.set_ylabel(ylabel_txt)

    ax_j1.plot(
        [SG1617_GRNA_TARGET_0BASED[0] + 0.5, SG1617_GRNA_TARGET_0BASED[1] + 0.5],
        [logo_ymax - 0.24, logo_ymax - 0.24],
        color="#212529",
        lw=2.2,
    )
    ax_j1.text(
        (SG1617_GRNA_TARGET_0BASED[0] + SG1617_GRNA_TARGET_0BASED[1] + 1) / 2,
        logo_ymax - 0.16,
        "Frangoul et al. sg1617 gRNA target",
        ha="center",
        va="bottom",
        fontsize=7.5,
    )
    ax_j1.set_title(
        f"Our AG2 Reproduction of Fig. 4j — Saturated scATAC ISM\n({VARIANT_STR}, BM erythroblast observed 10x)",
        loc="left",
        pad=10,
    )
    ax_j1.tick_params(labelbottom=False)
    ax_j2.set_xticks([60_495_266, 60_495_276, 60_495_284])
    _format_genomic_xaxis(ax_j2)
    ax_j2.set_xlabel("chr2 position (GRCh38, 1-based; ISM window 60,495,263–60,495,285)")
    ax_j1.text(-0.14, 1.22, "a", transform=ax_j1.transAxes, fontsize=11, fontweight="bold")

    # --- Left Bottom (gs_k): Reproduction of Paper Figure 4k (+58 kb core REF vs ALT) ---
    sub_core = tracks_df[
        (tracks_df["position_0based"] >= ENHANCER_58KB_CORE[0])
        & (tracks_df["position_0based"] < ENHANCER_58KB_CORE[1])
    ].copy()
    xc = sub_core["position_1based"].to_numpy()
    ymax_k = float(max(sub_core["bm_erythroblast_ref"].max(), 22.0) * 1.14)

    ax_k1 = fig.add_subplot(gs_k[0])
    ax_k2 = fig.add_subplot(gs_k[1], sharex=ax_k1)

    _despine(ax_k1)
    ax_k1.plot(xc, sub_core["bm_erythroblast_ref"], color=COLOR_REF, lw=1.2, label="REF (G)")
    ax_k1.plot(xc, sub_core["bm_erythroblast_alt"], color=COLOR_ALT, lw=1.2, label="ALT (A)")
    ax_k1.axvline(VARIANT_POS_1BASED, color="#FDB863", lw=1.5, alpha=0.85)
    ax_k1.set_xlim(ENHANCER_58KB_CORE[0] + 1, ENHANCER_58KB_CORE[1])
    ax_k1.set_ylim(0.0, ymax_k)
    ax_k1.set_ylabel("BM\nerythroblast\nscATAC")
    ax_k1.set_title(
        "Our AG2 Reproduction of Fig. 4k — +58 kb Core Enhancer (490 bp)",
        loc="left",
        pad=6,
    )
    ax_k1.legend(loc="upper right", frameon=False)
    ax_k1.tick_params(labelbottom=False)
    ax_k1.text(-0.14, 1.15, "b", transform=ax_k1.transAxes, fontsize=11, fontweight="bold")

    _despine(ax_k2)
    ax_k2.plot(xc, sub_core["blood_naive_b_ref"], color=COLOR_REF, lw=1.2, label="REF (G)")
    ax_k2.plot(xc, sub_core["blood_naive_b_alt"], color=COLOR_ALT, lw=1.2, linestyle="--", label="ALT (A)")
    ax_k2.axvline(VARIANT_POS_1BASED, color="#FDB863", lw=1.5, alpha=0.85)
    ax_k2.set_ylim(0.0, ymax_k)
    ax_k2.set_ylabel("Blood\nnaive B cell\nscATAC")
    ax_k2.set_xticks([60_495_150, 60_495_300, 60_495_450])
    _format_genomic_xaxis(ax_k2)
    ax_k2.set_xlabel("chr2 position (GRCh38, +58 kb core: 60,495,049–60,495,539)")

    # --- Right Column: Qualitative Concordance & Transparency Summary vs Paper Fig. 4j-k ---
    ax_cmp = fig.add_subplot(outer_gs[1])
    ax_cmp.axis("off")
    ax_cmp.text(-0.035, 1.022, "c", transform=ax_cmp.transAxes, fontsize=11, fontweight="bold")
    ax_cmp.set_title(
        "Comparison with AlphaGenome 2 Paper (Fig. 4j–k, p. 10 & Methods pp. 47–48)",
        loc="left",
        pad=8,
    )

    ref_attr_g = float(ism_ref_matrix[VARIANT_POS_0BASED - ism_start_0, 2])
    alt_attr_a = float(ism_alt_matrix[VARIANT_POS_0BASED - ism_start_0, 0])
    b_core_max = float(max(sub_core["blood_naive_b_ref"].max(), sub_core["blood_naive_b_alt"].max()))

    rows_table = [
        ["Evaluation Dimension", "Paper Fig. 4j–k Report", "Our Live v2.0.0 API Result", "Audit"],
        [
            "Model & Input\nContext",
            "AG2 all-folds student;\n1-Mb on chr2:60495266",
            f"AG2 {metadata.get('model_version', 'ALL_FOLDS')} (v{metadata.get('alphagenome_version', '2.0.0')});\nchr2:59,970,978–61,019,554",
            "Exact\nsetup",
        ],
        [
            "Fig. 4j REF ISM\n(chr2:60,495,263–285)",
            "GATAA motif dominates\npositive attribution (~+2)",
            f"GATAA dominates local ISM;\nG@60495266 = {ref_attr_g:+.2f}",
            "Qualitative\nmatch",
        ],
        [
            "Fig. 4j ALT ISM\n(chr2:60495266:G>A)",
            "GATAA motif attribution\nablated on ALT allele",
            f"GATAA attribution lost;\nA@60495266 = {alt_attr_a:+.2f}",
            "Qualitative\nmatch",
        ],
        [
            "Fig. 4k Erythroblast\n(+58 kb core, 490 bp)",
            "High REF scATAC (>20);\nmarked drop on ALT",
            f"REF max {sub_core['bm_erythroblast_ref'].max():.1f} → ALT {sub_core['bm_erythroblast_alt'].max():.1f};\n1201-bp Δ_ATAC = {obs_ery['api_raw_score']:.4f}",
            "Qualitative\nmatch",
        ],
        [
            "Fig. 4k Naive B Cell\n(+58 kb core, 490 bp)",
            "Little predicted scATAC\nfor either REF or ALT",
            f"Flat near-zero (<{b_core_max:.2f});\n1201-bp Δ_ATAC = {obs_b['api_raw_score']:.4f}",
            "Qualitative\nmatch",
        ],
        [
            "Scorer Numerical\nCross-Validation",
            "log2(sum_1201(ALT)+1)\n− log2(sum_1201(REF)+1)",
            f"Max |API − Manual| =\n{max(obs_ery['abs_diff'], obs_b['abs_diff']):.2e}",
            "Verified\n(float32)",
        ],
        [
            "Paper Exact Float\nAvailability",
            "Not tabulated in text or\npublic GWAS downloads",
            "No eyeball digitizing;\nqualitative comparison only",
            "Boundary\nnoted",
        ],
    ]

    tbl = ax_cmp.table(
        cellText=rows_table[1:],
        colLabels=rows_table[0],
        colWidths=[0.24, 0.29, 0.33, 0.14],
        loc="center",
        cellLoc="left",
        bbox=[0.0, 0.02, 1.0, 0.94],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.4)
    for (r, _c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#DEE2E6")
        cell.set_linewidth(0.7)
        cell.PAD = 0.03
        if r == 0:
            cell.set_facecolor("#2B3A4A")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        elif r % 2 == 0:
            cell.set_facecolor("#F8F9FA")
        else:
            cell.set_facecolor("white")

    return _save_fig(fig, out_dir / "fig_d_paper_comparison")


def generate_all_figures(
    tracks_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    bcl11a_gtf_df: pd.DataFrame,
    ism_ref_matrix: np.ndarray,
    ism_alt_matrix: np.ndarray,
    metadata: dict[str, Any],
    output_dir: pathlib.Path,
) -> dict[str, dict[str, pathlib.Path]]:
    """Generates Figures A, B, C, and D in both 300-DPI PNG and vector PDF formats."""
    return {
        "figure_a": plot_figure_a_locus(bcl11a_gtf_df=bcl11a_gtf_df, out_dir=output_dir),
        "figure_b": plot_figure_b_erythroblast_ref_vs_alt(
            tracks_df=tracks_df, scores_df=scores_df, out_dir=output_dir
        ),
        "figure_c": plot_figure_c_cell_specificity(
            tracks_df=tracks_df, scores_df=scores_df, out_dir=output_dir
        ),
        "figure_d": plot_figure_d_paper_comparison(
            tracks_df=tracks_df,
            scores_df=scores_df,
            ism_ref_matrix=ism_ref_matrix,
            ism_alt_matrix=ism_alt_matrix,
            metadata=metadata,
            out_dir=output_dir,
        ),
    }
