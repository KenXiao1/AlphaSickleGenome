"""Canonical genomic coordinates, conditions, and scoring parameters from AlphaGenome 2 (2026)."""

from __future__ import annotations

# Reference genome & target gene (GENCODE v46 / GRCh38; Paper Methods p. 47)
REFERENCE_ASSEMBLY = "GRCh38"
UCSC_GENOME_ID = "hg38"
BCL11A_GENE_SYMBOL = "BCL11A"
BCL11A_GENE_ID = "ENSG00000119866"
BCL11A_MANE_TRANSCRIPT_ID = "ENST00000642384"
BCL11A_STRAND = "-"
BCL11A_MANE_TSS_1BASED_OR_0POINT = 60_553_654

# Variant specification (1-based genomic coordinate on positive reference strand; Paper Fig. 4j-k)
VARIANT_CHROM = "chr2"
VARIANT_POS_1BASED = 60_495_266
VARIANT_POS_0BASED = VARIANT_POS_1BASED - 1  # 60_495_265
VARIANT_REF = "G"
VARIANT_ALT = "A"
VARIANT_STR = f"{VARIANT_CHROM}:{VARIANT_POS_1BASED}:{VARIANT_REF}>{VARIANT_ALT}"

# Model input & scoring window sizes (Paper Methods pp. 20, 35, 47-48)
INPUT_SEQUENCE_LENGTH_1MB = 2**20  # 1,048,576 bp
SC_ATAC_CENTER_MASK_WIDTH = 1_201
SC_ATAC_PSEUDOCOUNT = 1.0

# Enhancer & guide RNA intervals in 0-based half-open [start, end) coordinates
# Paper Methods p. 47 ("Locus analysis at BCL11A"):
ENHANCER_62KB_PAPER = (60_490_400, 60_491_800)
ENHANCER_58KB_CORE = (60_495_049, 60_495_539)  # Intersection of Bauer et al. & Javierre et al. (Fig. 4k window)
ENHANCER_55KB_PAPER = (60_497_800, 60_499_200)

# Official tutorial (`enhancer_to_gene_bcl11a.ipynb`) broader DHS coordinates for comparison/documentation:
ENHANCER_62KB_DHS_TUTORIAL = (60_490_409, 60_491_729)
ENHANCER_58KB_DHS_TUTORIAL = (60_494_920, 60_496_208)
ENHANCER_55KB_DHS_TUTORIAL = (60_497_676, 60_498_941)

# Therapeutic sg1617 gRNA target & GATAA motif (Frangoul et al. 2021; Wu et al. 2019; Paper Methods p. 48)
SG1617_GRNA_TARGET_0BASED = (60_495_263, 60_495_284)
GATAA_MOTIF_0BASED = (60_495_265, 60_495_270)  # 1-based inclusive 60,495,266-60,495,270 -> 'GATAA'
ISM_WINDOW_0BASED = (60_495_262, 60_495_285)  # 23-bp saturated ISM window in Paper Fig. 4j (Methods p. 48)

# Saved slice window covering all three erythroid enhancers (+62/+58/+55 kb) at 1-bp resolution
SAVED_LOCUS_SLICE_0BASED = (60_490_000, 60_500_000)

# Verified official AlphaGenome 2 reference condition CURIEs (from `dna_model.supported_conditions()`)
CONDITION_BM_ERYTHROBLAST = "UBERON:0002371_CL:0000765"  # bone marrow erythroblast (Condition A, main)
CONDITION_BLOOD_NAIVE_B = "UBERON:0000178_CL:0000788"  # blood naive B cell (Condition B, Paper Fig. 4k control)
CONDITION_BM_NAIVE_B_TUTORIAL = "UBERON:0002371_CL:0000788"  # bone marrow naive B cell (used in E2G tutorial)

# Observed single-cell ATAC channel shared across both conditions
OBSERVED_SC_ATAC_TRACK_SUFFIX = "scATAC-seq 10X-cellranger-nuclear"

# Upstream reference URLs for reproducible verification
UCSC_HG38_SEQUENCE_URL = (
    "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr2;start=60495000;end=60495600"
)
ALPHAGENOME_GENCODE_V46_FEATHER_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather"
)
