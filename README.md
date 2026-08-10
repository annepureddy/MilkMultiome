# Distinct chromatin accessibility programs are linked to metabolic regulation in human lactocytes

Code repository for the manuscript: **"Distinct chromatin accessibility programs are linked to metabolic regulation in human lactocytes"**

## Overview

Matched snRNA-seq and snATAC-seq data were generated from somatic cells isolated from breast milk of 11 donors across a postpartum window of 11.6--97.3 weeks using the 10x Genomics Chromium Next GEM Single Cell Multiome ATAC + Gene Expression kit. Libraries were sequenced on the Illumina NextSeq 2000.

## Repository structure

```
MilkMultiome_2026/
├── README.md
├── 01_multiome_processing.R
├── session_info.txt
├── figures/
│   ├── Figure1.R
│   ├── Figure2.R
│   ├── Figure3.R
│   ├── Figure4.R
│   └── figure5/
│       ├── Figure5A_COMPASS.ipynb
│       ├── Figure5B_SCENITH.R
│       ├── reactions.tsv
│       ├── Pseudo_df_metadata.tsv
│       └── reaction_metadata.csv
```

### Processing

| File | Description |
|------|-------------|
| `01_multiome_processing.R` | End-to-end preprocessing: CellBender ambient RNA removal, QC filtering, doublet detection (scDblFinder, dbr = 0.08), SCTransform normalization (regressing percent mitochondrial reads), Harmony batch correction, RNA and ATAC clustering (resolution 0.25, SLM algorithm), epithelial subclustering, WNN integration (RNA dims 1:12 + ATAC dims 2:10), MACS2 peak calling, LinkPeaks peak-to-gene correlation, DEG (Wilcoxon) and DAP (LR test) identification |
| `session_info.txt` | Consolidated R session information with all package versions |

### Figure scripts

| File | Description |
|------|-------------|
| `figures/Figure1.R` | Donor metadata lollipop plot; UMAP visualizations (RNA + ATAC); DEG dot plot; ATAC-RNA proportion heatmap; Pearson correlation heatmaps; LC1/LC2 module score violin plot; GO enrichment barplots |
| `figures/Figure2.R` | Lactocyte subcluster characterization and gene expression analysis |
| `figures/Figure3.R` | Chromatin accessibility analysis: DAP stacked barplots; genomic feature annotation; butterfly lollipop plots (Epithelial 0 vs 4); EHD4 and ELF5 coverage and expression plots |
| `figures/Figure4.R` | TF motif enrichment and gene regulatory analysis: TF expression heatmap (z-scored); ChromVAR, RNA, and footprinting plots for FOSL2 and NR3C2; GO enrichment heatmaps per TF per cluster; ACSL1, CD36, and ACSS1 coverage, expression, and motif plots |
| `figures/figure5/Figure5A_COMPASS.ipynb` | COMPASS metabolic flux analysis (Python): PCA of reaction consistency scores; Cohen's d dot-plot by RECON2 subsystem; volcano plots for citric acid cycle and fatty acid oxidation pathways |
| `figures/figure5/Figure5B_SCENITH.R` | SCENITH/CENITH metabolic profiling: CD36+ metabolic parameter heatmap; effect size dot plot; raw MFI barplots for CD36+ and CD36- populations |

## Data requirements

To reproduce the figures, the following files are needed:

| File | Source |
|------|--------|
| `Epithelial_subset.rds` | Processed Seurat object (GEO) |
| `TableS5_Links_whole_dataset.csv` | Peak-to-gene links (Supplemental Table) |
| `TableS4_DAPs_Epi_LR.csv` | Differentially accessible peaks (Supplemental Table) |
| `TableS2.csv` | Differentially expressed genes (Supplemental Table) |
| `TableS11_SCENITH.xlsx` | SCENITH MFI data (Supplemental Table S11) |
| `TableS9_reactions.tsv` | COMPASS penalty matrix (included in `figure5/`) |
| `Pseudo_df_metadata.tsv` | Pseudobulk metadata (included in `figure5/`) |
| `reaction_metadata.csv` | RECON2 annotations (included in `figure5/`) |

Place the Seurat object and CSV files in the same directory as the figure scripts. COMPASS input files are included in the `figure5/` subdirectory.

## Software and dependencies

All R analyses were performed in R 4.5.2. COMPASS analysis was performed in Python 3.

| Package | Version | Purpose |
|---------|---------|---------|
| Seurat | 5.4.0 | snRNA-seq processing and visualization |
| Signac | 1.16.0 | snATAC-seq processing and visualization |
| Harmony | 1.2.4 | Batch correction |
| SCTransform | 0.4.1 | Normalization |
| clusterProfiler | 4.8.14 | GO enrichment analysis |
| ComplexHeatmap | 2.22.0 | Heatmap visualization |
| chromVAR | 1.28.0 | TF motif activity scores |
| JASPAR2024 | -- | TF motif database |
| ChIPseeker | -- | Peak annotation |
| COMPASS | -- | Metabolic flux analysis (RECON2) |

Full session information is provided in `session_info.txt`.

### Additional tools

- **Alignment:** Cell Ranger ARC v2.0.2 (GRCh38 reference)
- **Ambient RNA removal:** CellBender v0.3.2
- **Peak calling:** MACS2 v2.2.9.1
- **Doublet detection:** scDblFinder (doublet rate = 0.08)

## Data availability

Raw and processed sequencing data are deposited in NCBI Gene Expression Omnibus (GEO) under accession GSE342361. Supplemental tables referenced in the scripts are available as Supplementary Data in the published manuscript.

## Contact

For questions regarding the code: Laasya Devi Annepureddy (laasya.devi.annepureddy.th@dartmouth.edu)
For questions regarding the study: Brittany A. Goods (britt.goods@unimelb.edu.au)
