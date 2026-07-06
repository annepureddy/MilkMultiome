# Distinct chromatin accessibility programs are linked to metabolic regulation in human lactocytes

**Laasya Devi Annepureddy<sup>1</sup>, Onyeka Ozuzu<sup>1</sup>, Olha Kholod<sup>1</sup>, Olivia Palmer<sup>1</sup>, Madison Spivak<sup>1</sup>, Ibrahim Ozgenc<sup>1,2</sup>, Irma Vlasac<sup>3</sup>, Michael C Rudolph<sup>4</sup>, Carol Ringeberg<sup>5</sup>, Elizabeth Sergison<sup>6</sup>, Fred Kolling<sup>5</sup>, Yu-Ying Chen<sup>7\*</sup>, Brock Christensen<sup>3,8+\*</sup>, Brittany A. Goods<sup>1,9+\*</sup>**

*In review*

---

## Abstract

The mammary epithelium undergoes extensive remodeling during pregnancy to support lactation, yet how specialized milk-producing epithelial (lactocyte) states are established in the human mammary gland remains poorly defined. Here, we use single-nucleus multiome profiling of human milk-derived cells to map the paired transcriptional and regulatory landscapes of lactocytes. We identify five lactocyte subclusters with distinct transcriptomic and epigenomic profiles. Secretory lactocytes exhibit increased chromatin accessibility at loci encoding metabolic enzymes involved in fatty acid handling and lipid secretion. Motif analysis implicates several transcription factors, including XBP1, NFIA, and NR3C2/1, in accessible regions near genes involved in lipid import and acyl-CoA biosynthesis, including CD36, ACSL1, and ACSS1. Functional metabolic profiling of these cells suggests that secretory epithelial cells, but not non-secretory cells, as marked by CD36 expression, are capable of engaging both glycolysis and fatty acid oxidation. These findings link lactocyte subtype identity to chromatin organization and metabolic phenotype, providing a framework for understanding epithelial regulation and functional specialization during human lactation.

---

## Repository structure

```
MilkMultiome_2026/
├── README.md
├── figures/
│   ├── Figure1.R       # Donor metadata lollipop plot; UMAP visualizations (RNA + ATAC);
│   │                   # DEG dot plot; ATAC-RNA proportion heatmap; Pearson correlation
│   │                   # heatmaps; LC1/LC2 module score violin plot; GO enrichment barplots
│   ├── Figure2.R       # [description — to be updated]
│   ├── Figure3.R       # [description — to be updated]
│   └── Figure4.R       # [description — to be updated]
└── session_info.txt    # R and package versions for all analyses
```

---

## Data availability

Raw and processed data are deposited in NCBI Gene Expression Omnibus (GEO) under accession **GSE[XXX]**.

Supplemental tables referenced in the scripts (including DEG tables and gene signature files) are available as Supplementary Data in the published manuscript.

> **To reproduce figures:** download the processed Seurat object (`Epithelial_final_links_footprint.rds`) from GEO, set `data_dir` at the top of each script to point to your local download location, and run the scripts in order.

---

## Software and dependencies

All analyses were performed in **R** (version — see `session_info.txt`).

Key packages:

| Package | Purpose |
|---|---|
| Seurat / Signac | Single-nucleus RNA-seq and ATAC-seq processing and visualization |
| clusterProfiler | GO enrichment analysis |
| rrvgo | GO term semantic similarity and simplification |
| pheatmap | Heatmap visualization |
| ggplot2 / patchwork | Figure assembly |

Full session information including all package versions is provided in `session_info.txt`.

---

## Contact

For questions regarding the code, please contact **Laasya Devi Annepureddy** (laasya.devi.annepureddy.th@dartmouth.edu).

For questions regarding the study, please contact **Brittany A. Goods** (britt.goods@unimelb.edu.au)
