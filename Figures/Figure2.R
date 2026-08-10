# ============================================================
# Figure 2: 
#
# Annepureddy et al., MilkMultiome manuscript (2026)
#
# Panels reproduced here:
#   A     Paired heatmaps of z-scored RNA expression (top) and promoter
#         accessibility (bottom) for lactose synthesis, MFG synthesis, and
#         HMO biosynthesis gene modules across epithelial subclusters (Epi 0–4).
#         Promoter accessibility = mean ATAC-seq signal within ±2 kb of TSS.
#         Columns ordered LC1-like (Epi 0, Epi 3) → intermediate (Epi 1)
#         → LC2-like (Epi 2, Epi 4).
#   B     Genome browser–style coverage plots (±5 kb from gene body) at:
#         MFG synthesis loci (left): BTN1A1, CIDEA, FASN
#         Lactose synthesis loci (right): UGP2, LALBA, CSN2
#         RNA expression overlaid as a separate track below each
#         accessibility profile. Track colors = epithelial subcluster identity.
# ============================================================

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(patchwork)
library(dplyr)

# ------------------------------------------------------------
# Paths
# Set data_dir to the location of your downloaded GEO files.
# All figure outputs are written to output_dir.
# ------------------------------------------------------------
data_dir   <- "data"
output_dir <- "output/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
# Epithelial subset: paired snRNA-seq / snATAC-seq, re-clustered
Epi_subset <- readRDS(file.path(data_dir, "Epithelial_subset.rds"))

# ------------------------------------------------------------
# Color palette (Epithelial subclusters 0–4)
# Columns ordered LC1-like → intermediate → LC2-like throughout
# ------------------------------------------------------------
my_colors <- c(
  "Epithelial 0" = "#E24A33",
  "Epithelial 1" = "#348ABD",
  "Epithelial 2" = "#988ED5",
  "Epithelial 3" = "#FBC15E",
  "Epithelial 4" = "#8EBA42"
)

cluster_order <- c("Epithelial 0", "Epithelial 3", "Epithelial 1",
                   "Epithelial 2", "Epithelial 4")

## ============================================================
## Panel A: Paired RNA / ATAC heatmaps for lactation gene modules
## ============================================================

# --- Define lactation gene modules ---
# Genes drawn from Nyquist & Annepureddy et al.
lactation_genes <- list(
  lactose_synthesis = c("PGM2", "UGP2", "GALE", "SLC35A2", "LALBA",
                        "B4GALT1", "CSN2"),
  mfg_genes         = c("FABP3", "BTN1A1", "XDH", "CIDEA", "PLIN2", "FASN"),
  hmo_biosynthesis  = c("ST3GAL6", "B4GALT6", "FUT9", "FUT4", "B3GALNT1",
                        "B3GNT4", "B3GNT2", "FUT10", "B3GAT2", "ST3GAL5",
                        "ST3GAL3", "FUT8", "ST3GAL1", "FUT11", "FUT6",
                        "ST6GALNAC2", "B4GALT3", "B3GNT3", "B3GNT7",
                        "B3GALT5", "ST6GALNAC4", "B4GALT5", "GCNT3", "FUT3",
                        "B3GALT6", "B3GALT4", "B3GALNT2", "B3GAT3",
                        "B4GALT1", "B3GNT9", "ST6GALNAC6", "B4GALT7",
                        "B4GALT4", "FUT2", "ST6GAL1", "B4GALT2", "GCNT2",
                        "GCNT1")
)

all_genes <- unique(unlist(lactation_genes, use.names = FALSE))

# Gene-to-module mapping for row annotation
gene_to_category <- data.frame(
  gene     = unlist(lactation_genes),
  category = rep(names(lactation_genes), lengths(lactation_genes)),
  stringsAsFactors = FALSE
) |> dplyr::distinct(gene, .keep_all = TRUE)

# --- Retrieve TSS coordinates from EnsDb (GRCh38 / Ensembl v86) ---
gene_coords <- ensembldb::genes(
  EnsDb.Hsapiens.v86,
  filter  = GeneNameFilter(all_genes),
  columns = c("gene_name", "seq_name", "gene_seq_start",
              "gene_seq_end", "strand")
)
gene_coords <- keepStandardChromosomes(gene_coords, pruning.mode = "coarse")

# Promoter regions: ±2 kb from TSS
promoter_regions <- promoters(gene_coords, upstream = 2000, downstream = 2000)
names(promoter_regions) <- promoter_regions$gene_name

found_genes   <- unique(promoter_regions$gene_name)
missing_genes <- setdiff(all_genes, found_genes)
if (length(missing_genes) > 0) {
  message("Genes not found in EnsDb: ", paste(missing_genes, collapse = ", "))
}

# --- Map promoter regions to ATAC peaks ---
DefaultAssay(Epi_subset) <- "peaksbyc"
peak_granges <- granges(Epi_subset)

seqlevelsStyle(promoter_regions) <- "UCSC"   # match UCSC-style peak coords
overlaps <- findOverlaps(promoter_regions, peak_granges)

overlap_df <- data.frame(
  gene = promoter_regions$gene_name[queryHits(overlaps)],
  peak = GRangesToString(peak_granges[subjectHits(overlaps)]),
  stringsAsFactors = FALSE
)

message(length(unique(overlap_df$gene)), " genes have promoter-proximal peaks")
message(length(unique(overlap_df$peak)), " unique peaks mapped")

# --- Average RNA expression per cluster ---
DefaultAssay(Epi_subset) <- "RNA"
Idents(Epi_subset) <- "Celltype"

rna_genes_present <- intersect(found_genes, rownames(Epi_subset[["RNA"]]))
rna_avg <- AverageExpression(
  Epi_subset,
  assays   = "RNA",
  features = rna_genes_present,
  slot     = "data"
)$RNA

# --- Average ATAC promoter accessibility per gene per cluster ---
# Multiple peaks overlapping a gene's promoter are mean-averaged
DefaultAssay(Epi_subset) <- "peaksbyc"

atac_peaks_present <- intersect(unique(overlap_df$peak),
                                rownames(Epi_subset[["peaksbyc"]]))
atac_avg <- AverageExpression(
  Epi_subset,
  assays   = "peaksbyc",
  features = atac_peaks_present,
  slot     = "data"
)$peaksbyc

genes_with_atac <- unique(overlap_df$gene[overlap_df$peak %in% atac_peaks_present])

atac_gene_avg <- do.call(rbind, lapply(genes_with_atac, function(g) {
  peaks <- intersect(overlap_df$peak[overlap_df$gene == g], rownames(atac_avg))
  if (length(peaks) == 1) {
    return(atac_avg[peaks, , drop = FALSE])
  } else if (length(peaks) > 1) {
    return(matrix(colMeans(atac_avg[peaks, , drop = FALSE]),
                  nrow = 1, dimnames = list(g, colnames(atac_avg))))
  }
}))
rownames(atac_gene_avg) <- genes_with_atac

# --- Align genes present in both modalities ---
shared_genes <- intersect(rownames(rna_avg), rownames(atac_gene_avg))
message(length(shared_genes), " genes with both RNA and promoter-proximal ATAC data")

# --- Module annotation and display labels ---
keep_categories <- c("hmo_biosynthesis", "lactose_synthesis", "mfg_genes")

category_labels <- c(
  hmo_biosynthesis  = "HMO Biosynthesis",
  lactose_synthesis = "Lactose Synthesis",
  mfg_genes         = "MFG Synthesis"
)

category_colors <- c(
  hmo_biosynthesis  = "#66C2A5",
  lactose_synthesis = "#FC8D62",
  mfg_genes         = "#8DA0CB"
)

gene_order <- gene_to_category |>
  dplyr::filter(gene %in% shared_genes, category %in% keep_categories) |>
  dplyr::arrange(factor(category, levels = keep_categories))

shared_genes_ordered <- gene_order$gene

row_split <- factor(
  gene_order$category,
  levels = keep_categories,
  labels = c("HMO\nBiosynthesis", "Lactose\nSynthesis", "MFG\nSynthesis")
)

# --- Subset matrices and z-score across clusters (row-wise) ---
rna_mat    <- as.matrix(rna_avg[shared_genes_ordered, cluster_order])
atac_mat   <- as.matrix(atac_gene_avg[shared_genes_ordered, cluster_order])
rna_scaled  <- t(scale(t(rna_mat)))
atac_scaled <- t(scale(t(atac_mat)))

# --- Row annotation (gene module) ---
row_anno <- rowAnnotation(
  Category = gene_order$category,
  col = list(Category = category_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Category = list(
      title  = "Category",
      at     = names(category_labels),
      labels = unname(category_labels)
    )
  )
)

# --- Color scales (viridis for both RNA and ATAC) ---
rna_col  <- colorRamp2(c(-2, 0, 2), viridis(3))
atac_col <- colorRamp2(c(-2, 0, 2), viridis(3))

# --- Heatmap: RNA expression ---
ht_rna <- Heatmap(
  rna_scaled,
  name              = "RNA\n(z-score)",
  col               = rna_col,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_title      = "Gene Expression (RNA)",
  row_names_side    = "left",
  row_names_gp      = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp   = gpar(fontsize = 10),
  width             = unit(5, "cm"),
  row_split         = row_split,
  row_title_rot     = 0,
  row_title_gp      = gpar(fontsize = 9, fontface = "bold"),
  row_gap           = unit(3, "mm")
)

# --- Heatmap: Promoter accessibility ---
ht_atac <- Heatmap(
  atac_scaled,
  name              = "ATAC\n(z-score)",
  col               = atac_col,
  cluster_rows      = FALSE,   # row order matches RNA heatmap
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_title      = "Promoter Accessibility (ATAC)",
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp   = gpar(fontsize = 10),
  width             = unit(5, "cm"),
  row_gap           = unit(3, "mm")
)

# --- Draw and save ---
pdf(file.path(output_dir, "Fig2A_lactation_RNA_ATAC_heatmap.pdf"),
    width  = 13,
    height = max(6, length(shared_genes_ordered) * 0.25 + 2))
draw(
  row_anno + ht_rna + ht_atac,
  column_title    = paste("Lactation Gene Expression and Promoter Accessibility",
                          "Across Epithelial Subclusters"),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  merge_legend    = TRUE,
  heatmap_legend_side = "right"
)
dev.off()

## ============================================================
## Panel B: Coverage plots (±5 kb from gene body)
## MFG synthesis (left): BTN1A1, CIDEA, FASN
## Lactose synthesis (right): UGP2, LALBA, CSN2
## RNA expression overlaid as track below each accessibility profile
## ============================================================
DefaultAssay(Epi_subset) <- "peaksbyc"
Idents(Epi_subset) <- "Celltype"

# Set consistent cluster display order
levels(Epi_subset) <- cluster_order

mfg_genes_cov    <- c("BTN1A1", "CIDEA", "FASN")
lactose_genes_cov <- c("UGP2", "LALBA", "CSN2")

make_cov_plot <- function(gene, obj, cols) {
  CoveragePlot(
    object             = obj,
    region             = gene,
    features           = gene,
    assay              = "peaksbyc",
    expression.assay   = "RNA",
    show.bulk          = TRUE,
    peaks              = TRUE,
    links              = FALSE,
    extend.upstream    = 5000,
    extend.downstream  = 5000
  ) & scale_fill_manual(values = cols)
}

mfg_plots    <- lapply(mfg_genes_cov,    make_cov_plot, obj = Epi_subset, cols = my_colors)
lactose_plots <- lapply(lactose_genes_cov, make_cov_plot, obj = Epi_subset, cols = my_colors)

# Individual pathway PDFs
ggsave(file.path(output_dir, "Fig2B_coverage_mfg_synthesis.pdf"),
       wrap_plots(mfg_plots, ncol = 1) + plot_annotation(title = "MFG Synthesis"),
       width = 8, height = 12)

ggsave(file.path(output_dir, "Fig2B_coverage_lactose_synthesis.pdf"),
       wrap_plots(lactose_plots, ncol = 1) + plot_annotation(title = "Lactose Synthesis"),
       width = 8, height = 12)

# Combined figure: MFG (left column) | Lactose (right column), 3 genes per column
p_combined <- (mfg_plots[[1]] | lactose_plots[[1]]) /
              (mfg_plots[[2]] | lactose_plots[[2]]) /
              (mfg_plots[[3]] | lactose_plots[[3]])

ggsave(file.path(output_dir, "Fig2B_coverage_combined.pdf"),
       p_combined, width = 14, height = 15)

## ============================================================
## Session info
## ============================================================
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
