# ============================================================
# Figure 1: 
# Annepureddy et al., MilkMultiome manuscript (2026)
#
# Panels reproduced here:
#   A(ii) Lollipop plot of weeks postpartum at sample collection (n = 11),
#         dot color = donor age, dashed line = cohort median
#   B(i)  UMAP: lactocyte subclusters, snRNA-seq
#   B(ii) UMAP: snATAC-seq, RNA-defined cell type labels projected onto
#         ATAC space
#   B(iii)UMAP: independent unsupervised clustering of snATAC-seq data
#   C     Dot plot of top 5 DEGs per lactocyte subcluster (ranked by padj)
#   D     Heatmap: proportional correspondence of each ATAC cluster by
#         RNA-defined lactocyte subtype (columns sum to 100%)
#   E     Pearson correlation heatmaps (RNA and ATAC), top 2000 HVGs/peaks,
#         hierarchical clustering
#   F     Violin plot: LC1 vs. LC2 gene signature module score across
#         lactocyte subclusters, ordered by median
#   G     Barplot: top 5 significant GO Biological Process terms per
#         epithelial subcluster (top 500 DEGs), simplified with rrvgo
# ============================================================

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(forcats)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(writexl)
library(ggthemes)

# ------------------------------------------------------------
# Paths
# Set these once for your local/HPC environment. All inputs are
# expected in `data_dir`; all figures/tables are written to `output_dir`.
# ------------------------------------------------------------
data_dir   <- "data"
output_dir <- "output/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
# Epithelial subset: paired snRNA-seq / snATAC-seq, re-clustered
Epi_subset <- readRDS(file.path(data_dir, "Epithelial_subset.rds"))

# Table S2 - significant DEGs (epithelial subclusters)
degs_RNA <- read.csv(file.path(data_dir, "TableS2.csv"))
degs_RNA <- degs_RNA[degs_RNA$p_val_adj < 0.05, ]

# Luminal cell gene signatures (Nyquist et al.)
gene_sig <- read.csv(file.path(data_dir, "Nyquist_LC1_LC2_signatures.csv"))

# ------------------------------------------------------------
# Color palette (Epithelial subclusters 0-4)
# ------------------------------------------------------------
my_colors <- c(
  "Epithelial 0" = "#E24A33",
  "Epithelial 1" = "#348ABD",
  "Epithelial 2" = "#988ED5",
  "Epithelial 3" = "#FBC15E",
  "Epithelial 4" = "#8EBA42"
)

# atac_colors used for independent ATAC clusters (Panel B-iii, D)
atac_colors <- c("#A8D8EA", "#F6C3B7", "#B5DABB", "#C9B1D0", "#F9E4AA")

## ============================================================
## Panel A(ii): Lollipop plot - weeks postpartum by donor, colored by age
## ============================================================
donors <- data.frame(
  donor_id = c("FC003", "FC004", "FC006", "FC007", "FC008",
               "FC009", "FC011", "FC012", "FC013", "FC014", "FC018"),
  age      = c(33.93, 36.32, 24.45, 33.53, 26.43,
               31.74, 25.75, 34.60, 36.24, 27.83, 34.27),
  days_pp  = c(449, 151, 229, 257, 218, 81, 681, 265, 142, 227.5, 135)
) %>%
  mutate(
    weeks_pp = days_pp / 7,
    donor_id = reorder(donor_id, weeks_pp)
  )

p_A_ii <- ggplot(donors, aes(x = weeks_pp, y = donor_id)) +
  geom_segment(aes(x = 0, xend = weeks_pp, yend = donor_id),
               colour = "grey80", linewidth = 0.5) +
  geom_point(aes(colour = age), size = 4) +
  geom_vline(xintercept = median(donors$weeks_pp),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  annotate("text",
           x = median(donors$weeks_pp) + 1.5, y = 0.65,
           label = paste0("median ", round(median(donors$weeks_pp), 1), " wks"),
           hjust = 0, size = 2.6, colour = "grey35") +
  scale_colour_gradientn(
    colours = c("#cce5f6", "#4a9cc9", "#1a3a6b"),
    name = "Age (years)",
    breaks = c(25, 30, 35),
    guide = guide_colourbar(barwidth = 5, barheight = 0.55,
                             title.position = "top", title.hjust = 0.5)
  ) +
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
                      expand = c(0.01, 0)) +
  labs(x = "Weeks postpartum", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),
    axis.title.x       = element_text(size = 10),
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3)
  )

ggsave(file.path(output_dir, "Fig1Aii.pdf"),
       plot = p_A_ii, width = 4.5, height = 3.5, units = "in", device = cairo_pdf)

## ============================================================
## Panel B: UMAPs
##  (i)   Lactocyte subclusters, snRNA-seq
##  (ii)  snATAC-seq, RNA-defined labels projected onto ATAC space
##  (iii) Independent unsupervised clustering of snATAC-seq
## ============================================================
DefaultAssay(Epi_subset) <- "RNA"
Idents(Epi_subset) <- Epi_subset$Celltype

p_B_i   <- DimPlot(Epi_subset, reduction = "umap.rna.integrated",
                    cols = my_colors) + NoLegend()
p_B_ii  <- DimPlot(Epi_subset, reduction = "umap.atac.integrated",
                    cols = my_colors) + NoLegend()
p_B_iii <- DimPlot(Epi_subset, reduction = "umap.atac.integrated",
                    group.by = "peaksbyc_snn_res.0.2",
                    cols = atac_colors, label = TRUE) + NoLegend()

ggsave(file.path(output_dir, "Fig1B_umaps.pdf"),
       plot = p_B_i + p_B_ii + p_B_iii + plot_layout(nrow = 1) &
         theme_void() & NoAxes() & NoLegend(),
       width = 20, height = 5, dpi = 600)

## ============================================================
## Panel C: Dot plot, top 5 DEGs per lactocyte subcluster (ranked by padj)
## ============================================================
p_C <- DotPlot(Epi_subset, assay = "RNA", features = degs_RNA$gene) +
  scale_color_viridis_c(option = "cividis") +
  RotatedAxis()

ggsave(file.path(output_dir, "Fig1C.pdf"), plot = p_C,
       width = 10, height = 3, units = "in", dpi = 300)

## ============================================================
## Panel D: ATAC cluster x RNA subtype proportion heatmap
## Column values sum to 100% (fraction of cells in each ATAC cluster
## corresponding to each RNA-defined transcriptional subtype)
## ============================================================
props <- as.data.frame(table(
  RNA  = Epi_subset$Celltype,
  ATAC = Epi_subset$peaksbyc_snn_res.0.2
)) %>%
  group_by(ATAC) %>%
  mutate(Pct = Freq / sum(Freq) * 100) %>%
  ungroup()

mat <- props %>%
  dplyr::select(RNA, ATAC, Pct) %>%
  pivot_wider(names_from = ATAC, values_from = Pct, values_fill = 0) %>%
  column_to_rownames("RNA") %>%
  as.matrix()

rna_colors <- my_colors
names(rna_colors) <- rownames(mat)
names(atac_colors) <- colnames(mat)

annotation_row <- data.frame(RNA_Celltype = rownames(mat), row.names = rownames(mat))
ann_colors_D <- list(RNA_Celltype = rna_colors, ATAC_Cluster = atac_colors)

pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         color = colorRampPalette(c("white", "#f7c6c0", "#c0392b", "#8B0000"))(100),
         annotation_row = annotation_row,
         annotation_colors = ann_colors_D,
         angle_col = 0,
         fontsize = 10,
         fontsize_number = 8,
         main = "Cell Type Proportion (%)",
         filename = file.path(output_dir, "Fig1D_proportion_heatmap.pdf"),
         width = 8, height = 5)

## ============================================================
## Panel E: Pearson correlation heatmaps, RNA and ATAC
## Top 2000 highly variable genes/peaks, hierarchical clustering
## ============================================================
ann_colors_E <- list(Celltype = setNames(my_colors, names(my_colors)))

# --- RNA correlation ---
DefaultAssay(Epi_subset) <- "RNA"
Idents(Epi_subset) <- "Celltype"

avg_expr  <- AverageExpression(Epi_subset, assays = "RNA", slot = "scale.data")$RNA
top_genes <- head(VariableFeatures(Epi_subset), 2000)
top_genes <- top_genes[top_genes %in% rownames(avg_expr)]
cor_rna   <- cor(as.matrix(avg_expr[top_genes, ]), use = "pairwise.complete.obs",
                  method = "pearson")
rng_rna   <- range(cor_rna, na.rm = TRUE)

annotation_df_rna <- data.frame(Celltype = colnames(cor_rna), row.names = colnames(cor_rna))

pdf(file.path(output_dir, "Fig1E_rna_correlation_heatmap.pdf"), width = 6, height = 5)
pheatmap(
  cor_rna,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  border_color = NA,
  color = colorRampPalette(c("lightgrey", "#6A040F"))(100),
  breaks = seq(rng_rna[1], rng_rna[2], length.out = 101),
  annotation_row = annotation_df_rna,
  annotation_col = annotation_df_rna,
  annotation_colors = ann_colors_E,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "RNA Pearson Correlation (Top 2000 HVGs)"
)
dev.off()

# --- ATAC correlation ---
DefaultAssay(Epi_subset) <- "peaksbyc"
Idents(Epi_subset) <- "peaksbyc_snn_res.0.2"

avg_access <- AverageExpression(Epi_subset, assays = "peaksbyc", slot = "data")$peaksbyc
top_peaks  <- head(VariableFeatures(Epi_subset[["peaksbyc"]]), 2000)
top_peaks  <- top_peaks[top_peaks %in% rownames(avg_access)]
cor_atac   <- cor(as.matrix(avg_access[top_peaks, ]), use = "pairwise.complete.obs",
                   method = "pearson")
rng_atac   <- range(cor_atac, na.rm = TRUE)

pdf(file.path(output_dir, "Fig1E_atac_correlation_heatmap.pdf"), width = 6, height = 5)
pheatmap(
  cor_atac,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  border_color = NA,
  color = colorRampPalette(c("lightgrey", "#6A040F"))(100),
  breaks = seq(rng_atac[1], rng_atac[2], length.out = 101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "ATAC Pearson Correlation (Top 2000 Variable Peaks)"
)
dev.off()

## ============================================================
## Panel F: LC1 vs. LC2 module score violin plot, ordered by median
## Gene signatures: top 250 LC1 / LC2 genes from Nyquist et al.
## (positive = LC1 enrichment, negative = LC2 enrichment)
## ============================================================
DefaultAssay(Epi_subset) <- "SCT"

n_top <- 250
gene_modules <- list(
  Nyquist_LC1 = gene_sig$Nyquist_LC1[1:n_top],
  Nyquist_LC2 = gene_sig$Nyquist_LC2[1:n_top]
)

# Drop NAs, restrict to genes present in the dataset
gene_modules <- lapply(gene_modules, function(g) {
  g <- g[!is.na(g)]
  intersect(g, rownames(Epi_subset))
})

for (module_name in names(gene_modules)) {
  Epi_subset <- AddModuleScore(
    object   = Epi_subset,
    features = list(gene_modules[[module_name]]),
    name     = module_name,
    seed     = 42
  )
}

# AddModuleScore appends "1" to each name; rename for clarity
score_cols  <- paste0(names(gene_modules), "1")
names_clean <- names(gene_modules)
colnames(Epi_subset@meta.data)[match(score_cols, colnames(Epi_subset@meta.data))] <- names_clean

Epi_subset$Nyquist_LC1_vs_LC2 <- Epi_subset$Nyquist_LC1 - Epi_subset$Nyquist_LC2

# Order subclusters by descending median LC1-vs-LC2 score
cluster_order <- Epi_subset@meta.data %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(med = median(Nyquist_LC1_vs_LC2, na.rm = TRUE)) %>%
  dplyr::arrange(desc(med)) %>%
  dplyr::pull(Celltype)

Epi_subset$Celltype <- factor(Epi_subset$Celltype, levels = cluster_order)

p_F <- VlnPlot(Epi_subset, features = "Nyquist_LC1_vs_LC2", group.by = "Celltype",
               pt.size = 0, cols = my_colors) +
  geom_boxplot(width = 0.1, outlier.size = 0, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  ggtitle("Nyquist: LC1 vs LC2 Score") +
  NoLegend()

ggsave(file.path(output_dir, "Fig1F_Nyquist_LC1_vs_LC2_violin.pdf"),
       plot = p_F, width = 8, height = 6, dpi = 300)

## ============================================================
## Panel G: GO Biological Process enrichment, top 5 terms per
## epithelial subcluster (top 500 DEGs by pct_diff, padj < 0.05)
## Terms simplified with rrvgo (Sayols et al., PMID: 37151216)
## ============================================================
top500 <- degs_RNA %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(pct_diff), .by_group = TRUE) %>%
  slice_head(n = 500) %>%
  summarise(genes = list(gene), .groups = "drop")

run_enrichGO <- function(genes, cluster_name) {
  ego_df <- as.data.frame(enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    minGSSize     = 5,
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  ))
  ego_df$cluster <- cluster_name
  return(ego_df)
}

all_go <- top500 %>%
  mutate(enrichGO_result = map2(genes, cluster, run_enrichGO)) %>%
  dplyr::select(-genes) %>%
  unnest(enrichGO_result, names_sep = "_") %>%
  filter(enrichGO_result_p.adjust < 0.05, enrichGO_result_Count >= 5) %>%
  dplyr::select(-enrichGO_result_cluster) %>%
  rename_with(~ str_remove(., "^enrichGO_result_"))

# Semantic similarity matrix across all significant GO BP terms
sim_matrix <- calculateSimMatrix(
  all_go$ID,
  orgdb  = "org.Hs.eg.db",
  ont    = "BP",
  method = "Rel"
)

# Reduce redundant terms per cluster (rrvgo), keep top 5 by padj
reduced_go <- all_go %>%
  group_by(cluster) %>%
  group_modify(function(df, grp) {
    ids <- df$ID
    sim_sub <- sim_matrix[ids[ids %in% rownames(sim_matrix)],
                           ids[ids %in% colnames(sim_matrix)]]
    scores <- setNames(-log10(df$p.adjust[df$ID %in% rownames(sim_sub)]),
                        df$ID[df$ID %in% rownames(sim_sub)])
    reduced <- reduceSimMatrix(sim_sub, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
    df %>%
      filter(ID %in% reduced$go[reduced$parent == reduced$go]) %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)   # top 5 terms per cluster, per Fig 1G legend
  }) %>%
  ungroup() %>%
  mutate(log_padj = -log10(p.adjust))

# Full and reduced results, for supplementary reference
write_xlsx(
  list(
    all_go     = all_go %>% filter(p.adjust < 0.05),
    reduced_go = reduced_go %>% filter(p.adjust < 0.05)
  ),
  file.path(output_dir, "Fig1G_GO_enrichment_results.xlsx")
)

# Bar plots: -log10(padj) for top 5 reduced terms, one per cluster
plots_G <- reduced_go %>%
  split(.$cluster) %>%
  imap(function(df, clust) {
    p <- ggplot(df, aes(x = log_padj, y = fct_reorder(Description, log_padj))) +
      geom_col() +
      labs(title = paste("Cluster", clust),
           x = "-log10(Adjusted p-value)",
           y = NULL) +
      theme_few()

    ggsave(file.path(output_dir, paste0("Fig1G_Cluster_", clust, "_GO_top5.pdf")),
           plot = p, width = 8, height = 5, dpi = 300)
    return(p)
  })

## ============================================================
## Session info
## Recorded for reproducibility (package/R versions), written
## alongside figures for GitHub submission.
## ============================================================
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
