## =============================================================================
## Figure 4: TF motif enrichment and gene regulatory analysis in lactocyte
## subclusters
## =============================================================================

library(Signac)
library(Seurat)
library(TFBSTools)
library(JASPAR2024)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(chromVAR)
library(motifmatchr)
library(scico)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggthemes)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

set.seed(1234)

## =============================================================================
## Load data
## =============================================================================

Epi_subset <- readRDS("Epithelial_subset.rds")
Idents(Epi_subset) <- "Celltype"

## Load whole-dataset peak-to-gene links and store in object
links_df <- read.csv("Links_whole_dataset.csv")
links_gr <- makeGRangesFromDataFrame(links_df, keep.extra.columns = TRUE)
Links(Epi_subset) <- links_gr

## Cluster colors (consistent across figures)
cols <- c(
  "Epithelial 0" = "#E24A33",
  "Epithelial 1" = "#348ABD",
  "Epithelial 2" = "#988ED5",
  "Epithelial 3" = "#FBC15E",
  "Epithelial 4" = "#8EBA42"
)

## =============================================================================
## SECTION 1: ADD MOTIFS TO OBJECT
## =============================================================================

DefaultAssay(Epi_subset) <- "peaksbyc"

## Uncomment below if motifs are not already stored in the object
# jaspar   <- JASPAR2024()
# sq24     <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
# motifs24 <- TFBSTools::getMatrixSet(
#   sq24,
#   list(species      = "Homo sapiens",
#        collection   = "CORE",
#        tax_group    = "vertebrates",
#        all_versions = TRUE)
# )
#
# Epi_subset <- AddMotifs(
#   object = Epi_subset,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   pfm    = motifs24
# )

## =============================================================================
## SECTION 2: LOAD DIFFERENTIALLY ACCESSIBLE PEAKS
## =============================================================================

daps <- read.csv("DAPs_Epi_LR.csv")
daps <- daps[daps$p_val_adj < 0.05, ]
colnames(daps)[colnames(daps) == "gene"] <- "peak"

## =============================================================================
## SECTION 3: FIND ENRICHED MOTIFS PER SUBCLUSTER
## =============================================================================

open.peaks   <- AccessiblePeaks(Epi_subset)
meta.feature <- GetAssayData(Epi_subset, assay = "peaksbyc", layer = "meta.features")

cluster_names <- c("Epithelial 0", "Epithelial 1", "Epithelial 2",
                   "Epithelial 3", "Epithelial 4")

peaks_per_cluster <- lapply(cluster_names, function(cl) {
  unique(daps$peak[daps$cluster == cl])
})
names(peaks_per_cluster) <- cluster_names

cat("Peaks used for motif enrichment:\n")
for (cl in cluster_names) {
  cat(sprintf("  %-30s : %d peaks\n", cl, length(peaks_per_cluster[[cl]])))
}

motif.list <- lapply(cluster_names, function(cl) {
  message("Finding motifs for: ", cl)
  peaks_cl <- peaks_per_cluster[[cl]]
  if (length(peaks_cl) == 0) {
    message("  No peaks for ", cl, ", skipping")
    return(NULL)
  }

  background <- MatchRegionStats(
    meta.feature  = meta.feature[open.peaks, ],
    query.feature = meta.feature[peaks_cl, ],
    n             = 50000
  )

  result <- FindMotifs(
    object     = Epi_subset,
    features   = peaks_cl,
    background = background
  )

  write.csv(result,
            file      = paste0(gsub(" ", "", cl), "_Motifs.csv"),
            row.names = FALSE)
  return(result)
})
names(motif.list) <- cluster_names

## =============================================================================
## SECTION 4: OVERLAP ENRICHED MOTIFS WITH DEGs -> TF RANKINGS
## =============================================================================

degs_RNA <- read.csv("degs_sig_epithelial_v2.csv")
epi.DEGs <- split(degs_RNA, degs_RNA$cluster)

tf.rankings <- list()

for (clust in cluster_names) {
  message("Building TF rankings for: ", clust)
  motifs.df <- motif.list[[clust]]
  if (is.null(motifs.df)) next
  sig.motifs <- motifs.df %>% filter(p.adjust < 0.05)
  if (nrow(sig.motifs) == 0) {
    message("  No significant motifs, skipping"); next
  }
  degs <- epi.DEGs[[clust]]
  if (is.null(degs) || nrow(degs) == 0) {
    message("  No DEGs for ", clust, ", skipping"); next
  }
  tf.overlap <- sig.motifs %>%
    filter(motif.name %in% degs$gene)

  if (nrow(tf.overlap) == 0) {
    message("  No TF motifs overlap with DEGs, skipping"); next
  }
  tf.rankings[[clust]] <- tf.overlap %>%
    arrange(desc(fold.enrichment)) %>%
    dplyr::select(motif.id = motif,
                  TF = motif.name,
                  fold.enrichment,
                  p.adjust)

  message("  ", nrow(tf.rankings[[clust]]), " TFs retained")
}

for (clust in names(tf.rankings)) {
  write.csv(tf.rankings[[clust]],
            file      = paste0(gsub(" ", "_", clust), "_tf_rankings.csv"),
            row.names = FALSE)
}

## =============================================================================
## SECTION 5: BUILD top_TFs_per_cluster (motif IDs per cluster)
## =============================================================================

motif.positions <- Epi_subset@assays[["peaksbyc"]]@motifs@positions
top_TFs_per_cluster <- lapply(tf.rankings, function(df) df$motif.id)

## =============================================================================
## SECTION 6: LOAD & FILTER PEAK-GENE LINKS FOR TF ANALYSIS
## =============================================================================

## Re-read links as data frame for TF-peak-gene overlap analysis
links <- read.csv("Links_whole_dataset.csv")
links <- links %>% dplyr::filter(p.adjust < 0.05, score > 0)

## Annotate links with DAP cluster identity
links_dap <- links %>%
  inner_join(
    daps %>% dplyr::select(peak, dap_cluster = cluster),
    by = "peak"
  )

## Apply z-score stringency filter
links_filt <- links_dap %>%
  filter(abs(zscore) > 0.5)

## =============================================================================
## SECTION 7: SELECT TOP 10 TFs PER CLUSTER
## =============================================================================

top10_TFs_per_cluster <- bind_rows(tf.rankings, .id = "cluster") %>%
  filter(!grepl("::", TF)) %>%
  group_by(cluster) %>%
  slice_max(order_by = fold.enrichment, n = 10, with_ties = FALSE) %>%
  ungroup()

## =============================================================================
## SECTION 8: GET TF -> PEAK -> GENE LINKS
## =============================================================================

get_tf_linked_genes <- function(motif_id, tf_name, cluster_name,
                                motif.positions, links_filt) {

  tf_gr <- motif.positions[[motif_id]]
  if (is.null(tf_gr) || length(tf_gr) == 0) return(NULL)

  cluster_links <- links_filt %>% filter(dap_cluster == cluster_name)
  if (nrow(cluster_links) == 0) return(NULL)

  cluster_gr <- GRanges(
    seqnames = cluster_links$seqnames,
    ranges   = IRanges(start = cluster_links$start,
                       end   = cluster_links$end),
    strand   = cluster_links$strand,
    peak     = cluster_links$peak,
    gene     = cluster_links$gene,
    zscore   = cluster_links$zscore,
    p.adjust = cluster_links$p.adjust
  )

  hits <- findOverlaps(tf_gr, cluster_gr)
  if (length(hits) == 0) return(NULL)

  unique_peaks <- unique(cluster_gr$peak[subjectHits(hits)])

  cluster_links %>%
    filter(peak %in% unique_peaks) %>%
    mutate(TF       = tf_name,
           motif.id = motif_id) %>%
    dplyr::select(TF, motif.id, peak, gene, dap_cluster, zscore, p.adjust)
}

## =============================================================================
## SECTION 9: RUN TF-PEAK-GENE ANALYSIS (TOP 10 TFs PER CLUSTER)
## =============================================================================

tf_peak_gene_list <- list()

for (clust in unique(top10_TFs_per_cluster$cluster)) {
  message("Processing: ", clust)
  clust_tfs <- top10_TFs_per_cluster %>% filter(cluster == clust)

  cluster_results <- map2_dfr(
    clust_tfs$motif.id,
    clust_tfs$TF,
    ~get_tf_linked_genes(
      motif_id        = .x,
      tf_name         = .y,
      cluster_name    = clust,
      motif.positions = motif.positions,
      links_filt      = links_filt
    )
  )

  tf_peak_gene_list[[clust]] <- cluster_results
}

tf_peak_gene_df <- bind_rows(tf_peak_gene_list)

## =============================================================================
## SECTION 10: ANNOTATE WITH DEG STATUS & FILTER CONCORDANT
## =============================================================================

tf_peak_gene_annotated <- tf_peak_gene_df %>%
  inner_join(
    degs_RNA %>% dplyr::select(gene, deg_cluster = cluster),
    by = "gene"
  ) %>%
  mutate(concordant = dap_cluster == deg_cluster)

tf_peak_gene_final <- tf_peak_gene_annotated %>%
  filter(concordant) %>%
  dplyr::select(TF, motif.id, peak, gene, cluster = dap_cluster,
                zscore, p.adjust) %>%
  distinct()

write.csv(tf_peak_gene_annotated, "TF_peak_gene_full_annotated.csv", row.names = FALSE)
write.csv(tf_peak_gene_final, "TF_peak_gene_concordant.csv", row.names = FALSE)

## =============================================================================
## PANEL 4A: TF EXPRESSION HEATMAP (all significant TFs, z-scored)
## =============================================================================

## Combine all TF rankings into single data frame
all_tf <- bind_rows(
  lapply(names(tf.rankings), function(cl) {
    df <- tf.rankings[[cl]]
    df$cluster <- gsub("Epithelial ", "Epi ", cl)
    return(df)
  })
)

sig_tfs     <- all_tf %>% dplyr::filter(p.adjust < 0.05)
tf_genes_all <- unique(sig_tfs$TF)

DefaultAssay(Epi_subset) <- "SCT"
tf_genes_all <- tf_genes_all[tf_genes_all %in% rownames(Epi_subset)]

## Average expression matrix
avg_exp    <- AverageExpression(Epi_subset, features = tf_genes_all, group.by = "Celltype")$SCT
mat        <- as.matrix(avg_exp)
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

## Column order
col_order  <- c("Epithelial 0", "Epithelial 3", "Epithelial 1", "Epithelial 2", "Epithelial 4")
mat_scaled <- mat_scaled[, col_order]

## Row clustering
row_clust <- hclust(dist(mat_scaled), method = "average")

## Rank TFs per cluster
ranked_tfs <- all_tf %>%
  dplyr::filter(p.adjust < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(fold.enrichment)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

## Key TFs to always label
key_tfs <- c("ELF5", "EHF", "MAFF", "CREB3L1", "ELF3")

label_tfs <- ranked_tfs %>%
  dplyr::filter(rank <= 10 | TF %in% key_tfs) %>%
  dplyr::filter(TF %in% rownames(mat_scaled))

cluster_colors <- c("Epi 0" = "#E24A33", "Epi 3" = "#FBC15E",
                    "Epi 1" = "#348ABD", "Epi 2" = "#988ED5", "Epi 4" = "#8EBA42")

## Label info with rank strings and colors
label_info <- label_tfs %>%
  mutate(rank_str = paste0(gsub("Epi ", "", cluster), ":#", rank)) %>%
  group_by(TF) %>%
  summarise(
    rank_label   = paste(rank_str, collapse = ", "),
    best_cluster = cluster[which.min(rank)],
    best_rank    = min(rank),
    .groups      = "drop"
  ) %>%
  mutate(full_label  = paste0(TF, " (", rank_label, ")"),
         label_color = cluster_colors[best_cluster])

## Get row indices in clustered order
row_order_idx     <- row_clust$order
row_names_ordered <- rownames(mat_scaled)[row_order_idx]
labeled_positions <- which(row_names_ordered %in% label_info$TF)

label_df <- data.frame(
  TF  = row_names_ordered[labeled_positions],
  pos = labeled_positions,
  stringsAsFactors = FALSE
) %>%
  left_join(label_info, by = "TF")

col_fun <- colorRamp2(
  seq(-2, 2, length.out = 100),
  scico(100, palette = "vikO")
)

ha_mark <- rowAnnotation(
  TF = anno_mark(
    at     = labeled_positions,
    labels = label_df$full_label,
    labels_gp = gpar(fontsize = 7, fontface = "italic", col = label_df$label_color),
    link_width = unit(8, "mm"),
    link_gp    = gpar(col = label_df$label_color, lwd = 0.5)
  )
)

ht <- Heatmap(mat_scaled,
              name             = "z-score",
              col              = col_fun,
              cluster_rows     = as.dendrogram(row_clust),
              cluster_columns  = FALSE,
              column_names_gp  = gpar(fontsize = 10),
              column_names_rot = 45,
              border           = FALSE,
              column_title     = "TF Expression (z-scored)",
              column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
              show_row_dend    = TRUE,
              row_dend_width   = unit(20, "mm"),
              width            = unit(5, "cm"))

pdf("Fig4A_TF_heatmap_ranked_annotated.pdf", width = 5, height = 7)
draw(ht + ha_mark, padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

## =============================================================================
## PANEL 4C: CHROMVAR + RNA + FOOTPRINTING FOR FOSL2 AND NR3C2
## =============================================================================

## --- ChromVAR feature plots ---
DefaultAssay(Epi_subset) <- "chromvar"

fosl2_motif <- all_tf %>%
  dplyr::filter(TF == "FOSL2") %>%
  slice_max(fold.enrichment, n = 1) %>%
  pull(motif.id)

nr3c2_motif <- all_tf %>%
  dplyr::filter(TF == "NR3C2") %>%
  slice_max(fold.enrichment, n = 1) %>%
  pull(motif.id)

p_chromvar_fosl2 <- FeaturePlot(Epi_subset,
                                reduction = "wnn.umap",
                                features  = fosl2_motif,
                                min.cutoff = 0,
                                order = TRUE) +
  scale_color_gradientn(colors = c("grey90", "darkred")) +
  theme_void() +
  theme(plot.title = element_text(face = "bold.italic"))

p_chromvar_nr3c2 <- FeaturePlot(Epi_subset,
                                reduction = "wnn.umap",
                                features  = nr3c2_motif,
                                min.cutoff = 0,
                                order = TRUE) +
  scale_color_gradientn(colors = c("grey90", "darkred")) +
  theme_void() +
  theme(plot.title = element_text(face = "bold.italic"))

## --- SCT RNA feature plots ---
DefaultAssay(Epi_subset) <- "SCT"

p_rna_fosl2 <- FeaturePlot(Epi_subset,
                           features   = "FOSL2",
                           reduction  = "wnn.umap",
                           min.cutoff = 0.7,
                           order = TRUE) +
  theme_void() +
  theme(plot.title = element_text(face = "bold.italic"))

p_rna_nr3c2 <- FeaturePlot(Epi_subset,
                           features   = "NR3C2",
                           reduction  = "wnn.umap",
                           min.cutoff = "q10",
                           order = TRUE) +
  theme_void() +
  theme(plot.title = element_text(face = "bold.italic"))

## Combined ChromVAR + RNA (no legends for assembly in Illustrator)
p1 <- p_chromvar_fosl2 + theme(legend.position = "none")
p2 <- p_rna_fosl2      + theme(legend.position = "none")
p3 <- p_chromvar_nr3c2  + theme(legend.position = "none")
p4 <- p_rna_nr3c2       + theme(legend.position = "none")

combined <- (p1 | p2) / (p3 | p4)
ggsave("Fig4C_chromvar_RNA_FOSL2_NR3C2.pdf", combined, width = 10, height = 10)

## Save individual legends for Illustrator
leg_chromvar_fosl2 <- get_legend(p_chromvar_fosl2)
leg_rna_fosl2      <- get_legend(p_rna_fosl2)
leg_chromvar_nr3c2 <- get_legend(p_chromvar_nr3c2)
leg_rna_nr3c2      <- get_legend(p_rna_nr3c2)

pdf("Fig4C_legend_chromvar_FOSL2.pdf", width = 2, height = 4)
grid.newpage(); grid.draw(leg_chromvar_fosl2)
dev.off()

pdf("Fig4C_legend_RNA_FOSL2.pdf", width = 2, height = 4)
grid.newpage(); grid.draw(leg_rna_fosl2)
dev.off()

pdf("Fig4C_legend_chromvar_NR3C2.pdf", width = 2, height = 4)
grid.newpage(); grid.draw(leg_chromvar_nr3c2)
dev.off()

pdf("Fig4C_legend_RNA_NR3C2.pdf", width = 2, height = 4)
grid.newpage(); grid.draw(leg_rna_nr3c2)
dev.off()

## --- Footprinting plots ---
DefaultAssay(Epi_subset) <- "peaksbyc"

p_foot_fosl2 <- PlotFootprint(Epi_subset, features = "FOSL2") +
  ggtitle("FOSL2 - Footprint") +
  theme_few() +
  theme(plot.title = element_text(face = "bold.italic"),
        legend.position = "right")
ggsave("Fig4C_footprint_FOSL2.pdf", p_foot_fosl2, width = 5, height = 5)

p_foot_nr3c2 <- PlotFootprint(Epi_subset, features = "NR3C2") +
  ggtitle("NR3C2 - Footprint") +
  theme_few() +
  theme(plot.title = element_text(face = "bold.italic"),
        legend.position = "right")
ggsave("Fig4C_footprint_NR3C2.pdf", p_foot_nr3c2, width = 5, height = 5)

## Combined: rows = TFs, columns = ChromVAR | RNA | Footprint
p_tf_combined <- (p_chromvar_fosl2 | p_rna_fosl2 | p_foot_fosl2) /
  (p_chromvar_nr3c2 | p_rna_nr3c2 | p_foot_nr3c2) +
  plot_annotation(title = "TF Activity: ChromVAR, RNA Expression & Footprinting",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("Fig4C_TF_feature_footprint_FOSL2_NR3C2.pdf",
       p_tf_combined, width = 18, height = 10)

## =============================================================================
## PANEL 4D: GO ENRICHMENT HEATMAPS PER TF PER CLUSTER
## =============================================================================

## Get unique TF-cluster pairs from concordant results
tf_cluster_pairs <- tf_peak_gene_final %>%
  distinct(TF, cluster)

go_results_list <- list()

for (i in seq_len(nrow(tf_cluster_pairs))) {
  tf    <- tf_cluster_pairs$TF[i]
  clust <- tf_cluster_pairs$cluster[i]

  target_genes <- tf_peak_gene_final %>%
    filter(TF == !!tf, cluster == !!clust) %>%
    pull(gene) %>%
    unique()

  if (length(target_genes) < 3) next

  ego <- tryCatch({
    enrichGO(gene          = target_genes,
             OrgDb         = org.Hs.eg.db,
             keyType       = "SYMBOL",
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.1,
             readable      = FALSE)
  }, error = function(e) NULL)

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next

  ego_df <- as.data.frame(ego) %>%
    filter(p.adjust < 0.05) %>%
    mutate(TF = tf,
           cluster = clust,
           neg_log10_padj = -log10(p.adjust))

  go_results_list[[paste0(clust, "_", tf)]] <- ego_df
}

go_all <- bind_rows(go_results_list)

## Filter pathways with fewer than 3 genes, add enrichment metrics
go_all <- go_all %>%
  filter(Count >= 3) %>%
  mutate(
    gene_ratio_num  = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio)),
    bg_ratio_num    = as.numeric(sub("/.*", "", BgRatio)) /
      as.numeric(sub(".*/", "", BgRatio)),
    fold_enrichment = gene_ratio_num / bg_ratio_num,
    combined_score  = fold_enrichment * neg_log10_padj
  )

write.csv(go_all, "GO_enrichment_per_TF_per_cluster.csv", row.names = FALSE)

## Generate per-cluster GO heatmaps
for (clust in unique(go_all$cluster)) {
  clust_go <- go_all %>% filter(cluster == clust)
  if (nrow(clust_go) == 0) next

  top_pathways <- clust_go %>%
    group_by(Description) %>%
    summarise(best_score = max(combined_score), .groups = "drop") %>%
    slice_max(order_by = best_score, n = 15, with_ties = FALSE) %>%
    pull(Description)

  mat <- clust_go %>%
    filter(Description %in% top_pathways) %>%
    dplyr::select(TF, Description, combined_score) %>%
    complete(TF, Description, fill = list(combined_score = 0)) %>%
    pivot_wider(names_from = TF, values_from = combined_score, values_fill = 0) %>%
    column_to_rownames("Description") %>%
    as.matrix()

  tf_order_clust <- top10_TFs_per_cluster %>%
    filter(cluster == clust) %>%
    arrange(desc(fold.enrichment)) %>%
    pull(TF)
  tf_order_clust <- tf_order_clust[tf_order_clust %in% colnames(mat)]
  mat <- mat[, tf_order_clust, drop = FALSE]

  col_fun_go <- colorRamp2(c(0, max(mat) / 2, max(mat)),
                           c("white", "#CC0000", "darkred"))

  ht_go <- Heatmap(
    mat,
    name                       = "Combined Score\n(FE x -log10 padj)",
    col                        = col_fun_go,
    cluster_rows               = TRUE,
    clustering_distance_rows   = "euclidean",
    clustering_method_rows     = "ward.D2",
    cluster_columns            = FALSE,
    row_names_side             = "left",
    row_names_gp               = gpar(fontsize = 9),
    column_names_gp            = gpar(fontsize = 10, fontface = "italic"),
    column_names_rot           = 45,
    column_title               = clust,
    column_title_gp            = gpar(fontsize = 12, fontface = "bold"),
    rect_gp                    = gpar(col = "white", lwd = 0.5),
    show_row_dend              = TRUE,
    show_column_dend           = FALSE
  )

  pdf(paste0("Fig4D_GO_heatmap_", gsub(" ", "_", clust), ".pdf"),
      width = 10, height = 6)
  draw(ht_go)
  dev.off()
}

## =============================================================================
## PANELS 4E-F: COVERAGE + EXPRESSION + MOTIFS FOR METABOLIC GENES
## =============================================================================

## Identify Epi 4 TF-linked metabolic genes
genes_of_interest <- c("CD36", "ACSL1", "ACSS1")

epi4_links <- tf_peak_gene_final %>%
  filter(cluster == "Epithelial 4",
         gene %in% genes_of_interest)

## Cross-reference with GO pathway hits
epi4_top15 <- go_all %>%
  filter(cluster == "Epithelial 4") %>%
  group_by(Description) %>%
  summarise(best_score = max(combined_score), .groups = "drop") %>%
  slice_max(order_by = best_score, n = 15, with_ties = FALSE) %>%
  pull(Description)

pathway_gene_hits <- go_all %>%
  filter(cluster == "Epithelial 4",
         Description %in% epi4_top15) %>%
  dplyr::select(TF, Description, geneID, combined_score) %>%
  rowwise() %>%
  mutate(
    genes_in_pathway = list(strsplit(geneID, "/")[[1]]),
    hit_genes        = list(intersect(genes_in_pathway, genes_of_interest)),
    n_hits           = length(hit_genes)
  ) %>%
  filter(n_hits > 0) %>%
  unnest(hit_genes) %>%
  dplyr::select(TF, Description, hit_genes, combined_score) %>%
  arrange(hit_genes, desc(combined_score))

verified_pairs <- pathway_gene_hits %>%
  distinct(TF, hit_genes) %>%
  rename(gene = hit_genes)

epi4_links_verified <- epi4_links %>%
  inner_join(verified_pairs, by = c("TF", "gene"))

## --- ACSL1 coverage plot ---
DefaultAssay(Epi_subset) <- "peaksbyc"
Epi_subset$Celltype <- factor(Epi_subset$Celltype, levels = names(cols))

acsl1_highlight <- GRanges(
  seqnames = c("chr4", "chr4", "chr4"),
  ranges   = IRanges(
    start = c(184805365, 184817702, 184898770),
    end   = c(184806057, 184818430, 184900319)
  )
)

p_acsl1 <- CoveragePlot(
  object     = Epi_subset,
  region     = "chr4-184800000-184905000",
  group.by   = "Celltype",
  annotation = TRUE,
  peaks      = TRUE,
  links      = TRUE,
  ymax       = 220
) & scale_fill_manual(values = cols)

pdf("Fig4_coverage_ACSL1.pdf", width = 8, height = 4)
print(p_acsl1)
dev.off()

## ACSL1 expression
DefaultAssay(Epi_subset) <- "SCT"
p_exp_acsl1 <- VlnPlot(
  object   = Epi_subset,
  features = "ACSL1",
  group.by = "Celltype",
  cols     = cols,
  pt.size  = 0
) +
  theme(text         = element_text(family = "Arial", size = 8),
        axis.text    = element_text(family = "Arial", size = 8),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(face = "italic", family = "Arial", size = 8))

pdf("Fig4_expr_ACSL1.pdf", width = 3, height = 4)
print(p_exp_acsl1)
dev.off()

## ACSL1 motifs
DefaultAssay(Epi_subset) <- "peaksbyc"
acsl1_motif_ids <- epi4_links_verified %>%
  filter(gene == "ACSL1") %>%
  distinct(TF, motif.id)

acsl1_motif_plots <- lapply(seq_len(nrow(acsl1_motif_ids)), function(i) {
  MotifPlot(object = Epi_subset, motifs = acsl1_motif_ids$motif.id[i], assay = "peaksbyc") +
    ggtitle(acsl1_motif_ids$TF[i]) +
    theme(text       = element_text(size = 8),
          plot.title = element_text(size = 8, face = "italic"))
})

pdf("Fig4_motifs_ACSL1.pdf", width = 10, height = 4)
print(wrap_plots(acsl1_motif_plots, ncol = 4))
dev.off()

## --- CD36 coverage plot ---
cd36_highlight <- GRanges(
  seqnames = "chr7",
  ranges   = IRanges(start = 80654702, end = 80655387)
)

DefaultAssay(Epi_subset) <- "peaksbyc"
p_cd36 <- CoveragePlot(
  object           = Epi_subset,
  region           = "chr7-80649000-80661000",
  annotation       = TRUE,
  peaks            = TRUE,
  links            = TRUE,
  region.highlight = cd36_highlight
) & scale_fill_manual(values = cols)

pdf("Fig4_coverage_CD36.pdf", width = 8, height = 4)
print(p_cd36)
dev.off()

## CD36 expression
DefaultAssay(Epi_subset) <- "SCT"
p_exp_cd36 <- VlnPlot(
  object   = Epi_subset,
  features = "CD36",
  group.by = "Celltype",
  cols     = cols,
  pt.size  = 0
) +
  NoLegend() +
  theme(text         = element_text(family = "Arial", size = 8),
        axis.text    = element_text(family = "Arial", size = 8),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(face = "italic", family = "Arial", size = 8))

pdf("Fig4_expr_CD36.pdf", width = 3, height = 4)
print(p_exp_cd36)
dev.off()

## CD36 motifs
cd36_motif_ids <- epi4_links_verified %>%
  filter(gene == "CD36") %>%
  distinct(TF, motif.id)

DefaultAssay(Epi_subset) <- "peaksbyc"
cd36_motif_plots <- lapply(seq_len(nrow(cd36_motif_ids)), function(i) {
  MotifPlot(object = Epi_subset, motifs = cd36_motif_ids$motif.id[i], assay = "peaksbyc") +
    ggtitle(cd36_motif_ids$TF[i]) +
    theme(text       = element_text(size = 8),
          plot.title = element_text(size = 8, face = "italic"))
})

pdf("Fig4_motifs_CD36.pdf", width = 10, height = 4)
print(wrap_plots(cd36_motif_plots, ncol = 4))
dev.off()

## --- ACSS1 coverage plot ---
acss1_highlight <- GRanges(
  seqnames = c("chr20", "chr20"),
  ranges   = IRanges(
    start = c(25039527, 25047845),
    end   = c(25040766, 25048905)
  )
)

DefaultAssay(Epi_subset) <- "peaksbyc"
p_acss1 <- CoveragePlot(
  object           = Epi_subset,
  region           = "chr20-25034000-25054000",
  annotation       = TRUE,
  peaks            = TRUE,
  links            = TRUE,
  region.highlight = acss1_highlight
) & scale_fill_manual(values = cols)

pdf("Fig4_coverage_ACSS1.pdf", width = 8, height = 4)
print(p_acss1)
dev.off()

## ACSS1 expression
DefaultAssay(Epi_subset) <- "SCT"
p_exp_acss1 <- VlnPlot(
  object   = Epi_subset,
  features = "ACSS1",
  group.by = "Celltype",
  cols     = cols,
  pt.size  = 0
) +
  NoLegend() +
  theme(text         = element_text(family = "Arial", size = 8),
        axis.text    = element_text(family = "Arial", size = 8),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(face = "italic", family = "Arial", size = 8))

pdf("Fig4_expr_ACSS1.pdf", width = 3, height = 4)
print(p_exp_acss1)
dev.off()

## ACSS1 motifs
acsl1_motif_ids_acss1 <- epi4_links_verified %>%
  filter(gene == "ACSS1") %>%
  distinct(TF, motif.id)

DefaultAssay(Epi_subset) <- "peaksbyc"
acss1_motif_plots <- lapply(seq_len(nrow(acsl1_motif_ids_acss1)), function(i) {
  MotifPlot(object = Epi_subset, motifs = acsl1_motif_ids_acss1$motif.id[i], assay = "peaksbyc") +
    ggtitle(acsl1_motif_ids_acss1$TF[i]) +
    theme(text       = element_text(size = 8),
          plot.title = element_text(size = 8, face = "italic"))
})

pdf("Fig4_motifs_ACSS1.pdf", width = 10, height = 4)
print(wrap_plots(acss1_motif_plots, ncol = 3))
dev.off()
