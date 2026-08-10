## =============================================================================
## Figure 3: Chromatin accessibility reveals distinct regulatory architectures
## in LC1- and LC2-like lactocytes
## =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(patchwork)
library(pheatmap)
library(grid)

set.seed(1234)

## =============================================================================
## Load data
## =============================================================================

Epi_subset <- readRDS("Epithelial_subset.rds")

my_colors <- c(
  "Epithelial 0" = "#E24A33",
  "Epithelial 1" = "#348ABD",
  "Epithelial 2" = "#988ED5",
  "Epithelial 3" = "#FBC15E",
  "Epithelial 4" = "#8EBA42"
)

# Load DAPs, DEGs, and peak-to-gene links
daps <- read.csv("DAPs_Epi_LR.csv")
daps <- daps[daps$p_val_adj < 0.05, ]
colnames(daps)[colnames(daps) == "gene"] <- "peak"

degs_RNA <- read.csv("degs_sig_epithelial_v2.csv")
degs_RNA <- degs_RNA[degs_RNA$p_val_adj < 0.05, ]

links <- read.csv("Links_whole_dataset.csv")
links <- links[links$p.adjust < 0.05 & links$score > 0, ]

## =============================================================================
## Helper: classify DAPs by linkage status
## =============================================================================

epi_clusters <- paste0("Epithelial ", 0:4)

classify_dap <- function(peak, cluster, links, degs_RNA) {
  linked <- links[links$peak == peak, ]
  if (nrow(linked) == 0) return("No gene linked")

  # Check if any linked gene is a DEG in the SAME cluster
  for (i in seq_len(nrow(linked))) {
    deg_match <- degs_RNA[degs_RNA$gene == linked$gene[i] & degs_RNA$cluster == cluster, ]
    if (nrow(deg_match) > 0) return("Linked to same-cluster DEG")
  }

  # Check if any linked gene is a DEG in a DIFFERENT epithelial cluster
  other_epi <- setdiff(epi_clusters, cluster)
  for (i in seq_len(nrow(linked))) {
    deg_match <- degs_RNA[degs_RNA$gene == linked$gene[i] & degs_RNA$cluster %in% other_epi, ]
    if (nrow(deg_match) > 0) return("Linked to other-cluster DEG")
  }

  return("No gene linked")
}

## =============================================================================
## Fig 3A: Stacked bar chart of DAPs by linkage status
## =============================================================================

daps_epi <- daps[daps$cluster %in% epi_clusters, ]
daps_epi$cluster_num <- as.integer(gsub("Epithelial ", "", daps_epi$cluster))

daps_epi$category <- mapply(
  classify_dap,
  peak = daps_epi$peak,
  cluster = daps_epi$cluster,
  MoreArgs = list(links = links, degs_RNA = degs_RNA)
)

# Summarize counts
plot_data <- daps_epi %>%
  group_by(cluster_num, category) %>%
  summarise(count = n(), .groups = "drop")

totals <- plot_data %>%
  group_by(cluster_num) %>%
  summarise(total = sum(count), .groups = "drop")

# Order: LC1-like (0, 3) -> Intermediate (1) -> LC2-like (4, 2)
plot_data$cluster_num <- factor(plot_data$cluster_num, levels = c(0, 3, 1, 4, 2))
totals$cluster_num    <- factor(totals$cluster_num, levels = c(0, 3, 1, 4, 2))

plot_data$category <- factor(plot_data$category, levels = c(
  "Linked to same-cluster DEG",
  "Linked to other-cluster DEG",
  "No gene linked"
))

p3a <- ggplot(plot_data, aes(x = cluster_num, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = totals,
            aes(x = cluster_num, y = total, label = total),
            inherit.aes = FALSE, vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c("No gene linked" = "#D9D9D9",
               "Linked to same-cluster DEG"  = "#D55E00",
               "Linked to other-cluster DEG" = "#0072B2"),
    name = NULL
  ) +
  scale_x_discrete(labels = paste0("Epi ", c(0, 3, 1, 4, 2))) +
  # Lineage group annotations
  annotate("segment", x = 0.6, xend = 2.4, y = -120, yend = -120,
           linewidth = 0.8, color = "black") +
  annotate("text", x = 1.5, y = -220, label = "LC1-like",
           size = 4, fontface = "italic") +
  annotate("segment", x = 2.7, xend = 3.3, y = -120, yend = -120,
           linewidth = 0.8, color = "black") +
  annotate("text", x = 3, y = -220, label = "Intermediate",
           size = 4, fontface = "italic") +
  annotate("segment", x = 3.6, xend = 5.4, y = -120, yend = -120,
           linewidth = 0.8, color = "black") +
  annotate("text", x = 4.5, y = -220, label = "LC2-like",
           size = 4, fontface = "italic") +
  coord_cartesian(clip = "off", ylim = c(0, NA)) +
  labs(x = NULL, y = "Number of DAPs") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        plot.margin = margin(t = 10, r = 10, b = 50, l = 10),
        axis.text.x = element_text(size = 12))

ggsave("Fig3A_DAPs_stacked_barplot.pdf", p3a, width = 8, height = 7)

## =============================================================================
## Fig 3B: Genomic feature annotation of DAPs (DEG-linked vs unlinked)
## =============================================================================

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_to_gr <- function(peak_strings) {
  parts <- strsplit(peak_strings, "-")
  GRanges(seqnames = sapply(parts, `[`, 1),
          ranges = IRanges(start = as.integer(sapply(parts, `[`, 2)),
                           end   = as.integer(sapply(parts, `[`, 3))))
}

# Annotate all DAPs with ChIPseeker
cluster_order <- c(0, 3, 1, 4, 2)
categories <- c("Linked to same-cluster DEG", "Linked to other-cluster DEG", "No gene linked")

all_anno_details <- list()
for (cat in categories) {
  for (cl in cluster_order) {
    subset_df <- daps_epi[daps_epi$cluster_num == cl & daps_epi$category == cat, ]
    if (nrow(subset_df) < 2) next
    gr <- peak_to_gr(subset_df$peak)
    peakAnno <- annotatePeak(gr, TxDb = txdb, level = "gene",
                             annoDb = "org.Hs.eg.db", verbose = FALSE)
    all_anno_details[[paste0(cat, "||Epi ", cl)]] <- peakAnno
  }
}

# Extract annotation table
anno <- do.call(rbind, lapply(names(all_anno_details), function(key) {
  parts <- strsplit(key, "\\|\\|")[[1]]
  det <- as.data.frame(all_anno_details[[key]])
  det$category <- parts[1]
  det$cluster <- parts[2]
  det
}))

write.csv(anno, "DAPs_annotation_for_stats.csv", row.names = FALSE)

# Simplify annotation categories
simplify_annotation <- function(x) {
  dplyr::case_when(
    grepl("Promoter (<=1kb)", x, fixed = TRUE)  ~ "Promoter (<=1kb)",
    grepl("Promoter (1-2kb)", x, fixed = TRUE)  ~ "Promoter (1-2kb)",
    grepl("Promoter (2-3kb)", x, fixed = TRUE)  ~ "Promoter (2-3kb)",
    grepl("5' UTR", x)                          ~ "5' UTR",
    grepl("3' UTR", x)                          ~ "3' UTR",
    grepl("Exon", x)                            ~ "Exon",
    grepl("Intron", x)                           ~ "Intron",
    grepl("Downstream", x)                       ~ "Downstream",
    grepl("Distal Intergenic", x)                ~ "Distal Intergenic",
    TRUE                                         ~ x
  )
}

# Exclude Epi 2 (insufficient DEG-linked DAPs)
anno <- anno[anno$cluster != "Epi 2", ]
anno$annotation_simple <- factor(
  simplify_annotation(anno$annotation),
  levels = c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
             "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
)

# Split into DEG-linked vs unlinked
linked_df   <- anno[anno$category == "Linked to same-cluster DEG", ]
linked_df$bar_group <- "DEG-linked DAPs"
unlinked_df <- anno[anno$category != "Linked to same-cluster DEG", ]
unlinked_df$bar_group <- "Unlinked DAPs"
plot_df <- rbind(linked_df, unlinked_df)

# Compute proportions
prop_df <- plot_df %>%
  count(cluster, bar_group, annotation_simple, .drop = FALSE) %>%
  group_by(cluster, bar_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_df$bar_group <- factor(prop_df$bar_group,
                            levels = c("DEG-linked DAPs", "Unlinked DAPs"))

count_df <- plot_df %>%
  count(cluster, bar_group, name = "total")
count_df$bar_group <- factor(count_df$bar_group,
                             levels = c("DEG-linked DAPs", "Unlinked DAPs"))

p3b <- ggplot(prop_df, aes(x = bar_group, y = prop, fill = annotation_simple)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_text(data = count_df,
            aes(x = bar_group, y = -0.05, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 2.5) +
  facet_wrap(~ cluster, nrow = 1, strip.position = "bottom") +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0.1, 0.02))) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(
    name = "Genomic Feature",
    values = c("Promoter (<=1kb)"  = "#E64B35", "Promoter (1-2kb)"  = "#F39B7F",
               "Promoter (2-3kb)"  = "#FDB863", "5' UTR"            = "#91D1C2",
               "3' UTR"            = "#8491B4", "Exon"              = "#3C5488",
               "Intron"            = "#00A087", "Downstream"        = "#B09C85",
               "Distal Intergenic" = "#7E6148")
  ) +
  labs(x = NULL, y = "Proportion of DAPs") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        plot.margin = margin(t = 5, r = 10, b = 25, l = 5))

ggsave("Fig3B_genomic_features.pdf", p3b, width = 6, height = 6, dpi = 300)

## =============================================================================
## Fig 3C: Butterfly lollipop plot — top 15 DEGs by linked DAP count
##         (Lactocyte 0 vs Lactocyte 4)
## =============================================================================

# Count DAPs linked to each DEG per cluster
epi_clusters_04 <- paste0("Epithelial ", c(0, 1, 3, 4))
gene_dap_counts <- list()
dap_gene_mapping <- list()

for (cl in epi_clusters_04) {
  cl_num <- gsub("Epithelial ", "", cl)
  cl_daps <- daps[daps$cluster == cl, ]
  cl_degs <- degs_RNA[degs_RNA$cluster == cl, ]

  cl_links <- links %>%
    dplyr::filter(peak %in% cl_daps$peak, gene %in% cl_degs$gene)

  if (nrow(cl_links) == 0) next

  daps_per_gene <- cl_links %>%
    group_by(gene) %>%
    summarise(n_daps = n_distinct(peak),
              peaks = paste(unique(peak), collapse = ";"),
              mean_link_score = mean(score), .groups = "drop") %>%
    arrange(desc(n_daps)) %>%
    mutate(cluster = paste0("Epi ", cl_num))

  gene_dap_counts[[cl]] <- daps_per_gene
  dap_gene_mapping[[cl]] <- cl_links %>% mutate(cluster = paste0("Epi ", cl_num))
}

all_gene_dap <- bind_rows(gene_dap_counts)
write.csv(all_gene_dap, "DAPs_per_gene_counts.csv", row.names = FALSE)

# Butterfly plot: Epi 0 (left) vs Epi 4 (right), top 15 each
top0 <- all_gene_dap %>%
  dplyr::filter(cluster == "Epi 0") %>%
  arrange(desc(n_daps)) %>%
  slice_head(n = 15) %>%
  mutate(n_daps_plot = -n_daps)

top4 <- all_gene_dap %>%
  dplyr::filter(cluster == "Epi 4") %>%
  arrange(desc(n_daps)) %>%
  slice_head(n = 15) %>%
  mutate(n_daps_plot = n_daps)

plot_df <- bind_rows(top0, top4) %>%
  mutate(gene = reorder(gene, abs(n_daps_plot), FUN = max))

max_val <- max(abs(plot_df$n_daps_plot))
cluster_colors <- c("Epi 0" = "#E24A33", "Epi 4" = "#8EBA42")

p3c <- ggplot(plot_df, aes(x = n_daps_plot, y = gene, color = cluster)) +
  geom_segment(aes(x = 0, xend = n_daps_plot, y = gene, yend = gene),
               linewidth = 0.5) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = cluster_colors) +
  scale_x_continuous(
    breaks = seq(-max_val, max_val, by = 1),
    labels = abs(seq(-max_val, max_val, by = 1)),
    limits = c(-max_val - 0.5, max_val + 0.5)
  ) +
  annotate("text", x = -max_val/2, y = Inf, label = "Epi 0",
           vjust = -0.5, fontface = "bold", color = "#E24A33", size = 4) +
  annotate("text", x = max_val/2, y = Inf, label = "Epi 4",
           vjust = -0.5, fontface = "bold", color = "#8EBA42", size = 4) +
  coord_cartesian(clip = "off") +
  labs(x = "Number of linked DAPs", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(t = 20, r = 10, b = 5, l = 5))

ggsave("Fig3C_butterfly_lollipop_Epi0_Epi4.pdf", p3c, width = 5, height = 5)

## =============================================================================
## Fig 3D: EHD4 coverage plot (LC1-like, Lactocyte 0)
## =============================================================================

# Store whole-dataset peak-to-gene links in the epithelial object
links_gr <- makeGRangesFromDataFrame(
  links[, c("seqnames", "start", "end", "strand", "score", "gene")],
  keep.extra.columns = TRUE
)
Links(Epi_subset) <- links_gr

DefaultAssay(Epi_subset) <- "peaksbyc"
Idents(Epi_subset) <- Epi_subset$Celltype

# Get EHD4 DAPs from Epi 0
ehd4_peaks <- all_gene_dap %>%
  dplyr::filter(gene == "EHD4", cluster == "Epi 0") %>%
  dplyr::pull(peaks) %>%
  paste(collapse = ";") %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

ehd4_parts <- strsplit(ehd4_peaks, "-")
ehd4_highlight <- GRanges(
  seqnames = sapply(ehd4_parts, `[`, 1),
  ranges = IRanges(start = as.numeric(sapply(ehd4_parts, `[`, 2)),
                   end   = as.numeric(sapply(ehd4_parts, `[`, 3)))
)

epi_idents <- c("Epithelial 0", "Epithelial 1", "Epithelial 2",
                "Epithelial 3", "Epithelial 4")

cov_ehd4 <- CoveragePlot(
  object = Epi_subset,
  idents = epi_idents,
  region = "chr15-42050000-42100000",
  extend.downstream = 1000,
  region.highlight = ehd4_highlight,
  links = FALSE,
  ymax = 230
) & scale_fill_manual(values = my_colors) &
  theme(text       = element_text(family = "Arial", size = 5),
        axis.text  = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5),
        strip.text = element_text(family = "Arial", size = 5))

ggsave("Fig3D_EHD4_coverage.pdf", cov_ehd4,
       width = 3.5, height = 2.5, device = cairo_pdf)

# EHD4 expression violin
DefaultAssay(Epi_subset) <- "RNA"

expr_ehd4 <- ExpressionPlot(
  object = Epi_subset,
  features = "EHD4",
  idents = epi_idents
) & scale_fill_manual(values = my_colors) &
  theme(text       = element_text(family = "Arial", size = 5),
        axis.text  = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5))

ggsave("Fig3D_EHD4_expression.pdf", expr_ehd4,
       width = 2, height = 2.5, device = cairo_pdf)

## =============================================================================
## Fig 3E: ELF5 coverage plot (LC2-like, Lactocyte 4 and Lactocyte 1)
##         with colored region highlights by cluster specificity
## =============================================================================

DefaultAssay(Epi_subset) <- "peaksbyc"
Idents(Epi_subset) <- Epi_subset$Celltype

# Get ELF5 DAPs from Epi 1 and Epi 4
epi1_peaks <- all_gene_dap %>%
  dplyr::filter(gene == "ELF5", cluster == "Epi 1") %>%
  dplyr::pull(peaks) %>%
  strsplit(";") %>%
  unlist()

epi4_peaks <- all_gene_dap %>%
  dplyr::filter(gene == "ELF5", cluster == "Epi 4") %>%
  dplyr::pull(peaks) %>%
  strsplit(";") %>%
  unlist()

shared_peaks   <- intersect(epi1_peaks, epi4_peaks)
epi1_only      <- setdiff(epi1_peaks, epi4_peaks)
epi4_only      <- setdiff(epi4_peaks, epi1_peaks)

# Convert to GRanges
peaks_to_gr <- function(peak_strings) {
  parts <- strsplit(peak_strings, "-")
  GRanges(seqnames = sapply(parts, `[`, 1),
          ranges = IRanges(start = as.numeric(sapply(parts, `[`, 2)),
                           end   = as.numeric(sapply(parts, `[`, 3))))
}

shared_gr   <- if (length(shared_peaks) > 0)  peaks_to_gr(shared_peaks)   else GRanges()
epi1_gr     <- if (length(epi1_only) > 0)     peaks_to_gr(epi1_only)      else GRanges()
epi4_gr     <- if (length(epi4_only) > 0)     peaks_to_gr(epi4_only)      else GRanges()

# Coverage plot with all ELF5-linked DAPs highlighted
# Note: CoveragePlot region.highlight takes a single GRanges; for the
# multi-color highlight (green=Epi4, blue=Epi1, black=shared), generate
# three separate plots or overlay manually in Illustrator.

# Option 1: Combined highlight (all peaks in one color)
all_elf5_gr <- c(shared_gr, epi1_gr, epi4_gr)

cov_elf5 <- CoveragePlot(
  object = Epi_subset,
  idents = epi_idents,
  region = "ELF5",
  extend.upstream = 5000,
  extend.downstream = 5000,
  region.highlight = all_elf5_gr,
  links = FALSE,
  ymax = 230
) & scale_fill_manual(values = my_colors) &
  theme(text       = element_text(family = "Arial", size = 5),
        axis.text  = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5),
        strip.text = element_text(family = "Arial", size = 5))

ggsave("Fig3E_ELF5_coverage_all_highlights.pdf", cov_elf5,
       width = 3.5, height = 2.5, device = cairo_pdf)

# Option 2: Separate plots per highlight color for Illustrator overlay
# Epi 4-specific DAPs (green boxes in figure)
cov_elf5_epi4 <- CoveragePlot(
  object = Epi_subset, idents = epi_idents, region = "ELF5",
  extend.upstream = 5000, extend.downstream = 5000,
  region.highlight = epi4_gr, links = FALSE, ymax = 230
) & scale_fill_manual(values = my_colors) &
  theme(text = element_text(family = "Arial", size = 5),
        axis.text = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5),
        strip.text = element_text(family = "Arial", size = 5))

ggsave("Fig3E_ELF5_coverage_Epi4_DAPs.pdf", cov_elf5_epi4,
       width = 3.5, height = 2.5, device = cairo_pdf)

# Epi 1-specific DAPs (blue boxes in figure)
cov_elf5_epi1 <- CoveragePlot(
  object = Epi_subset, idents = epi_idents, region = "ELF5",
  extend.upstream = 5000, extend.downstream = 5000,
  region.highlight = epi1_gr, links = FALSE, ymax = 230
) & scale_fill_manual(values = my_colors) &
  theme(text = element_text(family = "Arial", size = 5),
        axis.text = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5),
        strip.text = element_text(family = "Arial", size = 5))

ggsave("Fig3E_ELF5_coverage_Epi1_DAPs.pdf", cov_elf5_epi1,
       width = 3.5, height = 2.5, device = cairo_pdf)

# Shared DAPs (black boxes in figure)
cov_elf5_shared <- CoveragePlot(
  object = Epi_subset, idents = epi_idents, region = "ELF5",
  extend.upstream = 5000, extend.downstream = 5000,
  region.highlight = shared_gr, links = FALSE, ymax = 230
) & scale_fill_manual(values = my_colors) &
  theme(text = element_text(family = "Arial", size = 5),
        axis.text = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5),
        strip.text = element_text(family = "Arial", size = 5))

ggsave("Fig3E_ELF5_coverage_shared_DAPs.pdf", cov_elf5_shared,
       width = 3.5, height = 2.5, device = cairo_pdf)

# ELF5 expression violin
DefaultAssay(Epi_subset) <- "RNA"

expr_elf5 <- ExpressionPlot(
  object = Epi_subset,
  features = "ELF5",
  idents = epi_idents
) & scale_fill_manual(values = my_colors) &
  theme(text       = element_text(family = "Arial", size = 5),
        axis.text  = element_text(family = "Arial", size = 5),
        axis.title = element_text(family = "Arial", size = 5))

ggsave("Fig3E_ELF5_expression.pdf", expr_elf5,
       width = 2, height = 2.5, device = cairo_pdf)
