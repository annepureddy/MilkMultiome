## =============================================================================
## Figure 5B: SCENITH / CENITH Metabolic Profiling
##
## Input: Table S11 (SCENITH_MFI_concat.xlsx)
##   - One sheet per donor, columns: Drug, CD36_MFI, CD36neg_MFI,
##     Viability, DaysPP
##
## Panels:
##   5B-i:   CD36+ metabolic parameter heatmap
##   5B-ii:  CD36+ effect size dot plot (delta normalised MFI vs Control)
##   5B-iii: Raw MFI barplots — CD36+ (connected donor dots)
##   5B-iv:  Raw MFI barplots — CD36- (connected donor dots)
## =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggbeeswarm)
library(coin)
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(grid)

set.seed(1234)

## =============================================================================
## SHARED CONFIG
## =============================================================================

FILE    <- "SCENITH_MFI_concat.xlsx"
SHEETS  <- excel_sheets(FILE)
EXCLUDE <- c("FC041", "FC048")
DRUGS   <- c("Control", "2DG", "Oligomyocin", "Etomoxir", "2DG_Oligo")
LABELS  <- c("Control", "2DG", "Oligomycin", "Etomoxir", "2DG+Olig")

DONOR_COLORS <- c(
  FC035 = "#A63030", FC042 = "#A68430", FC047 = "#72A630",
  FC052 = "#30A641", FC053 = "#30A696", FC054 = "#3061A6",
  FC056 = "#5330A6", FC057 = "#A630A4", FC058 = "#A6304F"
)

COLOR_POS <- "#3A7DC9"
COLOR_NEG <- "#D94F3D"

## =============================================================================
## LOAD DATA (both CD36+ and CD36-)
## =============================================================================

raw <- map_dfr(setdiff(SHEETS, EXCLUDE), function(donor) {
  df <- read_excel(FILE, sheet = donor) |>
    rename_with(trimws) |>
    filter(!is.na(Drug)) |>
    mutate(donor = donor)

  viab <- df$Viability[!is.na(df$Viability)][1]
  days <- df$DaysPP[!is.na(df$DaysPP)][1]

  df |>
    select(donor, drug = Drug,
           CD36pos = CD36_MFI,
           CD36neg = CD36neg_MFI) |>         ## <-- adjust column name if needed
    mutate(
      Viability = if (!is.null(viab) && !is.na(viab)) round(as.numeric(viab) * 100, 1) else NA_real_,
      DaysPP    = if (!is.null(days) && !is.na(days)) as.integer(days) else NA_integer_
    )
})

## =============================================================================
## PANEL 5B-i: CD36+ METABOLIC PARAMETER HEATMAP
## =============================================================================

CAP <- 150

wide <- raw |>
  select(donor, drug, CD36pos) |>
  pivot_wider(names_from = drug, values_from = CD36pos)

params <- wide |>
  mutate(
    dDG      = Control - `2DG`,
    dOG      = Control - Oligomyocin,
    dETO     = Control - Etomoxir,
    dDGO     = Control - `2DG_Oligo`,
    Gluc_Dep  = 100 * dDG  / dDGO,
    Mitoc_Dep = 100 * dOG  / dDGO,
    Glyc_Cap  = 100 - Mitoc_Dep,
    FAO_Cap   = 100 - Gluc_Dep,
    LCFA_Dep  = 100 * dETO / dDGO
  ) |>
  select(donor, Gluc_Dep, Mitoc_Dep, Glyc_Cap, FAO_Cap, LCFA_Dep) |>
  mutate(across(where(is.numeric), \(x) pmax(pmin(x, CAP), -CAP)))

meta <- raw |> distinct(donor, Viability, DaysPP)

mat <- params |>
  column_to_rownames("donor") |>
  as.matrix()

colnames(mat) <- c("Glucose\nDependency", "Mitochondrial\nDependency",
                    "Glycolytic\nCapacity", "FAO_AAO\nCapacity", "FAO\nDependency")

## Clustering
dist_donors <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs", method = "pearson"))
clust_donors <- hclust(dist_donors, method = "average")

## Colour scales
col_main <- colorRamp2(c(-150, -75, 0, 75, 150),
                       c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))

viab_range <- range(meta$Viability, na.rm = TRUE)
col_viab <- colorRamp2(c(viab_range[1], mean(viab_range), viab_range[2]),
                       c("#edf8e9", "#66c2a4", "#006d2c"))

weeks <- meta$DaysPP / 7
weeks_range <- range(weeks, na.rm = TRUE)
days_range  <- range(meta$DaysPP, na.rm = TRUE)
col_days <- colorRamp2(c(days_range[1], mean(days_range), days_range[2]),
                       c("#F0F0F0", "#800080", "purple4"))

meta_ordered <- meta |>
  filter(donor %in% rownames(mat)) |>
  column_to_rownames("donor")

ra <- rowAnnotation(
  `Viability (%)` = meta_ordered[rownames(mat), "Viability"],
  `Days PP`       = meta_ordered[rownames(mat), "DaysPP"],
  col = list(`Viability (%)` = col_viab, `Days PP` = col_days),
  annotation_legend_param = list(
    `Viability (%)` = list(title = "Viability (%)", title_gp = gpar(fontsize = 9),
                           labels_gp = gpar(fontsize = 8)),
    `Days PP`       = list(title = "Days PP", title_gp = gpar(fontsize = 9),
                           labels_gp = gpar(fontsize = 8))
  ),
  width = unit(1.4, "cm"), gap = unit(2, "mm"), border = TRUE,
  annotation_name_gp = gpar(fontsize = 9), annotation_name_side = "top"
)

means <- colMeans(mat, na.rm = TRUE)
ba <- HeatmapAnnotation(
  Mean = anno_barplot(means,
                      gp = gpar(fill = col_main(means), col = "white", lwd = 0.5),
                      ylim = c(-150, 150), axis = TRUE,
                      axis_param = list(gp = gpar(fontsize = 7)),
                      height = unit(2.0, "cm")),
  annotation_name_gp = gpar(fontsize = 9), annotation_name_side = "left"
)

ht <- Heatmap(
  mat, name = "Score", col = col_main,
  cluster_rows = clust_donors, cluster_columns = FALSE,
  show_row_dend = TRUE, row_dend_side = "right", row_dend_width = unit(1.5, "cm"),
  show_column_dend = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9), column_names_rot = 0,
  column_names_centered = TRUE,
  left_annotation = ra, bottom_annotation = ba,
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "Score", at = c(-150, -75, 0, 75, 150),
    labels = c("≤−150", "−75", "0", "75", "≥150"),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8), legend_height = unit(4, "cm")
  ),
  width = unit(9, "cm"), height = unit(7, "cm")
)

pdf("Fig5B_heatmap_cd36pos.pdf", width = 10, height = 7)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(5, 5, 10, 5), "mm"))
dev.off()

## =============================================================================
## PANEL 5B-ii: EFFECT SIZE DOT PLOT (CD36+)
## =============================================================================

## Normalise to donor's own Control
ctrl <- raw |>
  filter(drug == "Control") |>
  select(donor, ctrl_mfi = CD36pos)

data_es <- raw |>
  left_join(ctrl, by = "donor") |>
  mutate(norm = CD36pos / ctrl_mfi) |>
  filter(drug != "Control") |>
  mutate(
    delta      = norm - 1,
    drug       = factor(drug, levels = DRUGS[-1]),
    drug_label = factor(drug, levels = DRUGS[-1], labels = LABELS[-1])
  )

## Statistics
stats_res <- data_es |>
  group_by(drug, drug_label) |>
  summarise(
    n         = n(),
    mean_d    = mean(delta, na.rm = TRUE),
    sem_d     = sd(delta, na.rm = TRUE) / sqrt(n()),
    ci95      = qt(0.975, df = n() - 1) * sem_d,
    shapiro_p = shapiro.test(delta)$p.value,
    p_raw     = {
      if (shapiro.test(delta)$p.value > 0.05) {
        t.test(delta, mu = 0)$p.value
      } else {
        wilcox.test(delta, mu = 0, exact = FALSE)$p.value
      }
    },
    test_used = if_else(shapiro.test(delta)$p.value > 0.05,
                        "paired t-test", "Wilcoxon"),
    .groups = "drop"
  ) |>
  mutate(
    p_adj = p.adjust(p_raw, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    sig   = p_adj < 0.05,
    label = paste0(stars, "\np[adj]=", formatC(p_adj, digits = 3, format = "f"))
  )

label_df <- data_es |>
  group_by(drug_label) |>
  summarise(y_max = max(delta, na.rm = TRUE), .groups = "drop") |>
  left_join(stats_res |> select(drug_label, mean_d, ci95, stars, p_adj, sig),
            by = "drug_label") |>
  mutate(
    y_label    = y_max + 0.06,
    label_text = paste0(stars, "\np[adj]=", sprintf("%.3f", p_adj)),
    star_color = if_else(sig, "#b01010", "#777777"),
    population = "CD36⁺ Enriched"
  )

p_es <- ggplot(data_es, aes(x = drug_label, y = delta)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.6) +
  geom_beeswarm(aes(colour = donor), size = 3.5, cex = 2.8, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.42,
               colour = COLOR_POS, linewidth = 0.9, middle.linewidth = 1.8,
               show.legend = FALSE) +
  geom_errorbar(data = stats_res,
                aes(x = drug_label, y = mean_d,
                    ymin = mean_d - ci95, ymax = mean_d + ci95),
                width = 0.13, colour = COLOR_POS, linewidth = 0.9,
                inherit.aes = FALSE) +
  geom_text(data = label_df,
            aes(x = drug_label, y = y_label, label = label_text, colour = star_color),
            size = 2.8, lineheight = 0.95, vjust = 0,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(
    name = "Donor",
    values = c(DONOR_COLORS, "#b01010" = "#b01010", "#777777" = "#777777"),
    breaks = names(DONOR_COLORS), labels = names(DONOR_COLORS)
  ) +
  guides(colour = guide_legend(title = "Donor",
                               override.aes = list(size = 3.5, alpha = 1, shape = 16),
                               ncol = 1)) +
  labs(x = NULL,
       y = expression(Delta ~ "Normalised Puromycin MFI"),
       subtitle = paste0("n = ", length(unique(data_es$donor)),
                         " donors | bar = mean, whiskers = 95% CI | BH-adjusted p")) +
  facet_wrap(~ population) +
  theme_classic(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "#E8EEF5", colour = "grey70", linewidth = 0.5),
    strip.text          = element_text(face = "bold", size = 12, colour = "#1A3A5C"),
    plot.subtitle       = element_text(size = 9, hjust = 0.5, colour = "grey45"),
    axis.text.x         = element_text(size = 11, colour = "black"),
    panel.grid.major.y  = element_line(colour = "grey92", linewidth = 0.4),
    legend.position     = "right"
  )

ggsave("Fig5B_effect_size_cd36pos.pdf", plot = p_es, width = 8.5, height = 5.5)

## =============================================================================
## PANELS 5B-iii & 5B-iv: RAW MFI BARPLOTS WITH CONNECTED DONOR DOTS
## =============================================================================

## Reshape to long format with both populations
mfi_long <- raw |>
  select(donor, drug, CD36pos, CD36neg) |>
  pivot_longer(cols = c(CD36pos, CD36neg),
               names_to = "population", values_to = "MFI") |>
  mutate(
    drug       = factor(drug, levels = DRUGS, labels = LABELS),
    population = factor(population,
                        levels = c("CD36pos", "CD36neg"),
                        labels = c("CD36+", "CD36-"))
  )

## Helper function for paired MFI barplot
plot_mfi_paired <- function(df, pop_label, pop_color) {

  ## Mean MFI per condition (for bars)
  mean_df <- df |>
    group_by(drug) |>
    summarise(mean_MFI = mean(MFI, na.rm = TRUE),
              sem      = sd(MFI, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  ggplot() +
    ## Bars: mean MFI per condition
    geom_col(data = mean_df,
             aes(x = drug, y = mean_MFI),
             fill = pop_color, alpha = 0.3, width = 0.6) +

    ## Error bars (SEM)
    geom_errorbar(data = mean_df,
                  aes(x = drug, ymin = mean_MFI - sem, ymax = mean_MFI + sem),
                  width = 0.2, linewidth = 0.5, colour = "grey30") +

    ## Connected donor lines
    geom_line(data = df,
              aes(x = drug, y = MFI, group = donor, colour = donor),
              alpha = 0.5, linewidth = 0.4) +

    ## Donor dots
    geom_point(data = df,
               aes(x = drug, y = MFI, colour = donor),
               size = 2.5, alpha = 0.85) +

    scale_colour_manual(name = "Donor", values = DONOR_COLORS) +

    labs(x = NULL, y = "Puromycin MFI",
         title = pop_label) +

    theme_classic(base_size = 12) +
    theme(
      plot.title          = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x         = element_text(size = 11, colour = "black", angle = 30, hjust = 1),
      panel.grid.major.y  = element_line(colour = "grey92", linewidth = 0.4),
      legend.position     = "right"
    )
}

## CD36+ barplot
p_mfi_pos <- plot_mfi_paired(
  df = mfi_long |> filter(population == "CD36+"),
  pop_label = "CD36+ Raw MFI",
  pop_color = COLOR_POS
)

ggsave("Fig5B_raw_MFI_CD36pos.pdf", plot = p_mfi_pos, width = 7, height = 5)

## CD36- barplot
p_mfi_neg <- plot_mfi_paired(
  df = mfi_long |> filter(population == "CD36-"),
  pop_label = "CD36- Raw MFI",
  pop_color = COLOR_NEG
)

ggsave("Fig5B_raw_MFI_CD36neg.pdf", plot = p_mfi_neg, width = 7, height = 5)
