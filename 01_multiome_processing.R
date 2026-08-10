## =============================================================================
## Multiome (snRNA-seq + snATAC-seq) Preprocessing, Clustering, and
## Differential Analysis
## =============================================================================
## Associated manuscript: Distinct chromatin accessibility programs are linked to 
## metabolic regulation in human lactocytes
##
## Packages and versions used:
##   Seurat v5.4.0, Signac v1.6, scCustomize v3.3.0, scDblFinder v1.14,
##   harmony v1.2.4, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38,
##   clustree v0.5.1, ChIPseeker v1.47.1, clusterProfiler v4.8.14,
##   org.Hs.eg.db, MACS2 v2.2.9.1 (called via Signac)
##   R v4.5.1
##
## Full session information is printed at the end of this script.
## =============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(clustree)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(scCustomize)
library(scDblFinder)
library(harmony)
library(pheatmap)

set.seed(1234)
setwd("")

## =============================================================================
## 1. Generate unified ATAC peak set across all samples
## =============================================================================

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"

# Read Cell Ranger ARC peak calls per sample
sample_dirs <- c("FC014-4", "FC018-2", "FC003-5", "FC004-3", "FC006-1",
                 "FC007", "FC008-2", "FC009-1", "FC011-1", "FC012-2", "FC013-1")

gr_list <- lapply(sample_dirs, function(d) {
  peaks <- read.table(file.path(d, "atac_peaks.bed"), col.names = c("chr", "start", "end"))
  makeGRangesFromDataFrame(peaks)
})

# Create unified peak set, retaining peaks between 20-10,000 bp
peaks.uni <- reduce(x = do.call(c, gr_list))
peaksUse <- peaks.uni[width(peaks.uni) < 10000 & width(peaks.uni) > 20]

## =============================================================================
## 2. Create Seurat objects with RNA (CellBender-corrected) and ATAC assays
## =============================================================================

clean_data_path <- "clean_data/"
h5_files <- list.files(clean_data_path, pattern = "_out_filtered.h5$", full.names = TRUE)
h5_files <- h5_files[1:11]

for (h5_file in h5_files) {
  donor_id <- gsub("_out_filtered.h5", "", basename(h5_file))

  # Read CellBender-corrected counts
  seurat_data <- Read_CellBender_h5_Mat(file_name = h5_file)
  rna_features <- rownames(seurat_data)[!grepl("chr", rownames(seurat_data))]
  rna_counts <- seurat_data[rna_features, ]

  seurat_obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA",
                                   min.cells = 10, names.field = 1, names.delim = "_")

  # Add ATAC assay from fragment files
  frag_file <- file.path(donor_id, "atac_fragments.tsv.gz")
  if (file.exists(frag_file)) {
    frags <- CreateFragmentObject(frag_file)
    fragcounts <- FeatureMatrix(fragments = frags, features = peaksUse, cells = colnames(seurat_obj))
    seurat_obj[["ATAC"]] <- CreateChromatinAssay(fragcounts, fragments = frags,
                                                  sep = c(":", "-"), annotation = annotation)
  }

  # Add per-barcode metrics from Cell Ranger ARC
  metadata_file <- file.path(donor_id, "per_barcode_metrics.csv")
  if (file.exists(metadata_file)) {
    metadata <- read.csv(metadata_file, row.names = 1)
    common_cells <- intersect(rownames(metadata), colnames(seurat_obj))
    seurat_obj <- subset(seurat_obj, cells = common_cells)
    metadata <- metadata[common_cells, , drop = FALSE]
    seurat_obj <- AddMetaData(seurat_obj, metadata)
  }

  saveRDS(seurat_obj, file.path(donor_id, paste0(donor_id, "_obj_cellbender.rds")))
}

## =============================================================================
## 3. Quality control and doublet detection
## =============================================================================

# Load individual sample objects
sample_ids <- c("FC014-4", "FC018-2", "FC003-5", "FC004-3", "FC006-1",
                "FC007", "FC008-2", "FC009-1", "FC011-1", "FC012-2", "FC013-1")

seurat_list <- lapply(sample_ids, function(sid) {
  readRDS(file.path(sid, paste0(sid, "_obj_cellbender.rds")))
})

# Add QC metrics
seurat_list <- lapply(seurat_list, function(obj) {
  Add_Cell_QC_Metrics(obj, species = "human")
})

# Add sample metadata
# Note: Sample FC018-2 was originally labeled FC019-1 in sequencing;
# corrected here to reflect the true donor ID (FC018)
donor_labels  <- c("FC014", "FC018", "FC003", "FC004", "FC006",
                   "FC007", "FC008", "FC009", "FC011", "FC012", "FC013")
batch_labels  <- c("Batch_3", "Batch_3", "Batch_1", "Batch_1", "Batch_1",
                   "Batch_1", "Batch_2", "Batch_2", "Batch_2", "Batch_3", "Batch_3")
days_pp       <- c(228, 135, 450, 151, 229, 257, 218, 81, 681, 265, 142)

for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$orig.ident <- donor_labels[i]
  seurat_list[[i]]$batch      <- batch_labels[i]
  seurat_list[[i]]$days_PP    <- days_pp[i]
  seurat_list[[i]]$milk_stage <- "Mature"
}

# Compute ATAC QC metrics and apply initial filtering
seurat_list <- lapply(seurat_list, function(obj) {
  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  obj$pct_reads_in_peaks <- obj$atac_peak_region_fragments / obj$atac_fragments * 100

  overlaps <- findOverlaps(query = obj[["ATAC"]], subject = blacklist_hg38)
  hit.regions <- queryHits(x = overlaps)
  data.matrix <- GetAssayData(object = obj, assay = "ATAC", slot = "counts")[hit.regions, , drop = FALSE]
  obj$blacklist_region_fragments <- colSums(data.matrix)
  obj$blacklist_fraction <- obj$blacklist_region_fragments / obj$atac_peak_region_fragments

  subset(obj, subset = nCount_RNA > 300 & nFeature_RNA > 200 &
           percent_mito < 35 & nCount_ATAC > 300)
})

# Doublet detection with scDblFinder (dbr = 0.08)
# RNA-based doublet detection
for (i in seq_along(seurat_list)) {
  rna_counts <- GetAssayData(seurat_list[[i]], assay = "RNA", slot = "counts")
  res <- scDblFinder(rna_counts, dbr = 0.08, clusters = TRUE)
  metadata_df <- data.frame(
    scDblFinder_score_rna = res$scDblFinder.score,
    scDblFinder_class_rna = res$scDblFinder.class,
    row.names = colnames(rna_counts)
  )
  seurat_list[[i]] <- AddMetaData(seurat_list[[i]], metadata = metadata_df)
}

# ATAC-based doublet detection
for (i in seq_along(seurat_list)) {
  atac_counts <- GetAssayData(seurat_list[[i]], assay = "ATAC", slot = "counts")
  res <- scDblFinder(atac_counts, dbr = 0.08, artificialDoublets = 1,
                     aggregateFeatures = TRUE, nfeatures = 25, processing = "normFeatures")
  metadata_df <- data.frame(
    scDblFinder_score_atac = res$scDblFinder.score,
    scDblFinder_class_atac = res$scDblFinder.class,
    row.names = colnames(atac_counts)
  )
  seurat_list[[i]] <- AddMetaData(seurat_list[[i]], metadata = metadata_df)
}

# Summarize doublet rates
doublet_summary <- do.call(rbind, lapply(seq_along(seurat_list), function(i) {
  md <- seurat_list[[i]]@meta.data
  total <- nrow(md)
  data.frame(
    Sample = donor_labels[i],
    Total_Cells = total,
    RNA_Doublets = sum(md$scDblFinder_class_rna == "doublet", na.rm = TRUE),
    RNA_Doublet_Pct = round(sum(md$scDblFinder_class_rna == "doublet", na.rm = TRUE) / total * 100, 2),
    ATAC_Doublets = sum(md$scDblFinder_class_atac == "doublet", na.rm = TRUE),
    ATAC_Doublet_Pct = round(sum(md$scDblFinder_class_atac == "doublet", na.rm = TRUE) / total * 100, 2),
    Both_Doublets = sum(md$scDblFinder_class_rna == "doublet" & md$scDblFinder_class_atac == "doublet", na.rm = TRUE)
  )
}))
write.csv(doublet_summary, "doublet_summary.csv", row.names = FALSE)

## =============================================================================
## 4. Merge samples, remove doublets, and apply stringent QC filters
## =============================================================================

Merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:11],
  add.cell.ids = c("FC014_4", "FC018_2", "FC003_5", "FC004_3",
                   "FC006_1", "FC007", "FC008_2", "FC009_1",
                   "FC011_1", "FC012_2", "FC013_1")
)

# Apply stringent per-cell QC filters (as described in Methods)
Filtered <- subset(x = Merged,
                   subset = nCount_RNA > 450 &
                     nFeature_RNA > 200 &
                     percent_mito < 35 &
                     percent_ribo < 15 &
                     nCount_ATAC > 450 &
                     nucleosome_signal < 2.5 &
                     pct_reads_in_peaks > 15 &
                     TSS.enrichment > 2)

saveRDS(Filtered, "obj_clean_filtered.rds")

## =============================================================================
## 5. RNA processing: SCTransform, PCA, Harmony integration, clustering
## =============================================================================

Filtered <- readRDS("obj_clean_filtered.rds")
DefaultAssay(Filtered) <- "RNA"

# SCTransform normalization with mitochondrial percentage regressed out
Filtered <- SCTransform(Filtered, vars.to.regress = c("percent_mito"),
                        min_cells = 0, return.only.var.genes = FALSE, verbose = TRUE)

DefaultAssay(Filtered) <- "SCT"
Filtered <- RunPCA(Filtered)

# Harmony batch correction on PCA embeddings
Filtered <- IntegrateLayers(
  object = Filtered, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony_rna",
  reduction.name = "rna.harmony",
  normalization.method = "SCT",
  verbose = TRUE
)

# RNA UMAP and clustering (12 PCs, resolution 0.25, SLM algorithm)
Filtered <- RunUMAP(Filtered, dims = 1:12, reduction = "harmony_rna",
                    reduction.name = "umap.rna.integrated")
Filtered <- FindNeighbors(Filtered, reduction = "harmony_rna", dims = 1:12)
Filtered <- FindClusters(Filtered, resolution = 0.25, algorithm = 3, verbose = FALSE)

# Assign broad cell type annotations based on marker gene expression
Idents(Filtered) <- Filtered$SCT_snn_res.0.25
new.cluster.ids <- c("Epithelial cells", "Epithelial cells", "Macrophage Cluster 1",
                     "Epithelial cells", "Epithelial cells", "Epithelial cells",
                     "Macrophage Cluster 2", "Epithelial cells", "Epithelial cells",
                     "T-cells", "Epithelial cells", "Epithelial cells")
names(new.cluster.ids) <- levels(Filtered)
Filtered <- RenameIdents(Filtered, new.cluster.ids)
Filtered$General_celltypes <- Idents(Filtered)

## =============================================================================
## 6. Epithelial (lactocyte) subclustering
## =============================================================================

Idents(Filtered) <- Filtered$General_celltypes
Epi_subset <- subset(Filtered, idents = "Epithelial cells")

# Re-normalize and subcluster epithelial cells
DefaultAssay(Epi_subset) <- "RNA"
Epi_subset <- SCTransform(Epi_subset, vars.to.regress = c("percent_mito"),
                          min_cells = 0, return.only.var.genes = FALSE, verbose = TRUE)
DefaultAssay(Epi_subset) <- "SCT"
Epi_subset <- RunPCA(Epi_subset)
Epi_subset <- IntegrateLayers(
  object = Epi_subset, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony_rna",
  reduction.name = "rna.harmony",
  normalization.method = "SCT",
  verbose = TRUE
)
Epi_subset <- RunUMAP(Epi_subset, dims = 1:12, reduction = "harmony_rna",
                      reduction.name = "umap.rna.integrated")
Epi_subset <- FindNeighbors(Epi_subset, reduction = "harmony_rna", dims = 1:12)
Epi_subset <- FindClusters(Epi_subset, resolution = 0.25, algorithm = 3, verbose = FALSE)

# Assign lactocyte subcluster identities (update labels as appropriate)
# Epi_subset$Celltype <- paste0("Lactocyte ", Epi_subset$seurat_clusters)
# -- or assign manually based on marker expression: --
# Idents(Epi_subset) <- Epi_subset$SCT_snn_res.0.25
# new.epi.ids <- c("Lactocyte 0", "Lactocyte 1", "Lactocyte 2",
#                  "Lactocyte 3", "Lactocyte 4")
# names(new.epi.ids) <- levels(Epi_subset)
# Epi_subset <- RenameIdents(Epi_subset, new.epi.ids)
# Epi_subset$Celltype <- Idents(Epi_subset)

# Merge granular annotations back into full dataset
Filtered$Granular_celltype <- as.character(Filtered$General_celltypes)
Filtered$Granular_celltype[Cells(Epi_subset)] <- as.character(Epi_subset$Celltype)
Idents(Filtered) <- "Granular_celltype"

## =============================================================================
## 7. ATAC processing: MACS2 peak calling, LSI, Harmony, clustering
## =============================================================================

# Call peaks per cell type using MACS2
DefaultAssay(Filtered) <- "ATAC"
# Sys.which("macs2") ## update macs path
peaks_all <- CallPeaks(Filtered, group.by = "Granular_celltype",
                       macs2.path = "")
peaks_all <- keepStandardChromosomes(peaks_all, pruning.mode = "coarse")
peaks_all <- subsetByOverlaps(x = peaks_all, ranges = blacklist_hg38, invert = TRUE)

# Create new ChromatinAssay from MACS2-called peaks
macs2_counts_all <- FeatureMatrix(fragments = Fragments(Filtered),
                                  features = peaks_all, cells = colnames(Filtered))
Filtered[["peaksbyc"]] <- CreateChromatinAssay(counts = macs2_counts_all,
                                                fragments = Fragments(Filtered),
                                                annotation = annotation)

# Latent semantic indexing (LSI) normalization
DefaultAssay(Filtered) <- "peaksbyc"
Filtered <- RunTFIDF(Filtered)
Filtered <- FindTopFeatures(Filtered, min.cutoff = "q0")
Filtered <- RunSVD(Filtered)

# Harmony batch correction on LSI embeddings (components 2-10, excluding component 1)
Filtered <- RunHarmony(Filtered, group.by.vars = "orig.ident",
                       reduction.use = "lsi", reduction.save = "harmony_atac",
                       project.dim = FALSE)

Filtered <- RunUMAP(Filtered, reduction = "harmony_atac", dims = 2:10,
                    reduction.name = "umap.atac.integrated")
Filtered <- FindNeighbors(Filtered, reduction = "harmony_atac", dims = 2:10)
Filtered <- FindClusters(Filtered, graph = "peaksbyc_snn", resolution = 0.25,
                         algorithm = 3, verbose = FALSE)

## =============================================================================
## 8. Weighted nearest neighbor (WNN) integration and joint UMAP
## =============================================================================

Filtered <- FindMultiModalNeighbors(Filtered,
                                    reduction.list = list("harmony_rna", "harmony_atac"),
                                    dims.list = list(1:12, 2:10))
Filtered <- RunUMAP(Filtered, nn.name = "weighted.nn",
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

DimPlot(Filtered, reduction = "wnn.umap", group.by = "Granular_celltype")

saveRDS(Filtered, "WNN_SCT_cellbender.rds")

## =============================================================================
## 9. Peak-to-gene linkage analysis (Spearman correlation)
## =============================================================================

DefaultAssay(Filtered) <- "peaksbyc"

# Link peaks to gene expression using Spearman correlation
# Retains only positive links with adjusted p-value < 0.05 (BH correction)
Filtered <- LinkPeaks(
  object = Filtered,
  peak.assay = "peaksbyc",
  expression.assay = "SCT"
)

## =============================================================================
## 10. Differentially expressed genes (DEGs) across lactocyte subclusters
## =============================================================================

# Subset to lactocytes only
Idents(Filtered) <- "Granular_celltype"
lactocyte_ids <- grep("^Lactocyte", levels(Idents(Filtered)), value = TRUE)
Lactocytes <- subset(Filtered, idents = lactocyte_ids)

# DEGs: Wilcoxon rank-sum test, min.pct = 0.05, adjusted p-value < 0.05
DefaultAssay(Lactocytes) <- "SCT"
Idents(Lactocytes) <- "Granular_celltype"

DEGs_lactocytes <- FindAllMarkers(
  Lactocytes,
  assay = "SCT",
  only.pos = FALSE,
  min.pct = 0.05,
  test.use = "wilcox",
  verbose = TRUE
)
DEGs_lactocytes_sig <- DEGs_lactocytes %>% filter(p_val_adj < 0.05)
write.csv(DEGs_lactocytes_sig, "DEGs_lactocyte_subclusters.csv", row.names = FALSE)

## =============================================================================
## 11. Differentially accessible peaks (DAPs) across lactocyte subclusters
## =============================================================================

# DAPs: Likelihood-ratio test with total peak count as latent variable,
# min.pct = 0.1, adjusted p-value < 0.05
DefaultAssay(Lactocytes) <- "peaksbyc"

DAPs_lactocytes <- FindAllMarkers(
  Lactocytes,
  assay = "peaksbyc",
  only.pos = FALSE,
  min.pct = 0.1,
  test.use = "LR",
  latent.vars = "nCount_peaksbyc",
  verbose = TRUE
)
DAPs_lactocytes_sig <- DAPs_lactocytes %>% filter(p_val_adj < 0.05)
write.csv(DAPs_lactocytes_sig, "DAPs_lactocyte_subclusters.csv", row.names = FALSE)

## =============================================================================
## Session information
## =============================================================================

sessionInfo()
