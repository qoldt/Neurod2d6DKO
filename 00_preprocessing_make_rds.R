# =============================================================================
# PRE-PROCESSING - how scNeuroD_DKO_2025_doublet_Removed_azimuth.rds was made
#
# This is the provenance for the object both analyses read. It is NOT part of
# either notebook and is NOT re-run as part of them: the .rds already exists and
# is treated as read-only input.
#
# Extracted verbatim from the `seurat_standard` chunk of
# scNeuroD_DKO_2026_revision.Rmd (the file that wrote this .rds). Kept as a
# separate script because it needs things the analysis does not:
#
#   * the Cell Ranger outputs (a network volume; see NEUROD_CELLRANGER_DIR)
#   * packages that are NOT installed here and are not needed downstream:
#       DoubletFinder  (github: chris-mcginnis-ucsf/DoubletFinder)
#       Azimuth        (github: satijalab/azimuth - not on CRAN/Bioconductor)
#       celldex, glmGamPoi
#
# What it does: reads the WT and NeuroD2/NeuroD6 DKO filtered matrices, merges
# them with a Genotype label, computes mitochondrial/ribosomal QC, filters to
# nFeature_RNA > 400 & percent.mt < 5, runs the standard Seurat workflow, removes
# doublets per genotype with DoubletFinder, re-runs the workflow on the cleaned
# object (this is the stored res-0.5 clustering the analyses use), annotates with
# Azimuth against mousecortexref, and saves the .rds.
#
# TO RE-RUN: delete or move the existing .rds first. The guard below refuses to
# overwrite it, so an accidental source() cannot destroy the analysis input.
# =============================================================================

# ---------------------------------------------------------------------------
# Paths. PROJ_DIR is the repository root - this script lives there - so nothing
# below is an absolute path. Override with NEUROD_PROJ_DIR.
#
# The Cell Ranger run directory is machine-specific and has no sensible default,
# so it must be supplied as NEUROD_CELLRANGER_DIR: the folder holding the
# per-sample output directories A4794_SP278_NeuroD2_D6DKO/ and A4794_SP278_wt/.
# ---------------------------------------------------------------------------
PROJ_DIR <- Sys.getenv("NEUROD_PROJ_DIR", unset = normalizePath(getwd()))
OUT_RDS  <- file.path(PROJ_DIR, "scNeuroD_DKO_2025_doublet_Removed_azimuth.rds")

if (file.exists(OUT_RDS)) {
  stop("Refusing to run - the analysis object already exists:\n  ", OUT_RDS,
       "\n  Both notebooks read this file. Move it aside first if you really ",
       "want to regenerate it.", call. = FALSE)
}

cellranger_dir <- Sys.getenv("NEUROD_CELLRANGER_DIR")
if (!nzchar(cellranger_dir)) {
  stop("Set NEUROD_CELLRANGER_DIR to the directory containing the Cell Ranger ",
       "outputs A4794_SP278_NeuroD2_D6DKO/ and A4794_SP278_wt/.", call. = FALSE)
}

library(Seurat)
library(Azimuth)
library(celldex)
library(glmGamPoi)
library(ggplot2)
library(cowplot)
library(future)
library(DoubletFinder)  
library(dplyr)          


plotdir <- file.path(PROJ_DIR, "plots", "preprocessing")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)



# Load data
# Each path points to a Cell Ranger sample directory; the FILTERED matrix is
# what is read (see NEUROD_CELLRANGER_DIR above, and the layout in the README).
# This Data is from E15 cortex


NeuroD2_D6DKO <- Read10X(data.dir = file.path(cellranger_dir,
                                              "A4794_SP278_NeuroD2_D6DKO/outs",
                                              "filtered_feature_bc_matrix"))

WT <- Read10X(data.dir = file.path(cellranger_dir,
                                   "A4794_SP278_wt/outs",
                                   "filtered_feature_bc_matrix"))


# Prefix cell barcodes to make them unique
colnames(NeuroD2_D6DKO) <- paste0("NeuroD2_D6DKO_", colnames(NeuroD2_D6DKO))
colnames(WT) <- paste0("WT_", colnames(WT))

# Combine matrices (rows are genes, columns are cells)

#combined_counts <- cbind(NeuroD1KO, NeuroD2_D6DKO, NeuroDTKO, WT)
combined_counts <- cbind(NeuroD2_D6DKO, WT)

# Add metadata
meta_data <- data.frame(
  Genotype = factor(
    c(rep("NeuroD2_D6DKO", ncol(NeuroD2_D6DKO)),
      rep("WT", ncol(WT))),
    levels = c("WT", "NeuroD2_D6DKO")
  ),
  row.names = colnames(combined_counts)
)

# Create a Seurat object

options(Seurat.object.assay.calcn = TRUE)
sc_seurat <- CreateSeuratObject(
  counts = combined_counts,
  min.cells = 3,
  min.features = 200,
  meta.data = meta_data
)

#Remove sparse matrices to save RAM
rm(NeuroD2_D6DKO,  WT)

saveRDS(sc_seurat, file = file.path(PROJ_DIR, "seurat_pre_filter.rds"))
# Calculate QC metrics: percent of mitochondrial and ribosomal genes

sc_seurat[["percent.mt"]] <- PercentageFeatureSet(sc_seurat, pattern = "^mt-")
sc_seurat[["percent.rbp"]] <- PercentageFeatureSet(sc_seurat, pattern = "^Rp[sl]")

# Pre-filter QC plots

VlnPlot(sc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4)
ggsave(filename = file.path(plotdir, "QC_combined_pre-filter.pdf"), width = 20, height = 18, units = "cm")

# Filter cells
# Adjust thresholds as appropriate 
sc_seurat <- subset(sc_seurat, subset = nFeature_RNA > 400 & percent.mt < 5)

# Post-filter QC plots
VlnPlot(sc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4)
ggsave(filename = file.path(plotdir, "QC_combined_post-filter.pdf"), width = 20, height = 18, units = "cm")

# allow up to 40 GB
options(future.globals.maxSize = 40 * 1024^3)
plan("multicore")

# Standard workflow: normalization, feature selection, scaling, PCA, UMAP
sc_seurat <- sc_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

DimPlot(sc_seurat, split.by = "Genotype")
ggsave(filename = file.path(plotdir, "UMAP_post-filter_pre_doublet_removal.pdf"), width = 30, height = 15, units = "cm")

# Doublet detection with DoubletFinder
# Split by genotype, run DoubletFinder on each subset, then merge back

seurat_list <- SplitObject(sc_seurat, split.by = "Genotype")

doublet_dir <- file.path(plotdir, "DoubletFinder")
dir.create(doublet_dir, showWarnings = FALSE)

doublet_metrics <- list()

for (i in seq_along(seurat_list)) {
  sample_name <- names(seurat_list)[i]
  message("Processing sample: ", sample_name)
  
  DefaultAssay(seurat_list[[i]]) <- "RNA"
  
  # Find neighbors/clusters
  seurat_list[[i]] <- FindNeighbors(seurat_list[[i]], dims = 1:30)
  seurat_list[[i]] <- FindClusters(seurat_list[[i]], resolution = 0.4)
  
  # Parameter sweep
  sweep.res <- paramSweep(seurat_list[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK
  optimal.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  
  # Homotypic doublet proportion
  annotations <- seurat_list[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # Number of expected doublets (default 5%; adjust as needed)
  nExp_poi <- round(ncol(seurat_list[[i]]) * 0.05)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_list[[i]] <- doubletFinder(
    seurat_list[[i]],
    PCs = 1:30,
    pN = 0.25,
    pK = optimal.pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  # Identify the classification column from DoubletFinder
  classification_col <- grep("^DF.classifications_", 
                             colnames(seurat_list[[i]]@meta.data), 
                             value = TRUE)
  
  # Visualize doublets vs. singlets
  DimPlot(seurat_list[[i]], reduction = "umap", group.by = classification_col) +
    ggtitle(paste0("DoubletFinder Results: ", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = file.path(doublet_dir, paste0(sample_name, "_DoubletFinder_UMAP.pdf")),
         width = 10, height = 8)
  
  # Collect metrics
  total_cells <- ncol(seurat_list[[i]])
  doublet_count <- sum(seurat_list[[i]]@meta.data[[classification_col]] == "Doublet")
  singlet_count <- sum(seurat_list[[i]]@meta.data[[classification_col]] == "Singlet")
  
  doublet_metrics[[sample_name]] <- data.frame(
    Sample = sample_name,
    Total_Cells = total_cells,
    Expected_Doublets = nExp_poi,
    Adjusted_Expected_Doublets = nExp_poi.adj,
    Doublet_Count = doublet_count,
    Singlet_Count = singlet_count
  )
  
  # Keep only singlets
  singlet_cells <- rownames(
    seurat_list[[i]]@meta.data[
      seurat_list[[i]]@meta.data[[classification_col]] == "Singlet", ]
  )
  seurat_list[[i]] <- subset(seurat_list[[i]], cells = singlet_cells)
  
  # Optionally remove other assays to reduce object size
  seurat_list[[i]] <- DietSeurat(seurat_list[[i]], assays = "RNA")
}

# -------------------------------------------------------------------------
# Merge cleaned samples back together
# -------------------------------------------------------------------------
sc_seurat <- merge(seurat_list[[1]], y = seurat_list[-1])
sc_seurat <- JoinLayers(sc_seurat)  # if needed for Seurat layering

# Combine doublet metrics into a single data frame
metrics_df <- do.call(rbind, doublet_metrics)

# Save metrics to CSV
write.csv(metrics_df,
          file = file.path(doublet_dir, "DoubletFinder_Metrics_All_Samples.csv"),
          row.names = FALSE)

# Re-run standard workflow after doublet removal

sc_seurat <- sc_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5)

DimPlot(sc_seurat, split.by = "Genotype")
ggsave(filename = file.path(plotdir, "UMAP_post-filter_post_doublet_removal.pdf"), width = 30, height = 15, units = "cm")

# Run Azimuth
DefaultAssay(sc_seurat) <- "RNA"
sc_seurat <- NormalizeData(sc_seurat) 
sc_seurat <- Azimuth::RunAzimuth(sc_seurat, reference = "mousecortexref")

saveRDS(sc_seurat, file = OUT_RDS)
