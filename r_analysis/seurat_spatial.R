rm(list = ls())
system("type R")

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)

# Params
project.name <- "2022_03_novaseq_PM22082"
sample.name <- "211125_20"
#sample.name <- "211013_03"

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
matrix_dir <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge")
project_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/", sample.name)
rctd_path <- paste0(project_output_dir, "/", sample.name, "_rctd.rds")
coord.path <- paste0(project_dir, "/", sample.name, ".csv")
dir.create(project_output_dir)

# Read RCTD
rctd <- readRDS(rctd_path)
results <- rctd@results
spatialRNA <- rctd@spatialRNA

# Load expression
exp_matrix <- Read10X(data.dir = matrix_dir)
slide.seq = CreateSeuratObject(counts = exp_matrix, project = sample.name, assay="Spatial")
exp_ids = unlist(slide.seq@assays$Spatial@counts@Dimnames[2])

# Load bead
bead <- read.csv(file = coord.path)
rownames(x = bead) <- bead$SeqBarcode
bead <- bead[, 3:4]
bead_ids <- rownames(x = bead)

# Check for mismatch and filter out of data table if needed
bead_mistmatch <- bead_ids[!(bead_ids %in% exp_ids)]
if( length(bead_mistmatch) > 0) {
  bead <- bead[!(row.names(bead) %in% bead_mistmatch),]
  #slide.seq <- slide.seq[,!colnames(slide.seq) %in% bead_mistmatch]
}

# Integrate data
slide.seq[['image']] <- new(Class = 'SlideSeq', assay = "Spatial", coordinates = bead)
slide.seq <- AddMetaData(slide.seq, metadata = results$results_df)
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
slide.seq <- FilterSlideSeq(object = slide.seq, radius = 2450, do.plot = FALSE)

# Initial puck plots
# https://satijalab.org/seurat/reference/spatialplot
png(
  file = paste0(project_output_dir, "/", "seurat_spatial", "_ncount.png"),
  width = 1800,
  height = 1200,
  res = 200
)
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial", crop = TRUE) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()

# Pre-process
#slide.seq <- NormalizeData(slide.seq, normalization.method = "LogNormalize", scale.factor = 10000)
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
var_features <- VariableFeatures(slide.seq)

# Spatial dim plots
png(
  file = paste0(project_output_dir, "/", "seurat_spatial", "_spatialdim.png"),
  width = 1800,
  height = 1200,
  res = 200
)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE, group.by = "first_type", label.size = 2, repel = TRUE, label.box = TRUE, label.color = "white")
plot2 <- SpatialDimPlot(slide.seq, stroke = 0, group.by = "first_type")
plot1 + plot2
dev.off()

# Spatial dim plots with rctd
png(
  file = paste0(project_output_dir, "/", "seurat_spatial", "_spatialdim_rctd.png"),
  width = 1800,
  height = 1200,
  res = 200
)
p1 <- SpatialDimPlot(slide.seq, group.by = "first_type")
#p2 <- SpatialDimPlot(slide.seq, group.by = "second_type")
#p1 | p2
p1
dev.off()

# Find spatially variable features
DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(
  slide.seq, 
  assay = "SCT", 
  slot = "scale.data", 
  features = var_features[1:1000],
  selection.method = "moransi",
  x.cuts = 100,
  y.cuts = 100
)

# Plot top spatially var features
png(
  file = paste0(project_output_dir, "/", "seurat_spatial", "_spatialvar_top6.png"),
  width = 1800,
  height = 1200,
  res = 200
)
SpatialFeaturePlot(slide.seq, features = 
  head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"), 6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
dev.off()

png(
  file = paste0(project_output_dir, "/", "seurat_spatial", "_spatialvar_targets.png"),
  width = 1800,
  height = 1200,
  res = 200
)
SpatialFeaturePlot(slide.seq, features = c("Igkc", "Ighm", "Igha", "Sgk1", "Alas1", "Angpt2", "Krt18", "Cry1"), ncol = 4, alpha = c(0.1, 1), max.cutoff = "q95")
dev.off()

