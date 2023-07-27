rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnsDb.Mmusculus.v79)

# Params
project.name <- "2022_03_novaseq_PM22082"
sample.name <- "211013_17"

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
project_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/", sample.name)
rctd_path <- paste0(project_output_dir, "/", sample.name, "_rctd.rds")
seurat_path <- paste0(project_output_dir, "/", sample.name, "_seurat.rds")
dir.create(project_output_dir)

# Read RCTD results
rctd <- readRDS(rctd_path)
results <- rctd@results
spatialRNA <- rctd@spatialRNA

sk <- CreateSeuratObject(counts = spatialRNA@counts, project = "space_kidney", min.cells = 3, min.features = 100)
sk <- NormalizeData(sk, normalization.method = "LogNormalize", scale.factor = 10000)
sk <- FindVariableFeatures(sk, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sk), 10)
all.genes <- rownames(sk)
sk <- ScaleData(sk, features = all.genes)
sk <- RunPCA(sk, features = VariableFeatures(object = sk))
sk <- RunUMAP(sk, dims = 1:10)

# Filter labels on cells kept after seurat processing
matched_barcodes <- match(colnames(spatialRNA@counts), colnames(sk@assays$RNA@counts))
drop <- which(is.na(matched_barcodes))
spatialRNA_filt <- spatialRNA@counts[,-drop]
cell_labels <- results$results_df$first_type[-drop]
sk[["rctd_label"]] <- cell_labels

#cell_labels <- results$results_df$first_type
#sk[["rctd_label"]] <- cell_labels

# Set idents from meta
Idents(sk) <- 'rctd_label'

pdf(paste0(project_output_dir, "/count_qc_1.pdf"))
plot <- VlnPlot(sk, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
print(plot)
dev.off()

pdf(paste0(project_output_dir, "/count_qc_2.pdf"))
plot <- FeatureScatter(sk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot)
dev.off()

# plot variable features with and without labels
pdf(paste0(project_output_dir, "/high_var_genes.pdf"))
plot1 <- VariableFeaturePlot(sk)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot = plot1 + plot2
print(combined_plot)
dev.off()

pdf(paste0(project_output_dir, "/pca.pdf"))
plot <- DimPlot(sk, reduction = "pca")
print(plot)
dev.off()

pdf(paste0(project_output_dir, "/umap.pdf"))
#plot <- DimPlot(sk, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot <- DimPlot(sk, reduction = "umap")
print(plot)
dev.off()

sk.markers <- FindAllMarkers(sk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find DE markers for PCT
#cluster.markers.PCT <- FindMarkers(sk, ident.1 = "PCT", min.pct = 0.25)
#rownames(cluster.markers.PCT) <- gene_symbols$SYMBOL
#gene_symbols <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(cluster.markers.PCT), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
#head(cluster.markers.PCT, n = 20)

# Find DE markers for PST
#cluster.markers.PST <- FindMarkers(sk, ident.1 = "PST", min.pct = 0.25)
#gene_symbols <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(cluster.markers.PST), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
#rownames(cluster.markers.PST) <- gene_symbols$SYMBOL
#head(cluster.markers.PST, n = 20)

# Feature plot target genes
target_genes <- c("Lrp2", "Aqp11", "Slc34a1", "Hnf4a", "Cldn2", "Fbp1")
target_genes_ids <- ensembldb::select(EnsDb.Mmusculus.v79, keys=target_genes, keytype="SYMBOL", columns=c("SYMBOL","GENEID"))

pdf(paste0(project_output_dir, "/features.pdf"))
plot <- FeaturePlot(sk, features = target_genes_ids$GENEID)
print(plot)
dev.off()

#target_genes <- c("Scx", "Slc04a1", "Avpr2", "Hk1", "Egr1", "Fos")
#target_genes_ids <- ensembldb::select(EnsDb.Mmusculus.v79, keys=target_genes, keytype="SYMBOL", columns=c("SYMBOL","GENEID"))

#pdf(paste0(project_output_dir, "/features_neg.pdf"))
#plot <- FeaturePlot(sk, features = target_genes_ids$GENEID)
#print(plot)
#dev.off()

sk.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf(paste0(project_output_dir, "/heatmap.pdf"))
plot <- DoHeatmap(sk, features = top10$gene) + NoLegend()
print(plot)
dev.off()




# Save seurat object
saveRDS(sk, file = seurat_path)