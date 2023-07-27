rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(plyr)
library(dplyr)
library(EnsDb.Mmusculus.v79)

# Params
project.name <- "2022_03_novaseq_PM22082"
sample.name <- "211013_17"

# RCTD Params
rctd_gene_cutoff <- 0.0001
rctd_gene_cutoff_reg <- 0.0002
rctd_fc_cutoff <- 0.5
rctd_fc_cutoff_reg <- 0.75
rctd_UMI_min <- 100
rctd_UMI_max <- 2e+07

rctd_cell_min_distance <- 25 # minimum number of cells required per cell type. Default 25, can be lowered if desired.
rctd_confidence_thresh <- 3 # the minimum change in likelihood (compared to other cell types) necessary to determine a cell type identity with confidence
rctd_doublet_thresh <- 20 #  the penalty weight of predicting a doublet instead of a singlet for a pixel

rctd_N_epoch <- 8
rctd_N_epoch_bulk <- 30
rctd_min_change_bulk <- .0001 
rctd_N_X <- 50000
rctd_N_fit <- 500
rctd_K_val <- 100

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
project_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/", sample.name)
ref_matrix.path <- paste0(data_dir, "/GSE157079_P0_adult_counts.rds")
ref_clusters.path <- paste0(data_dir, "/GSE157079_P0_adult_clusters.txt")
barcodes.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/barcodes.tsv.gz")
features.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/features.tsv.gz")
matrix.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/matrix.mtx.gz")
coord.path <- paste0(project_dir, "/", sample.name, ".csv")
rctd_path <- paste0(project_output_dir, "/", sample.name, "_rctd.rds")
dir.create(project_output_dir)

# Load reference table
ref_counts <- readRDS(ref_matrix.path)
ref_meta <- read.table(ref_clusters.path,sep="\t", header=TRUE)

# Get gene symbols and convert to gene ids and match for found
ref_gene_names <- rownames(ref_counts)
ref_gene_ids <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ref_gene_names, keytype="SYMBOL", columns=c("SYMBOL","GENEID"))
matched_names <- match(ref_gene_names, ref_gene_ids$SYMBOL)

# Drop rows that arnt found
rownames(ref_counts) <- ref_gene_ids$GENEID[matched_names]
drop <- which(is.na(rownames(ref_counts)))
ref_counts <- ref_counts[-drop,]

# Get cell type mappings and create reference object
cell_types <- ref_meta$clusters; names(cell_types) <- ref_meta$barcodes
cell_types <- as.factor(cell_types)
reference <- Reference(ref_counts, cell_types)

# Load count matrix and annotate with gene ids
count_mat <- readMM(file = matrix.path)
features.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcodes.names = read.delim(barcodes.path, header = FALSE, stringsAsFactors = FALSE)
colnames(count_mat) = barcodes.names$V1
rownames(count_mat) = features.names$V1

# Create coordinate and nUMI inputs for puck
coords <- read.csv(coord.path)
coords_filt <- dplyr::select(coords, c(1,3,4))
rownames(coords_filt) <- coords_filt[,1]; coords_filt[,1] <- NULL
nUMI <- colSums(count_mat)

# Create puck and RCTD objects and run RCTD
puck <- SpatialRNA(coords_filt, count_mat, nUMI)

#gene_cutoff, fc_cutoff, gene_cutoff_reg, fc_cutoff_reg: are used for differentially expressed gene selection, 
#with gene_cutoff filtering for average expression and fc_cutoff filtering for log-fold-change across cell types.
#UMI_min, UMI_max: are the minimum and maximum read depth for pixels in the SpatialRNA dataset.

rctd <- create.RCTD(puck, 
                    reference, 
                    max_cores = 8, 
                    gene_cutoff=rctd_gene_cutoff, 
                    gene_cutoff_reg=rctd_gene_cutoff_reg, 
                    fc_cutoff=rctd_fc_cutoff, 
                    fc_cutoff_reg=rctd_fc_cutoff_reg, 
                    UMI_min=rctd_UMI_min, 
                    UMI_max=rctd_UMI_max,
                    CELL_MIN_INSTANCE=rctd_cell_min_distance,
                    CONFIDENCE_THRESHOLD=rctd_confidence_thresh)
rctd <- run.RCTD(rctd, doublet_mode = 'doublet')

# Save rctd object
saveRDS(rctd, file = rctd_path)