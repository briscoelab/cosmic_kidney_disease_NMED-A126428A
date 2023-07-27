rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(dplyr)
library(ggplot2)

# Params
project.name <- "2022_03_novaseq_PM22082"
sample.name <- "211004_12"

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
project_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/", sample.name)
barcodes.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/barcodes.tsv.gz")
features.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/features.tsv.gz")
matrix.path <- paste0(data_dir, "/", project.name, "/", sample.name, "_dge/matrix.mtx.gz")
coord.path <- paste0(project_dir, "/", sample.name, ".csv")
rctd_path <- paste0(project_output_dir, "/", sample.name, "_rctd.rds")
dir.create(project_output_dir)

# Read RCTD results
rctd <- readRDS(rctd_path)
results <- rctd@results
spatialRNA <- rctd@spatialRNA

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

# Create reporting objects and slices
norm_weights = normalize_weights(results$weights) 
cell_type_names <- rctd@cell_type_info$info[[2]] #list of cell type names
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",]
doub_occur <- table(doublets$second_type, doublets$first_type)

# Create puck-based reports

pdf(paste0(project_output_dir, "/nUMI.pdf"))
plot <- hist(log(puck@nUMI,2))
print(plot)
dev.off()

pdf(paste0(project_output_dir, "/puck.pdf"))
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
plot <- plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
print(plot)
dev.off()


# Create RCTD reports
plot <- plot_weights(cell_type_names, spatialRNA, project_output_dir, norm_weights) 
plot <- plot_weights_unthreshold(cell_type_names, spatialRNA, project_output_dir, norm_weights)
plot <- plot_cond_occur(cell_type_names, project_output_dir, norm_weights, spatialRNA)
plot <- plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, project_output_dir) 
plot <- plot_doublets(spatialRNA, doublets, project_output_dir, cell_type_names)
plot <- plot_doublets_type(spatialRNA, doublets, project_output_dir, cell_type_names) 
plot <- plot_doub_occur_stack(doub_occur, project_output_dir, cell_type_names)

pdf(paste0(project_output_dir, "/first_type.pdf"))
plot <- barplot(summary(results$results_df$first_type))
print(plot)
dev.off()

pdf(paste0(project_output_dir, "/spot_class.pdf"))
plot <- barplot(summary(results$results_df$spot_class))
print(plot)
dev.off()