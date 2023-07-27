rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(tidyverse)

# Params
project.name <- "2022_03_novaseq_PM22082"
output.name <- "sham_vs_1gy_SPE"

#                           1-F41-CTRL   2-F50-CTRL   3-F44-CTRL   
samples_ids_sham_ctrl <- c("211013_23", "211125_11", "211013_03")

#                     1-F61-GRB    2-F68-GRB    3-F65-GRB   
samples_ids_grb <- c("211125_22", "211125_04", "211125_21")

#                     1-F01-GCR    2-F08-GCR    3-F03-GCR   
samples_ids_gcr <- c("211013_21", "211125_20", "211004_12")

#                     1-F21-SPE    2-F27-SPE    3-F22-SPE   
samples_ids_spe <- c("211004_02", "211125_18", "211013_17")

target_samples <- c(samples_ids_sham_ctrl, samples_ids_spe)

# Tissue types
tissue_types <- c("DCT", "Early PT", "Endo", "IC", "LOH", "Macro", "Neutro", "NP", "PC", "PCT", "Podo", "Proliferating", "PST", "Stroma 1")

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
deseq_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/deseq")
deseq_output_path <- paste0(deseq_output_dir, "/", output.name, ".rds")
dir.create(deseq_output_dir)

# Init
gene_sums <- list()
combined = NULL

# Loop tissue type
for(j in 1:length(tissue_types)) {
  t_type <- tissue_types[j]

  #tissue_output_dir <- paste0(deseq_output_dir, t_type)
  name = paste0(output.name, "_", t_type)
  deseq_output_path <- paste0(deseq_output_dir, "/", name, ".rds")
  table_output_path <- paste0(deseq_output_dir, "/", name, ".tsv")
  dir.create(deseq_output_dir)

  # Loop samples
  for (i in 1:length(target_samples)) {
    sample <- target_samples[i]
    project_output_dir <- paste0(data_dir, "/", project.name, "_gen/", sample)
    rctd_path <- paste0(project_output_dir, "/", sample, "_rctd.rds")
    rctd <- readRDS(rctd_path)
    spatialRNA <- rctd@originalSpatialRNA
    #spatialRNA <- rctd@spatialRNA
    results_df <- rctd@results$results_df
    
    barcodes <- colnames(spatialRNA@counts)
    filt_ind <- results_df$first_type == t_type
    target_barcodes <- barcodes[filt_ind]
    target_counts <- spatialRNA@counts[,target_barcodes]
    #target_counts <- spatialRNA@counts
    
    print(length(rownames(target_counts)))
    
    row_sums <- rowSums(target_counts)
    row_names <- rownames(target_counts)
    df_gene_sums <- data.frame(row_sums)
    rownames(df_gene_sums) <- row_names
    colnames(df_gene_sums) <- c(sample)
    df_gene_sums <- tibble::rownames_to_column(df_gene_sums, "gene_id")

    gene_sums[[i]] <- df_gene_sums
  }

  combined <- gene_sums %>% reduce(inner_join, by = "gene_id")
  head(combined)
  saveRDS(combined, file = deseq_output_path)
  
  write.table(
    combined,
    file = table_output_path,
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
  )
}

