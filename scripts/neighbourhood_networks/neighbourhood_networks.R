## CRC visiumHD - neighbourhood networks

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(ComplexHeatmap)
library(parallel)
library(future)
library(future.apply)
library(doParallel) 
library(osfr)
library(circlize)
library(RColorBrewer)
library(openxlsx)
library(scico)
library(stats)
library(rstatix)
library(igraph)
library(ggraph)
library(viridis)

#--------------------------------------
# Load functions
#--------------------------------------
source("CRC_program_neighbourhood_network.R")
source("calculate_outgoing_and_incoming_connections.R")


#--------------------------------------
# CRC spatial data
#--------------------------------------
## Prepare input
usage_matrix = "topics_per_spot_P5_CRC_GRCh38_50_24um.csv"
annot_file = "ss_annotation_all_corr_Level2_pearson_genes_pearson.csv"

usage_norm <- read.csv(usage_matrix, row.names = 1)
annot <- read.csv(annot_file, row.names = 1)
annot$celltype <- gsub("T cell","Tcell",annot$celltype)
annot$celltype <- gsub(" ","_", annot$celltype)
annot$celltype <- gsub("\\(","_",annot$celltype)
annot$celltype <- gsub("\\)","",annot$celltype)

for(i in 1:ncol(usage_norm)){
    if(colnames(usage_norm)[i] %in% annot$program){
        colnames(usage_norm)[i] <- paste0( gsub("ot_","ot", colnames(usage_norm)[i]), "_", annot$celltype[which(annot$program == colnames(usage_norm)[i])] )
    }else{
        colnames(usage_norm)[i] <- gsub("ot_","ot", colnames(usage_norm)[i])
    }
}

## Set dynamic thresholds in usage matrix
usage_norm_dynamic_thresh <- usage_norm
for(i in 1:ncol(usage_norm_dynamic_thresh)){
    thresh = quantile(usage_norm_dynamic_thresh[,i], 0.90)
    usage_norm_dynamic_thresh[,i][which(usage_norm_dynamic_thresh[,i] < thresh)] <- 0 
}


spotmetadata_all <- data.frame(Bin = rownames(usage_norm), sample_id = "CRC", row.names = rownames(usage_norm))


## Get neighbourhood networks

cc_interaction_networks(usage_matrix = usage_norm_dynamic_thresh,
                                        metadata = spotmetadata_all, 
                                        metadata_column = "sample_id",
                                        fileprefix = "CRC_",
                                        cond1 = "CRC",
                                        cond2 = NULL, 
                                        save_folder = "/work/morrissy_lab/vthoppey/CRC",
                                        GLOBAL_THRESHOLD = 0,
                                        edge_threshold = "gt02",
                                        n_bins = 1000, 
                                        connections_heatmap = TRUE,
                                        plot_pval_heatmap = FALSE)



#--------------------------------------
# PDX spatial data
#--------------------------------------
## Prepare input
usage_matrix = "/work/morrissy_lab/vthoppey/SpatialData/neighbourhood_networks/topics_per_spot_pdx_merge_all_mm10_90.csv"
annot_file = "/work/morrissy_lab/vthoppey/SpatialData/neighbourhood_networks/programs_annotations.csv"
TMEprog <- paste0("ot_", c(76, 19, 16, 89, 5, 29, 21, 2, 87, 61, 1, 64, 75, 27, 20, 84, 36, 3, 23, 56, 47, 8))

usage_norm <- read.csv(usage_matrix, row.names = 1)
usage_norm <- usage_norm[,TMEprog]

annot <- read.csv(annot_file)
annot$Manual.Annotation <- gsub("T cell","Tcell",annot$Manual.Annotation)
annot$Manual.Annotation <- gsub(" ","_", annot$Manual.Annotation)
annot$Manual.Annotation <- gsub("\\(","_",annot$Manual.Annotation)
annot$Manual.Annotation <- gsub("\\)","",annot$Manual.Annotation)
annot <- annot %>% filter(Manual.Annotation != "")

for(i in 1:ncol(usage_norm)){
    if(colnames(usage_norm)[i] %in% annot$Program){
        colnames(usage_norm)[i] <- paste0( gsub("ot_","ot", colnames(usage_norm)[i]), "_", annot$Manual.Annotation[which(annot$Program == colnames(usage_norm)[i])] )
    }else{
        colnames(usage_norm)[i] <- gsub("ot_","ot", colnames(usage_norm)[i])
    }
}

## Set dynamic thresholds in usage matrix
usage_norm_dynamic_thresh <- usage_norm
for(i in 1:ncol(usage_norm_dynamic_thresh)){
    thresh = quantile(usage_norm_dynamic_thresh[,i], 0.90)
    usage_norm_dynamic_thresh[,i][which(usage_norm_dynamic_thresh[,i] < thresh)] <- 0 
}

spotmetadata_all <- data.frame(Bin = rownames(usage_norm), sample_id = "GBMxeno", row.names = rownames(usage_norm))

## Get neighbourhood networks
source("CRC_program_neighbourhood_network.R")
source("calculate_outgoing_and_incoming_connections.R")
cc_interaction_networks(usage_matrix = usage_norm_dynamic_thresh,
                                        metadata = spotmetadata_all, 
                                        metadata_column = "sample_id",
                                        fileprefix = "GBMxeno_",
                                        cond1 = "GBMxeno",
                                        cond2 = NULL, 
                                        save_folder = "/work/morrissy_lab/vthoppey/SpatialData/neighbourhood_networks",
                                        GLOBAL_THRESHOLD = 0,
                                        edge_threshold = "gt02",
                                        n_bins = 1000,
                                        plot_pval_heatmap = FALSE,
                                        connections_heatmap = FALSE)
