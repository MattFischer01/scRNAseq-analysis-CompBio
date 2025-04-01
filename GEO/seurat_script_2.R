# install and load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(future)

#### read in data and make seurat object ####
setwd('/home/2025/hjensen2/scRNAseq-analysis-CompBio/GEO')

#### Function to perform everything before integration (QC and normalization) ####

process_dge_file <- function(file_path, min_cells = 3, min_features = 200, mt_threshold = 5, nFeature_min = 200, nFeature_max = 2500) {
  
  # extract GSM ID from the filename (everything before the first underscore)
  gsm_id <- sub("_.*", "", basename(file_path))
  
  # read file without hardcoding column count by reading number of rows 
  dge_file <- read.table(gzfile(file_path), header = TRUE, row.names = 1, 
                         colClasses = c("character", rep("numeric", ncol(read.table(gzfile(file_path), header = TRUE, nrows = 1)) - 1)))
  
  # create seurat object with GSM ID as project name so i can track it later
  dge <- CreateSeuratObject(counts = dge_file, project = gsm_id, min.cells = min_cells, min.features = min_features)
  
  # calculate mitochondrial percentage and add it to the metadata
  dge[["percent.mt"]] <- PercentageFeatureSet(dge, pattern = "^mt-")
  
  # filter out doublets and unwanted cells (anything under and over a certain RNA feature number, and over 5% mt%)
  dge <- subset(dge, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < mt_threshold)
  
  # normalize data with defaults (log normalization)
  dge <- NormalizeData(dge)
  
  return(dge)
}

# read text file to get file names
file_list <- readLines("geo_list.txt")

# process files and store them in a named list
dge_list <- lapply(file_list, process_dge_file)

# dont think we need below command anymore
# names(dge_list) <- paste0("dge_", tools::file_path_sans_ext(basename(file_list)))

# change original indent to everything before .dge.txt.gz instead of full file name
dge_list <- lapply(dge_list, function(dge) {
  dge@meta.data$orig.ident <- gsub("\\.dge.*", "", dge@meta.data$orig.ident)
  return(dge)
})

# read in metadata file with sample ID and condition (user will need to provide this based on their needs)
# must be comma delineated
metadata <- read.table("metadata_seurat.txt", sep = ',', header = TRUE)

# add an empty condition column to populate later, not working to just add it directly
dge_list <- lapply(dge_list, function(dge) {
  dge@meta.data$condition <- NA  # initialize empty column with NAs
  return(dge)
})

# match the condition to the right id based on the metadata text file
dge_list <- lapply(dge_list, function(dge) {
  sample_id <- unique(dge@meta.data$orig.ident)  # get sample ID from meta data
  
  # find matching condition from metadata
  condition_value <- metadata$condition[match(sample_id, metadata$sample)]
  
  # if a match is found, update condition column in the seurat object metadata
  if (!is.na(condition_value)) {
    dge@meta.data$condition <- condition_value  
  }
  
  return(dge)
})

#### end of part before canonical correlation analysis and integration ####
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# above might be good tutorial

# finding variable features for each dataset
dge_list <- lapply(dge_list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dge_list)

# find integration anchors, this step runs CCA to control for batch effects
# this takes a long time to run, almost two hrs so going to put it in background using future package
# this will need to be taken out though for final wrapper script
parallel::detectCores() # check cores 

# increasing size limit for future package, for some reason cant do jobs with objects above .5GB
options(future.globals.maxSize = 1000 * 1024^2)

plan(multisession, workers = 2)

# rename cells before integration so samples wont be lost, sample names need to be changed eventually
for (i in seq_along(dge_list)) {
  dge_list[[i]] <- RenameCells(dge_list[[i]], add.cell.id = paste0("sample", i))
}

# use the features found earlier as the anchor.features
future_anchors <- future({
  FindIntegrationAnchors(object.list = dge_list, anchor.features = features)
})

# retrieve results
anchors <- value(future_anchors)

# back to regular R sequential running
plan(sequential)

# integrate the data
mice.combined <- IntegrateData(anchorset = anchors)

# saving everything from previous analysis because it will take forever to run everything if envs is lost
save(anchors, dge_list, mice.combined, file = "seurat_post_int.RData")

# specify that the new assay is integrated instead of RNA so we are working at higher level
DefaultAssay(mice.combined) <- "integrated"

#### standard workflow from this point on ####

# scale data with defaults
mice.combined <- ScaleData(object = mice.combined)

# run PCA with features set as VariableFeatures(integrated_data)
mice.combined <- RunPCA(mice.combined, features = VariableFeatures(object = mice.combined))

# can look at cells and features that define the PCA
VizDimLoadings(mice.combined, dims = 1:2, reduction = "pca")
DimPlot(mice.combined, reduction = "pca") + NoLegend()
DimHeatmap(mice.combined, dims = 1:15, cells = 500, balanced = TRUE)

# determining the dimensionality of the dataset
ElbowPlot(mice.combined)
# shoot for something btwn 15 and 20?

#saving post pca
save(anchors, dge_list, mice.combined, file = "post_pca_seurat_geo_data.RData")

# run tsne based on first 15 PCs
mice.combined <- RunTSNE(mice.combined, reduction = "pca", dims = 1:15)

# find neighbors with same first 15 PCs
mice.combined <- FindNeighbors(mice.combined, reduction = "pca", dims = 1:15)

# find clusters, can play with resolution depending on how our results look
mice.combined <- FindClusters(mice.combined, resolution = 0.5)

# visualization
p1 <- DimPlot(mice.combined, reduction = "tsne", label = TRUE, repel = TRUE)
p1




