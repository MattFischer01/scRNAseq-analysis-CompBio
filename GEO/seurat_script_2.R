# install and load packages
library(dplyr)
library(MAST)
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
load("seurat_post_int.RData")


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
save(anchors, dge_list, mice.combined, file = "post_pca_seurat_metadata.RData")
load("post_pca_seurat_metadata.RData")

# run tsne based on first 15 PCs
mice.combined <- RunTSNE(mice.combined, reduction = "pca", dims = 1:15)

# find neighbors with same first 15 PCs
mice.combined <- FindNeighbors(mice.combined, reduction = "pca", dims = 1:15)

# find clusters, can play with resolution depending on how our results look
mice.combined <- FindClusters(mice.combined, resolution = 0.5) # the paper used 0.5

# visualization, this looks at only the mice with EAE in the condition name
p1 <- DimPlot(subset(mice.combined, subset = condition %in% c("EAE", "EAE priming", "EAE remission")), 
              reduction = "tsne", label = TRUE, repel = TRUE)
# this plot groups by condition
p2 <- DimPlot(mice.combined, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "condition")
# this plot is all of the samples
p3 <- DimPlot(mice.combined, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "seurat_clusters")

p1|p2|p3 # plot all of them side by side

save(dge_list, anchors, mice.combined, file = "post_clusters_metadata.RData")

#### MAST ####

# need to find a tutorial on how to do DE analysis for each cluster compared to all other cells
# looks like in the FindMarkers() function, there is a test.use parameter that we can do MAST

# set the active identity of object to the clusters
Idents(mice.combined) <- "seurat_clusters"

# set default assay to "RNA"
DefaultAssay(mice.combined) <- "RNA"

# join layers, find markers wont run without doing this and it suggested it
mice.combined_joined <- JoinLayers(mice.combined)


# make a list to store the DE results for each cluster
mast_results <- list()

# this takes awhile, so save rdata after this step
# for each cluster, compare it to all other cells and return DE results to mast_results 
for (cluster in levels(Idents(mice.combined_joined))) {
  mast_results[[cluster]] <- FindMarkers(mice.combined_joined,
                                         ident.1 = cluster, # comparing current cluster to all other cells 
                                         test.use = "MAST", # MAST algorithm
                                         min.pct = 0.25,  # only test genes expressed in at least 25% of cells in a cluster
                                         logfc.threshold = 0.25) # limits testing on genes, speeds up process (defualt is 0.1)
}


# save mast_results
save(dge_list, mice.combined, mice.combined_joined, mast_results, file = "mast_results.RData")

# look at top DE genes for a specific cluster
head(mast_results[["0"]]) # replace 0 with cluster of interest
head(mast_results[["1"]])
# we will need to somehow compare markers to cell type, not sure how we will do this in our pipeline

# feature plots from the paper, they look diff but the cluster coloring and # of clusters colored is similar
FeaturePlot(mice.combined, c('S100b'))
FeaturePlot(mice.combined, c('Gja1'))
FeaturePlot(mice.combined, c('Aldh1l1'))
FeaturePlot(mice.combined, c('Gfap')) # this had some cells removed bc it didnt have this feature
FeaturePlot(mice.combined, c('Aqp4'))


