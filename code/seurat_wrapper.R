# install and load packages
library(dplyr)
library(MAST)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please include path to DGE files!")
  stop("Requires command line argument.")
}

#### read in data and make seurat object ####
setwd(args[1])

#### Function to perform everything before integration (QC and normalization) ####

process_dge_file <- function(file_path, min_cells = 3, min_features = 200, mt_threshold = 5, nFeature_min = 200, nFeature_max = 2500) {
  
  # extract GSM ID from the filename
  gsm_id <- sub("_.*", "", basename(file_path))
  
  # read .gz file using gzfile() 
  dge_file <- read.table(gzfile(file_path), header = TRUE, row.names = 1, 
                         colClasses = c("character", rep("numeric", ncol(read.table(gzfile(file_path), header = TRUE, nrows = 1)) - 1)))
  
  # output the number of features and cells
  cat("Number of features (genes):", nrow(dge_file), "\n")
  cat("Number of cells:", ncol(dge_file), "\n")
  
  # create Seurat object
  dge <- CreateSeuratObject(counts = dge_file, project = gsm_id, min.cells = min_cells, min.features = min_features)
  
  # calculate mitochondrial percentage
  dge[["percent.mt"]] <- PercentageFeatureSet(dge, pattern = "^mt-")
  
  # log the number of cells before filtering
  cat("Cells before filtering:", ncol(dge), "\n")
  
  # filter out doublets and unwanted cells
  dge <- subset(dge, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < mt_threshold)
  
  # log the number of cells after filtering
  cat("Cells after filtering:", ncol(dge), "\n")
  
  # normalize data
  dge <- NormalizeData(dge)
  
  return(dge)
}

# read text file to get paths, 
lines <- readLines("DGE-paths.txt")

# parse out diff elements based on white space
parsed <- strsplit(lines, "\\s+")

labels <- sapply(parsed, function(x) x[1]) # gets labels as first part 
paths <- sapply(parsed, function(x) x[3]) # gets paths as third part
file_list <- paths

# process files and store in a named list with error handling
dge_list <- list()
clean_labels <- c()

# loops through all file paths
for (i in seq_along(paths)) {
  file_path <- paths[i]
  label <- labels[i]
  
  message("Processing: ", file_path) # this says which file is being processed

  # if there is an error processing DGE files, still store the dge files that do run, but send an error message
  result <- tryCatch({
    dge <- process_dge_file(file_path) # processes the DGE file with function above
    dge
  }, error = function(e) {
    message("Error processing file: ", file_path)
    message("   → ", conditionMessage(e)) # gives the error for the file that errored out
    return(NULL) # sets results to NULL instead of crashing process
  })

# this just filters out any results that errored, adds the result to dge_list, and adds the associated label to clean_labels
  if (!is.null(result)) {
    dge_list[[length(dge_list) + 1]] <- result
    clean_labels <- c(clean_labels, label) 
  }
}

# update orig.ident to the dge files
dge_list <- lapply(dge_list, function(dge) {
  dge@meta.data$orig.ident <- gsub("\\.dge.*", "", dge@meta.data$orig.ident)
  return(dge)
})

# assign conditions to the metdata based on labels
dge_list <- Map(function(dge, label) {
  dge@meta.data$condition <- label
  return(dge)
}, dge_list, clean_labels)


save(dge_list, clean_labels, file = "../pre_integration.RData")

#### end of part before canonical correlation analysis and integration ####
# https://satijalab.org/seurat/archive/v4.3/integration_introduction
# above might be good tutorial

# finding variable features for each dataset
dge_list <- lapply(dge_list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dge_list)

# find integration anchors, this step runs CCA to control for batch effects
unique_labels <- make.unique(clean_labels, sep = "_")

# rename cells before integration so samples wont be lost
for (i in seq_along(dge_list)) {
  sample_name <- unique_labels[i]  # use cleaned, consistent sample name
  dge_list[[i]] <- RenameCells(dge_list[[i]], add.cell.id = sample_name) # rename cells by adding the sample name as prefix
  dge_list[[i]]@meta.data$orig.ident <- sample_name # update the 'orig.ident' metadata field to reflect the sample name
}

# use the features found earlier as the anchor.features
anchors <- FindIntegrationAnchors(object.list = dge_list, anchor.features = features)

# integrate the data with the anchors found earlier
mice.combined <- IntegrateData(anchorset = anchors)

# saving everything from previous analysis so it can be loaded in after wrapper
save(anchors, dge_list, mice.combined, file = "../seurat_post_int.RData")

