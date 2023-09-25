library(Seurat)
library(tidyverse)

obj_list = readRDS("/afs/crc.nd.edu/user/m/mzarodn2/Private/scrnaseq/data/Brain_atlas/obj_list.rds")

integrate_seurat <- function(obj_list) {
  
  options(future.globals.maxSize = 4000 * 1024^2)
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  
  integ_features <- SelectIntegrationFeatures(object.list = obj_list) 
 
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- ScaleData(x, features = integ_features, verbose = FALSE)
    x <- RunPCA(x, features = integ_features, verbose = FALSE)
  })
  
  # Find best buddies - can take a while to run
  integ_anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                          anchor.features = integ_features,
                                          reduction = "rpca")
  
  # Integrate across conditions
  obj <- IntegrateData(anchorset = integ_anchors, k.weight = 50)
  return(obj)
}

integrated_seurat <- integrate_seurat(obj_list)

saveRDS(integrated_seurat, "/afs/crc.nd.edu/user/m/mzarodn2/Private/scrnaseq/data/Brain_atlas/integrated_seurat.rds"")