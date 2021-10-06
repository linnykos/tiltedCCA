form_metacells <- function(mat_1, mat_2 = NA,
                           clustering_resolution = 1,
                           dims_1 = NA, dims_2 = NA,
                           center_1 = T, center_2 = T,
                           scale_1 = T, scale_2 = T,
                           verbose = F){
  if(all(is.na(mat_2))){
    mat_1 <- consensus_pca(mat_1, mat_2,
                           dims_1 = dims_1, dims_2 = dims_2,
                           center_1 = center_1, center_2 = center_2,
                           scale_1 = scale_1, scale_2 = scale_2)
  }
  
  if(length(rownames(mat_1)) == 0) rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
  if(length(colnames(mat_1)) == 0) colnames(mat_1) <- paste0("p", 1:ncol(mat_1))
  seurat_obj <- Seurat::CreateAssayObject(data = t(mat_1))
  
  snn_res <- Seurat::FindNeighbors(mat_1,
                                   compute.SNN = T)
  
  snn_res$snn@assay.used <- "RNA"
  seurat_obj[["newgraph"]] <- snn_res$snn
  seurat_obj <- Seurat::FindClusters(seurat_obj, 
                                     graph.name = "newgraph",
                                     resolution = clustering_resolution)
  
  as.factor(seurat_obj@meta.data[,"seurat_clusters"])
}

##################3

.compute_metamatrix <- function(mat, 
                                meta_clustering,
                                verbose){
  stopifnot(is.factor(meta_clustering),
            length(meta_clustering) == nrow(mat))
  num_meta <- length(levels(meta_clustering))
  
  mat_meta <- t(sapply(1:num_meta, function(x){
    if(verbose && num_meta > 10 && x %% floor(num_meta/10) == 0) cat('*')
    idx <- which(meta_clustering == levels(meta_clustering)[x])
    
    if(inherits(mat, "dgCMatrix")){
      sparseMatrixStats::colMeans2(mat[idx,,drop = F])
    } else {
      matrixStats::colMeans2(mat[idx,,drop = F])
    }
  }))
  
  mat_meta
}