#' Create Seurat dimension reduction objects of Tilted-CCA 
#'
#' @param input_obj a \code{multiSVD} object that is the output of \code{tiltedCCA_decomposition}
#' @param what a character, either \code{"common"}, \code{"distinct_1"}, or \code{"distinct_2"} for 
#' which embedding to construct the visualization for
#' @param aligned_umap_assay either \code{NULL} or a UMAP assay in \code{seurat_assay} 
#' (in which case the resulting UMAP will be rotated to best mimic the relative orientation
#' of cells in \code{seurat_obj[[aligned_umap_assay]]})
#' @param scale_max_1 numeric or \code{NULL}, to threshold Modality 1 in magnitude prior to computing latent dimensions
#' @param scale_max_2 numeric or \code{NULL}, to threshold Modality 2 in magnitude prior to computing latent dimensions
#' @param seurat_obj the \code{Seurat} object that was used to compute \code{input_obj}, the \code{multiSVD_obj}
#' @param seurat_assay the assay in \code{seurat_obj} to assign the resulting embedding to
#' @param verbose non-negative integer
#' @param ... 
#'
#' @return a \code{DimReduc} \code{SeuratObject}
#' @export
create_SeuratDim <- function(input_obj,
                             what,
                             aligned_umap_assay = NULL,
                             scale_max_1 = NULL,
                             scale_max_2 = NULL,
                             seurat_obj = NULL,
                             seurat_assay = "RNA",
                             verbose = 0,
                             ...){
  stopifnot(what %in% c("common", "distinct_1", "distinct_2"))
  
  if(!is.null(scale_max_1)) input_obj$param$svd_scale_max_1 <- scale_max_1
  if(!is.null(scale_max_2)) input_obj$param$svd_scale_max_2 <- scale_max_2
  
  if(what == "common"){
    param <- .get_param(input_obj)
    if(!param$svd_normalize_singular_value){
      warning("The normalize_singular_value is FALSE in input_obj, arguably nonsensical.")
    }
    
    input_obj <- .set_defaultAssay(input_obj, assay = 1)
    if("common_mat_1" %in% names(input_obj)){
      dimred_1 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_mat")
    } else if("common_dimred_1" %in% names(input_obj)){
      dimred_1 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
    } else {
      stop("Cannot find the appropriate common matrix for modality 1")
    }
    input_obj <- .set_defaultAssay(input_obj, assay = 2)
    if("common_mat_2" %in% names(input_obj)){
      dimred_2 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_mat")
    } else if("common_dimred_2" %in% names(input_obj)){
      dimred_2 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
    } else {
      stop("Cannot find the appropriate common matrix for modality 2")
    }
    dimred <- cbind(dimred_1, dimred_2)
    
  } else if(what == "distinct_1"){
    input_obj <- .set_defaultAssay(input_obj, assay = 1)
    if("distinct_mat_1" %in% names(input_obj)){
      dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
    } else if("distinct_dimred_1" %in% names(input_obj)){
      dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_dimred")
    } else {
      stop("Cannot find the appropriate distinct matrix for modality 1")
    }
    
  } else {
    input_obj <- .set_defaultAssay(input_obj, assay = 2)
    if("distinct_mat_2" %in% names(input_obj)){
      dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
    } else if("distinct_dimred_2" %in% names(input_obj)){
      dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_dimred")
    } else {
      stop("Cannot find the appropriate distinct matrix for modality 2")
    }
  } 
  
  seurat_umap <- Seurat::RunUMAP(dimred, 
                                 assay = seurat_assay,
                                 verbose = (verbose != 0),
                                 ...)
  rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
  
  if(!all(is.null(seurat_obj)) & !is.null(aligned_umap_assay)){
    seurat_umap@cell.embeddings <- .rotate_matrix(source_mat = seurat_obj[[aligned_umap_assay]]@cell.embeddings,
                                                  target_mat = seurat_umap@cell.embeddings)
  }
  
  Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                               assay = seurat_assay)
}

create_reducedSeuratObj <- function(input_obj,
                                    what,
                                    aligned_umap_assay = NULL,
                                    seurat_celltype = NULL,
                                    seurat_obj = NULL,
                                    verbose = 0,
                                    ...){
  stopifnot(what %in% c("common_basis", "laplacian_1",
                        "laplacian_2", "laplacian_common"))
  
  if(what == "common_basis"){
    dimred <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_basis")
    
  } else if(what == "laplacian_1"){
    input_obj <- .set_defaultAssay(input_obj, assay = 1)
    dimred <- .get_Laplacian(input_obj, bool_common = F)
    
  } else if(what == "laplacian_2"){
    input_obj <- .set_defaultAssay(input_obj, assay = 2)
    dimred <- .get_Laplacian(input_obj, bool_common = F)
    
  } else {
    dimred <- .get_Laplacian(input_obj, bool_common = T)
  }
  
  colnames(dimred) <- paste0("tmp", 1:ncol(dimred))
  new_obj <- Seurat::CreateSeuratObject(counts = t(dimred))
  seurat_umap <- Seurat::RunUMAP(dimred, assay = "RNA", 
                                 verbose = (verbose != 0),
                                 ...)
  rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
  
  
  if(!all(is.null(seurat_obj)) & !is.null(seurat_celltype) & !is.null(aligned_umap_assay)){
    metacell_list <- .get_metacell(input_obj, 
                                   resolution = "cell",
                                   type = "list",
                                   what = "metacell_clustering")
    if(!all(is.null(metacell_list))){
      stopifnot(all(names(metacell_list) == rownames(dimred)))
      reduced_umap <- t(sapply(metacell_list, function(idx_vec){
        colMeans(seurat_obj[[aligned_umap_assay]]@cell.embeddings[idx_vec,,drop=F])
      }))
    } else {
      reduced_umap <- seurat_obj[[aligned_umap_assay]]@cell.embeddings
    }
    
    new_obj$celltype <- .translate_celltype(input_obj = input_obj,
                                            celltype_vec = seurat_obj@meta.data[,seurat_celltype],
                                            metacell_list = metacell_list)
    
    seurat_umap@cell.embeddings <- .rotate_matrix(source_mat = reduced_umap,
                                                  target_mat = seurat_umap@cell.embeddings)
    new_obj[["umap"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                                      assay = "RNA")
    return(new_obj)
    
  } else {
    return(Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                        assay = "RNA"))
  }
}

#################################

.rotate_matrix <- function(source_mat, 
                           target_mat){
  stopifnot(all(dim(source_mat) == dim(target_mat)))
  
  tmp <- svd(t(source_mat) %*% target_mat)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- target_mat %*% t(rotation_mat)
  rownames(tmp) <- rownames(target_mat)
  colnames(tmp) <- colnames(target_mat)
  tmp
}

.translate_celltype <- function(input_obj,
                                celltype_vec,
                                metacell_list){
  if(all(is.null(metacell_list))) return(celltype_vec)
  
  vec <- sapply(metacell_list, function(idx_vec){
    tmp <- sort(table(celltype_vec[idx_vec]), decreasing = T)
    names(tmp)[1]
  })
  
  as.factor(vec)
}

###################

#' Rotate embedding in a Seurat object
#' 
#' This method rotates the embedding in \code{target_embedding}
#' to look (visually) most similar to \code{source_embedding}
#'
#' @param seurat_obj object of class \code{Seurat}
#' @param source_embedding character vector, so that \code{seurat_obj[[source_embedding]]} 
#' is a object of class \code{DimReduc}
#' @param target_embedding character vector, so that \code{seurat_obj[[target_embedding]]} 
#' is a object of class \code{DimReduc}
#'
#' @return an updated \code{seurat_obj}
#' @export
rotate_seurat_embeddings <- function(seurat_obj,
                                     source_embedding,
                                     target_embedding){
  source_mat <- seurat_obj[[source_embedding]]@cell.embeddings
  target_mat <- seurat_obj[[target_embedding]]@cell.embeddings
  
  res <- .rotate_matrix(source_mat = source_mat,
                        target_mat = target_mat)
  
  seurat_obj[[target_embedding]]@cell.embeddings <- res
  seurat_obj
}