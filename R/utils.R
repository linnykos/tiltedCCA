#' Construct cell-type subsamples
#' 
#' Sample at most \code{min_subsample_cell} cells from each
#' unique level of \code{membership_vec}
#'
#' @param membership_vec factor vector
#' @param min_subsample_cell positive integer
#'
#' @return indices
#' @export
construct_celltype_subsample <- function(membership_vec, min_subsample_cell){
  stopifnot(is.factor(membership_vec))
  
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    stopifnot(length(idx) > 2)
    
    if(length(idx) <= min_subsample_cell) return(idx)
    
    sample(idx, min_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}


form_seurat_obj <- function(mat_1, mat_2){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1); p_1 <- ncol(mat_1); p_2 <- ncol(mat_2)
  
  obj <- Seurat::CreateSeuratObject(counts = t(mat_1), assay = "mode1")
  obj[["mode2"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
  
  obj
}

.l2norm <- function(x){sqrt(sum(x^2))}

#######################

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

###########################

.combine_two_named_lists <- function(list1, list2){
  idx <- which(!names(list2) %in% names(list1))
  for(i in idx){
    list1[[names(list2)[i]]] <- list2[[i]]
  }
  
  list1
}

# return idx such that vec1[idx] == vec2
.matching_idx <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  
  ord_1 <- order(vec1, decreasing = F)
  rank_2 <- rank(vec2)
  
  ord_1[rank_2]
}

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= nrow(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

## [[note to self: allow target_obj to be a Seurat obj]]
.append_rowcolnames <- function(bool_colnames,
                                bool_rownames,
                                source_obj,
                                target_obj){
  stopifnot(bool_colnames | bool_rownames)
  
  if(inherits(target_obj, "svd")){
    if(inherits(source_obj, "svd")){
      if(bool_rownames) {
        stopifnot(nrow(source_obj$u) == nrow(target_obj$u))
        rownames(target_obj$u) <- rownames(source_obj$u)
      }
      if(bool_colnames){
        stopifnot(nrow(source_obj$v) == nrow(target_obj$v))
        rownames(target_obj$v) <- rownames(source_obj$v)
      }
      
    } else if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      if(bool_rownames) {
        stopifnot(nrow(source_obj) == nrow(target_obj$u))
        rownames(target_obj$u) <- rownames(source_obj)
      }
      if(bool_colnames){
        stopifnot(ncol(source_obj) == nrow(target_obj$v))
        rownames(target_obj$v) <- colnames(source_obj)
      }
    } else {
      stop("Cannot identify class of source_obj")
    }
    
  } else if(inherits(target_obj, c("matrix", "dgCMatrix"))){
    if(inherits(source_obj, "svd")){
      if(bool_rownames) {
        stopifnot(nrow(source_obj$u) == nrow(target_obj))
        rownames(target_obj) <- rownames(source_obj$u)
      }
      if(bool_colnames){
        stopifnot(nrow(source_obj$v) == ncol(target_obj))
        colnames(target_obj) <- rownames(source_obj$v)
      }
      
    } else if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      if(bool_rownames) {
        stopifnot(nrow(source_obj) == nrow(target_obj))
        rownames(target_obj) <- rownames(source_obj)
      }
      if(bool_colnames){
        stopifnot(ncol(source_obj) == ncol(target_obj))
        colnames(target_obj) <- colnames(source_obj)
      }
    } else {
      stop("Cannot identify class of source_obj")
    }
    
  } else {
    stop("Cannot identify class of target_obj")
  }
  
  target_obj
}
