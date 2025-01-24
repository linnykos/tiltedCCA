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
  membership_vec <- droplevels(membership_vec)
  
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    stopifnot(length(idx) > 2)
    
    if(length(idx) <= min_subsample_cell) return(idx)
    
    sample(idx, min_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}


form_seurat_obj <- function(mat_1, mat_2, suppress_warnings = TRUE){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1); p_1 <- ncol(mat_1); p_2 <- ncol(mat_2)
  
  if(suppress_warnings){
    suppressWarnings(obj <- Seurat::CreateSeuratObject(counts = t(mat_1), assay = "mode1"))
  } else {
    obj <- Seurat::CreateSeuratObject(counts = t(mat_1), assay = "mode1")
  }
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

# see https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
.combine_two_named_lists <- function(list1, list2){
  idx <- which(!names(list2) %in% names(list1))
  for(i in idx){
    if(all(is.null(list2[[i]]))){
      list1 <- c(list1, list(TEMP_NAME = NULL))
      names(list1)[which(names(list1) == "TEMP_NAME")] <- names(list2)[i]
    } else {
      list1[[names(list2)[i]]] <- list2[[i]]
    }
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
            col_idx > 0, col_idx <= ncol(mat))
  
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

## if the source_obj is a Seurat object, and bool_colnames is TRUE, this grabs ALL the features
.append_rowcolnames <- function(bool_colnames,
                                bool_rownames,
                                source_obj,
                                target_obj){
  stopifnot(bool_colnames | bool_rownames)
  
  if(bool_rownames){
    if(inherits(source_obj, "svd")){
      rowname_vec <- rownames(source_obj$u)
    } else if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      rowname_vec <- rownames(source_obj)
    } else if(inherits(source_obj, "multiSVD")){
      rowname_vec <- rownames(.get_SVD(source_obj)$u)
    } else if(inherits(source_obj, "Seurat")){
      # this refers to the "cells" in this TCCA package
      rowname_vec <- Seurat::Cells(source_obj)
    } else {
      stop("Not valid class for source_obj")
    }
    
    if(inherits(target_obj, "svd") & length(rowname_vec) > 0){
      rownames(target_obj$u) <- rowname_vec
    } else if(inherits(target_obj, c("matrix", "dgCMatrix")) & length(rowname_vec) > 0){
      rownames(target_obj) <- rowname_vec
    } 
  }
  
  if(bool_colnames){
    if(inherits(source_obj, "svd")){
      colname_vec <- rownames(source_obj$v)
    } else if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      colname_vec <- colnames(source_obj)
    } else if(inherits(source_obj, "multiSVD")){
      colname_vec <- rownames(.get_SVD(source_obj)$v)
    } else if(inherits(source_obj, "Seurat")){
      colname_vec <- SeuratObject::Features(source_obj)
    } else {
      stop("Cannot identify class for source_obj")
    }
    
    if(inherits(target_obj, "svd") & length(colname_vec) > 0){
      rownames(target_obj$v) <- colname_vec
    } else if(inherits(target_obj, c("matrix", "dgCMatrix")) & length(colname_vec) > 0){
      colnames(target_obj) <- colname_vec
    } 
  }
  
  target_obj
}

.convert_factor2list <- function(vec){
  stopifnot(is.factor(vec))
  vec <- droplevels(vec)
  level_vec <- levels(vec)
  lis <- lapply(level_vec, function(level_val){
    which(vec == level_val)
  })
  names(lis) <- level_vec
  lis
}

.convert_list2factor <- function(lis, n){
  stopifnot(is.list(lis), n >= max(unlist(lis)))
  vec <- rep(NA, n)
  if(length(names(lis)) != length(lis)){
    name_vec <- as.character(1:length(lis))
  } else {
    name_vec <- names(lis)
  }
  
  for(i in 1:length(lis)){
    vec[lis[[i]]] <- name_vec[i]
  }
  vec <- as.factor(vec)
  
  vec
}

.linear_regression <- function(bool_include_intercept,
                               bool_center_x,
                               bool_center_y,
                               bool_scale_x,
                               bool_scale_y,
                               return_type, 
                               x_mat,
                               y_vec,
                               tol = 1e-6){
  if(!is.matrix(x_mat)) x_mat <- matrix(x_mat, nrow = length(x_mat), ncol = 1)
  stopifnot(nrow(x_mat) == length(y_vec), return_type %in% c("r_squared"))
  
  if(stats::sd(y_vec) <= tol) {
    if(return_type == "r_squared") return(NA)
  }
  
  if(bool_center_x | bool_scale_x) x_mat <- scale(x_mat, 
                                                  center = bool_center_x,
                                                  scale = bool_scale_x)
  if(bool_center_y | bool_scale_y) y_vec <- scale(y_vec, 
                                                  center = bool_center_y,
                                                  scale = bool_scale_y)
  
  df <- data.frame(cbind(y_vec, x_mat))
  colnames(df)[1] <- "y"
  
  if(bool_include_intercept){
    lm_res <- stats::lm(y ~ ., data = df)
  } else {
    lm_res <- stats::lm(y ~ . - 1, data = df)
  }
  
  if(return_type == "r_squared"){
    summary(lm_res)$r.squared
  }
}