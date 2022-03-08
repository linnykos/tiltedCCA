create_multiSVD <- function(mat_1, mat_2,
                            dims_1, dims_2,
                            center_1 = T, center_2 = T,
                            normalize_row = T,
                            normalize_singular_value = T,
                            recenter_1 = F, recenter_2 = F,
                            rescale_1 = F, rescale_2 = F,
                            scale_1 = T, scale_2 = T){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1)
  svd_1 <- .get_SVD(center = center_1, input_obj = mat_1,
                    dims = dims_1, scale = scale_1)
  svd_2 <- .get_SVD(center = center_2, input_obj = mat_2,
                    dims = dims_2, scale = scale_2)
  
  param <- .form_multiSVD_param(center_1 = center_1, center_2 = center_2,
                                n = n,
                                normalize_row = normalize_row,
                                normalize_singular_value = normalize_singular_value,
                                recenter_1 = recenter_1, recenter_2 = recenter_2,
                                rescale_1 = rescale_1, rescale_2 = rescale_2,
                                scale_1 = scale_1, scale_2 = scale_2)
  
  structure(list(svd_1 = svd_1, svd_2 = svd_2,
                 default_assay = 1,
                 param = param),
            class = "multiSVD")
}

###############
#' @export
.normalize_svd <- function(input_obj, ...) UseMethod(".normalize_svd")

#' @export
.normalize_svd.svd <- function(input_obj,
                               averaging_mat,
                               normalize_row,
                               normalize_singular_value,
                               recenter,
                               rescale, ...){
  stopifnot(inherits(input_obj, "svd"))
  
  n <- nrow(input_obj$u)
  dimred <- .get_Dimred(input_obj = input_obj, 
                        normalize_singular_value = normalize_singular_value)
  
  .normalize_svd(input_obj = dimred,
                 averaging_mat = averaging_mat,
                 normalize_row = normalize_row,
                 recenter = recenter,
                 rescale = rescale)
}

#' @export
.normalize_svd.matrix <- function(input_obj,
                                  averaging_mat,
                                  normalize_row,
                                  recenter,
                                  rescale, ...){
  if(recenter | rescale) {
    input_obj <- sapply(1:ncol(input_obj), function(k){scale(input_obj, 
                                                             center = recenter,
                                                             scale = rescale)})
  }
  
  if(!all(is.null(averaging_mat))){
    input_obj <- averaging_mat %*% input_obj
  }
  
  if(normalize_row){
    l2_vec <- apply(input_obj, 1, function(x){.l2norm(x)})
    .mult_vec_mat(1/l2_vec, input_obj)
  }
  
  input_obj
}

####################################

.svd_truncated <- function(mat, K, symmetric, rescale,
                           mean_vec, sd_vec,
                           K_full_rank){
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  if(K == min(dim(mat))) K_full_rank <- T
  
  mean_vec <- .compute_matrix_mean(mat, mean_vec)
  sd_vec <- .compute_matrix_sd(mat, sd_vec)
  
  if(min(dim(mat)) > 2*(K+2)){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      if(symmetric){
        tmp <- irlba::partial_eigen(mat, n = ifelse(K_full_rank, K, K+2),
                                    center = mean_vec, scale = sd_vec)
        list(u = tmp$vectors, d = tmp$values, v = tmp$vectors)
      } else {
        irlba::irlba(mat, nv = ifelse(K_full_rank, K, K+2),
                     center = mean_vec, scale = sd_vec)
      }
    }, warning = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    }, error = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    })
  } else {
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  class(res) <- "svd"
  
  # pass row-names and column-names
  res <- .append_rowcolnames(bool_colnames = T, bool_rownames = T,
                             source_obj = mat, target_obj = res)
  
  # useful only if your application requires only the singular vectors
  # if the number of rows or columns is too large, the singular vectors themselves
  # are often a bit too small numerically
  if(rescale){
    n <- nrow(mat); p <- ncol(mat)
    res$u <- res$u * sqrt(n)
    res$v <- res$v * sqrt(p)
    res$d <- res$d / (sqrt(n)*sqrt(p))
  }
  
  res
}

.check_svd <- function(svd_res, dims, tol = 1e-6){
  idx <- sort(intersect(which(svd_res$d > tol), dims), decreasing = F)
  if(length(idx) == length(svd_res$d)) return(svd_res)
  
  svd_res$u <- svd_res$u[,idx, drop = F]
  svd_res$v <- svd_res$v[,idx, drop = F]
  svd_res$d <- svd_res$d[idx]
  
  svd_res
}

.compute_matrix_mean <- function(mat, mean_vec){
  if(length(mean_vec) == 1 && !is.null(mean_vec)){
    if(mean_vec){
      mean_vec <- Matrix::colMeans(mat)
    } else{
      mean_vec <- NULL
    }
  }
  
  mean_vec
}

.compute_matrix_sd <- function(mat, sd_vec){
  if(length(sd_vec) == 1 && !is.null(sd_vec)){
    if(sd_vec){
      if(inherits(x = mat, what = c('dgCMatrix', 'dgTMatrix'))){
        sd_vec <- sparseMatrixStats::colSds(mat)
      } else {
        sd_vec <- matrixStats::colSds(mat)
      }
    } else{
      sd_vec <- NULL
    }
  }
  
  sd_vec
}

#########################################

.form_multiSVD_param <- function(center_1, center_2,
                                 n,
                                 normalize_row,
                                 normalize_singular_value,
                                 recenter_1, recenter_2,
                                 rescale_1, rescale_2,
                                 scale_1, scale_2){
  list(svd_center_1 = center_1, svd_center_2 = center_2,
       svd_n = n,
       svd_normalize_row = normalize_row,
       svd_normalize_singular_value = normalize_singular_value,
       svd_recenter_1 = recenter_1, svd_recenter_2 = recenter_2,
       svd_rescale_1 = rescale_1, svd_rescale_2 = rescale_2,
       svd_scale_1 = scale_1, svd_scale_2 = scale_2)
}