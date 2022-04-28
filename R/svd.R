create_multiSVD <- function(mat_1, mat_2,
                            dims_1, dims_2,
                            center_1 = T, center_2 = T,
                            normalize_row = T,
                            normalize_singular_value = T,
                            recenter_1 = F, recenter_2 = F,
                            rescale_1 = F, rescale_2 = F,
                            scale_max_1 = NULL, scale_max_2 = NULL, 
                            scale_1 = T, scale_2 = T,
                            verbose = 0){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1)
  svd_1 <- .get_SVD(center = center_1, input_obj = mat_1,
                    dims = dims_1, scale = scale_1, scale_max = scale_max_1)
  svd_2 <- .get_SVD(center = center_2, input_obj = mat_2,
                    dims = dims_2, scale = scale_2, scale_max = scale_max_2)
  
  param <- .form_multiSVD_param(center_1 = center_1, center_2 = center_2,
                                dims_1 = dims_1, dims_2 = dims_2,
                                n = n,
                                normalize_row = normalize_row,
                                normalize_singular_value = normalize_singular_value,
                                recenter_1 = recenter_1, recenter_2 = recenter_2,
                                rescale_1 = rescale_1, rescale_2 = rescale_2,
                                scale_1 = scale_1, scale_2 = scale_2,
                                scale_max_1 = scale_max_1, scale_max_2 = scale_max_2)
  
  structure(list(svd_1 = svd_1, svd_2 = svd_2,
                 default_assay = 1,
                 param = param),
            class = "multiSVD")
}

###############

.normalize_svd <- function(input_obj,
                           averaging_mat,
                           normalize_row,
                           normalize_singular_value,
                           recenter,
                           rescale, 
                           tol = 1e-4,
                           ...){
  
  dimred <- .get_Dimred(input_obj = input_obj, 
                        normalize_singular_value = normalize_singular_value,
                        ...)
  n <- nrow(dimred)
  
  if(recenter | rescale) {
    dimred <- scale(dimred, center = recenter, scale = rescale)
  }
  
  if(!all(is.null(averaging_mat))){
    dimred <- as.matrix(averaging_mat %*% dimred)
  }
  
  if(normalize_row){
    l2_vec <- apply(dimred, 1, function(x){.l2norm(x)})
    l2_vec[l2_vec <= tol] <- tol
    dimred <- .mult_vec_mat(1/l2_vec, dimred)
  }
  
  dimred
}

####################################

.svd_safe <- function(mat,
                      check_stability, # boolean
                      K, # positive integer
                      mean_vec, # boolean, NULL or vector
                      rescale, # boolean
                      scale_max, # NULL or positive integer
                      sd_vec){ # boolean, NULL or vector
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  
  mean_vec <- .compute_matrix_mean(mat, mean_vec)
  sd_vec <- .compute_matrix_sd(mat, sd_vec)
  
  res <- .svd_in_sequence(check_stability = check_stability,
                          K = K,
                          mat = mat, 
                          mean_vec = mean_vec,
                          scale_max = scale_max,
                          sd_vec = sd_vec)
  res <- list(d = res$d, u = res$u, v = res$v, method = res$method)
  class(res) <- "svd"
  
  # pass row-names and column-names
  res <- .append_rowcolnames(bool_colnames = T, 
                             bool_rownames = T,
                             source_obj = mat, 
                             target_obj = res)
  
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

.svd_in_sequence <- function(check_stability,
                             K,
                             mat, 
                             mean_vec,
                             scale_max,
                             sd_vec){
  res <- tryCatch({
    .irlba_custom(check_stability = check_stability,
                  K = K,
                  mat = mat, 
                  mean_vec = mean_vec,
                  scale_max = scale_max,
                  sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "irlba"; return(res)}
  
  ##
  
  res <- tryCatch({
    .rpsectra_custom(check_stability = check_stability,
                     K = K,
                     mat = mat, 
                     mean_vec = mean_vec,
                     scale_max = scale_max,
                     sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "RSpectra"; return(res)}
  
  ##
  
  if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
  if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
  if(!is.null(scale_max)){
    mat[mat > abs(scale_max)] <- abs(scale_max)
    mat[mat < -abs(scale_max)] <- -abs(scale_max)
  }
  res <- svd(mat)
  res$method <- "base"
  res
}

.irlba_custom <- function(check_stability,
                          K,
                          mat, 
                          mean_vec,
                          scale_max,
                          sd_vec){
  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using irlba")
    tmp <- irlba::irlba(A = mat,
                        nv = K,
                        work = min(c(K + 10, dim(mat))),
                        scale = sd_vec,
                        center = mean_vec)
    
    if(check_stability & K > 5) {
      tmp2 <- irlba::irlba(A = mat,
                           nv = 5,
                           scale = sd_vec,
                           center = mean_vec)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }
    
    return(tmp)
    
  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }
    
    tmp <- irlba::irlba(A = mat, nv = K)
    
    if(check_stability & K > 5) {
      tmp2 <- irlba::irlba(A = mat, nv = 5)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }
    
    return(tmp)
  }
}

.rpsectra_custom <- function(check_stability,
                             K,
                             mat,
                             mean_vec,
                             scale_max,
                             sd_vec){
  
  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(mean_vec))) warning("mean_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(sd_vec))) warning("sd_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using RSpectra")
  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }
  }
  
  tmp <- RSpectra::svds(A = mat, k = K)
  
  if(check_stability & K > 5) {
    tmp2 <- RSpectra::svds(A = mat, k = 5)
    ratio_vec <- tmp2$d/tmp$d[1:5]
    if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("RSpectra is potentially unstable")
  }
  
  tmp
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
                                 dims_1, dims_2,
                                 n,
                                 normalize_row,
                                 normalize_singular_value,
                                 recenter_1, recenter_2,
                                 rescale_1, rescale_2,
                                 scale_1, scale_2,
                                 scale_max_1, scale_max_2){
  list(svd_center_1 = center_1, svd_center_2 = center_2,
       svd_dims_1 = dims_1, svd_dims_2 = dims_2,
       svd_n = n,
       svd_normalize_row = normalize_row,
       svd_normalize_singular_value = normalize_singular_value,
       svd_recenter_1 = recenter_1, svd_recenter_2 = recenter_2,
       svd_rescale_1 = rescale_1, svd_rescale_2 = rescale_2,
       svd_scale_1 = scale_1, svd_scale_2 = scale_2,
       svd_scale_max_1 = scale_max_1, svd_scale_max_2 = scale_max_2)
}