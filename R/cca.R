#' Perform CCA
#' 
#' This function takes either two matrices or two SVDs. Both \code{input_1}
#' and \code{input_2} must be "of the same type." Calls the \code{.compute_cca_aggregate_matrix} function.
#'
#' @param input_1 first input
#' @param input_2 second input
#' @param dims_1 desired latent dimensions of data matrix 1. 
#' Only used if \code{input_1} is a matrix, not if it's  a list representing the SVD
#' @param dims_2 desired latent dimensions of data matrix 2.
#'  Only used if \code{input_2} is a matrix, not if it's  a list representing the SVD
#' @param return_scores boolean. If \code{TRUE}, return the scores (i.e., matrices where the rows are the cells).
#' If \code{FALSE}, return the loadings (i.e., matrices where the rows are the variables).
#' Either way, one of the output matrices will have \code{rank_1} columns and another
#' will have \code{rank_2} columns
#' @param tol small numeric
#'
#' @return list
.cca <- function(input_1, input_2, dims_1, dims_2, 
                 return_scores, tol = 1e-6){
  if(is.matrix(input_1) & is.matrix(input_2)){
    rank_1 <- max(dims_1); rank_2 <- max(dims_2)
    stopifnot(nrow(input_1) == nrow(input_2), 
              rank_1 <= ncol(input_1), rank_2 <= ncol(input_2))
    svd_1 <- .svd_safe(mat = input_1,
                       check_stability = T,
                       K = rank_1,
                       mean_vec = NULL,
                       rescale = F,
                       scale_max = NULL,
                       sd_vec = NULL)
    svd_2 <- .svd_safe(mat = input_2,
                       check_stability = T,
                       K = rank_2,
                       mean_vec = NULL,
                       rescale = F,
                       scale_max = NULL,
                       sd_vec = NULL)
    
    svd_1 <- .check_svd(svd_1, dims_1); svd_2 <- .check_svd(svd_2, dims_2)
  } else {
    stopifnot(inherits(input_1, "svd"), inherits(input_2, "svd"))
    
    rank_1 <- length(which(input_1$d >= 1e-6)); rank_2 <- length(which(input_2$d >= 1e-6))
    svd_1 <- input_1; svd_2 <- input_2
  }
  
  rank_1 <- ncol(svd_1$u); rank_2 <- ncol(svd_2$u)
  n <- nrow(svd_1$u)
  
  # perform CCA
  cov_1_invhalf <- .mult_mat_vec(svd_1$v, sqrt(n)/svd_1$d)
  cov_2_invhalf <- .mult_mat_vec(svd_2$v, sqrt(n)/svd_2$d)
  
  agg_mat <-  .compute_cca_aggregate_matrix(svd_1, svd_2, augment = T)
  svd_res <- svd(agg_mat)
  full_rank <- min(c(rank_1, rank_2))
  
  loading_1 <- cov_1_invhalf %*% svd_res$u[1:ncol(cov_1_invhalf), 1:rank_1, drop = F]
  loading_2 <- cov_2_invhalf %*% svd_res$v[1:ncol(cov_2_invhalf), 1:rank_2, drop = F]
  
  # return
  if(!return_scores){
    list(loading_1 = loading_1, loading_2 = loading_2, obj_vec = svd_res$d[1:full_rank])
  } else {
    if(!is.matrix(input_1) | !is.matrix(input_2)){
      input_1 <- tcrossprod(.mult_mat_vec(input_1$u, input_1$d), input_1$v) 
      input_2 <- tcrossprod(.mult_mat_vec(input_2$u, input_2$d), input_2$v) 
    }
    
    list(score_1 = input_1 %*% loading_1, score_2 = input_2 %*% loading_2, obj_vec = svd_res$d[1:full_rank])
  }
  
}

#' Helper function with the CCA function
#' 
#' Called by the \code{.cca} function. Recall when computing CCA,
#' the main matrix we need compute is, roughly speaking, 
#' half-inverse of the first covariance times the cross-covariance matrix
#' times the half-inverse of the second covariance. If we had the SVD
#' of the two original matrices, this is actually equivalent to the product
#' of the left singular vectors.
#'
#' @param svd_1 SVD of the denoised variant of \code{mat_1} from \code{dcca_factor}
#' @param svd_2 SVD of the denoised variant of \code{mat_2} from \code{dcca_factor}
#' @param augment boolean. If \code{TRUE}, augment the matrix with either rows or columns
#' with 0's so the dimension of the output matrix matches those in \code{svd_1} and \code{svd_2}
#'
#' @return matrix
.compute_cca_aggregate_matrix <- function(svd_1, svd_2, augment){
  res <- crossprod(svd_1$u, svd_2$u)
  
  if(augment){
    if(nrow(res) < ncol(res)) res <- rbind(res, matrix(0, ncol(res)-nrow(res), ncol(res)))
    if(ncol(res) < nrow(res)) res <- cbind(res, matrix(0, nrow(res), nrow(res)-ncol(res)))
  }
  
  res
}

#' Using the CCA solution, compute the score matrices.
#' 
#' This is called by \code{.tiltedCCA_common_score}. It's called unnormalized
#' scores since if there are \code{n} rows (i.e., \code{nrow(svd_1$u)}),
#' then the return matrices will be orthogonal matrices where the
#' matrix multiplied by itself is a diagonal matrix with all values \code{n}.
#'
#' @param svd_1 SVD of the denoised variant of \code{mat_1} from \code{dcca_factor}
#' @param svd_2 SVD of the denoised variant of \code{mat_2} from \code{dcca_factor}
#' @param cca_res returned object from \code{.cca}
#'
#' @return list of two matrices
.compute_unnormalized_scores <- function(svd_1, svd_2, cca_res){
  score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, cca_res$loading_1) 
  score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, cca_res$loading_2)
  
  list(score_1 = score_1, score_2 = score_2)
}