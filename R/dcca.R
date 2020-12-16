# from Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.

#' D-CCA
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param rank_1 desired rank of data matrix 1
#' @param rank_2 desired rank of data matrix 1
#' @param rank_12 desired rank of cross-covariance matrix
#' @param enforce_rank boolean 
#' @param verbose boolean
#'
#' @return list
#' @export
dcca <- function(mat_1, mat_2, rank_1, rank_2, rank_12, enforce_rank = T, verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2), rank_12 <= min(c(rank_1, rank_2)), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
  n <- nrow(mat_1)
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  if(verbose) print("D-CCA: Starting matrix shrinkage")
  if(enforce_rank | nrow(mat_1) < 2*ncol(mat_1)) mat_1 <- .spoet(mat_1, rank_1)
  if(enforce_rank | nrow(mat_2) < 2*ncol(mat_2)) mat_2 <- .spoet(mat_2, rank_2)
  
  if(verbose) print("D-CCA: Computing covariance matrices")
  cov_1 <- stats::cov(mat_1) * (n-1)/n
  cov_2 <- stats::cov(mat_2) * (n-1)/n
  cov_12 <- crossprod(mat_1, mat_2)/n
  full_rank <- Matrix::rankMatrix(cov_12)
  
  if(verbose) print("D-CCA: Computing CCA")
  cca_res <- .cca(cov_1, cov_2, cov_12, K = full_rank)
  score_1 <- mat_1 %*% cca_res$factor_1 #note: this is unnormalized
  score_2 <- mat_2 %*% cca_res$factor_2
  
  R_vec <- sapply(cca_res$obj_vec, function(x){1-sqrt((1-x)/(1+x))})
  
  if(verbose) print("D-CCA: Computing common matrices")
  common_factors <- .mult_mat_vec((score_1+score_2)/2, R_vec)
  common_mat_1 <- common_factors[,1:rank_12, drop = F] %*% crossprod(score_1[,1:rank_12, drop = F], mat_1)/n
  common_mat_2 <- common_factors[,1:rank_12, drop = F] %*% crossprod(score_2[,1:rank_12, drop = F], mat_2)/n
  
  if(verbose) print("D-CCA: Computing distinctive matrices")
  if(full_rank > rank_12){
    common_mat_1_rem <- common_factors[,(rank_12+1):full_rank, drop = F] %*% crossprod(score_1[,(rank_12+1):full_rank, drop = F], mat_1)/n
    common_mat_2_rem <- common_factors[,(rank_12+1):full_rank, drop = F] %*% crossprod(score_2[,(rank_12+1):full_rank, drop = F], mat_2)/n
    distinct_mat_1 <- mat_1 - common_mat_1 - common_mat_1_rem
    distinct_mat_2 <- mat_2 - common_mat_2 - common_mat_2_rem
  } else {
    distinct_mat_1 <- mat_1 - common_mat_1
    distinct_mat_2 <- mat_2 - common_mat_2
  }
  
  list(common_factors = common_factors, cca_obj = cca_res$obj_vec,
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2,
       distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2)
}

#################################

.spoet <- function(mat, K){
  n <- nrow(mat); p <- ncol(mat); m <- min(n, p)
  target_full_dim <- min(c(nrow(mat), ncol(mat), K*10))
  svd_res <- .svd_truncated(mat, target_full_dim)
  tau <- sum((svd_res$d[(K+1):length(svd_res$d)])^2)/(n*p - n*K - p*K)
  sing_vec <- sapply(svd_res$d[1:K], function(x){
    sqrt(max(c(x^2-tau*p, 0)))
  })
  
  tcrossprod(svd_res$u[,1:K] %*% .diag_matrix(sing_vec), svd_res$v[,1:K])
}

.cca <- function(cov_1, cov_2, cov_12, K){
  stopifnot(nrow(cov_1) == ncol(cov_1), nrow(cov_2) == ncol(cov_2),
            nrow(cov_1) == nrow(cov_12), nrow(cov_2) == ncol(cov_12))
  
  svd_1 <- .svd_truncated(cov_1, nrow(cov_1))
  svd_2 <- .svd_truncated(cov_2, nrow(cov_2))
  stopifnot(all(svd_1$d >= 0), all(svd_2$d >= 0))
  
  cov_1_invhalf <- .inverse_onehalf(cov_1)
  cov_2_invhalf <- .inverse_onehalf(cov_2)
  
  agg_mat <-  crossprod(cov_1_invhalf, cov_12) %*% cov_2_invhalf
  svd_res <- .svd_truncated(agg_mat, K)
  
  list(factor_1 = cov_1_invhalf %*% svd_res$u,
       factor_2 = cov_2_invhalf %*% svd_res$v,
       obj_vec = svd_res$d)
}
