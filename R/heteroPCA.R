# from Zhang, Anru, T. Tony Cai, and Yihong Wu. "Heteroskedastic PCA: Algorithm, optimality, and applications." arXiv preprint arXiv:1810.08316 (2018).

#' HeteroPCA
#'
#' @param mat data matrix
#' @param K desired rank
#' @param max_iter numeric
#' @param tol numeric
#' @param verbose boolean
#'
#' @return covariance matrix
#' @export
heteroPCA <- function(mat, K, max_iter = 25, tol = 1e-4, verbose = T){
  stopifnot(min(dim(mat)) >= K, K >= 1)
  cov_mat <- stats::cov(mat)
  diag(cov_mat) <- 0
  sing_vec <- rep(Inf, K)
  
  iter <- 1
  while(iter < max_iter){
    if(verbose & iter %% floor(max_iter/10) == 0) cat('*')
    svd_res <- .svd_truncated(cov_mat, K)
    if(all(abs(sing_vec - svd_res$d)/svd_res$d < tol)) break()
    
    tmp_mat <- tcrossprod(.mult_mat_vec(svd_res$u, svd_res$d), svd_res$v)
    diag(cov_mat) <- diag(tmp_mat)
    
    sing_vec <- svd_res$d
    iter <- iter + 1
  }
  
  if(verbose) print("heteroPCA: Final iteration")
  res <- .svd_truncated(cov_mat, K, symmetric = T)
  tcrossprod(.mult_mat_vec(res$u, res$d), res$v)
}