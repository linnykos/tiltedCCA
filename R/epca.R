# from Liu, Lydia T., Edgar Dobriban, and Amit Singer. "$ e $ PCA: High dimensional exponential family PCA." The Annals of Applied Statistics 12.4 (2018): 2121-2150.

epca <- function(mat, K, mean_var_func = function(x){x}){
  stopifnot(K >= 1, all(dim(mat) > K+2))
  
  # compute mean and cov
  p <- ncol(mat); n <- nrow(mat)
  mean_vec <- colMeans(mat); cov_mat <- stats::cov(mat)
  
  # whiten and diagonally debias
  diag_vec <- sapply(mean_vec, mean_var_func)
  stopifnot(all(diag_vec > 0))
  tmp <- (1/diag_vec)^(1/2)
  cov_mat <- .mult_mat_vec(.mult_vec_mat(tmp, cov_mat), tmp)
  diag(cov_mat) <- diag(cov_mat) - 1
  
  # shrink eigen
  gamma <- p/n
  svd_res <- .svd_truncated(cov_mat, K)
  shrunk_eig <- sapply(svd_res$d, function(x){
    if(x > (1+gamma)^2){
      .5*(x+1-gamma+sqrt((x+1-gamma)^2-4*x))
    } else {
      x + sqrt(gamma)
    }
  })
  cov_mat <- svd_res$u %*% diag(shrunk_eig) %*% t(svd_res$v)
  
  # recolor
  tmp <- diag_vec^(1/2)
  eigen(.mult_mat_vec(.mult_vec_mat(tmp, cov_mat), tmp))$vectors[,1:K]
}