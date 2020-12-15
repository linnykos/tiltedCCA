# from Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.

dcca <- function(mat_1, mat_2, rank_1, rank_2, rank_12){
  mat_1 <- scale(mat_1, center = T)
  mat_2 <- scale(mat_2, center = T)
  
  mat_1 <- .spoet(mat_1, rank_1)
  mat_2 <- .spoet(mat_2, rank_2)
  
  cov_1 <- stats::cov(mat_1)
  cov_2 <- stats::cov(mat_2)
  cov_12 <- crossprod(mat_1, mat_2)
  
  cca_res <- .cca(cov_1, cov_2, cov_12, rank_12)
  score_1 <- mat_1 %*% cca_res$factor_1
  score_2 <- mat_2 %*% cca_res$factor_2
  
  R_diag <- diag(sapply(cca_res$obj_vec, function(x){1-sqrt((1-x)/(1+x))}))
  
  common_factors <- (score_1+score_2)/2 %*% R_diag
  common_mat_1 <- common_factors %*% crossprod(score_1, mat_1)
  common_mat_2 <- common_factors %*% crossprod(score_2, mat_2)
  distinct_mat_1 <- mat_1 - common_mat_1
  distinct_mat_2 <- mat_2 - common_mat_2
  
  list(common_factors = common_factors, 
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
    sqrt(max(x^2-tau*p), 0)
  })
  
  svd_res$u %*% diag(sing_vec) %*% t(svd_res$v)
}

.cca <- function(cov_1, cov_2, cov_12, K){
  stopifnot(nrow(cov_1) == ncol(cov_1), nrow(cov_2) == ncol(cov_2),
            nrow(cov_1) == nrow(cov_12), nrow(cov_2) == ncol(cov_12))
  
  svd_1 <- .svd_truncated(cov_1, nrow(cov_1))
  svd_2 <- .svd_truncated(cov_2, nrow(cov_2))
  stopifnot(all(svd_1$d >= 0), all(svd_2$d >= 0))
  
  cov_1_invhalf <- svd_1$u %*% diag(svd_1$d^(-1/2))
  cov_2_invhalf <- svd_2$u %*% diag(svd_2$d^(-1/2))
  
  agg_mat <- t(cov_1_invhalf) %*% cov_12 %*% cov_2_invhalf
  svd_res <- .svd_truncated(agg_mat, K)
  
  list(factor_1 = cov_1_invhalf %*% svd_res$u,
       factor_2 = cov_2_invhalf %*% svd_res$v,
       obj_vec <- svd_res$d)
}
