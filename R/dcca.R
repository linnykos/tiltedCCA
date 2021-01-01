# from Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.

#' D-CCA Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param rank_1 desired rank of data matrix 1
#' @param rank_2 desired rank of data matrix 1
#' @param enforce_rank boolean 
#' @param verbose boolean
#'
#' @return list
#' @export
dcca_factor <- function(mat_1, mat_2, rank_1, rank_2,
                 enforce_rank = T, verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
  n <- nrow(mat_1)
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  if(verbose) print("D-CCA: Starting matrix shrinkage")
  if(enforce_rank | nrow(mat_1) < 2*ncol(mat_1)) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
  if(enforce_rank | nrow(mat_2) < 2*ncol(mat_2)) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  if(verbose) print("D-CCA: Computing CCA")
  cca_res <- .cca(svd_1, svd_2)
  
  if(verbose) print("D-CCA: Computing unnormalized scores")
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  if(verbose) print("D-CCA: Computing common factors")
  R_vec <- sapply(cca_res$obj_vec, function(x){1-sqrt((1-x)/(1+x))})
  common_factors <- .mult_mat_vec((score_1+score_2)/2, R_vec)
  
  if(verbose) print("D-CCA: Done")
  structure(list(common_factors = common_factors, cca_obj = cca_res$obj_vec,
       rank_1 = rank_1, rank_2 = rank_2, 
       score_1 = score_1, score_2 = score_2, 
       svd_1 = svd_1, svd_2 = svd_2), class = "dcca")
}

#' D-CCA Decomposition
#'
#' @param dcca_res output from \code{dcca_factor}
#' @param rank_12 desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return list
#' @export
dcca_decomposition <- function(dcca_res, rank_12, verbose = T){
  stopifnot(class(dcca_res) == "dcca")
  n <- nrow(dcca_res$svd_1$u)
  full_rank <- length(dcca_res$cca_obj)
  
  if(verbose) print("D-CCA: Form denoised observation matrices")
  mat_1 <- tcrossprod(.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d) , dcca_res$svd_1$v)
  mat_2 <- tcrossprod(.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d) , dcca_res$svd_2$v)
  
  if(verbose) print("D-CCA: Computing common matrices")
  common_mat_1 <- dcca_res$common_factors[,1:rank_12, drop = F] %*% crossprod(dcca_res$score_1[,1:rank_12, drop = F], mat_1)/n
  common_mat_2 <- dcca_res$common_factors[,1:rank_12, drop = F] %*% crossprod(dcca_res$score_2[,1:rank_12, drop = F], mat_2)/n
  
  if(verbose) print("D-CCA: Computing distinctive matrices")
  if(full_rank > rank_12){
    common_mat_1_rem <- dcca_res$common_factors[,(rank_12+1):full_rank, drop = F] %*% 
      crossprod(dcca_res$score_1[,(rank_12+1):full_rank, drop = F], mat_1)/n
    common_mat_2_rem <- dcca_res$common_factors[,(rank_12+1):full_rank, drop = F] %*% 
      crossprod(dcca_res$score_2[,(rank_12+1):full_rank, drop = F], mat_2)/n
    distinct_mat_1 <- mat_1 - common_mat_1 - common_mat_1_rem
    distinct_mat_2 <- mat_2 - common_mat_2 - common_mat_2_rem
  } else {
    distinct_mat_1 <- mat_1 - common_mat_1
    distinct_mat_2 <- mat_2 - common_mat_2
  }
  
  if(verbose) print("D-CCA: Done")
  list(common_factors = dcca_res$common_factors[,1:rank_12, drop = F],
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
  
  svd_res$d <- sing_vec
  svd_res$u <- svd_res$u[,1:K,drop = F]
  svd_res$v <- svd_res$v[,1:K,drop = F]
  
  svd_res
}

.cca <- function(svd_1, svd_2){
  stopifnot(all(svd_1$d >= 0), all(svd_2$d >= 0), nrow(svd_1$u) == nrow(svd_2$u))
  n <- nrow(svd_1$u)
  
  cov_1_invhalf <- .mult_mat_vec(svd_1$v, sqrt(n)/svd_1$d)
  cov_2_invhalf <- .mult_mat_vec(svd_2$v, sqrt(n)/svd_2$d)
                                    
  agg_mat <-  .compute_cca_aggregate_matrix(svd_1, svd_2)
  svd_res <- svd(agg_mat)
  
  list(factor_1 = cov_1_invhalf %*% svd_res$u,
       factor_2 = cov_2_invhalf %*% svd_res$v,
       obj_vec = svd_res$d)
}

.compute_cca_aggregate_matrix <- function(svd_1, svd_2){
  crossprod(svd_1$u, svd_2$u)
}

# unnoramlized means that while score_1 is orthogonal, it is not orthonormal
.compute_unnormalized_scores <- function(svd_1, svd_2, cca_res){
  score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, cca_res$factor_1) #note: this is unnormalized
  score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, cca_res$factor_2)
  
  list(score_1 = score_1, score_2 = score_2)
}
