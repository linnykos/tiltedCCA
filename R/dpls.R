#' D-PLS Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param rank_1 desired rank of data matrix 1
#' @param rank_2 desired rank of data matrix 2
#' @param apply_shrinkage boolean 
#' @param verbose boolean
#'
#' @return list of class \code{dcca}
#' @export
dpls_factor <- function(mat_1, mat_2, rank_1, rank_2, apply_shrinkage = T,
                        verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)), rank_1 == rank_2)
  n <- nrow(mat_1)
  if(verbose) print("D-PLS: Rescaling matrices")
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  if(verbose) print(paste0("D-PLS: Starting matrix shrinkage"))
  if(apply_shrinkage) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
  if(apply_shrinkage) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)

  # alternatively, apply D-CCA to all cells
  msg <- " (all cells)"
  if(verbose) print(paste0("D-PLS", msg, ": Computing PLS"))
  pls_res <- .pls(svd_1, svd_2, rank = min(rank_1, rank_2))
  
  res <- .pls_common_score(svd_1, svd_2, pls_res,
                           verbose = verbose, msg = msg)
  
  class(res) <- "dpls"
  res
}

#' D-PLS Decomposition
#'
#' @param dpls_res output from \code{dpls_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return list of class \code{dpls_decomp}
#' @export
dcca_decomposition <- function(dpls_res, verbose = T){
  stopifnot(class(dcca_res) == "dpls")
  n <- nrow(dpls_res$svd_1$u)
  full_rank <- ncol(dpls_res$pls_res$u)
  
  if(verbose) print("D-PLS: Form coefficient matrices")
  # WARNING: this seems to current assume that the ranks are equal
  coef_1 <- tcrossprod(solve(crossprod(svd_1$v, pls_res$pls_res$u)), svd_1$v)
  coef_2 <- tcrossprod(solve(crossprod(svd_2$v, pls_res$pls_res$v)), svd_2$v)
  
  if(verbose) print("D-PLS: Computing common matrices")
  common_mat_1 <- pls_res$common_mat_1 %*% coef_1
  common_mat_2 <- pls_res$common_mat_2 %*% coef_2
  
  if(verbose) print("D-PLS: Computing distinctive matrices")
  distinct_mat_1 <- (pls_res$score_1 - pls_res$common_mat_1) %*% coef_1
  distinct_mat_2 <- (pls_res$score_2 - pls_res$common_mat_2) %*% coef_1
  
  if(verbose) print("D-PLS: Done")
  structure(list(common_mat_1 = common_mat_1, common_mat_2 = common_mat_2, 
                 distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
                 pls_obj = pls_res$pls_res$d), class = "dpls_decomp")
}

#############################

.pls <- function(svd_1, svd_2, rank){
  .svd_truncated(crossprod(.mult_mat_vec(svd_1$v, svd_1$d), svd_1$u) %*% tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v), K = rank)
}

.pls_common_score <- function(svd_1, svd_2, pls_res,
                              verbose = T, msg = ""){
  n <- nrow(svd_1$u); K <- ncol(pls_res$u)
  score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, pls_res$u)
  score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, pls_res$v)
  
  common_mat_1 <- matrix(NA, nrow = n, ncol = K)
  common_mat_2 <- matrix(NA, nrow = n, ncol = K)
  
  for(j in 1:K){
    basis_res <- .representation_2d(score_1[,j], score_2[,j])
    tmp_decomp <- .decomposition_2d(basis_res$rep1, basis_res$rep2, plotting = F)
    common_mat_1[,j] <- basis_res$basis_mat %*% tmp_decomp$common_vec1
    common_mat_2[,j] <- basis_res$basis_mat %*% tmp_decomp$common_vec2
  }
  
  list(svd_1 = svd_1, svd_2 = svd_2, pls_res = pls_res,
       score_1 = score_2, score_2 = score_2,
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2)
}
