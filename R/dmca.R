#' D-MCA Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param rank_1 desired rank of data matrix 1
#' @param rank_2 desired rank of data matrix 2
#' @param apply_shrinkage boolean 
#' @param verbose boolean
#'
#' @return list of class \code{dmca}
#' @export
dmca_factor <- function(mat_1, mat_2, rank_1, rank_2, apply_shrinkage = T,
                        verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)), rank_1 == rank_2)
  n <- nrow(mat_1)
  if(verbose) print("D-MCA: Rescaling matrices")
  mat_1 <- scale(mat_1, center = T, scale = T)
  mat_2 <- scale(mat_2, center = T, scale = T)
  
  if(verbose) print(paste0("D-MCA: Starting matrix shrinkage"))
  if(apply_shrinkage) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
  if(apply_shrinkage) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)

  # alternatively, apply D-CCA to all cells
  msg <- " (all cells)"
  if(verbose) print(paste0("D-MCA", msg, ": Computing MCA"))
  mca_res <- .mca(svd_1, svd_2, rank = min(rank_1, rank_2))
  
  res <- .mca_common_score(svd_1, svd_2, mca_res,
                           verbose = verbose, msg = msg)
  
  class(res) <- "dmca"
  res
}

#' D-MCA Decomposition
#'
#' @param dmca_res output from \code{dmca_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dmca_factor}
#' @param verbose boolean
#'
#' @return list of class \code{dmca_decomp}
#' @export
dmca_decomposition <- function(dmca_res, verbose = T){
  stopifnot(class(dmca_res) == "dmca")
  n <- nrow(dmca_res$svd_1$u)
  full_rank <- ncol(dmca_res$mca_res$u)
  
  if(verbose) print("D-MCA: Form coefficient matrices")
  # WARNING: this seems to current assume that the ranks are equal
  coef_1 <- .mca_coef_mat(dmca_res$svd_1$v, dmca_res$mca_res$u)
  coef_2 <- .mca_coef_mat(dmca_res$svd_2$v, dmca_res$mca_res$v)
  
  if(verbose) print("D-MCA: Computing common matrices")
  common_mat_1 <- dmca_res$common_mat_1 %*% coef_1
  common_mat_2 <- dmca_res$common_mat_2 %*% coef_2
  
  if(verbose) print("D-MCA: Computing distinctive matrices")
  distinct_score_1 <- dmca_res$score_1 - dmca_res$common_mat_1
  distinct_score_2 <- dmca_res$score_2 - dmca_res$common_mat_2
  distinct_mat_1 <- distinct_score_1 %*% coef_1
  distinct_mat_2 <- distinct_score_2 %*% coef_2
  
  if(verbose) print("D-MCA: Done")
  structure(list(common_mat_1 = common_mat_1, common_mat_2 = common_mat_2, 
                 distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
                 common_score_1 = dmca_res$common_mat_1, common_score_2 = dmca_res$common_mat_2,
                 distinct_score_1 = distinct_score_1, distinct_score_2 = distinct_score_2,
                 mca_obj = dmca_res$mca_res$d), class = "dmca_decomp")
}

#############################

.mca <- function(svd_1, svd_2, rank){
  .svd_truncated(tcrossprod(svd_1$v %*% crossprod(.mult_mat_vec(svd_1$u, svd_1$d), .mult_mat_vec(svd_2$u, svd_2$d)), svd_2$v), K = rank)
}

.mca_common_score <- function(svd_1, svd_2, mca_res,
                              verbose = T, msg = ""){
  n <- nrow(svd_1$u); K <- ncol(mca_res$u)
  score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, mca_res$u)
  score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, mca_res$v)
  
  common_mat_1 <- matrix(NA, nrow = n, ncol = K)
  common_mat_2 <- matrix(NA, nrow = n, ncol = K)
  
  for(j in 1:K){
    basis_res <- .representation_2d(score_1[,j], score_2[,j])
    stopifnot(all(basis_res$rep1 >= 0), all(basis_res$rep2 >= 0))
    
    tmp_decomp <- .decomposition_2d(basis_res$rep1, basis_res$rep2, plotting = F)
    common_mat_1[,j] <- basis_res$basis_mat %*% tmp_decomp$common_vec1
    common_mat_2[,j] <- basis_res$basis_mat %*% tmp_decomp$common_vec2
  }
  
  if(length(rownames(svd_1$u)) != 0) rownames(common_mat_1) <- rownames(svd_1$u)
  if(length(colnames(svd_2$u)) != 0) rownames(common_mat_2) <- colnames(svd_2$u)
  
  list(svd_1 = svd_1, svd_2 = svd_2, mca_res = mca_res,
       score_1 = score_1, score_2 = score_2,
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2)
}

.mca_coef_mat <- function(singular_mat, loading_mat){
  stopifnot(nrow(singular_mat) == nrow(loading_mat),
            ncol(singular_mat) == ncol(loading_mat))
  
  tcrossprod(solve(crossprod(singular_mat, loading_mat)), singular_mat)
}
