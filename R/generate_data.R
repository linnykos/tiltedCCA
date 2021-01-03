#' Generate data
#'
#' @param common_loading \code{n} by \code{rank_12} orthogonal matrix
#' @param distinct_loading_1 \code{n} by \code{rank_1} orthogonal matrix where \code{rank_1} is larger than \code{rank_12}
#' @param distinct_loading_2 \code{n} by \code{rank_2} orthogonal matrix where \code{rank_2} is larger than \code{rank_12}
#' @param coef_mat_1 \code{rank_1} by \code{p_1} matrix
#' @param coef_mat_2 \code{rank_2} by \code{p_2} matrix
#' @param noise_func function that takes in a matrix and outputs a matrix of the same dimension
#'
#' @return list of two matrices, one of dimension \code{n} by \code{p_1} and another of dimension \code{n} by \code{p_2}
#' @export
generate_data <- function(common_loading, distinct_loading_1, distinct_loading_2,
                          coef_mat_1, coef_mat_2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}){
  stopifnot(nrow(common_loading) == nrow(distinct_loading_1), nrow(common_loading) == nrow(distinct_loading_2),
            ncol(common_loading) <= ncol(distinct_loading_1), ncol(common_loading) <= ncol(distinct_loading_2),
            ncol(distinct_loading_1) == nrow(coef_mat_1), ncol(distinct_loading_2) == nrow(coef_mat_2))
  tmp1 <- crossprod(common_loading); tmp2 <- crossprod(distinct_loading_1); tmp3 <- crossprod(distinct_loading_2)
  stopifnot(abs(sum(abs(tmp1)) - sum(abs(diag(tmp1)))) <= 1e-6,  
            abs(sum(abs(tmp2)) - sum(abs(diag(tmp2)))) <= 1e-6,
            abs(sum(abs(tmp3)) - sum(abs(diag(tmp3)))) <= 1e-6)
  
  rank_12 <- ncol(common_loading)
  mat_1 <- common_loading %*% coef_mat_1[1:rank_12,,drop = F] + distinct_loading_1 %*% coef_mat_1
  mat_2 <- common_loading %*% coef_mat_2[1:rank_12,,drop = F] + distinct_loading_2 %*% coef_mat_2
  
  mat_1 <- noise_func(mat_1); mat_2 <- noise_func(mat_2)
  
  list(mat_1 = mat_1, mat_2 = mat_2)
}