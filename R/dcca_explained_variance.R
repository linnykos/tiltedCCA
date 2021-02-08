#' D-CCA Variance decomposition
#'
#' @param dcca_res Output of \code{dcca_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return object of class \code{dcca_decomp2}
#' @export
dcca_variance_decomposition <- function(dcca_res, rank_c, verbose = T){
  stopifnot(class(dcca_res) == "dcca")
  n <- nrow(dcca_res$svd_1$u)
  full_rank <- length(dcca_res$cca_obj)
  
  if(verbose) print("D-CCA: Form denoised observation matrices")
  mat_1 <- tcrossprod(.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d) , dcca_res$svd_1$v)
  mat_2 <- tcrossprod(.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d) , dcca_res$svd_2$v)
  
  if(verbose) print("D-CCA: Computing common percentage")
  common_perc_1 <- .common_percentage(dcca_res$score_1, dcca_res$common_score)
  common_perc_2 <- .common_percentage(dcca_res$score_1, dcca_res$common_score)
  
  if(verbose) print("D-CCA: Computing matrices")
  coef_mat_1 <- crossprod(dcca_res$score_1, mat_1)/n
  coef_mat_2 <- crossprod(dcca_res$score_2, mat_2)/n
  
  common_mat_1 <- common_perc_1 %*% coef_mat_1
  common_mat_2 <- common_perc_2 %*% coef_mat_2
 
  if(verbose) print("D-CCA: Done")
  structure(list(common_perc_1 = common_perc_1, common_perc_2 = common_perc_1,
                 common_mat_1 = common_mat_1, common_mat_2 = common_mat_2,
                 mat_1 = mat_1, mat_2 = mat_2), class = "dcca_decomp2")
}

#' Weights for cells and variables
#'
#' @param dcca_decomp2 Output of \code{dcca_variance_decomposition}
#' @param weight_func function for scalars
#' @param verbose boolean
#'
#' @return 4 vectors as a list
#' @export
explained_variance <- function(dcca_decomp2, weight_func = function(x){x},
                               verbose = T){
  stopifnot(class(dcca_decomp2) == "dcca_decomp2")
  n <- nrow(dcca_decomp2$mat_1)
  p1 <- ncol(dcca_decomp2$mat_1); p2 <- ncol(dcca_decomp2$mat_2)
  
  if(verbose) print("D-CCA: Computing sample weights")
  cell_weight_vec1 <- sapply(1:n, function(i){
    val1 <- weight_func(.l2norm(dcca_decomp2$common_mat_1[i,]))
    val2 <- weight_func(.l2norm(dcca_decomp2$mat_1[i,] - dcca_decomp2$common_mat_1[i,]))
    val2/(val1+val2)
  })
  
  cell_weight_vec2 <- sapply(1:n, function(i){
    val1 <- weight_func(.l2norm(dcca_decomp2$common_mat_2[i,]))
    val2 <- weight_func(.l2norm(dcca_decomp2$mat_2[i,] - dcca_decomp2$common_mat_2[i,]))
    val2/(val1+val2)
  })
  
  if(verbose) print("D-CCA: Computing variable weights")
  var_weight_vec1 <- sapply(1:p1, function(j){
    val1 <- weight_func(.l2norm(dcca_decomp2$common_mat_1[,j]))
    val2 <- weight_func(.l2norm(dcca_decomp2$mat_1[,j] - dcca_decomp2$common_mat_1[,j]))
    val2/(val1+val2)
  })
  
  var_weight_vec2 <- sapply(1:p2, function(j){
    val1 <- weight_func(.l2norm(dcca_decomp2$common_mat_2[,j]))
    val2 <- weight_func(.l2norm(dcca_decomp2$mat_2[,j] - dcca_decomp2$common_mat_2[,j]))
    val2/(val1+val2)
  })
  
  if(verbose) print("D-CCA: Done")
  structure(list(cell_weight_vec1 = cell_weight_vec1,
                 cell_weight_vec2 = cell_weight_vec2,
                 var_weight_vec1 = var_weight_vec1,
                 var_weight_vec2 = var_weight_vec2), class = "dcca_variance")
}

#####################################

.common_percentage <- function(score_mat, common_score_mat){
  stopifnot(nrow(score_mat) == nrow(common_score_mat), ncol(score_mat) >= ncol(common_score_mat))
  n <- nrow(score_mat)
  
  rank_1 <- ncol(score_mat)
  rank_c <- ncol(common_score_mat)
  mat <- matrix(0, n, rank_1)
  for(j in 1:rank_c){
    tmp <- .orthogonal_vec2vec(common_score_mat[,j], score_mat[,j])
    mat[,j] <- common_score_mat[,j] - tmp
  }
  
  mat
}
