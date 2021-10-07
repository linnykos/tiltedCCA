# from Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.

#' D-CCA Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param dims_1 desired latent dimensions of data matrix 1
#' @param dims_2 desired latent dimensions of data matrix 2
#' @param center_1 boolean, on whether or not to center \code{mat_1} prior to SVD
#' @param center_2 boolean, on whether or not to center \code{mat_2} prior to SVD
#' @param scale_1 boolean, on whether or not to rescale \code{mat_1} prior to SVD
#' @param scale_2 boolean, on whether or not to rescale \code{mat_2} prior to SVD
#' @param metacell_clustering optional clustering
#' @param cell_max number of cells used to compute the distinct percentaage
#' @param fix_distinct_perc boolean. If \code{TRUE}, the output \code{distinct_perc_2} will be fixed to at 0.5,
#' meaning the common scores will be the "middle" of \code{score_1} and \code{score_2}.
#' If \code{FALSE}, \code{distinct_perc_2} will be adaptively estimated via the
#' \code{.common_decomposition} function.
#' @param verbose boolean
#'
#' @return list of class \code{dcca}
#' @export
dcca_factor <- function(mat_1, mat_2, dims_1, dims_2, 
                        center_1 = T, center_2 = T,
                        scale_1 = T, scale_2 = T,
                        cell_max = nrow(mat_1),
                        clustering_resolution = 1,
                        discretization_gridsize = 9, 
                        fix_tilt_perc = F, 
                        form_meta_matrix = F,
                        metacell_clustering = NA,
                        num_neigh = min(30, round(nrow(mat_1)/20)),
                        verbose = T){
  rank_1 <- max(dims_1); rank_2 <- max(dims_2)
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
  
  n <- nrow(mat_1)
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: Starting matrix shrinkage"))
  svd_1 <- .svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                          mean_vec = center_1, sd_vec = scale_1, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                          mean_vec = center_2, sd_vec = scale_2, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = dims_1)
  svd_2 <- .check_svd(svd_2, dims = dims_2)
  
  if(all(is.na(metacell_clustering))){
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Constructing meta-cells"))
    metacell_clustering <- form_metacells(svd_1, svd_2, 
                                      clustering_resolution = clustering_resolution,
                                      dims_1 = NA, dims_2 = NA,
                                      center_1 = center_1, center_2 = center_2,
                                      scale_1 = scale_1, scale_2 = scale_2,
                                      verbose = verbose)
  }
  
  if(form_meta_matrix){
    # apply D-CCA to meta-cells
    msg <- " (meta-cells)"
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Constructing meta-cells matrix 1"))
    mat_1_meta <- .compute_metamatrix(mat_1, metacell_clustering, verbose)
    
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Constructing meta-cells matrix 2"))
    mat_2_meta <- .compute_metamatrix(mat_2, metacell_clustering, verbose)
    
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing CCA"))
    cca_res <- .cca(mat_1_meta, mat_2_meta, dims_1 = dims_1, dims_2 = dims_2,
                    return_scores = F)
  } else {
    # alternatively, apply D-CCA to all cells
    msg <- " (all cells)"
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing CCA"))
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  }
  
  res <- .dcca_common_score(cca_res = cca_res, 
                            cell_max = cell_max,
                            check_alignment = all(!is.na(metacell_clustering)), 
                            discretization_gridsize = discretization_gridsize,
                            fix_tilt_perc = fix_tilt_perc, 
                            metacell_clustering = metacell_clustering,
                            num_neigh = num_neigh,
                            svd_1 = svd_1, 
                            svd_2 = svd_2, 
                            verbose = verbose, msg = msg)
  
  class(res) <- "dcca"
  res
}

#' D-CCA Decomposition
#'
#' @param dcca_res output from \code{dcca_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return list of class \code{dcca_decomp}
#' @export
dcca_decomposition <- function(dcca_res, rank_c = NA, verbose = T){
  stopifnot(class(dcca_res) == "dcca")
  
  if(is.na(rank_c)) rank_c <- min(c(ncol(dcca_res$distinct_score_1), ncol(dcca_res$distinct_score_2)))
  stopifnot( rank_c <= min(c(ncol(dcca_res$distinct_score_1), ncol(dcca_res$distinct_score_2))))
  n <- nrow(dcca_res$svd_1$u)
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: Form denoised observation matrices"))
  mat_1 <- tcrossprod(.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d) , dcca_res$svd_1$v)
  mat_2 <- tcrossprod(.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d) , dcca_res$svd_2$v)
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: Computing common matrices"))
  coef_mat_1 <- crossprod(dcca_res$score_1, mat_1)/n
  coef_mat_2 <- crossprod(dcca_res$score_2, mat_2)/n
  
  common_mat_1 <- dcca_res$common_score[,1:rank_c, drop = F] %*% coef_mat_1[1:rank_c,,drop = F]
  common_mat_2 <- dcca_res$common_score[,1:rank_c, drop = F] %*% coef_mat_2[1:rank_c,,drop = F]
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: Computing distinctive matrices"))
  distinct_mat_1 <- dcca_res$distinct_score_1 %*% coef_mat_1
  distinct_mat_2 <- dcca_res$distinct_score_2 %*% coef_mat_2
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: Done"))
  structure(list(common_score = dcca_res$common_score[,1:rank_c, drop = F],
                 distinct_score_1 = dcca_res$distinct_score_1,
                 distinct_score_2 = dcca_res$distinct_score_2,
                 score_1 = dcca_res$score_1, score_2 = dcca_res$score_2,
                 svd_1 = dcca_res$svd_1, svd_2 = dcca_res$svd_2, 
                 common_mat_1 = common_mat_1, common_mat_2 = common_mat_2, 
                 distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
                 cca_obj = dcca_res$cca_obj, distinct_perc_2 = dcca_res$distinct_perc_2), class = "dcca_decomp")
}


#################################

#' Compute the distinct scores
#' 
#' Given \code{score_1} and \code{score_2}, and having already
#' computed \code{common_score}, compute the distinct scores.
#' This is more-or-less a simple subtraction, but we need to handle situations
#' where we might need to "fill-in extra dimensions"
#'
#' @param score_1 matrix
#' @param score_2 matrix
#' @param common_score matrix
#'
#' @return list of two matrices
.compute_distinct_score <- function(score_1, score_2, common_score){
  full_rank <- ncol(common_score)
  distinct_score_1 <- score_1[,1:full_rank, drop = F] - common_score
  distinct_score_2 <- score_2[,1:full_rank, drop = F] - common_score
  
  # handle different ranks
  distinct_score_1 <- cbind(distinct_score_1, score_1[,-(1:full_rank), drop = F])
  distinct_score_2 <- cbind(distinct_score_2, score_2[,-(1:full_rank), drop = F])
  
  list(distinct_score_1 = distinct_score_1, distinct_score_2 = distinct_score_2)
}

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
    svd_1 <- .svd_truncated(input_1, rank_1, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    svd_2 <- .svd_truncated(input_2, rank_2, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    
    svd_1 <- .check_svd(svd_1, dims_1); svd_2 <- .check_svd(svd_2, dims_2)
  } else {
    stopifnot(is.list(input_1), is.list(input_2), 
              all(c("u","d","v") %in% names(input_1)),
              all(c("u","d","v") %in% names(input_2)),
              all(input_1$d >= 0), all(input_2$d >= 0), 
              all(diff(input_1$d) <= 1e-6), all(diff(input_2$d) <= 1e-6),
              nrow(input_1$u) == nrow(input_2$u), all(is.na(dims_1)), all(is.na(dims_2)))
    
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
#' This is called by \code{.dcca_common_score}. It's called unnormalized
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
