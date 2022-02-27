#' Tilted-CCA Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param dims_1 desired latent dimensions of data matrix 1 (increasing integer vector of length 2)
#' @param dims_2 desired latent dimensions of data matrix 2 (increasing integer vector of length 2)
#' @param center_1 boolean, on whether or not to center \code{mat_1} prior to SVD
#' @param center_2 boolean, on whether or not to center \code{mat_2} prior to SVD
#' @param scale_1 boolean, on whether or not to rescale \code{mat_1} prior to SVD
#' @param scale_2 boolean, on whether or not to rescale \code{mat_2} prior to SVD
#' @param cell_max positive integer for how many cells to subsample (useful if  
#'                 \code{nrow(mat_1) is too large})
#' @param discretization_gridsize positive integer for how many values between 0 and 1 (inclusive) to search the 
#'                                appropriate amount of tilt over
#' @param enforce_boundary boolean, on whether or not the tilt is required to stay between
#'                         the two canonical score vectors                               
#' @param fix_tilt_perc boolean or a numeric. If \code{FALSE}, then the tilt is adaptively
#'                     determined, and if \code{TRUE}, then the tilt is set to be equal to 
#'                     \code{0.5}. If numeric, the value should be between \code{0} and \code{1},
#'                     which the tilt will be set to.
#' @param num_neigh positive integer for how many NNs are used to construct the relevant
#'                  NN graphs. 
#' @param verbose boolean                
#'                              
#' The \code{cell_max} parameter is used specifically to limit the 
#' number of cells involved in \code{.common_decomposition} (for 
#' less cells to be involved in \code{.determine_cluster})
#' 
#' For the tilt values (possibly set in \code{fix_tilt_perc}),
#' values close to 0 or 1 means the common space resembles the 
#' canonical scores of \code{mat_2} or \code{mat_1} respectively.
#' 
#' \code{metacell_clustering_1} and \code{metacell_clustering_2} need to be
#' both either \code{NA} or factor vectors. If the former (i.e.,
#' \code{metacell_clustering_1=NA} and \code{metacell_clustering_2=NA}),
#' then the appropriate amount of tilt is determined by Jaccard dissimilarity
#' based on the NN structure. If the latter (i.e.,
#' \code{is.factor(metacell_clustering_1)=TRUE} and \code{is.factor(metacell_clustering_2)=TRUE}),
#' then the appropriate amount of tilt is determined by KL-divergences of
#' each cell's NNs' factor proportions. These are all done in \code{.determine_cluster}.
#'
#' \code{num_neigh} is used to construct the Shared NN graphs (in \code{.form_snns})
#' for \code{mat_1} and \code{mat_2}
#' as well as in \code{.determine_cluster} (to call \code{.form_snns} to construct
#' the analogous Shared NN graph for the common space).
#'
#' @return list of class \code{dcca}
#' @export
tiltedCCA <- function(mat_1, mat_2, dims_1, dims_2, 
                      target_dimred,
                      center_1 = T, center_2 = T,
                      scale_1 = T, scale_2 = T,
                      discretization_gridsize = 21, 
                      enforce_boundary = F,
                      fix_tilt_perc = F, 
                      metacell_clustering = lapply(1:nrow(mat_1), function(i){i}),
                      snn_bool_intersect = F,
                      snn_k = 20,
                      snn_min_deg = 0,
                      snn_num_neigh = 30,
                      verbose = T){
  rank_1 <- max(dims_1); rank_2 <- max(dims_2)
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
  
  n <- nrow(mat_1)
  
  if(verbose) print(paste0(Sys.time(),": Tilted-CCA: Starting matrix shrinkage"))
  svd_1 <- .svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                          mean_vec = center_1, sd_vec = scale_1, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                          mean_vec = center_2, sd_vec = scale_2, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = dims_1)
  svd_2 <- .check_svd(svd_2, dims = dims_2)
  
  msg <- " (all cells)"
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing CCA"))
  cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  averaging_mat <- .generate_averaging_matrix(n, metacell_clustering)
  
  res <- .tiltedCCA_common_score(averaging_mat = averaging_mat,
                                 cca_res = cca_res, 
                                 discretization_gridsize = discretization_gridsize,
                                 enforce_boundary = enforce_boundary,
                                 fix_tilt_perc = fix_tilt_perc, 
                                 snn_bool_intersect = snn_bool_intersect,
                                 snn_k = snn_k,
                                 snn_min_deg = snn_min_deg,
                                 snn_num_neigh = snn_num_neigh,
                                 svd_1 = svd_1, 
                                 svd_2 = svd_2, 
                                 target_dimred = target_dimred,
                                 verbose = verbose, msg = msg)
  
  param_list <- list(center_1 = center_1, center_2 = center_2,
                     scale_1 = scale_1, scale_2 = scale_2,
                     discretization_gridsize = discretization_gridsize, 
                     enforce_boundary = enforce_boundary,
                     fix_tilt_perc = fix_tilt_perc,
                     snn_bool_intersect = snn_bool_intersect,
                     snn_k = snn_k,
                     snn_min_deg = snn_min_deg,
                     snn_num_neigh = snn_num_neigh)
  res$param_list <- param_list
  res$metacell_clustering <- metacell_clustering
  
  class(res) <- "tiltedCCA"
  res
}

#' Tilted-CCA Decomposition
#'
#' @param tiltedCCA_res output from \code{dcca_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return list of class \code{dcca_decomp}
#' @export
tiltedCCA_decomposition <- function(tiltedCCA_res, rank_c = NA, verbose = T){
  stopifnot(class(tiltedCCA_res) == "tiltedCCA")
  
  if(is.na(rank_c)) rank_c <- min(c(ncol(tiltedCCA_res$distinct_score_1), ncol(tiltedCCA_res$distinct_score_2)))
  stopifnot( rank_c <= min(c(ncol(tiltedCCA_res$distinct_score_1), ncol(tiltedCCA_res$distinct_score_2))))
  n <- nrow(tiltedCCA_res$svd_1$u)
  
  if(verbose) print(paste0(Sys.time(),": Tilted-CCA: Form denoised observation matrices"))
  mat_1 <- tcrossprod(.mult_mat_vec(tiltedCCA_res$svd_1$u, tiltedCCA_res$svd_1$d) , tiltedCCA_res$svd_1$v)
  mat_2 <- tcrossprod(.mult_mat_vec(tiltedCCA_res$svd_2$u, tiltedCCA_res$svd_2$d) , tiltedCCA_res$svd_2$v)
  
  if(verbose) print(paste0(Sys.time(),": Tilted-CCA: Computing common matrices"))
  coef_mat_1 <- crossprod(tiltedCCA_res$score_1, mat_1)/n
  coef_mat_2 <- crossprod(tiltedCCA_res$score_2, mat_2)/n
  
  common_mat_1 <- tiltedCCA_res$common_score[,1:rank_c, drop = F] %*% coef_mat_1[1:rank_c,,drop = F]
  common_mat_2 <- tiltedCCA_res$common_score[,1:rank_c, drop = F] %*% coef_mat_2[1:rank_c,,drop = F]
  
  if(verbose) print(paste0(Sys.time(),": Tilted-CCA: Computing distinctive matrices"))
  distinct_mat_1 <- tiltedCCA_res$distinct_score_1 %*% coef_mat_1
  distinct_mat_2 <- tiltedCCA_res$distinct_score_2 %*% coef_mat_2
  
  if(verbose) print(paste0(Sys.time(),": Tilted-CCA: Done"))
  structure(list(common_score = tiltedCCA_res$common_score[,1:rank_c, drop = F],
                 distinct_score_1 = tiltedCCA_res$distinct_score_1,
                 distinct_score_2 = tiltedCCA_res$distinct_score_2,
                 score_1 = tiltedCCA_res$score_1, 
                 score_2 = tiltedCCA_res$score_2,
                 svd_1 = tiltedCCA_res$svd_1, 
                 svd_2 = tiltedCCA_res$svd_2, 
                 common_mat_1 = common_mat_1, 
                 common_mat_2 = common_mat_2, 
                 distinct_mat_1 = distinct_mat_1, 
                 distinct_mat_2 = distinct_mat_2,
                 cca_obj = tiltedCCA_res$cca_obj, 
                 distinct_perc_2 = tiltedCCA_res$distinct_perc_2,
                 param_list = tiltedCCA_res$param_list), class = "dcca_decomp")
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

.generate_averaging_matrix <- function(n, subcluster_list){
  stopifnot(max(unlist(subcluster_list)) <= n)
  
  k <- length(subcluster_list)
  i_vec <- unlist(lapply(1:k, function(i){
    len <- length(subcluster_list[[i]]); rep(i, len)
  }))
  j_vec <- unlist(subcluster_list)
  x_vec <- unlist(lapply(subcluster_list, function(vec){
    len <- length(vec); rep(1/len, len)
  }))
  
  Matrix::sparseMatrix(i = i_vec,
                       j = j_vec,
                       x = x_vec,
                       dims = c(k,n),
                       repr = "C")
}
