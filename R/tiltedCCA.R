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
tiltedCCA <- function(input_obj,
                      discretization_gridsize = 21, 
                      enforce_boundary = F,
                      fix_tilt_perc = F, 
                      verbose = 0){
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Gathering relevant objects"))
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  svd_1 <- .get_SVD(input_obj)
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  svd_2 <- .get_SVD(input_obj)
  
  metacell_clustering <- .get_metacell(input_obj,
                                       resolution = "cell", 
                                       type = "list", 
                                       what = "metacell_clustering")
  if(!all(is.null(metacell_clustering))){
    averaging_mat <- .generate_averaging_matrix(n, metacell_clustering)
  } else {
    averaging_mat <- NULL
  }
  
  target_dimred <- .get_laplacian(input_obj, bool_common = T)
  param <- .get_param(input_obj)
  snn_k <- param$snn_latent_k
  snn_min_deg <- param$snn_min_deg
  snn_num_neigh <- param$snn_num_neigh

  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing CCA"))
  cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  
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
  res$averaging_mat <- averaging_mat
  
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
  structure(list(cca_obj = tiltedCCA_res$cca_obj, 
                 common_basis = tiltedCCA_res$common_basis,
                 common_mat_1 = common_mat_1, 
                 common_mat_2 = common_mat_2, 
                 common_score = tiltedCCA_res$common_score[,1:rank_c, drop = F],
                 df_percentage = tiltedCCA_res$df_percentage,
                 distinct_mat_1 = distinct_mat_1, 
                 distinct_mat_2 = distinct_mat_2,
                 distinct_score_1 = tiltedCCA_res$distinct_score_1,
                 distinct_score_2 = tiltedCCA_res$distinct_score_2,
                 param_list = tiltedCCA_res$param_list,
                 score_1 = tiltedCCA_res$score_1, 
                 score_2 = tiltedCCA_res$score_2, 
                 svd_1 = tiltedCCA_res$svd_1, 
                 svd_2 = tiltedCCA_res$svd_2, 
                 target_dimred = tiltedCCA_res$target_dimred,
                 tilt_perc = tiltedCCA_res$tilt_perc), class = "tiltedCCA_decomp")
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
