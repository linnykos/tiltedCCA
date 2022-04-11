#' Tilted-CCA Factorization
#'
#' @param input_obj \code{multiSVD} object
#' @param discretization_gridsize positive integer for how many values between 0 and 1 (inclusive) to search the 
#'                                appropriate amount of tilt over
#' @param enforce_boundary boolean, on whether or not the tilt is required to stay between
#'                         the two canonical score vectors                               
#' @param fix_tilt_perc boolean or a numeric. If \code{FALSE}, then the tilt is adaptively
#'                     determined, and if \code{TRUE}, then the tilt is set to be equal to 
#'                     \code{0.5}. If numeric, the value should be between \code{0} and \code{1},
#'                     which the tilt will be set to.
#' @param verbose non-negative integer             
#' 
#' For the tilt values (possibly set in \code{fix_tilt_perc}),
#' values close to 0 or 1 means the common space resembles the 
#' canonical scores of \code{mat_2} or \code{mat_1} respectively.
#' 
#' @return object of class \code{multiSVD}
#' @export
tiltedCCA <- function(input_obj,
                      discretization_gridsize = 21, 
                      enforce_boundary = F,
                      fix_tilt_perc = F, 
                      verbose = 0){
  stopifnot(inherits(input_obj, "multiSVD"))
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Gathering relevant objects"))
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  svd_1 <- .get_SVD(input_obj)
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  svd_2 <- .get_SVD(input_obj)
  
  n <- nrow(svd_1$u)
  metacell_clustering_list <- .get_metacell(input_obj,
                                       resolution = "cell", 
                                       type = "list", 
                                       what = "metacell_clustering")
  if(!all(is.null(metacell_clustering_list))){
    averaging_mat <- .generate_averaging_matrix(metacell_clustering_list = metacell_clustering_list,
                                                n = n)
  } else {
    averaging_mat <- NULL
  }
  
  target_dimred <- .get_Laplacian(input_obj, bool_common = T)
  param <- .get_param(input_obj)
  snn_bool_cosine <- param$snn_bool_cosine
  snn_bool_intersect <- param$snn_bool_intersect
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
                                 snn_bool_cosine = snn_bool_cosine,
                                 snn_bool_intersect = snn_bool_intersect,
                                 snn_k = snn_k,
                                 snn_min_deg = snn_min_deg,
                                 snn_num_neigh = snn_num_neigh,
                                 svd_1 = svd_1, 
                                 svd_2 = svd_2, 
                                 target_dimred = target_dimred,
                                 verbose = verbose)
  
  cca_obj <- .create_cca_obj(cca_obj = res$cca_obj, 
                             score_1 = res$score_1, 
                             score_2 = res$score_2)
  tcca_obj <- .create_tcca_obj(common_basis = res$common_basis,
                               common_score = res$common_score, 
                               distinct_score_1 = res$distinct_score_1,
                               distinct_score_2 = res$distinct_score_2,
                               df_percentage = res$df_percentage,
                               tilt_perc = res$tilt_perc)
  input_obj$cca_obj <- cca_obj
  input_obj$tcca_obj <- tcca_obj
  
  param <- .form_tiltedcca_param(discretization_gridsize = discretization_gridsize, 
                                 enforce_boundary = enforce_boundary,
                                 fix_tilt_perc = fix_tilt_perc)
  param <- .combine_two_named_lists(.get_param(input_obj), param)
  input_obj$param <- param
  
  input_obj
}

#' Tilted-CCA Decomposition
#'
#' @param input_obj output from \code{tiltedCCA}
#' @param bool_modality_1_full boolean
#' @param bool_modality_2_full boolean
#' @param verbose non-negative integer
#'
#' @return object of class \code{multiSVD}
#' @export
tiltedCCA_decomposition <- function(input_obj, 
                                    bool_modality_1_full = T,
                                    bool_modality_2_full = T,
                                    verbose = 0){
  stopifnot(inherits(input_obj, "multiSVD"))
  
  common_score <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_score")
  rank_c <- ncol(common_score)
  n <- nrow(common_score)
  
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  svd_1 <- .get_SVD(input_obj)
  score_1 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "score")
  distinct_score_1 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_score")
 
   if(bool_modality_1_full){
     mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_1$v)
     coef_mat_1 <- crossprod(score_1, mat_1)/n
   } else {
     coef_mat_1 <- crossprod(score_1, .mult_mat_vec(svd_1$u, svd_1$d))/n
   }
  common_mat_1 <- common_score %*% coef_mat_1[1:rank_c,,drop = F]
  distinct_mat_1 <- distinct_score_1 %*% coef_mat_1
     
  if(bool_modality_1_full){
    input_obj$common_mat_1 <- common_mat_1
    input_obj$distinct_mat_1 <- distinct_mat_1
  } else {
    input_obj$common_dimred_1 <- common_mat_1
    input_obj$distinct_dimred_1 <- distinct_mat_1
  }
  
  ####
  
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  svd_2 <- .get_SVD(input_obj)
  score_2 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "score")
  distinct_score_2 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_score")
  if(bool_modality_2_full){
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v)
    coef_mat_2 <- crossprod(score_2, mat_2)/n
  } else {
    coef_mat_2 <- crossprod(score_2, .mult_mat_vec(svd_2$u, svd_2$d))/n
  }
  common_mat_2 <- common_score %*% coef_mat_2[1:rank_c,,drop = F]
  distinct_mat_2 <- distinct_score_2 %*% coef_mat_2
  
  if(bool_modality_2_full){
    input_obj$common_mat_2 <- common_mat_2
    input_obj$distinct_mat_2 <- distinct_mat_2
  } else {
    input_obj$common_dimred_2 <- common_mat_2
    input_obj$distinct_dimred_2 <- distinct_mat_2
  }
  
  input_obj
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

##################################

.create_cca_obj <- function(cca_obj,
                            score_1,
                            score_2){
  structure(list(cca_obj = cca_obj, 
                 score_1 = score_1, 
                 score_2 = score_2),
            class = "cca")
}

.create_tcca_obj <- function(common_basis,
                             common_score,
                             distinct_score_1,
                             distinct_score_2,
                             df_percentage,
                             tilt_perc){
  structure(list(common_basis = common_basis,
                 common_score = common_score, 
                 distinct_score_1 = distinct_score_1,
                 distinct_score_2 = distinct_score_2,
                 df_percentage = df_percentage,
                 tilt_perc = tilt_perc),
            class = "tcca")
}

.form_tiltedcca_param <- function(discretization_gridsize, 
                                  enforce_boundary,
                                  fix_tilt_perc){
  list(tcca_discretization_gridsize = discretization_gridsize, 
       tcca_enforce_boundary = enforce_boundary,
       tcca_fix_tilt_perc = fix_tilt_perc)
}
