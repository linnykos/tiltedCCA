#' Main workhorse of dcca_factor
#' 
#' Given the two matrices (given by \code{svd_1} and \code{svd_2}) and the
#' CCA solution in \code{cca_res}, compute the common scores.
#' This calls the functions
#' \code{.common_decomposition} and \code{.compute_distinct_score}. 
#'
#' @param averaging_mat sparse matrix
#' @param cca_res returned object from \code{.cca}
#' @param discretization_gridsize positive integer for how many values between 0 and 1 (inclusive) to search the 
#'                                appropriate amount of tilt over
#' @param enforce_boundary boolean, on whether or not the tilt is required to stay between
#'                         the two canonical score vectors                               
#' @param fix_tilt_perc boolean or a numeric. If \code{FALSE}, then the tilt is adaptively
#'                     determined, and if \code{TRUE}, then the tilt is set to be equal to 
#'                     \code{0.5}. If numeric, the value should be between \code{0} and \code{1},
#'                     which the tilt will be set to.
#' @param snn_bool_cosine boolean
#' @param snn_bool_intersect boolean
#' @param snn_k integer
#' @param snn_min_deg integer
#' @param snn_num_neigh integer
#' @param svd_1 SVD of the denoised variant of \code{mat_1} from \code{dcca_factor}
#' @param svd_2 SVD of the denoised variant of \code{mat_2} from \code{dcca_factor}
#' @param target_dimred matrix
#' @param verbose non-negative integer
#'
#' @return list 
.tiltedCCA_common_score <- function(averaging_mat,
                                    cca_res, 
                                    discretization_gridsize,
                                    enforce_boundary,
                                    fix_tilt_perc, 
                                    snn_bool_cosine,
                                    snn_bool_intersect,
                                    snn_k,
                                    snn_min_deg,
                                    snn_num_neigh,
                                    svd_1, 
                                    svd_2, 
                                    target_dimred,
                                    verbose = 0){
  full_rank <- length(cca_res$obj_vec)
  n <- nrow(svd_1$u)
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing unnormalized scores"))
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
            nrow(score_1) == nrow(score_2))
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing common factors"))
  obj_vec <- cca_res$obj_vec
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing discrete tilt"))
  tmp <- .common_decomposition(averaging_mat = averaging_mat,
                               discretization_gridsize = discretization_gridsize,
                               enforce_boundary = enforce_boundary,
                               fix_tilt_perc = fix_tilt_perc,
                               score_1 = score_1,
                               score_2 = score_2,
                               snn_bool_cosine = snn_bool_cosine,
                               snn_bool_intersect = snn_bool_intersect,
                               snn_k = snn_k,
                               snn_min_deg = snn_min_deg,
                               snn_num_neigh = snn_num_neigh,
                               svd_1 = svd_1, 
                               svd_2 = svd_2,
                               target_dimred = target_dimred,
                               verbose = verbose)
  
  common_basis <- tmp$common_basis
  common_score <- tmp$common_score
  tilt_perc <- tmp$tilt_perc
  df_percentage <- tmp$df_percentage
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Done"))
  list(common_basis = common_basis,
       common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2,
       score_1 = score_1, score_2 = score_2, 
       svd_1 = svd_1, svd_2 = svd_2, 
       cca_obj = obj_vec, 
       df_percentage = df_percentage,
       target_dimred = target_dimred,
       tilt_perc = tilt_perc
  )
}
