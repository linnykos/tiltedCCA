#' Main workhorse of dcca_factor
#' 
#' Given the two matrices (given by \code{svd_1} and \code{svd_2}) and the
#' CCA solution in \code{cca_res}, compute the common scores.
#' This calls the functions
#' \code{.common_decomposition} and \code{.compute_distinct_score}. 
#'
#' @param svd_1 SVD of the denoised variant of \code{mat_1} from \code{dcca_factor}
#' @param svd_2 SVD of the denoised variant of \code{mat_2} from \code{dcca_factor}
#' @param cca_res returned object from \code{.cca}
#' @param num_neigh number of neighbors to consider to computed the common percentage
#' @param fix_tilt_perc boolean or numeric between 0 and 1 (inclusive). If \code{TRUE}, the output \code{distinct_perc_2} will be fixed to at 0.5,
#' meaning the common scores will be the "middle" of \code{score_1} and \code{score_2}.
#' If \code{FALSE}, \code{distinct_perc_2} will be adaptively estimated via the
#' \code{.common_decomposition} function.
#' @param cell_max number of cells used to compute the distinct percentaage
#' @param check_alignment boolean. If \code{TRUE}, recompute \code{score_1} and \code{score_2}
#' after using \code{.compute_unnormalized_scores}. This might be needed if the \code{.cca} solution
#' was not computed from exactly \code{svd_1} and \code{svd_2}
#' @param verbose boolean
#' @param msg character
#'
#' @return list 
.dcca_common_score <- function(cca_res, 
                               cell_max,
                               check_alignment, 
                               discretization_gridsize,
                               fix_tilt_perc, 
                               metacell_clustering_1,
                               metacell_clustering_2,
                               num_neigh,
                               svd_1, 
                               svd_2, 
                               verbose = T, msg = ""){
  stopifnot(cell_max > 10)
  full_rank <- length(cca_res$obj_vec)
  n <- nrow(svd_1$u)
  
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing unnormalized scores"))
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
            nrow(score_1) == nrow(score_2))
  
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing common factors"))
  if(check_alignment){
    # reparameterize the scores
    tmp <- .cca(score_1, score_2, dims_1 = 1:ncol(score_1), dims_2 = 1:ncol(score_2), 
                return_scores = T)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    stopifnot(is.matrix(score_1), is.matrix(score_2))
    obj_vec <- diag(crossprod(score_1, score_2))/n
  } else {
    obj_vec <- cca_res$obj_vec
  }
  
  # compute the common scores
  n <- nrow(score_1)
  if(cell_max < n){
    n_idx <- sample(1:n, size = cell_max)
  } else {
    n_idx <- 1:n
  }
  
  # [[note to self: use these n_idx somehow]]
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing discrete tilt"))
  tmp <- .common_decomposition(discretization_gridsize = discretization_gridsize,
                               fix_tilt_perc = fix_tilt_perc,
                               metacell_clustering_1 = metacell_clustering_1,
                               metacell_clustering_2 = metacell_clustering_2,
                               n_idx = n_idx,
                               num_neigh = num_neigh,
                               score_1 = score_1,
                               score_2 = score_2,
                               svd_1 = svd_1, 
                               svd_2 = svd_2,
                               verbose = verbose)
  
  common_score <- tmp$common_score
  tilt_perc <- tmp$tilt_perc
  df_percentage <- tmp$df_percentage
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Done"))
  list(common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2,
       score_1 = score_1, score_2 = score_2, 
       svd_1 = svd_1, svd_2 = svd_2, 
       cca_obj = obj_vec, 
       df_percentage = df_percentage,
       metacell_clustering_1 = metacell_clustering_1,
       metacell_clustering_2 = metacell_clustering_2,
       tilt_perc = tilt_perc
       )
}
