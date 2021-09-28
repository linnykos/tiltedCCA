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
#' @param fix_distinct_perc boolean. If \code{TRUE}, the output \code{distinct_perc_2} will be fixed to at 0.5,
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
.dcca_common_score <- function(svd_1, svd_2, 
                               cca_res, 
                               num_neigh, 
                               fix_distinct_perc, 
                               cell_max,
                               check_alignment, 
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
  if(fix_distinct_perc){
    tmp <- .common_decomposition(score_1, score_2, 
                                 snn_union = NA, 
                                 fix_distinct_perc = T,
                                 verbose = verbose)
  } else {
    if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing kNN"))
    n <- nrow(score_1)
    if(cell_max < n){
      n_idx <- sample(1:n, size = cell_max)
    } else {
      n_idx <- 1:n
    }
    
    rescaling_factor <- max(c(svd_1$d, svd_2$d))
    dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*rescaling_factor)
    snn_1 <- .form_snn(dimred_1[n_idx,,drop = F], num_neigh = num_neigh)
    dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*rescaling_factor)
    snn_2 <- .form_snn(dimred_2[n_idx,,drop = F], num_neigh = num_neigh)
    
    snn_union <- snn_1 + snn_2
    
    tmp <- .common_decomposition(score_1, score_2, 
                                 snn_union = snn_union,
                                 fix_distinct_perc = F, 
                                 verbose = verbose)
  }
  
  common_score <- tmp$common_score; distinct_perc_2 <- tmp$distinct_perc_2
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Done"))
  list(common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2,
       score_1 = score_1, score_2 = score_2, 
       svd_1 = svd_1, svd_2 = svd_2, 
       cca_obj = obj_vec, distinct_perc_2 = distinct_perc_2)
}

###################

.form_snn <- function(mat, num_neigh){
  stopifnot(num_neigh > 1)
  n <- nrow(mat)
  nn_mat <- RANN::nn2(mat, k = num_neigh)$nn.idx
  if(all(nn_mat[,1] == 1:n)){
    nn_mat <- nn_mat[,-1,drop = F]
  }
  
  i_vec <- rep(1:n, times = ncol(nn_mat))
  j_vec <- as.numeric(nn_mat)
  
  sparseMat <- Matrix::sparseMatrix(i = i_vec,
                                    j = j_vec,
                                    x = rep(1, length(i_vec)),
                                    repr = "C")
  sparseMat <- sparseMat * Matrix::t(sparseMat)
  sparseMat
}