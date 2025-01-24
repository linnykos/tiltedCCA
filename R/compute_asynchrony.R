#' Compute the asynchrony based on Tilted-CCA
#'
#' @param input_obj \code{multiSVD} class, after applying \code{tiltedCCA_decomposition()}
#' @param anchor_modality numeric of \code{1} or \code{2}, denoting which modality the synchrony is computed on
#'
#' @returns a matrix with \code{n} cells and 2 columns, which one column 
#' (named \code{asynchrony}) denotes the synchrony score (between 0 and 1, where a number closer to 1 means
#' the cell is synchronous between both modalities) and 
#' (named \code{asynchrony_rescaled}), which is simply a monotone transformation
#' of the first column for easier visualization
#' @export
compute_asynchrony <- function(input_obj,
                               anchor_modality = 1){
  stopifnot(anchor_modality %in% c(1,2),
            length(anchor_modality) == 1)
  
  # grab the base modality's common component
  input_obj <- .set_defaultAssay(input_obj, assay = anchor_modality)
  base_common <- .get_tCCAobj(input_obj, 
                              apply_postDimred = FALSE,
                              what = "common_mat")
  
  input_obj <- .set_defaultAssay(input_obj, assay = anchor_modality)
  svd_1 <- .get_SVD(input_obj)
  input_obj <- .set_defaultAssay(input_obj, assay = setdiff(c(1,2), anchor_modality))
  svd_2 <- .get_SVD(input_obj)
  tmp <- crossprod(svd_2$u, svd_1$u)
  svd_tmp <- svd(tmp)
  
  rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
  other_pred <- tcrossprod(.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
  n <- nrow(base_common)
  asynchrony_vec <- sapply(1:n, function(i){
    df <- data.frame(base = base_common[i,],
                     other = other_pred[i,])
    suppressWarnings(lm_res <- stats::lm(base ~ other, data = df))
    suppressWarnings(summary(lm_res)$r.squared)
  })
  
  # rescale the alignment for better visualization
  scaling_grid <- seq(0.1, 10, length.out = 100)
  rank_vec <- rank(asynchrony_vec)
  scaling_quality <- sapply(scaling_grid, function(val){
    suppressWarnings(stats::cor(asynchrony_vec^val, rank_vec, use = "complete.obs"))
  })
  if(all(is.na(scaling_quality))){
    asynchrony_vec_rescaled <- asynchrony_vec
  } else {
    asynchrony_vec_rescaled <- asynchrony_vec^(scaling_grid[which.max(scaling_quality)])
  }
  
  # add the names and then output
  res_df <- cbind(asynchrony = asynchrony_vec,
                  asynchrony_rescaled = asynchrony_vec_rescaled)
  res_df <- .append_rowcolnames(bool_colnames = FALSE,
                                bool_rownames = TRUE,
                                source_obj = base_common,
                                target_obj = res_df)
  
  return(res_df)
}