#' Computing Consensus PCA
#'
#' @param mat_1 data matrix of \code{n} cells and \code{p1} features for Modality 1
#' @param mat_2 data matrix of \code{n} cells and \code{p2} features for Modality 2
#' @param dims_1 vector of latent dimensions for \code{mat_1} used for analysis
#' @param dims_2 vector of latent dimensions for \code{mat_2} used for analysis
#' @param dims_consensus vector of latent dimensions for the Consensus PCA
#' @param apply_pca boolean, where PCA is combined set of latent dimensions (from both modalities) if \code{T}
#' or not if \code{F}
#' @param center_1 boolean, to center each the feature in Modality 1 prior to computing latent dimensions
#' @param center_2 boolean, to center each the feature in Modality 2 prior to computing latent dimensions
#' @param center_consensus boolean, to center the combined set of latent dimensions (from both modalities)
#' @param normalize_row boolean, to normalize each cell's latent vector after dimension-reduction for both modalities
#' @param normalize_singular_value boolean, to normalize each modality by its largest singular value
#' @param recenter_1 boolean, to center each latent dimension in Modality 1 after computing latent dimensions
#' @param recenter_2 boolean, to center each latent dimension in Modality 2 after computing latent dimensions
#' @param recenter_consensus boolean, to center the latent dimension after Consensus PCA
#' @param rescale_1 boolean, to rescale each latent dimension in Modality 1 after computing latent dimensions
#' @param rescale_2 boolean, to rescale each latent dimension in Modality 2 after computing latent dimensions
#' @param rescale_consensus boolean, to rescale the latent dimension after Consensus PCA
#' @param scale_1 boolean, to rescale each the feature in Modality 1 prior to computing latent dimensions
#' @param scale_2 boolean, to rescale each the feature in Modality 2 prior to computing latent dimensions
#' @param scale_consensus boolean, to rescale the combined set of latent dimensions (from both modalities)
#' @param scale_max_1 numeric or \code{NULL}, to threshold Modality 1 in magnitude prior to computing latent dimensions
#' @param scale_max_2 numeric or \code{NULL}, to threshold Modality 2 in magnitude prior to computing latent dimensions
#' @param scale_max_consensus  numeric or \code{NULL}, to threshold the combined set of latent dimensions in magnitude prior to computing latent dimensions
#' @param svd_1 list of \code{u}, \code{d}, \code{v} for the SVD of Modality 1 if it's already computed
#' @param svd_2 list of \code{u}, \code{d}, \code{v} for the SVD of Modality 2 if it's already computed
#' @param tol small positive number
#' @param verbose non-negative integer             
#'
#' @return \code{consensusPCA} object
#' @export
consensus_pca <- function(mat_1, mat_2,
                          dims_1, dims_2,
                          dims_consensus,
                          apply_pca = T,
                          center_1 = T, center_2 = T,
                          center_consensus = T,
                          normalize_row = F,
                          normalize_singular_value = T,
                          recenter_1 = F, recenter_2 = F,
                          recenter_consensus = F,
                          rescale_1 = F, rescale_2 = F,
                          rescale_consensus = F,
                          scale_1 = T, scale_2 = T,
                          scale_consensus = T,
                          scale_max_1 = NULL, scale_max_2 = NULL, 
                          scale_max_consensus = NULL,
                          svd_1 = NULL, svd_2 = NULL,
                          tol = 1e-3,
                          verbose = 0){
  
  if(verbose > 0) print("Extracting SVDs")
  if(all(is.null(svd_1))) svd_1 <- .get_SVD(center = center_1, input_obj = mat_1,
                                            dims = dims_1, scale = scale_1, scale_max = scale_max_1)
  if(all(is.null(svd_2))) svd_2 <- .get_SVD(center = center_2, input_obj = mat_2,
                                            dims = dims_2, scale = scale_2, scale_max = scale_max_2)
  
  n <- nrow(svd_1$u)
  stopifnot(n == nrow(svd_2$u))
  param <- .form_consensusPCA_param(apply_pca = apply_pca,
                                    center_1 = center_1, center_2 = center_2,
                                    center_consensus = center_consensus,
                                    dims_1 = dims_1, dims_2 = dims_2,
                                    dims_consensus = dims_consensus,
                                    n = n,
                                    normalize_row = normalize_row,
                                    normalize_singular_value = normalize_singular_value,
                                    recenter_1 = recenter_1, recenter_2 = recenter_2,
                                    recenter_consensus = recenter_consensus,
                                    rescale_1 = rescale_1, rescale_2 = rescale_2,
                                    rescale_consensus = rescale_consensus,
                                    scale_1 = scale_1, scale_2 = scale_2,
                                    scale_consensus = scale_consensus,
                                    scale_max_1 = scale_max_1, scale_max_2 = scale_max_2,
                                    scale_max_consensus = scale_max_consensus)
  
  if(verbose > 0) print("Computing PCAs")
  dimred_1 <- .normalize_svd(input_obj = svd_1,
                             averaging_mat = NULL,
                             normalize_row = normalize_row,
                             normalize_singular_value = normalize_singular_value,
                             recenter = recenter_1,
                             rescale = rescale_1)
  dimred_2 <- .normalize_svd(input_obj = svd_2,
                             averaging_mat = NULL,
                             normalize_row = normalize_row,
                             normalize_singular_value = normalize_singular_value,
                             recenter = recenter_2,
                             rescale = rescale_2,
                             tol = 1e-4)
  
  if(verbose > 0) print("Computing Consensus PCA")
  stopifnot(nrow(dimred_1) == nrow(dimred_2))
  dimred_combined <- cbind(dimred_1, dimred_2)
  
  if(apply_pca){
    svd_consensus <- .get_SVD(center = center_consensus, 
                              input_obj = dimred_combined,
                              dims = dims_consensus, 
                              scale = scale_consensus, 
                              scale_max = scale_max_consensus)
    dimred_consensus <- .normalize_svd(input_obj = svd_consensus,
                                       averaging_mat = NULL,
                                       normalize_row = normalize_row,
                                       normalize_singular_value = normalize_singular_value,
                                       recenter = recenter_consensus,
                                       rescale = rescale_consensus)
  } else {
    if(recenter_consensus | rescale_consensus) {
      dimred_combined <- scale(dimred_combined, center = recenter_consensus, scale = rescale_consensus)
    }
    
    if(normalize_row){
      l2_vec <- apply(dimred_combined, 1, function(x){.l2norm(x)})
      l2_vec[l2_vec <= tol] <- tol
      dimred_combined <- .mult_vec_mat(1/l2_vec, dimred_combined)
    }
  }
  
  structure(list(dimred_consensus = dimred_consensus,
                 dimred_1 = dimred_1,
                 dimred_2 = dimred_2,
                 param = param),
            class = "consensusPCA")
}

.form_consensusPCA_param <- function(apply_pca,
                                     center_1, center_2,
                                     center_consensus,
                                     dims_1, dims_2,
                                     dims_consensus,
                                     n,
                                     normalize_row,
                                     normalize_singular_value,
                                     recenter_1, recenter_2,
                                     recenter_consensus,
                                     rescale_1, rescale_2,
                                     rescale_consensus,
                                     scale_1, scale_2,
                                     scale_consensus,
                                     scale_max_1, scale_max_2,
                                     scale_max_consensus){
  list(svd_apply_pca = apply_pca,
       svd_center_1 = center_1, svd_center_2 = center_2,
       svd_center_consensus = center_consensus,
       svd_dims_1 = dims_1, svd_dims_2 = dims_2,
       svd_dims_consensus = dims_consensus,
       svd_n = n,
       svd_normalize_row = normalize_row,
       svd_normalize_singular_value = normalize_singular_value,
       svd_recenter_1 = recenter_1, svd_recenter_2 = recenter_2,
       svd_recenter_consensus = recenter_consensus,
       svd_rescale_1 = rescale_1, svd_rescale_2 = rescale_2,
       svd_rescale_consensus = rescale_consensus,
       svd_scale_1 = scale_1, svd_scale_2 = scale_2,
       svd_scale_consensus = scale_consensus,
       svd_scale_max_1 = scale_max_1, svd_scale_max_2 = scale_max_2,
       svd_scale_max_consensus = scale_max_consensus)
}