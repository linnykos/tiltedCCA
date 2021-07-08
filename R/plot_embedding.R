#' Side-by-side UMAPs of the common, distinct and everything matrix
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec factor vector
#' @param data_1 boolean
#' @param data_2 boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct "on the same scale" by adding appropriately-scaled Gaussian noise
#' @param col_vec vector of colors
#' @param pca boolean. If \code{TRUE}, plot the PCA embedding with the leading 2 components. 
#' If \code{FALSE}, plot the UMAP embedding.
#' @param only_embedding boolean
#' @param main_addition additional string to append to main of each plot
#' @param verbose boolean
#'
#' @return depends on \code{only_embedding}. 
#' If \code{TRUE}, returns three matrices as a list.
#' If \code{FALSE}, shows a plot but returns nothing
#' @export
plot_embeddings <- function(obj, membership_vec = NA, data_1 = T, data_2 = F, 
                            add_noise = T, 
                            col_vec = scales::hue_pal()(length(levels(membership_vec))),
                            pca = F, only_embedding = F,
                            main_addition = "", verbose = F){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(obj$common_score),
            class(obj) %in% c("dcca", "dcca_decomp"))
  stopifnot(data_1 | data_2)
  
  n <- nrow(obj$common_score)
  embedding <- .prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                   add_noise = add_noise)
  if(pca) {
    embedding <- .extract_pca_embedding(embedding)
  } else {
    embedding <- .extract_umap_embedding(embedding, only_embedding = T)
  }
  
  if(!only_embedding) {
    if(all(is.na(membership_vec))){
      col_cells <- rep(col_vec[1], n)
    } else{
      col_cells <- col_vec[as.numeric(membership_vec)]
    }
    
    main_vec <- c("Common view", "Distinct view", "Everything view")
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
    xlim <- range(sapply(embedding, function(x){x[,1]}))
    ylim <- range(sapply(embedding, function(x){x[,2]}))
    n_idx <- sample(1:nrow(obj$common_score))
    
    graphics::par(mfrow = c(1,3))
    for(i in 1:3){
      graphics::plot(embedding[[i]][n_idx,1], embedding[[i]][n_idx,2],
                     asp = T, pch = 16, col = col_cells[n_idx], 
                     main = paste0(main_vec[i], main_addition),
                     xlab = label1, ylab = label2, xlim = xlim, ylim = ylim)
    }
    invisible()
  } else {
    return(embedding)
  }
}

#############################

#' Extract UMAP embedding
#'
#' @param embedding return object from \code{.prepare_embeddings}
#' @param only_embedding boolean
#' @param reduction_key string for \code{Seurat::RunUMAP}
#'
#' @return list of three 2-column matrix or \code{Seurat} object, one for the
#' common, distinct and everything matrix
.extract_umap_embedding <- function(embedding, only_embedding, reduction_key = "UMAP"){
  stopifnot(length(embedding) == 3)
  for(i in 1:3){
    if(only_embedding){
      embedding[[i]] <- Seurat::RunUMAP(embedding[[i]], metric = "euclidean", 
                                        verbose = F)@cell.embeddings
    } else {
      embedding[[i]] <- Seurat::RunUMAP(embedding[[i]], metric = "euclidean", 
                                        reduction.key = reduction_key, verbose = F)
    }
  }
  
  embedding
}

.extract_pca_embedding <- function(embedding){
  stopifnot(length(embedding) == 3)
  for(i in 1:3){
    tmp <- .svd_truncated(embedding[[i]], K = 2, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    embedding[[i]] <- .mult_mat_vec(tmp$u, tmp$d)
  }
  
  embedding
}

###################################

#' Prepare the appropriate inputs
#'
#' @param obj output from either \code{generate_data} or \code{dcca_decomposition}
#' @param data_1 boolean, for computing the embedding for data 1
#' @param data_2 boolean, for computing the embedding for data 2
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct scores "on the same scale" by adding appropriately-scaled Gaussian noise
#'
#' @return a list of 3 elements, representing the low-dimensional
#' representation for the common, distinct and everything matrices
#' for the intended dataset(s)
.prepare_embeddings <- function(obj, data_1, data_2, add_noise){
  if(data_1){
    res_1 <- .prepare_embeddings_singleton(obj$common_score, obj$distinct_score_1,
                                           obj$svd_1, add_noise = add_noise)
  } else res_1 <- NA
  if(data_2){
    res_2 <- .prepare_embeddings_singleton(obj$common_score, obj$distinct_score_2,
                                           obj$svd_2, add_noise = add_noise)
  } else res_2 <- NA
  
  if(!all(is.na(res_1)) & all(is.na(res_2))){
    return(res_1)
  } else if(all(is.na(res_1)) & !all(is.na(res_2))){
    return(res_2)
  } else {
    res <- res_1
    for(i in 1:3){
      res[[i]] <- cbind(res[[i]], res_2[[i]])
    }
    return(res)
  }
  
}

.prepare_embeddings_singleton <- function(common_score, distinct_score, svd_res, add_noise){
  embedding <- vector("list", 3)
  names(embedding) <- c("common", "distinct", "everything")
  
  embedding[[1]] <- .extract_matrix_helper(common_score, distinct_score,
                                           svd_res, common_bool = T, distinct_bool = F,
                                           center = F, renormalize = F, add_noise = add_noise)
  embedding[[2]] <- .extract_matrix_helper(common_score, distinct_score,
                                           svd_res, common_bool = F, distinct_bool = T,
                                           center = F, renormalize = F, add_noise = add_noise)
  embedding[[3]] <- .extract_matrix_helper(common_score, distinct_score,
                                           svd_res, common_bool = T, distinct_bool = T,
                                           center = F, renormalize = F, add_noise = add_noise)
  
  for(i in 1:3){
    if(length(rownames(common_score)) != 0) rownames(embedding[[i]]) <- rownames(common_score)
  }
  
  embedding
}

#' Matrix extraction helper for plotting
#' 
#' This function is designed to be called for plotting based on PCA
#' or UMAP.
#' One of \code{common_bool} and \code{distinct_bool} must be \code{TRUE}.
#' The options of \code{center} and \code{renormalize} are there
#' to emulate the typical UMAP where the cosine distance is used
#' to measure the distance between samples.
#'
#' @param common_score matrix of common scores
#' @param distinct_score matrix of distinct scores, with same number of rows as \code{common_score}
#' @param svd_e list containing the SVD of the full matrix 
#' @param common_bool boolean
#' @param distinct_bool boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct scores "on the same scale" by adding appropriately-scaled Gaussian noise.
#' Which matrices get the added noise depends on \code{common_bool} or \code{distinct_bool})
#' @param center boolean. If \code{TRUE}, center each canonical variable (i.e., column)
#' where the centering is based on \code{common_score+distinct_score}
#' @param renormalize boolean. If \code{TRUE}, normalize each sample (i.e., row)
#' to have unit norm, where the rescaling factor is based on \code{common_score+distinct_score}.
#' This is suggested to be only used if \code{center=TRUE}, in which case the centering
#' happens before the renormalization
#'
#' @return a matrix
.extract_matrix_helper <- function(common_score, distinct_score,
                                   svd_e, common_bool, distinct_bool, add_noise,
                                   center, renormalize){
  stopifnot(nrow(common_score) == nrow(distinct_score),
            nrow(common_score) == nrow(svd_e$u), nrow(distinct_score) == nrow(svd_e$u),
            ncol(common_score) <= length(svd_e$d), ncol(distinct_score) == length(svd_e$d),
            common_bool | distinct_bool)
  
  n <- nrow(common_score)
  canonical_score <- .add_two_matrices(common_score, distinct_score)
  full_mat <- .mult_mat_vec(svd_e$u, svd_e$d/max(svd_e$d))
  full_mat <- canonical_score %*% crossprod(canonical_score, full_mat)/n # reorient for consistency for the rest of the pipeline
  center_vec <- apply(full_mat, 2, mean)
  tmp <- full_mat
  if(center) tmp <- sapply(1:ncol(full_mat), function(k){full_mat[,k] - center_vec[k]})
  if(renormalize) l2_vec <- apply(tmp, 1, function(x){.l2norm(x)})
  
  if(common_bool == distinct_bool){
    tmp <- full_mat
    
  } else {
    if(common_bool){ 
      if(add_noise){
        common_score <- .add_columnwise_noise(common_score, distinct_score)
      } else{
        if(ncol(common_score) < ncol(distinct_score)) {
          common_score <- cbind(common_score, matrix(0, nrow = n, ncol = ncol(distinct_score)-ncol(common_score)))
        }
      }
      tmp <- tcrossprod(common_score, canonical_score) %*% full_mat/n
      
    } else { 
      if(add_noise){
        distinct_score <- .add_columnwise_noise(distinct_score, common_score)
      } else{
        if(ncol(distinct_score) < ncol(common_score)) {
          distinct_score <- cbind(distinct_score, matrix(0, nrow = n, ncol = ncol(common_score)-ncol(distinct_score)))
        }
      }
      tmp <- tcrossprod(distinct_score, canonical_score) %*% full_mat/n
    }
    
    # center variables
    if(center){
      tmp <- sapply(1:ncol(tmp), function(k){tmp[,k] - center_vec[k]})
    }
  }
  
  # normalize cells
  if(renormalize){.mult_vec_mat(1/l2_vec, tmp)} else {tmp}
}

# helper function for when mat1 has less columns than mat2
.add_two_matrices <- function(mat1, mat2){
  stopifnot(nrow(mat1) == nrow(mat2), ncol(mat1) <= ncol(mat2))
  if(ncol(mat1) == ncol(mat2)) return(mat1+mat2)
  mat2[,1:ncol(mat1)] <- mat1 + mat2[,1:ncol(mat1)]
  mat2
}

# helper function to add noise
.add_columnwise_noise <- function(target_mat, second_mat){
  stopifnot(nrow(target_mat) == nrow(second_mat))
  n <- nrow(target_mat)
  r1 <- ncol(target_mat); r2 <- ncol(second_mat)
  
  if(r1 < r2){
    target_mat <- cbind(target_mat, matrix(0, n, r2-r1))
  }
  
  for(j in 1:min(r1,r2)){
    sd1 <- stats::sd(target_mat[,j])
    sd2 <- stats::sd(second_mat[,j])
    if(sd1 < sd2){
      sd_add <- min(sd1, sqrt(sd2^2-sd1^2))
      target_mat[,j] <- target_mat[,j] + stats::rnorm(n, sd = sd_add)
    } 
  }
  
  target_mat
}

