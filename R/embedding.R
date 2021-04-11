#' Extract UMAP embedding
#'
#' @param prep_list return object from \code{.prepare_umap_embedding}
#' @param common_1 boolean
#' @param common_2 boolean
#' @param distinct_1 boolean
#' @param distinct_2 boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct scores "on the same scale" by adding appropriately-scaled Gaussian noise
#' @param only_embedding boolean
#' @param reduction_key string for \code{Seurat::RunUMAP}
#'
#' @return 2-column matrix or \code{Seurat} object
.extract_umap_embedding <- function(prep_list, common_1 = T, common_2 = T,
                              distinct_1 = T, distinct_2 = T,
                              add_noise = T,
                              only_embedding = T, reduction_key = "UMAP"){
  stopifnot(common_1 | common_2 | distinct_1 | distinct_2)
  
  prep_list$svd_list$e1$d <- prep_list$svd_list$e1$d/prep_list$svd_list$e1$d[1]
  prep_list$svd_list$e2$d <- prep_list$svd_list$e2$d/prep_list$svd_list$e2$d[1]
 
  ###########
  
  tmp_list <- vector("list", 1); l <- 1
  tmp <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_1, 
                                svd_e = prep_list$svd_list$e1, 
                                common_bool = common_1, distinct_bool = distinct_1, 
                                add_noise = add_noise, center = T, renormalize = T)
  if(!all(is.na(tmp))) { tmp_list[[l]] <- tmp; l <- l+1 }
  tmp <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_2, 
                                svd_e = prep_list$svd_list$e2, 
                                common_bool = common_2, distinct_bool = distinct_2, 
                                add_noise = add_noise, center = T, renormalize = T)
  if(!all(is.na(tmp))) { tmp_list[[l]] <- tmp; l <- l+1 }
  
  #####
  
  tmp <- do.call(cbind, tmp_list)
  if(length(rownames(prep_list$common_score)) != 0){
    rownames(tmp) <- rownames(prep_list$common_score)
  }
  
  if(only_embedding){
    Seurat::RunUMAP(tmp, metric = "euclidean", verbose = F)@cell.embeddings
  } else {
    Seurat::RunUMAP(tmp, metric = "euclidean", reduction.key = reduction_key, verbose = F)
  }
}

#' Prepare the appropriate inputs for \code{.extract_umap_embedding}
#'
#' @param obj output from either \code{generate_data} or \code{dcca_decomposition}
#'
#' @return a list
.prepare_umap_embedding <- function(obj){
  stopifnot(class(obj) %in% c("dcca_data", "dcca_decomp"))
  
  rank_1 <- ncol(obj$distinct_score_1)
  rank_2 <- ncol(obj$distinct_score_2)
  
  svd_list <- vector("list", 2)
  svd_list[[1]] <- .svd_truncated(obj$common_mat_1 + obj$distinct_mat_1, rank_1, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[2]] <- .svd_truncated(obj$common_mat_2 + obj$distinct_mat_2, rank_2, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  names(svd_list) <- c("e1", "e2")
  
  if(length(rownames(obj$common_score)) != 0){
    for(i in 1:nrow(svd_list)){
      rownames(svd_list[[i]]$u) <- rownames(obj$common_score)
    }
  }
  
  list(common_score = obj$common_score,
       distinct_score_1 = obj$distinct_score_1, distinct_score_2 = obj$distinct_score_2,
       svd_list = svd_list)
}

#' Extract SVD embedding
#'
#' @param obj output from either \code{generate_data} or \code{dcca_decomposition}
#'
#' @return list
.extract_svd_embedding <- function(obj){
  stopifnot(class(obj) %in% c("dcca_data", "dcca_decomp"))
  
  rank_c <- ncol(obj$common_score)
  
  n <- nrow(obj$common_mat_1)
  rank_1 <- ncol(obj$distinct_score_1)
  rank_2 <- ncol(obj$distinct_score_2)
  
  svd_list <- vector("list", 6)
  
  svd_list[[1]] <- .svd_truncated(obj$common_mat_1, rank_c, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[2]] <- .svd_truncated(obj$common_mat_2, rank_c, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[3]] <- .svd_truncated(obj$distinct_mat_1, rank_1, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[4]] <- .svd_truncated(obj$distinct_mat_2, rank_2, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[5]] <- .svd_truncated(obj$common_mat_1 + obj$distinct_mat_1, rank_1, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  svd_list[[6]] <- .svd_truncated(obj$common_mat_2 + obj$distinct_mat_2, rank_2, 
                                  symmetric = F, rescale = F, K_full_rank = F)
  
  if(length(rownames(obj$common_mat_1)) != 0){
    for(i in 1:nrow(svd_list)){
      rownames(svd_list[[i]]$u) <- rownames(obj$common_mat_1)
    }
  }
  
  names(svd_list) <- c("c1", "c2", "d1", "d2", "e1", "e2")

  svd_list
}

###################################

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
            common_bool | distinct_bool)
  
  n <- nrow(common_score)
  full_mat <- .mult_mat_vec(svd_e$u, svd_e$d)
  full_mat <- canonical_score %*% crossprod(canonical_score, full_mat) # reorient for consistency for the rest of the pipeline
  center_vec <- apply(full_mat, 2, mean)
  tmp <- full_mat
  if(center) tmp <- sapply(1:ncol(full_mat), function(k){full_mat[,k] - center_vec[k]})
  if(renormalize) l2_vec <- apply(tmp, 1, function(x){.l2norm(x)})
  canonical_score <- .add_two_matrices(common_score, distinct_score)
  
  if(common_bool == distinct_bool){
    tmp <- full_mat
    
  } else {
    if(common_bool){ 
      if(add_noise){
        common_score <- .add_columnwise_noise(common_score, distinct_score)
      }
      tmp <- tcrossprod(common_score, canonical_score) %*% full_mat
      
    } else { 
      if(add_noise){
        distinct_score <- .add_columnwise_noise(distinct_score, common_score)
      }
      tmp <- tcrossprod(distinct_score, canonical_score) %*% full_mat
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

