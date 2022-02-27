#' Prepare the appropriate inputs
#'
#' @param obj output from either \code{generate_data} or \code{dcca_decomposition}
#' @param data_1 boolean, for computing the embedding for data 1
#' @param data_2 boolean, for computing the embedding for data 2
#' @param center boolean
#' @param renormalize boolean
#'
#' @return a list of 3 elements, representing the low-dimensional
#' representation for the common, distinct and everything matrices
#' for the intended dataset(s)
.prepare_embeddings <- function(obj, data_1, data_2, center, renormalize){
  if(data_1){
    res_1 <- .prepare_embeddings_singleton(obj$common_score, 
                                           obj$distinct_score_1,
                                           obj$svd_1,
                                           center = center, 
                                           renormalize = renormalize)
  } else res_1 <- NA
  if(data_2){
    res_2 <- .prepare_embeddings_singleton(obj$common_score, 
                                           obj$distinct_score_2,
                                           obj$svd_2, 
                                           center = center, 
                                           renormalize = renormalize)
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

.prepare_embeddings_singleton <- function(common_score, 
                                          distinct_score, 
                                          svd_res, 
                                          center, 
                                          renormalize){
  embedding <- vector("list", 3)
  names(embedding) <- c("common", "distinct", "everything")
  
  embedding[[1]] <- .extract_matrix_helper(common_score, 
                                           distinct_score,
                                           svd_res, 
                                           common_bool = T, 
                                           distinct_bool = F,
                                           center = center, 
                                           renormalize = renormalize)
  embedding[[2]] <- .extract_matrix_helper(common_score, 
                                           distinct_score,
                                           svd_res, 
                                           common_bool = F, 
                                           distinct_bool = T,
                                           center = center, 
                                           renormalize = renormalize)
  embedding[[3]] <- .extract_matrix_helper(common_score, 
                                           distinct_score,
                                           svd_res, 
                                           common_bool = T, 
                                           distinct_bool = T,
                                           center = center, 
                                           renormalize = renormalize)
  
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
#' @param center boolean. If \code{TRUE}, center each canonical variable (i.e., column)
#' where the centering is based on \code{common_score+distinct_score}
#' @param renormalize boolean. If \code{TRUE}, normalize each sample (i.e., row)
#' to have unit norm, where the rescaling factor is based on \code{common_score+distinct_score}.
#' This is suggested to be only used if \code{center=TRUE}, in which case the centering
#' happens before the renormalization
#'
#' @return a matrix
.extract_matrix_helper <- function(common_score, 
                                   distinct_score,
                                   svd_e, 
                                   common_bool, 
                                   distinct_bool,
                                   center, 
                                   renormalize){
  stopifnot(nrow(common_score) == nrow(distinct_score),
            nrow(common_score) == nrow(svd_e$u), nrow(distinct_score) == nrow(svd_e$u),
            ncol(common_score) <= length(svd_e$d), ncol(distinct_score) == length(svd_e$d),
            common_bool | distinct_bool)
  
  n <- nrow(common_score)
  
  # normalize so crossprod(common_score) and crossprod(distinct_score) have sqrt(n) along diagonal
  common_score <- common_score/n^(1/4)
  distinct_score <- distinct_score/n^(1/4)
  
  canonical_score <- .add_two_matrices(common_score, distinct_score)
  full_mat <- .mult_mat_vec(svd_e$u, svd_e$d/svd_e$d[1])
  tmp <- canonical_score %*% crossprod(canonical_score, full_mat) # reorient for consistency for the rest of the pipeline
  center_vec <- apply(tmp, 2, mean)
  if(center) tmp <- sapply(1:ncol(tmp), function(k){tmp[,k] - center_vec[k]})
  
  if(common_bool != distinct_bool){
    if(common_bool){ 
      if(ncol(common_score) < ncol(distinct_score)) {
        common_score <- cbind(common_score, matrix(0, nrow = n, ncol = ncol(distinct_score)-ncol(common_score)))
      }
      tmp <- common_score %*% crossprod(canonical_score, full_mat)
      
    } else { 
      if(ncol(distinct_score) < ncol(common_score)) {
        distinct_score <- cbind(distinct_score, matrix(0, nrow = n, ncol = ncol(common_score)-ncol(distinct_score)))
      }
      tmp <- distinct_score %*% crossprod(canonical_score, full_mat)
    }
    
    # center variables
    if(center){
      tmp <- sapply(1:ncol(tmp), function(k){tmp[,k] - center_vec[k]})
    }
  }
  
  # normalize cells
  if(renormalize) {
    l2_vec <- apply(tmp, 1, function(x){.l2norm(x)})
    .mult_vec_mat(1/l2_vec, tmp)
  } else {
    tmp
  }
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

