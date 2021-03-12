#' Extract UMAP embedding
#'
#' @param prep_list return object from \code{.prepare_umap_embedding}
#' @param common_1 boolean
#' @param common_2 boolean
#' @param distinct_1 boolean
#' @param distinct_2 boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct "on the same scale" by adding appropriately-scaled Gaussian noise
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
  tmp <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_1, prep_list$svd_list$e1, 
                                common_1, distinct_1, add_noise = add_noise)
  if(!all(is.na(tmp))) { tmp_list[[l]] <- tmp; l <- l+1 }
  tmp <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_2, prep_list$svd_list$e2, 
                                common_2, distinct_2, add_noise = add_noise)
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

.prepare_umap_embedding <- function(obj){
  rank_1 <- ncol(obj$distinct_score_1)
  rank_2 <- ncol(obj$distinct_score_2)
  
  svd_list <- vector("list", 2)
  svd_list[[1]] <- .svd_truncated(obj$common_mat_1 + obj$distinct_mat_1, rank_1)
  svd_list[[2]] <- .svd_truncated(obj$common_mat_2 + obj$distinct_mat_2, rank_2)
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
#' @param obj return object from \code{dcca_decomposition}
#'
#' @return list
.extract_svd_embedding <- function(obj){
  rank_c <- ncol(obj$common_score)
  
  n <- nrow(obj$common_mat_1)
  rank_1 <- ncol(obj$distinct_score_1)
  rank_2 <- ncol(obj$distinct_score_2)
  
  svd_list <- vector("list", 6)
  
  svd_list[[1]] <- .svd_truncated(obj$common_mat_1, rank_c)
  svd_list[[2]] <- .svd_truncated(obj$common_mat_2, rank_c)
  svd_list[[3]] <- .svd_truncated(obj$distinct_mat_1, rank_1)
  svd_list[[4]] <- .svd_truncated(obj$distinct_mat_2, rank_2)
  svd_list[[5]] <- .svd_truncated(obj$common_mat_1 + obj$distinct_mat_1, rank_1)
  svd_list[[6]] <- .svd_truncated(obj$common_mat_2 + obj$distinct_mat_2, rank_2)
  
  if(length(rownames(obj$common_mat_1)) != 0){
    for(i in 1:nrow(svd_list)){
      rownames(svd_list[[i]]$u) <- rownames(obj$common_mat_1)
    }
  }
  
  names(svd_list) <- c("c1", "c2", "d1", "d2", "e1", "e2")

  svd_list
}

###################################

.extract_matrix_helper <- function(common_score, distinct_score,
                                   svd_e, common_bool, distinct_bool, add_noise = T){
  stopifnot(nrow(common_score) == nrow(distinct_score))
  if(!common_bool & !distinct_bool) return(NA)
  
  n <- nrow(common_score)
  full_mat <- .mult_mat_vec(svd_e$u, svd_e$d)
  center_vec <- apply(full_mat, 2, mean)
  tmp <- sapply(1:ncol(full_mat), function(k){
    full_mat[,k] - center_vec[k]
  })
  l2_vec <- apply(tmp, 1, function(x){.l2norm(x)})
  
  if(common_bool == distinct_bool){
    tmp <- full_mat
  } else {
    if(ncol(common_score) < ncol(distinct_score)){
      common_score <- rbind(common_score, matrix(0, nrow = n, ncol = ncol(distinct_score) - ncol(common_score)))
    }
    canonical_score <- common_score+distinct_score
    
    if(common_bool){ 
      if(add_noise){
        for(i in 1:ncol(common_score)){
          sd1 <- stats::sd(common_score[,i]); sd2 <- stats::sd(distinct_score[,i])
          if(sd1 < sd2){common_score[,i] <- common_score[,i] + stats::rnorm(n, sd = max(sd2 - sd1, sd1/3))} 
        }
      }
      tmp <- tcrossprod(common_score, canonical_score) %*% full_mat
    } else { 
      if(add_noise){
        for(i in 1:ncol(common_score)){
          sd1 <- stats::sd(common_score[,i]); sd2 <- stats::sd(distinct_score[,i])
          if(sd2 < sd1){distinct_score[,i] <- distinct_score[,i] + stats::rnorm(n, sd = max(sd1 - sd2, sd2/3))} 
        }
      }
      tmp <- tcrossprod(distinct_score, canonical_score) %*% full_mat
    }
    
    # center variables
    tmp <- sapply(1:ncol(tmp), function(k){
      tmp[,k] - center_vec[k]
    })
  }
  
  # normalize cells
  if(common_bool | distinct_bool){ .mult_vec_mat(1/l2_vec, tmp)
  } else{
    NA
  }
}