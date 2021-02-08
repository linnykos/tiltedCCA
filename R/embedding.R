#' Extract UMAP embedding
#'
#' @param svd_list return object from \code{extract_svd_embedding}
#' @param common_1 boolean
#' @param common_2 boolean
#' @param distinct_1 boolean
#' @param distinct_2 boolean
#' @param noise_val non-negative numeric (set to 0 to essentially ignore)
#' @param vis_param list containing 4 scalars: \code{cmin_1}, \code{cmin_2}, \code{dmin_1}, \code{dmin_2}
#' @param only_embedding boolean
#' @param reduction_key string for \code{Seurat::RunUMAP}
#'
#' @return 2-column matrix or \code{Seurat} object
#' @export
extract_umap_embedding <- function(svd_list, common_1 = T, common_2 = T,
                              distinct_1 = T, distinct_2 = T,
                              noise_val = 0.05,
                              vis_param = NA,
                              only_embedding = T, reduction_key = "UMAP"){
  n <- nrow(svd_list[[1]]$u)
  
  c1 <- svd_list[[1]]$d[1]
  c2 <- svd_list[[2]]$d[1]
  d1 <- svd_list[[3]]$d[1]
  d2 <- svd_list[[4]]$d[1]
  
  if(all(!is.na(vis_param))){
    cmin_1 <- vis_param$cmin_1; cmin_2 <- vis_param$cmin_2
    dmin_1 <- vis_param$dmin_1; dmin_2 <- vis_param$dmin_2
  } else {
    cmin_1 <- mean(svd_list[[1]]$d); cmin_2 <- mean(svd_list[[2]]$d)
    dmin_1 <- mean(svd_list[[3]]$d); dmin_2 <- mean(svd_list[[4]]$d)
  }
  
  if(common_1){ 
    svd_list[[1]]$d <- svd_list[[1]]$d/(c1 + d1) 
  } else { 
    svd_list[[1]]$u <- matrix(stats::rnorm(n), nrow = n, ncol = 1)
    svd_list[[1]]$d <- cmin_1/(c1+d1)*ifelse(distinct_1, noise_val, 0)
  }
  
  if(common_2){ 
    svd_list[[2]]$d <- svd_list[[2]]$d/(c2 + d2) 
  } else { 
    svd_list[[2]]$u <- matrix(stats::rnorm(n), nrow = n, ncol = 1) 
    svd_list[[2]]$d <- cmin_2/(c2+d2)*ifelse(distinct_2, noise_val, 0)
  }
 
  if(distinct_1){ 
    svd_list[[3]]$d <- svd_list[[3]]$d/(c1 + d1) 
  } else { 
    svd_list[[3]]$u <- matrix(stats::rnorm(n), nrow = n, ncol = 1)
    svd_list[[3]]$d <- dmin_1/(c1+d1)*ifelse(common_1, noise_val, 0)
  }

  if(distinct_2){ 
    svd_list[[4]]$d <- svd_list[[4]]$d/(c2 + d2) 
  } else { 
    svd_list[[4]]$u <- matrix(stats::rnorm(n), nrow = n, ncol = 1) 
    svd_list[[4]]$d <- dmin_2/(c2+d2)*ifelse(common_2, noise_val, 0)
  }
  
  tmp <- do.call(cbind, lapply(svd_list, function(res){
    .mult_mat_vec(res$u, res$d)
  }))
  
  if(length(rownames(svd_list[[1]]$u)) != 0){
    rownames(tmp) <- rownames(svd_list[[1]]$u)
  }
  
  if(only_embedding){
    Seurat::RunUMAP(tmp, verbose = F)@cell.embeddings
  } else {
    Seurat::RunUMAP(tmp, reduction.key = reduction_key, verbose = F)
  }
}

#' Extract SVD embedding
#'
#' @param obj return object from \code{dcca_decomposition}
#' @param mode \code{"dcca"} or \code{"dmca"}
#'
#' @return list
#' @export
extract_svd_embedding <- function(obj, mode = "dcca"){
  if(class(obj) == "dcca_decomp" | mode == "dcca") {
    rank_c <- ncol(obj$common_score)
  } else {
    rank_c <- ncol(obj$common_score_1)
  }
  n <- nrow(obj$common_mat_1)
  rank_1 <- ncol(obj$distinct_score_1)
  rank_2 <- ncol(obj$distinct_score_2)
  
  svd_list <- vector("list", 4)
  
  svd_list[[1]] <- .svd_truncated(obj$common_mat_1, rank_c)
  svd_list[[2]] <- .svd_truncated(obj$common_mat_2, rank_c)
  svd_list[[3]] <- .svd_truncated(obj$distinct_mat_1, rank_1)
  svd_list[[4]] <- .svd_truncated(obj$distinct_mat_2, rank_2)
  
  if(length(rownames(obj$common_mat_1)) != 0){
    for(i in 1:4){
      rownames(svd_list[[i]]$u) <- rownames(obj$common_mat_1)
    }
  }

  svd_list
}