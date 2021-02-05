#' Extract UMAP embedding
#'
#' @param obj return object from \code{dcca_decomposition}
#' @param common boolean
#' @param distinct_1 boolean
#' @param distinct_2 boolean
#' @param only_embedding boolean
#' @param reduction_key string for \code{Seurat::RunUMAP}
#'
#' @return 2-column matrix
#' @export
extract_embedding <- function(obj, common_1 = T, common_2 = T,
                              distinct_1 = T, distinct_2 = T, mode = "dcca",
                              only_embedding = T, reduction_key = "UMAP"){
  if(class(obj) == "dcca_decomp" | mode == "dcca") {
    rank_c <- ifelse((common_1 | common_2), ncol(obj$common_score), 1)
  } else {
    rank_c <- ifelse((common_1 | common_2), ncol(obj$common_score_1), 1)
  }
  n <- nrow(obj$common_mat_1)
  rank_1 <- ifelse(distinct_1, ncol(obj$distinct_score_1), 1)
  rank_2 <- ifelse(distinct_2, ncol(obj$distinct_score_2), 1)
  
  svd_list <- vector("list", 0)
  len <- length(svd_list)
  
  tmp1 <- .svd_truncated(obj$common_mat_1, rank_c); c1 <- tmp1$d[1]
  tmp2 <- .svd_truncated(obj$common_mat_2, rank_c); c2 <- tmp2$d[1]
  tmp3 <- .svd_truncated(obj$distinct_mat_1, rank_1); d1 <- tmp3$d[1]
  tmp4 <- .svd_truncated(obj$distinct_mat_2, rank_2); d2 <- tmp4$d[1]
  
  if(common_1){ tmp1$d <- tmp1$d/(c1 + d1) } else { tmp1$u <- matrix(1, nrow = n, ncol = 1) ; tmp1$d <- c1/(c1+d1)}
  svd_list[[len+1]] <- tmp1; len <- len + 1
  if(common_2){ tmp2$d <- tmp2$d/(c2 + d2) } else { tmp2$u <- matrix(1, nrow = n, ncol = 1) ; tmp2$d <- c2/(c2+d2)}
  svd_list[[len+1]] <- tmp2; len <- len + 1
  if(distinct_1){ tmp3$d <- tmp3$d/(c1 + d1) } else { tmp3$u <- matrix(1, nrow = n, ncol = 1) ; tmp3$d <- d1/(c1+d1)}
  svd_list[[len+1]] <- tmp3; len <- len + 1
  if(distinct_2){ tmp4$d <- tmp4$d/(c2 + d2) } else { tmp4$u <- matrix(1, nrow = n, ncol = 1) ; tmp4$d <- d2/(c2+d2)}
  svd_list[[len+1]] <- tmp4; len <- len + 1
  
  # if(common_1){ tmp1$d <- tmp1$d/(c1 + d1); svd_list[[len+1]] <- tmp1; len <- len + 1 }
  # if(common_2){ tmp2$d <- tmp2$d/(c2 + d2); svd_list[[len+1]] <- tmp2; len <- len + 1 }
  # if(distinct_1){ tmp3$d <- tmp3$d/(c1 + d1); svd_list[[len+1]] <- tmp3; len <- len + 1 }
  # if(distinct_2){ tmp4$d <- tmp4$d/(c2 + d2); svd_list[[len+1]] <- tmp4; len <- len + 1 }
  
  tmp <- do.call(cbind, lapply(svd_list, function(res){
    .mult_mat_vec(res$u, res$d)
  }))
  
  if(length(rownames(obj$common_mat_1)) != 0){
    rownames(tmp) <- rownames(obj$common_mat_1)
  }
  
  if(only_embedding){
    Seurat::RunUMAP(tmp, verbose = F)@cell.embeddings
  } else {
    Seurat::RunUMAP(tmp, reduction.key = reduction_key, verbose = F)
  }
}