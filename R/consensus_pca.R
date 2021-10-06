consensus_pca <- function(input_1, input_2,
                          dims_1 = NA, dims_2 = NA,
                          center_1 = T, center_2 = T,
                          scale_1 = T, scale_2 = T){
  if(!is.list(input_1) & !is.list(input_2)){
    stopifnot(all(!is.na(dims_1)) & all(!is.na(dims_2)))
    
    svd_1 <- .svd_truncated(input_1, K = rank_1, symmetric = F, rescale = F, 
                            mean_vec = center_1, sd_vec = scale_1, K_full_rank = F)
    svd_2 <- .svd_truncated(input_2, K = rank_2, symmetric = F, rescale = F, 
                            mean_vec = center_2, sd_vec = scale_2, K_full_rank = F)
    
    svd_1 <- .check_svd(svd_1, dims = dims_1)
    svd_2 <- .check_svd(svd_2, dims = dims_2)
  } else {
    svd_1 <- input_1; svd_2 <- input_2
  }
  
  stopifnot(nrow(svd_1$u) == nrow(svd_2$u))
  n <- nrow(svd_1$u)
  svd_1$d <- svd_1$d/svd_1$d[1]*sqrt(n)
  svd_2$d <- svd_2$d/svd_2$d[1]*sqrt(n)
  
  pca_mat <- cbind(.mult_mat_vec(svd_1$u, svd_1$d),
                   .mult_mat_vec(svd_2$u, svd_2$d))
  
  if(length(rownames(svd_1$u)) != 0) rownames(pca_mat) <- rownames(svd_1$u)
  
  pca_mat
}