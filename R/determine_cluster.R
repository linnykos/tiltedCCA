.determine_cluster <- function(mat,
                               metacell_clustering,
                               num_neigh,
                               tol = 1e-3){
  stopifnot(is.factor(metacell_clustering))
  
  n <- nrow(mat)
  snn_mat <- .form_snn_mat(mat, num_neigh = num_neigh)
  density_mat <- .compute_density_matrix(as.matrix(snn_mat),
                                         metacell_clustering = metacell_clustering)
  K <- length(levels(metacell_clustering))
  
  quality_vec <- sapply(1:K, function(k){
    offdiag_mean <- mean(density_mat[k,-k])
    offdiag_sd <- sd(density_mat[k,-k])
    # print(paste0(k, ": ", round(density_mat[k,k],2)))
    # print(paste0(k, ": ", round(offdiag_mean+offdiag_sd,2)))
    
    density_mat[k,k]/max(c(offdiag_mean-offdiag_sd, tol))
  })
  
  mean(quality_vec) * sum(snn_mat)/n
}

.form_snn_mat <- function(mat, num_neigh){
  stopifnot(num_neigh > 1)
  n <- nrow(mat)
  nn_mat <- RANN::nn2(mat, k = num_neigh)$nn.idx
  if(all(nn_mat[,1] == 1:n)){
    nn_mat <- nn_mat[,-1,drop = F]
  }
  
  i_vec <- rep(1:n, times = ncol(nn_mat))
  j_vec <- as.numeric(nn_mat)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                    j = j_vec,
                                    x = rep(1, length(i_vec)),
                                    repr = "C")
  sparse_mat <- sparse_mat * Matrix::t(sparse_mat)
  sparse_mat
}

.compute_density_matrix <- function(mat, 
                                    metacell_clustering){
  diag(mat) <- 0
  level_vec <- levels(metacell_clustering)
  idx_list <- lapply(level_vec, function(clust){
    which(metacell_clustering == clust)
  })
  
  col_sum_mat <- sapply(1:length(level_vec), function(k){
    matrixStats::rowSums2(mat[,idx_list[[k]],drop = F])
  })
  
  sum_mat <- sapply(1:length(level_vec), function(k){
    matrixStats::colSums2(col_sum_mat[idx_list[[k]],,drop = F])
  })
  
  len_vec <- sapply(idx_list, length)
  len_mat <- tcrossprod(len_vec)
  
  density_mat <- sum_mat/len_mat
  
  rownames(density_mat) <- level_vec
  colnames(density_mat) <- level_vec
  density_mat
}
