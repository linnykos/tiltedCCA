.determine_cluster <- function(mat,
                               metacell_clustering,
                               num_neigh){
  stopifnot(is.factor(metacell_clustering))
  
  n <- nrow(mat)
  snn_mat <- .form_snn_mat(mat, num_neigh = num_neigh)
  level_vec <- levels(metacell_clustering)
  K <- length(level_vec)
  block_density <- sum(sapply(1:K, function(k){
    idx <- which(metacell_clustering == k)
    sum(snn_mat[idx,idx])
  }))
  overall_density <- sum(snn_mat)
  
  block_density/(overall_density-block_density) * overall_density/n
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
