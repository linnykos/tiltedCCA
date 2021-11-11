.form_snns <- function(num_neigh,
                       svd_1,
                       svd_2){
  tmp_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  snn_mat_1 <- .form_snn_mat(bool_intersect = T,
                                            mat = tmp_1, 
                                            num_neigh = num_neigh)
  metacell_clustering_1 <- lapply(1:nrow(snn_mat_1), function(i){
    .nonzero_col(snn_mat_1, i, bool_value = F)
  })
  
  tmp_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  snn_mat_2 <- .form_snn_mat(bool_intersect = T,
                                            mat = tmp_2, 
                                            num_neigh = num_neigh)
  metacell_clustering_2 <- lapply(1:nrow(snn_mat_2), function(i){
    .nonzero_col(snn_mat_2, i, bool_value = F)
  })
  
  list(metacell_clustering_1 = metacell_clustering_1,
       metacell_clustering_2 = metacell_clustering_2)
}

.form_snn_mat <- function(bool_intersect, mat, num_neigh){
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
                                     dims = c(n,n),
                                     repr = "C")
  
  if(bool_intersect) sparse_mat <- sparse_mat * Matrix::t(sparse_mat)
  sparse_mat
}


