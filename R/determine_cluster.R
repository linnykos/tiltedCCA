.determine_cluster <- function(mat,
                               metacell_clustering_1,
                               metacell_clustering_2,
                               n_idx,
                               num_neigh,
                               tol = 1e-3){
  stopifnot(is.list(metacell_clustering_1), 
            is.list(metacell_clustering_2), 
            length(metacell_clustering_1) == nrow(mat),
            length(metacell_clustering_2) == nrow(mat))
  
  if(length(n_idx) < nrow(mat)){
    mat <- mat[n_idx,,drop = F]
    metacell_clustering_1 <- metacell_clustering_1[n_idx]
    metacell_clustering_2 <- metacell_clustering_2[n_idx]
  }
  
  n <- nrow(mat)
  snn_mat <- .form_snn_mat(bool_intersect = T,
                           mat = mat, 
                           num_neigh = num_neigh)
  nn_idx_list <- lapply(1:nrow(snn_mat), function(i){
    .nonzero_col(snn_mat, col_idx = i, bool_value = F)
  })
  overlap_1 <- .compute_set_inclusion(list_1 = metacell_clustering_1,
                                      list_2 = nn_idx_list)
  overlap_2 <- .compute_set_inclusion(list_1 = metacell_clustering_2,
                                      list_2 = nn_idx_list)
  min(overlap_1, overlap_2)
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

.compute_set_inclusion <- function(list_1, list_2){
  stopifnot(is.list(list_1), is.list(list_2), length(list_1) == length(list_2))
  
  n <- length(list_1)
  overlap_vec <- sapply(1:n, function(i){
    # [[note to self: add a test for this corner case]]
    if(length(list_1[[i]]) == 0 && length(list_2[[i]]) == 0) return(0)
    
    length(intersect(list_1[[i]], list_2[[i]]))/length(unique(c(list_1[[i]], list_2[[i]])))
  })
  
  mean(overlap_vec)
}
