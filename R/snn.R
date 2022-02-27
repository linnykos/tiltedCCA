.form_snn_mat <- function(mat, 
                          num_neigh,
                          bool_cosine, # suggested: TRUE
                          bool_intersect, # suggested: TRUE
                          min_deg, #suggested: 1
                          verbose = T,
                          tol = 1e-4){
  stopifnot(num_neigh >= min_deg, min_deg >= 0)
  
  if(bool_cosine) {
    l2_vec <- apply(mat, 1, .l2norm)
    l2_vec[l2_vec <= tol] <- tol
    mat <- .mult_vec_mat(1/l2_vec, mat)
  }
  
  if(verbose) print("Compute NNs")
  n <- nrow(mat)
  nn_mat <- RANN::nn2(mat, k = num_neigh)$nn.idx
  if(all(nn_mat[,1] == 1:n)){
    nn_mat <- nn_mat[,-1,drop = F]
  }
  
  if(verbose) print("Forming NN matrix")
  i_vec <- rep(1:n, times = ncol(nn_mat))
  j_vec <- as.numeric(nn_mat)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                     j = j_vec,
                                     x = rep(1, length(i_vec)),
                                     dims = c(n,n),
                                     repr = "C")
  
  if(bool_intersect) {
    sparse_mat <- sparse_mat * Matrix::t(sparse_mat)
  } else {
    sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
    sparse_mat@x <- rep(1, length(sparse_mat@x))
  }
  
  if(min_deg > 0){
    deg_vec <- sparseMatrixStats::rowSums2(sparse_mat)
    if(min(deg_vec) < min_deg) {
      idx <- which(deg_vec < min_deg)
      if(verbose) print(paste0("Joining the ", length(idx), " nodes with too few neighbors"))
      
      for(i in idx) sparse_mat[i,nn_mat[i,1:min_deg]] <- 1
      
      sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
      sparse_mat@x <- rep(1, length(sparse_mat@x))
    }
  }
  
  sparse_mat
}

.compute_laplacian_basis <- function(sparse_mat,
                                     k, # suggested: 50
                                     verbose = T){
  if(verbose) print("Computing symmetrized Laplacians")
  deg_vec <- sparseMatrixStats::rowSums2(sparse_mat)
  deg_vec[deg_vec == 0] <- 1
  diag_mat <- Matrix::Diagonal(x = 1/sqrt(deg_vec))
  lap_mat <-  diag_mat %*% sparse_mat %*% diag_mat
  
  if(verbose) print("Converting to random walk Laplacian")
  deg_vec <- sparseMatrixStats::rowSums2(lap_mat)
  deg_vec[deg_vec == 0] <- 1
  diag_mat <- Matrix::Diagonal(x = 1/deg_vec)
  lap_mat <-  diag_mat %*% lap_mat
  
  if(verbose) print("Extracting basis")
  eigen_res <- irlba::partial_eigen(lap_mat, n = k, symmetric = F)
  dimred <- .mult_mat_vec(eigen_res$vectors, eigen_res$values)
  
  dimred
}


