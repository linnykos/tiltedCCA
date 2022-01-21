form_snns <- function(num_neigh,
                      svd_1,
                      svd_2,
                      bool_intersect = T,
                      distance_func = "cosine",
                      min_deg = 1,
                      verbose = T){
  stopifnot(distance_func %in% c("cosine", "euclidean"),
            length(distance_func) == 1)
  
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  
  if(distance_func == "cosine"){
    norm_vec <- apply(dimred_1, 1, .l2norm)
    dimred_1 <- .mult_vec_mat(1/norm_vec, dimred_1)
    
    norm_vec <- apply(dimred_2, 1, .l2norm)
    dimred_2 <- .mult_vec_mat(1/norm_vec, dimred_2)
  }
  
  snn_mat_1 <- .form_snn_mat(bool_intersect = bool_intersect,
                             mat = dimred_1, 
                             min_deg = min_deg,
                             num_neigh = num_neigh,
                             verbose = verbose)
  snn_mat_2 <- .form_snn_mat(bool_intersect = bool_intersect,
                             mat = dimred_2, 
                             min_deg = min_deg,
                             num_neigh = num_neigh,
                             verbose = verbose)
  
  rownames(snn_mat_1) <- rownames(svd_1$u)
  colnames(snn_mat_1) <- rownames(svd_1$u)
  rownames(snn_mat_2) <- rownames(svd_1$u)
  rownames(snn_mat_2) <- rownames(svd_1$u)
  
  list(snn_mat_1 = snn_mat_1,
       snn_mat_2 = snn_mat_2)
}

compute_laplacian_basis <- function(sparse_mat,
                                    k = 20,
                                    verbose = T){
  stopifnot(k <= ncol(sparse_mat), nrow(sparse_mat) == ncol(sparse_mat))
  
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
  eigen_res <- irlba::partial_eigen(lap_mat, n = k+1, symmetric = F)
  dimred <- .mult_mat_vec(eigen_res$vectors[,-1,drop = F], eigen_res$values[-1])
  
  rownames(dimred) <- rownames(sparse_mat)
  
  dimred
}

################################

.form_snn_mat <- function(bool_intersect,
                          mat, 
                          min_deg,
                          num_neigh,
                          verbose = T){
  stopifnot(num_neigh >= min_deg, min_deg >= 0)
  
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
      
      i_vec2 <- rep(1:n, times = ncol(nn_mat[,1:min_deg]))
      j_vec2 <- as.numeric(nn_mat[,1:min_deg])
      sparse_mat2 <- Matrix::sparseMatrix(i = i_vec2,
                                          j = j_vec2,
                                          x = rep(1, length(i_vec2)),
                                          dims = c(n,n),
                                          repr = "C")
      
      sparse_mat <- sparse_mat + sparse_mat2 + Matrix::t(sparse_mat2)
      sparse_mat@x <- rep(1, length(sparse_mat@x))
    }
  }
  
  sparse_mat
}



