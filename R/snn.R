compute_snns <- function(input_obj,
                         latent_k,
                         num_neigh,
                         bool_cosine = T,
                         bool_intersect = T,
                         min_deg = 1,
                         tol = 1e-4,
                         verbose = 0){
  stopifnot(inherits(input_obj), "multiSVD")
  
  metacell_clustering <- .get_metacell(input_obj,
                                       resolution = "cell", 
                                       type = "list", 
                                       what = "metacell_clustering")
  if(!all(is.null(metacell_clustering))){
    averaging_mat <- .generate_averaging_matrix(n, metacell_clustering)
  } else {
    averaging_mat <- NULL
  }
  
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  dimred_1 <- .get_postDimred(input_obj, averaging_mat = averaging_mat)
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  dimred_2 <- .get_postDimred(input_obj, averaging_mat = averaging_mat)
  
  snn_1 <- .form_snn_mat(mat = dimred_1,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg,
                         tol = tol,
                         verbose = verbose)
  snn_2 <- .form_snn_mat(mat = dimred_2,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg,
                         tol = tol,
                         verbose = verbose)
  
  clustering_1 <- .get_metacell(input_obj,
                                resolution = "metacell", 
                                type = "factor",
                                what = "large_clustering_1")
  clustering_2 <- .get_metacell(input_obj,
                                resolution = "metacell", 
                                type = "factor",
                                what = "large_clustering_2")
  common_snn <- .compute_common_snn(snn_1 = snn_1, 
                                    snn_2 = snn_2,
                                    clustering_1 = clustering_1, 
                                    clustering_2 = clustering_2,
                                    num_neigh = num_neigh,
                                    verbose = verbose)
  snn_list <- list(snn_1 = snn_1, snn_2 = snn_2, common_snn = common_snn)
  laplacian_list <- lapply(snn_list, function(snn_mat){
    .compute_laplacian_basis(latent_k = latent_k,
                             sparse_mat = snn_mat,
                             verbose = verbose)
  })
  names(laplacian_list) <- c("laplacian_1", "laplacian_2", "common_laplacian")
  
  input_obj$snn_list <- snn_list
  input_obj$laplacian_list <- laplacian_list
  param <- .form_snn_param(bool_cosine = bool_cosine,
                           bool_intersect = bool_intersect,
                           latent_k = latent_k,
                           min_deg = min_deg,
                           num_neigh = num_neigh)
  param <- .combine_two_named_lists(.get_param(input_obj), param)
  input_obj$param <- param
  class(input_obj) <- c("snn", "multiSVD")
  
  input_obj
}

##################################

.form_snn_mat <- function(mat, 
                          num_neigh,
                          bool_cosine, # suggested: TRUE
                          bool_intersect, # suggested: TRUE
                          min_deg, #suggested: 1
                          tol = 1e-4,
                          verbose = 0){
  stopifnot(num_neigh >= min_deg, min_deg >= 0)
  
  if(bool_cosine) {
    l2_vec <- apply(mat, 1, .l2norm)
    l2_vec[l2_vec <= tol] <- tol
    mat <- .mult_vec_mat(1/l2_vec, mat)
  }
  
  if(verbose >= 1) print("Compute NNs")
  n <- nrow(mat)
  nn_mat <- RANN::nn2(mat, k = num_neigh)$nn.idx
  if(all(nn_mat[,1] == 1:n)){
    nn_mat <- nn_mat[,-1,drop = F]
  }
  
  if(verbose >= 1) print("Forming NN matrix")
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
      if(verbose >= 1) print(paste0("Joining the ", length(idx), " nodes with too few neighbors"))
      
      for(i in idx) sparse_mat[i,nn_mat[i,1:min_deg]] <- 1
      
      sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
      sparse_mat@x <- rep(1, length(sparse_mat@x))
    }
  }
  
  if(length(rownames(mat)) > 0){
    rownames(sparse_mat) <- rownames(mat)
    colnames(sparse_mat) <- rownames(mat)
  }
  
  sparse_mat
}

.compute_laplacian_basis <- function(latent_k, # suggested: 50
                                     sparse_mat,
                                     verbose = 0){
  if(verbose >= 1) print("Computing symmetrized Laplacians")
  deg_vec <- sparseMatrixStats::rowSums2(sparse_mat)
  deg_vec[deg_vec == 0] <- 1
  diag_mat <- Matrix::Diagonal(x = 1/sqrt(deg_vec))
  lap_mat <-  diag_mat %*% sparse_mat %*% diag_mat
  
  if(verbose >= 1) print("Converting to random walk Laplacian")
  deg_vec <- sparseMatrixStats::rowSums2(lap_mat)
  deg_vec[deg_vec == 0] <- 1
  diag_mat <- Matrix::Diagonal(x = 1/deg_vec)
  lap_mat <-  diag_mat %*% lap_mat
  
  # [[note to self: Do I need to remove the first eigenvector?]]
  if(verbose >= 1) print("Extracting basis")
  eigen_res <- irlba::partial_eigen(lap_mat, n = latent_k, symmetric = F)
  dimred <- .mult_mat_vec(eigen_res$vectors, eigen_res$values)
  
  if(length(rownames(sparse_mat)) > 0){
    rownames(dimred) <- rownames(sparse_mat)
  }
  
  dimred
}

#######################

.form_snn_param <- function(bool_cosine,
                            bool_intersect,
                            latent_k,
                            min_deg,
                            num_neigh){
  list(snn_bool_cosine = bool_cosine,
       snn_bool_intersect = bool_intersect,
       snn_latent_k = latent_k,
       snn_min_deg = min_deg,
       snn_num_neigh = num_neigh)
}

