compute_min_subspace <- function(dimred_1, dimred_2,
                                 metacell_clustering_1,
                                 metacell_clustering_2,
                                 k = min(c(ncol(dimred_1), ncol(dimred_2))),
                                 num_neigh = 30,
                                 verbose = F){
  min_mat <- .min_embedding(dimred_1 = dimred_1, dimred_2 = dimred_2,
                            metacell_clustering_1 = metacell_clustering_1,
                            metacell_clustering_2 = metacell_clustering_2,
                            num_neigh = num_neigh,
                            verbose = verbose)
  
  ## [[TODO: If we're doing this, then we don't need to kernelize...]]
  min_mat@x <- rep(1, length(min_mat@x))
  subspace_mat <- .svd_truncated(min_mat, K = k, symmetric = F, rescale = F, 
                                 mean_vec = T, sd_vec = F, K_full_rank = F)$u
  
  if(length(rownames(dimred_1)) != 0){
    rownames(subspace_mat) <- rownames(dimred_1)
  }
  
  list(min_mat = min_mat, subspace_mat = subspace_mat)
}

.min_embedding <- function(dimred_1, dimred_2,
                           metacell_clustering_1,
                           metacell_clustering_2,
                           num_neigh, 
                           verbose = F){
  n <- nrow(dimred_1)
  
  if(verbose) print("Computing NNs")
  nn_res_1 <- RANN::nn2(dimred_1, k = num_neigh+1)$nn.idx[,-1]
  nn_res_2 <- RANN::nn2(dimred_2, k = num_neigh+1)$nn.idx[,-1]
  
  if(verbose) print("Joining NNs")
  if(!is.factor(metacell_clustering_1) & !is.factor(metacell_clustering_2)){
    nn_list <- lapply(1:n, function(i){
      if(verbose && i %% floor(n/10) == 0) cat('*')
      sort(unique(c(nn_res_1[i,], nn_res_2[i,])))
    })
  } else {
    nn_list <- lapply(1:n, function(i){
      if(verbose && i %% floor(n/10) == 0) cat('*')
      .compute_metacluster_subsample(metacell_clustering_1 = metacell_clustering_1,
                                     metacell_clustering_2 = metacell_clustering_2,
                                     nn_vec_1 = nn_res_1[i,],
                                     nn_vec_2 = nn_res_2[i,],
                                     num_neigh = num_neigh)
    })
  }
  
  if(verbose) print("Computing max value")
  max_dist_1 <- sapply(1:n, function(i){
    if(verbose && i %% floor(n/10) == 0) cat('*')
    max(sapply(nn_res_1[i,], function(j){
      .l2norm(dimred_1[i,] - dimred_1[j,])
    }))
  })
  max_dist_2 <- sapply(1:n, function(i){
    if(verbose && i %% floor(n/10) == 0) cat('*')
    max(sapply(nn_res_2[i,], function(j){
      .l2norm(dimred_2[i,] - dimred_2[j,])
    }))
  })
  
  if(verbose) print("Computing kernels for Modality 1")
  kernel_list_1 <- lapply(1:n, function(i){
    if(verbose && i %% floor(n/10) == 0) cat('*')
    dist_vec <- sapply(nn_list[[i]], function(j){
      .l2norm(dimred_1[i,] - dimred_1[j,])
    })
    
    exp(-dist_vec/max_dist_1[i])
  })
  
  if(verbose) print("Computing kernels for Modality 2")
  kernel_list_2 <- lapply(1:n, function(i){
    if(verbose && i %% floor(n/10) == 0) cat('*')
    dist_vec <- sapply(nn_list[[i]], function(j){
      .l2norm(dimred_2[i,] - dimred_2[j,])
    })
    
    exp(-dist_vec/max_dist_2[i])
  })
  
  if(verbose) print("Computing  kernel")
  dist_list <- lapply(1:n, function(i){
    if(verbose && i %% floor(n/10) == 0) cat('*')
    pmax(kernel_list_1[[i]], kernel_list_2[[i]])
  })
  
  if(verbose) print("Forming matrix")
  i_vec <- rep(1:n, times = sapply(nn_list, length))
  j_vec <- unlist(nn_list)
  x_vec <- unlist(dist_list)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                     j = j_vec,
                                     x = x_vec,
                                     dims = c(n,n),
                                     repr = "C")
  if(length(rownames(dimred_1)) != 0){
    colnames(sparse_mat) <- rownames(dimred_1)
    rownames(sparse_mat) <- rownames(dimred_1)
  }
  
  sparse_mat
}

.compute_metacluster_subsample <- function(metacell_clustering_1,
                                           metacell_clustering_2,
                                           nn_vec_1,
                                           nn_vec_2,
                                           num_neigh,
                                           tol = 0.1){
  prop_1 <- table(metacell_clustering_1[nn_vec_2])/length(nn_vec_2)
  prop_2 <- table(metacell_clustering_2[nn_vec_1])/length(nn_vec_1)
  desired_prop <- tcrossprod(prop_1, prop_2)
  desired_mat <- pmax(num_neigh*desired_prop, tol)
  all_nn <- unique(c(nn_vec_1, nn_vec_2))
  empirical_count <- table(metacell_clustering_1[all_nn], metacell_clustering_2[all_nn])

  idx_mat <- cbind(as.numeric(metacell_clustering_1[all_nn]),
                   as.numeric(metacell_clustering_2[all_nn]))
  prob_vec <- sapply(1:nrow(idx_mat), function(x){
    i <- idx_mat[x,1]; j <- idx_mat[x,2]
    desired_mat[i,j]/empirical_count[i,j]
  })
  
  sample(all_nn, size = num_neigh, replace = F, prob = prob_vec)
}