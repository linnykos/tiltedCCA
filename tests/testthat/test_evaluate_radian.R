context("Test .evaluate_radian")

compute_evaluate_radian_ingredients <- function(setting = 1){
  # setting 1 has modality 2 having no distinct information
  if(setting == 1){
    n_clust <- 100
    high <- 0.9; low <- 0.05
    B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                       0.1, 0.9, 0.1,
                       0.1, 0.1, 0.9), 3, 3, byrow = T)
    K <- ncol(B_mat1)
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
    n <- length(membership_vec); true_membership_vec <- membership_vec
    svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
    svd_u_2 <- multiomicCCA::generate_random_orthogonal(n, 2, centered = T)
    
    metacell_clustering_1 <- factor(membership_vec)
    metacell_clustering_2 <- factor(rep("A", n))
    
    p_1 <- 20; p_2 <- 40
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
    
  } else if(setting == 2){
    # setting 2 is where both modalities are the same
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:3, function(i){
      if(i == 1){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else if(i == 2){
        MASS::mvrnorm(n = n_each, mu = c(0,12), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:3, function(i){
      if(i == 1){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else if(i == 2){
        MASS::mvrnorm(n = n_each, mu = c(0,12), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    metacell_clustering_1 <- factor(true_membership_vec)
    metacell_clustering_2 <- factor(true_membership_vec)
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else if(setting == 3){
    # setting 3 is where modality 2 has information, but not as much
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:3, function(i){
      if(i %in% c(1,2)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(9,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:3, function(i){
      if(i %in% c(1,3)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(3,0), Sigma = diag(2)) 
      }
    }))
    
    metacell_clustering_1 <- factor(c(rep("A", 2*n_each), rep("B", n_each)))
    metacell_clustering_2 <- factor(c(rep("A", n_each), rep("B", n_each), rep("A", n_each)))
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else {
    # setting 4 is two modalities with high distinct information
    n_each <- 100
    true_membership_vec <- rep(1:4, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:4, function(i){
      if(i %in% c(1,2)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:4, function(i){
      if(i %in% c(1,3)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    metacell_clustering_1 <- factor(c(rep("A", 2*n_each), rep("B", 2*n_each)))
    metacell_clustering_2 <- factor(c(rep("A", n_each), rep("B", n_each), 
                                      rep("A", n_each), rep("B", n_each)))
    
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
  
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  }
  
  svd_1 <- .svd_truncated(mat_1, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = c(1:2))
  svd_2 <- .check_svd(svd_2, dims = c(1:2))
  
  cca_res <- .cca(svd_1, svd_2, 
                  dims_1 = NA, dims_2 = NA, 
                  return_scores = F)
  
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  r <- min(ncol(score_1), ncol(score_2))
  basis_list <- lapply(1:r, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  circle_list <- lapply(1:r, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  min_subspace <- compute_min_subspace(dimred_1 = .mult_mat_vec(svd_1$u, svd_1$d),
                                       dimred_2 = .mult_mat_vec(svd_2$u, svd_2$d),
                                       metacell_clustering_1 = metacell_clustering_1,
                                       metacell_clustering_2 = metacell_clustering_2,
                                       binarize = F,
                                       num_neigh = 30,
                                       verbose = F)
  
  list(basis_list = basis_list, 
       score_1 = score_1,
       score_2 = score_2,
       metacell_clustering_1 = metacell_clustering_1,
       metacell_clustering_2 = metacell_clustering_2,
       svd_1 = svd_1, 
       svd_2 = svd_2,
       circle_list = circle_list,
       min_mat = min_subspace$min_mat,
       target_subspace = min_subspace$subspace_mat)
}

## .evaluate_radian is correct

test_that(".evaluate_radian works", {
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients(setting = 1)
  basis_list <- tmp$basis_list; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  target_subspace <- tmp$target_subspace
  
  res <- .evaluate_radian(basis_list = basis_list, 
                          circle_list = circle_list,
                          enforce_boundary = T,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          n_idx = 1:nrow(score_1),
                          num_neigh = 30,
                          percentage = 0.5,
                          return_common_score = F,
                          score_1 = score_1,
                          score_2 = score_2,
                          svd_1 = svd_1, 
                          svd_2 = svd_2,
                          target_subspace = target_subspace)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  
  res <- .evaluate_radian(basis_list = basis_list, 
                          circle_list = circle_list,
                          enforce_boundary = T,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          n_idx = 1:nrow(score_1),
                          num_neigh = 30,
                          percentage = 0.5,
                          return_common_score = T,
                          score_1 = score_1,
                          score_2 = score_2,
                          svd_1 = svd_1, 
                          svd_2 = svd_2,
                          target_subspace = target_subspace)
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(nrow(score_1), 2)))
})

test_that(".evaluate_radian is maximized at appropriate values", {
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients(setting = 1)
  basis_list <- tmp$basis_list; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  min_mat <- as.matrix(tmp$min_mat)
  # image(min_mat, asp = T)
  target_subspace <- tmp$target_subspace
  # plot(target_subspace[,1], target_subspace[,2], asp = T, col = rep(1:3, each = 100))
  res <- sapply(seq(0, 1, length.out = 9), function(percentage){
    .evaluate_radian(basis_list = basis_list, 
                     circle_list = circle_list,
                     enforce_boundary = T,
                     metacell_clustering_1 = metacell_clustering_1,
                     metacell_clustering_2 = metacell_clustering_2,
                     n_idx = 1:nrow(score_1),
                     num_neigh = 30,
                     percentage = percentage,
                     return_common_score = F,
                     score_1 = score_1,
                     score_2 = score_2,
                     svd_1 = svd_1, 
                     svd_2 = svd_2,
                     target_subspace = target_subspace)
  })
  expect_true(length(res) == 9)
  expect_true(res[1] <= min(res[5:9]))
  
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients(setting = 3)
  basis_list <- tmp$basis_list; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  min_mat <- as.matrix(tmp$min_mat)
  target_subspace <- tmp$target_subspace
  # plot(target_subspace[,1], target_subspace[,2], asp = T, col = rep(1:3, each = 100))
  res <- sapply(seq(0, 1, length.out = 9), function(percentage){
    .evaluate_radian(basis_list = basis_list, 
                     circle_list = circle_list,
                     enforce_boundary = T,
                     metacell_clustering_1 = metacell_clustering_1,
                     metacell_clustering_2 = metacell_clustering_2,
                     n_idx = 1:nrow(score_1),
                     num_neigh = 30,
                     percentage = percentage,
                     return_common_score = F,
                     score_1 = score_1,
                     score_2 = score_2,
                     svd_1 = svd_1, 
                     svd_2 = svd_2,
                     target_subspace = target_subspace)
  })
  expect_true(max(res[1:4]) <= min(res[5:9]))
  
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients(setting = 4)
  basis_list <- tmp$basis_list; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  min_mat <- as.matrix(tmp$min_mat)
  target_subspace <- tmp$target_subspace
  res <- sapply(seq(0, 1, length.out = 9), function(percentage){
    .evaluate_radian(basis_list = basis_list, 
                     circle_list = circle_list,
                     enforce_boundary = T,
                     metacell_clustering_1 = metacell_clustering_1,
                     metacell_clustering_2 = metacell_clustering_2,
                     n_idx = 1:nrow(score_1),
                     num_neigh = 30,
                     percentage = percentage,
                     return_common_score = F,
                     score_1 = score_1,
                     score_2 = score_2,
                     svd_1 = svd_1, 
                     svd_2 = svd_2,
                     target_subspace = target_subspace)
  })
  expect_true(min(res[4:6]) <= min(res[1:3]))
  expect_true(min(res[4:6]) <= min(res[7:9]))
})

