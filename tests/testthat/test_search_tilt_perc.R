context("Test .search_tilt_perc")

compute_search_tilt_perc_ingredients <- function(setting = 1){
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
    
    p_1 <- 20; p_2 <- 40
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
    
    metacell_clustering_1 <- factor(membership_vec)
    metacell_clustering_2 <- factor(rep("A", n))
    
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
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
    
    metacell_clustering_1 <- factor(true_membership_vec)
    metacell_clustering_2 <- factor(true_membership_vec)
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
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
    
    metacell_clustering_1 <- factor(c(rep("A", 2*n_each), rep("B", n_each)))
    metacell_clustering_2 <- factor(c(rep("A", n_each), rep("B", n_each), rep("A", n_each)))
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
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
    
    metacell_clustering_1 <- factor(c(rep("A", 2*n_each), rep("B", 2*n_each)))
    metacell_clustering_2 <- factor(c(rep("A", n_each), rep("B", n_each), 
                                      rep("A", n_each), rep("B", n_each)))
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
  
  basis_list <- lapply(1:ncol(score_1), function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  circle_list <- lapply(1:ncol(score_1), function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  list(metacell_clustering_1 = metacell_clustering_1,
       metacell_clustering_2 = metacell_clustering_2,
       score_1 = score_1,
       score_2 = score_2,
       svd_1 = svd_1,
       svd_2 = svd_2,
       circle_list = circle_list,
       basis_list = basis_list)
}

test_that(".search_tilt_perc works", {
  set.seed(10)
  tmp <- compute_search_tilt_perc_ingredients(setting = 1)
  score_1 <- tmp$score_1
  score_2 <- tmp$score_2
  circle_list <- tmp$circle_list
  basis_list <- tmp$basis_list
  svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  
  discretization_gridsize <- 9
  res <- .search_tilt_perc(basis_list = basis_list,
                           circle_list = circle_list,
                           discretization_gridsize = discretization_gridsize,
                           enforce_boundary = T,
                           metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n_idx = 1:nrow(score_1),
                           num_neigh = 30,
                           score_1 = score_1,
                           score_2 = score_2,
                           svd_1 = svd_1,
                           svd_2 = svd_2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("df", "percentage"))))
  expect_true(all(colnames(res$df) == c("percentage", "ratio_val")))
  expect_true(nrow(res$df) == discretization_gridsize)
  expect_true(length(res$percentage) == 1)
  expect_true(is.numeric(res$percentage))
})

test_that(".search_tilt_perc gives reasonable values across different settings", {
  set.seed(10)
  tmp <- compute_search_tilt_perc_ingredients(setting = 1)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  circle_list <- tmp$circle_list; basis_list <- tmp$basis_list
  svd_1 <- tmp$svd_1; svd_2 <- tmp$svd_2
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  res <- .search_tilt_perc(basis_list = basis_list,
                           circle_list = circle_list,
                           discretization_gridsize = 9,
                           enforce_boundary = T,
                           metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n_idx = 1:nrow(score_1),
                           num_neigh = 30,
                           score_1 = score_1,
                           score_2 = score_2,
                           svd_1 = svd_1,
                           svd_2 = svd_2)
  expect_true(res$percentage <= 0.25)
  
  set.seed(10)
  tmp <- compute_search_tilt_perc_ingredients(setting = 3)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  circle_list <- tmp$circle_list; basis_list <- tmp$basis_list
  svd_1 <- tmp$svd_1; svd_2 <- tmp$svd_2
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  res <- .search_tilt_perc(basis_list = basis_list,
                           circle_list = circle_list,
                           discretization_gridsize = 9,
                           enforce_boundary = T,
                           metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n_idx = 1:nrow(score_1),
                           num_neigh = 30,
                           score_1 = score_1,
                           score_2 = score_2,
                           svd_1 = svd_1,
                           svd_2 = svd_2)
  expect_true(res$percentage < 0.5)
  
  set.seed(10)
  tmp <- compute_search_tilt_perc_ingredients(setting = 4)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  circle_list <- tmp$circle_list; basis_list <- tmp$basis_list
  svd_1 <- tmp$svd_1; svd_2 <- tmp$svd_2
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  res <- .search_tilt_perc(basis_list = basis_list,
                           circle_list = circle_list,
                           discretization_gridsize = 9,
                           enforce_boundary = T,
                           metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n_idx = 1:nrow(score_1),
                           num_neigh = 30,
                           score_1 = score_1,
                           score_2 = score_2,
                           svd_1 = svd_1,
                           svd_2 = svd_2)
  expect_true(res$percentage >= 0.25)
  expect_true(res$percentage <= 0.75)
})

