compute_tiltedCCA_ingredients <- function(setting = 1){
  set.seed(10)
  # setting 1 has modality 2 having no distinct information
  if(setting == 1){
    n_clust <- 100
    high <- 0.9; low <- 0.05
    B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                       0.1, 0.9, 0.1,
                       0.1, 0.1, 0.9), 3, 3, byrow = T)
    K <- ncol(B_mat1)
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
    clustering_1 <- factor(membership_vec)
    clustering_2 <- factor(rep(1, length(membership_vec)))
    n <- length(membership_vec); true_membership_vec <- membership_vec
    svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
    svd_u_2 <- generate_random_orthogonal(n, 2, centered = T)
  
    p_1 <- 20; p_2 <- 40
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- generate_random_orthogonal(p_1, 2)
    svd_v_2 <- generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
  } else if(setting == 2){
    # setting 2 is where both modalities are the same
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    clustering_1 <- factor(rep(1:3, n_each))
    clustering_2 <- factor(rep(1:3, n_each))
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
    svd_v_1 <- generate_random_orthogonal(p_1, 2)
    svd_v_2 <- generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else if(setting == 3){
    # setting 3 is where modality 2 has information, but not as much
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    clustering_1 <- factor(c(rep(1, 2*n_each), rep(2, n_each)))
    clustering_2 <- factor(c(rep(1, n_each), rep(2, n_each), rep(1, n_each)))
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
    svd_v_1 <- generate_random_orthogonal(p_1, 2)
    svd_v_2 <- generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else {
    # setting 4 is two modalities with high distinct information
    n_each <- 100
    true_membership_vec <- rep(1:4, each = n_each)
    clustering_1 <- factor(c(rep(1, 2*n_each), rep(2, 2*n_each)))
    clustering_2 <- factor(c(rep(1, n_each), rep(2, n_each), 
                             rep(1, n_each), rep(2, n_each)))
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
    svd_v_1 <- generate_random_orthogonal(p_1, 2)
    svd_v_2 <- generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  }
  
  ############################
  rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
  rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
  colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
  colnames(mat_2) <- paste0("p", 1:ncol(mat_2))
  
  # compute relevant ingredients based on mat_1 and mat_2
  n <- nrow(mat_1)
  svd_1 <- tiltedCCA:::.svd_truncated(mat_1, K = 2, symmetric = F, rescale = F, 
                                      mean_vec = F, sd_vec = F, K_full_rank = F)
  svd_2 <- tiltedCCA:::.svd_truncated(mat_2, K = 2, symmetric = F, rescale = F, 
                                      mean_vec = F, sd_vec = F, K_full_rank = F)
  
  set.seed(10)
  snn_1 <- tiltedCCA:::.form_snn_mat(mat = tiltedCCA:::.mult_mat_vec(svd_1$u, svd_1$d),
                                     num_neigh = 10,
                                     bool_cosine = T,
                                     bool_intersect = T,
                                     min_deg = 1,
                                     verbose = F)
  snn_2 <- tiltedCCA:::.form_snn_mat(mat = tiltedCCA:::.mult_mat_vec(svd_2$u, svd_2$d),
                                     num_neigh = 10,
                                     bool_cosine = T,
                                     bool_intersect = T,
                                     min_deg = 1,
                                     verbose = F)
  
  set.seed(10)
  common_mat <- tiltedCCA:::.compute_common_snn(snn_mat_1 = snn_1, 
                                                snn_mat_2 = snn_2,
                                                clustering_1 = clustering_1, 
                                                clustering_2 = clustering_2,
                                                num_neigh = 10,
                                                verbose = F)
  
  target_dimred <- tiltedCCA:::.compute_laplacian_basis(common_mat, 
                                                        k = 2, 
                                                        verbose = F)
  
  svd_1 <- tiltedCCA:::.check_svd(svd_1, dims = c(1:2))
  svd_2 <- tiltedCCA:::.check_svd(svd_2, dims = c(1:2))
  
  cca_res <- tiltedCCA:::.cca(svd_1, svd_2, 
                  dims_1 = NA, dims_2 = NA, 
                  return_scores = F)
  
  tmp <- tiltedCCA:::.compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  metacell_clustering <- lapply(1:nrow(mat_1), function(i){i})
  averaging_mat <- tiltedCCA:::.generate_averaging_matrix(n, metacell_clustering)
  
  res <- .common_decomposition(averaging_mat = averaging_mat,
                               discretization_gridsize = 9,
                               enforce_boundary = T,
                               fix_tilt_perc = 0.5,
                               score_1 = score_1,
                               score_2 = score_2,
                               snn_bool_intersect = T,
                               snn_k = 2,
                               snn_min_deg = 1,
                               snn_num_neigh = 10,
                               svd_1 = svd_1, 
                               svd_2 = svd_2,
                               target_dimred = target_dimred)
  common_score <- res$common_score
  
  list(averaging_mat = averaging_mat,
       cca_res_obj = cca_res$obj_vec,
       common_score = common_score,
       K = 2,
       mat_1 = mat_1,
       mat_2 = mat_2,
       score_1 = score_1,
       score_2 = score_2,
       svd_1 = svd_1,
       svd_2 = svd_2,
       target_dimred = target_dimred,
       true_membership_vec = as.factor(true_membership_vec))
}

for(setting in 1:4){
  test_data <- compute_tiltedCCA_ingredients(setting = setting)
  save(test_data, file = paste0("tests/assets/test_data", setting, ".RData"))
}
