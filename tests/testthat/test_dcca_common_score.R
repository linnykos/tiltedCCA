context("Test .dcca_common_score")

compute_dcca_common_score_ingredients <- function(setting = 1){
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
  
  list(cca_res = cca_res,
       svd_1 = svd_1, 
       svd_2 = svd_2,
       mat_1 = mat_1,
       mat_2 = mat_2)
}


## .dcca_common_score is correct

test_that("(Basic) .dcca_common_score works", {
  set.seed(5)
  tmp <- compute_dcca_common_score_ingredients()
  cca_res <- tmp$cca_res; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; mat_1 <- tmp$mat_1
  mat_2 <- tmp$mat_2
  
  tmp <- .form_snns(num_neigh = 30, svd_1 = svd_1, svd_2 = svd_2)
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  
  res <- .dcca_common_score(cca_res = cca_res, 
                            cell_max = nrow(mat_1),
                            discretization_gridsize = 9,
                            enforce_boundary = T,
                            fix_tilt_perc = F, 
                            metacell_clustering_1 = metacell_clustering_1,
                            metacell_clustering_2 = metacell_clustering_2,
                            num_neigh = 30,
                            svd_1 = svd_1, 
                            svd_2 = svd_2, 
                            verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "df_percentage",
                                             "metacell_clustering_1", 
                                             "metacell_clustering_2", 
                                             "tilt_perc"))))
  expect_true(all(dim(res$common_score) == c(nrow(mat_1), 2)))
})

test_that("(Coding) .dcca_common_score preserves rownames and colnames", {
  set.seed(5)
  tmp <- compute_dcca_common_score_ingredients()
  cca_res <- tmp$cca_res; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; mat_1 <- tmp$mat_1
  mat_2 <- tmp$mat_2
  
  n <- nrow(mat_1)
  p1 <- ncol(mat_1); p2 <- ncol(mat_2)
  
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  rownames(svd_1$u) <- rownames(mat_1)
  rownames(svd_2$u) <- rownames(mat_2)
  rownames(svd_1$v) <- colnames(mat_1)
  rownames(svd_2$v) <- colnames(mat_2)
  
  tmp <- .form_snns(num_neigh = 30, svd_1 = svd_1, svd_2 = svd_2)
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  
  res <- .dcca_common_score(cca_res = cca_res, 
                            cell_max = nrow(mat_1),
                            discretization_gridsize = 9,
                            enforce_boundary = T,
                            fix_tilt_perc = F, 
                            metacell_clustering_1 = metacell_clustering_1,
                            metacell_clustering_2 = metacell_clustering_2,
                            num_neigh = 30,
                            svd_1 = svd_1, 
                            svd_2 = svd_2, 
                            verbose = F)
  
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
})

test_that("(Math) .dcca_common_score yields uncorrelated residuals", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    tmp <- compute_dcca_common_score_ingredients()
    cca_res <- tmp$cca_res; svd_1 <- tmp$svd_1
    svd_2 <- tmp$svd_2; mat_1 <- tmp$mat_1
    mat_2 <- tmp$mat_2
    
    tmp <- .form_snns(num_neigh = 30, svd_1 = svd_1, svd_2 = svd_2)
    metacell_clustering_1 <- tmp$metacell_clustering_1
    metacell_clustering_2 <- tmp$metacell_clustering_2
    
    res <- .dcca_common_score(cca_res = cca_res, 
                              cell_max = nrow(mat_1),
                              discretization_gridsize = 9,
                              enforce_boundary = T,
                              fix_tilt_perc = F, 
                              metacell_clustering_1 = metacell_clustering_1,
                              metacell_clustering_2 = metacell_clustering_2,
                              num_neigh = 30,
                              svd_1 = svd_1, 
                              svd_2 = svd_2, 
                              verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})
