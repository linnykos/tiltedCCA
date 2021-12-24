context("Test clisi information")

compute_clisi_ingredients <- function(setting = 1){
  # setting 1 has modality 2 having no distinct information
  if(setting == 1){
    n_clust <- 100
    high <- 0.9; low <- 0.05
    B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                       0.1, 0.9, 0.1,
                       0.1, 0.1, 0.9), 3, 3, byrow = T)
    true_membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
    n <- length(true_membership_vec)
    svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, true_membership_vec, centered = T)[,1:2]
    svd_u_2 <- multiomicCCA::generate_random_orthogonal(n, 2, centered = T)
    
    p_1 <- 20; p_2 <- 20
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
    
    metacell_clustering_1 <- as.factor(true_membership_vec)
    metacell_clustering_2 <- as.factor(rep(1, n))
    
  } else if(setting == 2){
    # setting 1 is two modalities with high distinct information
    n_each <- 400
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
    
    metacell_clustering_1 <- as.factor(c(rep(1, each = 2*n_each), rep(2, each = 2*n_each)))
    metacell_clustering_2 <- as.factor(c(rep(1, each = n_each), 
                                         rep(2, each = n_each),
                                         rep(1, each = n_each), 
                                         rep(2, each = n_each)))
  }
  
  K <- 2
  svd_1 <- .svd_truncated(mat_1, K = K, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = K, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = c(1:K))
  svd_2 <- .check_svd(svd_2, dims = c(1:K))
  
  mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
  mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  
  list(metacell_clustering_1 = metacell_clustering_1,
       metacell_clustering_2 = metacell_clustering_2,
       true_membership_vec = as.factor(true_membership_vec),
       mat_1 = mat_1, mat_2 = mat_2, K = K)
}


## .clisi is correct

test_that(".clisi works", {
  set.seed(5)
  tmp <- compute_clisi_ingredients(setting = 1)
  mat_1 <- tmp$mat_1; mat_2 <- tmp$mat_2; K <- tmp$K
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  true_membership_vec <- tmp$true_membership_vec
  
  dcca_obj <- dcca_factor(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          fix_tilt_perc = F,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          verbose = F)
  list_g <- construct_frnn(dcca_obj, nn = 25, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
  res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  
  expect_true(length(res1) == 3)
  expect_true(all(sort(names(res1)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
  
  expect_true(all(sort(colnames(res1$df)) == sort(c("celltype", "value", "sd"))))
  expect_true(all(dim(res1$df) == c(3,3)))
  expect_true(all(dim(res1$clisi_mat) == c(3,3)))
  expect_true(all(dim(res1$clisi_cell_mat) == c(3,nrow(mat_1))))
  expect_true(all(res1$clisi_cell_mat >= 0))
  expect_true(all(res1$clisi_cell_mat <= 1))
  
  expect_true(all(res1$df[,"value"] <= 0.1))
  expect_true(all(res2$df[,"value"] >= 0.9))
})


test_that(".clisi works for a different setting", {
  set.seed(5)
  tmp <- compute_clisi_ingredients(setting = 2)
  mat_1 <- tmp$mat_1; mat_2 <- tmp$mat_2; K <- tmp$K
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  true_membership_vec <- tmp$true_membership_vec
  
  dcca_obj <- dcca_factor(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          fix_tilt_perc = F,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          verbose = F)
  list_g <- construct_frnn(dcca_obj, nn = 25, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
  res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  
  expect_true(all(res1$df[,"value"] <= 0.1))
  expect_true(all(res2$df[,"value"] >= 0.1))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  set.seed(5)
  tmp <- compute_clisi_ingredients(setting = 1)
  mat_1 <- tmp$mat_1; mat_2 <- tmp$mat_2; K <- tmp$K
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  true_membership_vec <- tmp$true_membership_vec
  
  dcca_obj <- dcca_factor(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          fix_tilt_perc = F,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          verbose = F)
  list_g <- construct_frnn(dcca_obj, nn = 25, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  res1 <- clisi_information(list_g$c_g, list_g$d_g, 
                            membership_vec = true_membership_vec,
                            verbose = F)
  
  expect_true(class(res1) == "clisi")
  expect_true(is.list(res1))
  expect_true(all(sort(names(res1)) == sort(c("common_clisi", "distinct_clisi"))))
  expect_true(all(sort(names(res1$common_clisi)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
})

#########################

## .clisi_cell is correct

test_that(".clisi_cell works", {
  set.seed(5)
  tmp <- compute_clisi_ingredients(setting = 1)
  mat_1 <- tmp$mat_1; mat_2 <- tmp$mat_2; K <- tmp$K
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  true_membership_vec <- tmp$true_membership_vec
  
  dcca_obj <- dcca_factor(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          fix_tilt_perc = F,
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          verbose = F)
  list_g <- construct_frnn(dcca_obj, nn = 25, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res <- .clisi_cell(g, 
                     membership_vec = true_membership_vec, 
                     position = 1)
  expect_true(length(res) == 3)
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
})
