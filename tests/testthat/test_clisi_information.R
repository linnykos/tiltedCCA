context("Test clisi information")

## .compute_radius is correct

test_that(".compute_radius works", {
  set.seed(10)
  mat <- matrix(rnorm(100), 20, 5)
  res <- .compute_radius(mat, 2, radius_quantile = 0.95, cell_subidx = 1:20)
  
  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
  
  res <- .compute_radius(mat, 2, radius_quantile = 0.95, cell_subidx = sample(1:20,10))
  
  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
})

test_that(".compute_radius increases as neighbor increases", {
  set.seed(10)
  n <- 100; p <- 5
  mat <- matrix(rnorm(n*p), n, p)
  k_vec <- 5:50
  rad_vec <- sapply(k_vec, function(k){
    .compute_radius(mat, k, radius_quantile = 0.95, cell_subidx = 1:n)
  })
  
  expect_true(all(diff(rad_vec) >= 0))
})

############################

## .construct_frnn is correct

test_that(".construct_frnn works", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  rad <- .compute_radius(mat, 5, radius_quantile = 0.5, cell_subidx = 1:n)
  res <- .construct_frnn(mat, radius = rad, nn = 5, frnn_approx = 0)
  
  expect_true(is.list(res))
  expect_true(length(res) == n)
  for(i in 1:n){
    expect_true(length(res[[i]]) > 0)
    expect_true(all(res[[i]] %% 1 == 0))
    expect_true(!i %in% res[[i]])
  }
})

######################

## .clisi is correct

test_that(".clisi works", {
  set.seed(10)
  n <- 100; p <- 2
  mat1 <- matrix(rnorm(n*p), n, p)
  mat2 <- matrix(rnorm(n*p), n, p)+50
  mat <- rbind(mat1, mat2)
  rad <- .compute_radius(mat, 10, radius_quantile = 0.95)
  g <- .construct_frnn(mat, rad, nn = 10, frnn_approx = 0)
  membership_vec1 <- as.factor(rep(c("a","b"), each = n))
  membership_vec2 <- as.factor(rep(c("a","b"), times = n))
  cell_subidx1 <- .construct_celltype_subsample(membership_vec1, 50)
  cell_subidx2 <- .construct_celltype_subsample(membership_vec2, 50)
  
  res1 <- .clisi(g, membership_vec1, cell_subidx1, full_idx = 1:length(g))
  res2 <- .clisi(g, membership_vec2, cell_subidx2, full_idx = 1:length(g))
  
  expect_true(length(res1) == 2)
  expect_true(all(sort(names(res1)) == sort(c("cell_info", "membership_info"))))
  expect_true(all(sort(names(res1$cell_info)) == sort(c("idx", "celltype", "len", "in_ratio", "clisi_score"))))
  expect_true(all(sapply(1:ncol(res1$cell_info), function(j){class(res1$cell_info[,j])}) == c("integer", "factor", "numeric", "numeric", "numeric")))
  expect_true(all(sapply(1:ncol(res1$membership_info), function(j){class(res1$membership_info[,j])}) == c("character", rep("numeric", 6))))
  expect_true(all(sort(names(res1$membership_info)) == sort(c("celltype", "mean_len", "mean_ratio", "mean_clisi", 
                                                              "sd_len", "sd_ratio", "sd_clisi"))))
  expect_true(all(dim(res1$cell_info) == c(length(cell_subidx1), 5)))
  expect_true(all(res1$membership_info$mean_ratio >= res2$membership_info$mean_ratio))
  expect_true(all(res1$membership_info$mean_clisi >= res2$membership_info$mean_clisi))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  set.seed(10)
  n_clust <- 100
  B_mat1 <- matrix(c(0.9, 0, 0, 
                     0, 0.9, 0,
                     0, 0, 0.9), 3, 3, byrow = T)
  B_mat2 <- matrix(c(0.9, 0.85, 0, 
                     0.85, 0.9, 0,
                     0, 0, 1), 3, 3, byrow = T)
  K <- ncol(B_mat1)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2, 
                                     noise_val = 0.1)
  K <- 2
  dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, verbose = F)
  dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  res <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                                  membership_vec = as.factor(membership_vec),
                                  rank_c = K, rank_d = K, nn = 50, frnn_approx = 0, 
                                  verbose = F)
  
  expect_true(class(res) == "clisi")
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_clisi", "distinct_clisi", "everything_clisi"))))
  expect_true(all(dim(res$common_clisi$cell_info) == dim(res$distinct_clisi$cell_info)))
  expect_true(all(dim(res$common_clisi$cell_info) == dim(res$everything_clisi$cell_info)))
  
  expect_true(all(dim(res$common_clisi$membership_info) == dim(res$distinct_clisi$membership_info)))
  expect_true(all(dim(res$common_clisi$membership_info) == dim(res$everything_clisi$membership_info)))
})

test_that("clisi_information works with max_subsample_frnn and max_subsample_clisi", {
  set.seed(10)
  n_clust <- 100
  B_mat1 <- matrix(c(0.9, 0, 0, 
                     0, 0.9, 0,
                     0, 0, 0.9), 3, 3, byrow = T)
  B_mat2 <- matrix(c(0.9, 0.85, 0, 
                     0.85, 0.9, 0,
                     0, 0, 1), 3, 3, byrow = T)
  K <- ncol(B_mat1)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2, 
                                     noise_val = 0.1)
  K <- 2
  dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, verbose = F)
  dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  res <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                           membership_vec = as.factor(membership_vec),
                           rank_c = K, rank_d = K, nn = 50, frnn_approx = 0, 
                           max_subsample_clisi = 50,
                           verbose = F)
  expect_true(class(res) == "clisi")
  
  set.seed(10)
  res1 <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                           membership_vec = as.factor(membership_vec),
                           rank_c = K, rank_d = K, nn = 50, frnn_approx = 0, 
                           max_subsample_frnn = 50, max_subsample_clisi = 20,
                           verbose = F)
  expect_true(class(res1) == "clisi")
  
  set.seed(10)
  res2 <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                           membership_vec = as.factor(membership_vec),
                           rank_c = K, rank_d = K, nn = 50, frnn_approx = 0, 
                           max_subsample_frnn = 50, max_subsample_clisi = 50,
                           verbose = F)
  expect_true(class(res2) == "clisi")
  expect_true(all(res1$common_clisi$cell_info$idx %in% res2$common_clisi$cell_info$idx))
})
