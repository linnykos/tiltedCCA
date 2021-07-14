context("Test frNN")

## .compute_radius is correct

test_that(".compute_radius works", {
  set.seed(10)
  mat <- matrix(rnorm(100), 20, 5)
  res <- .compute_radius(mat, 2, radius_quantile = 0.95)
  
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
    .compute_radius(mat, k, radius_quantile = 0.95)
  })
  
  expect_true(all(diff(rad_vec) >= 0))
})

############################

## .construct_frnn is correct

test_that(".construct_frnn works", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  rad <- .compute_radius(mat, 5, radius_quantile = 0.5)
  res <- .construct_frnn(mat, radius = rad, nn = 5, frnn_approx = 0)
  
  expect_true(is.list(res))
  expect_true(length(res$id) == n)
  expect_true(length(res$dist) == n)
  for(i in 1:n){
    expect_true(length(res$id[[i]]) > 0)
    expect_true(all(res$id[[i]] %% 1 == 0))
    expect_true(!i %in% res$id[[i]])
  }
})

test_that(".construct_frnn works in the presence of outliers", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  mat <- rbind(mat, matrix(c(10,-10,10,-10), 2, 2))
  n <- nrow(mat)
  rad <- .compute_radius(mat, 5, radius_quantile = 0.5)
  res <- .construct_frnn(mat, radius = rad, nn = 5, frnn_approx = 0)
  
  expect_true(is.list(res))
  expect_true(length(res$id) == n)
  expect_true(length(res$dist) == n)
  for(i in 1:n){
    expect_true(length(res$id[[i]]) > 0)
    expect_true(all(res$id[[i]] %% 1 == 0))
    expect_true(!i %in% res$id[[i]])
  }
})

#####################

## construct_frnn is correct

test_that("construct_frnn works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space1 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))
  common_space2 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))+20
  common_space <- rbind(common_space1, common_space2)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_obj <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  
  membership_vec1 <- as.factor(rep(c("a","b"), each = n/2))
  res <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec1, 
                        verbose = F, bool_matrix = T)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("c_g", "d_g", "e_g", "membership_vec", "original_radius"))))
  expect_true(inherits(res$c_g, "dgCMatrix"))
  expect_true(inherits(res$d_g, "dgCMatrix"))
  expect_true(inherits(res$e_g, "dgCMatrix"))
  
  res <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec1, 
                        verbose = F, bool_matrix = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("c_g", "d_g", "e_g", "membership_vec", "original_radius"))))
  expect_true(inherits(res$c_g, "frNN"))
  expect_true(inherits(res$d_g, "frNN"))
  expect_true(inherits(res$e_g, "frNN"))
})

################################

## combine_frnn is correct

test_that("combine_frnn works", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.2, 0.1, 
                    0.2, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3, byrow = T)
  K <- ncol(B_mat); rho <- 1
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  
  set.seed(10)
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- generate_random_orthogonal(p_2, K-1)
  
  set.seed(10)
  dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  dcca_obj <- dcca_factor(dat$mat_1, dat$mat_2, dims_1 = 1:(K-1), dims_2 = 1:(K-1), 
                          verbose = F)
  membership_vec <- as.factor(membership_vec)
  set.seed(10)
  list_g_1 <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
                           data_1 = T, data_2 = F,
                           verbose = F, bool_matrix = T)
  set.seed(10)
  list_g_2 <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
                             data_1 = T, data_2 = F,
                             verbose = F, bool_matrix = T)
  
  set.seed(10)
  res <- combine_frnn(dcca_obj, list_g_1$c_g, list_g_2$c_g, nn = 5)
  
  expect_true(all(dim(res) == n))
  expect_true(inherits(res, "dgCMatrix"))
})
