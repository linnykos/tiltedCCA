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
