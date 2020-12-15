context("Test DCCA")

## .spoet is correct

test_that(".spoet works", {
  set.seed(10)
  p <- 100; n <- 20; K <- 2
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  mat <- abs(MASS::mvrnorm(n = n, mu = c(200,rep(30, p/2-1), rep(5, p/2)), Sigma = cov_mat))
  
  res <- .spoet(mat, K = 2)
  
  expect_true(all(dim(mat) == dim(res)))
  expect_true(Matrix::rankMatrix(res) == 2)
})

################################

## .cca is correct

test_that(".cca works", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  cov_1 <- stats::cov(mat_1); cov_2 <- stats::cov(mat_2); cov_12 <- crossprod(mat_1, mat_2)/n
  
  res <- .cca(cov_1, cov_2, cov_12, K = 2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("factor_1", "factor_2", "obj_vec"))))
  expect_true(all(dim(res$factor_1) == c(p1, 2)))
  expect_true(all(dim(res$factor_2) == c(p2, 2)))
  expect_true(length(res$obj_vec) == 2)
})

test_that(".cca yields empirically uncorrelated canoncical scores", {
  set.seed(10)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  cov_1 <- stats::cov(mat_1) * (n-1)/n; cov_2 <- stats::cov(mat_2) * (n-1)/n
  cov_12 <- crossprod(mat_1, mat_2)/n
  
  cca_res <- .cca(cov_1, cov_2, cov_12, K = K)
  
  score_1 <- mat_1 %*% cca_res$factor_1/sqrt(n)
  score_2 <- mat_2 %*% cca_res$factor_2/sqrt(n)
  
  expect_true(sum(abs(crossprod(score_1) - diag(K))) <= 1e-4)
  expect_true(sum(abs(crossprod(score_2) - diag(K))) <= 1e-4)
})

test_that(".cca obtains the correlation that is the same value as the objective", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  cov_1 <- stats::cov(mat_1) * (n-1)/n; cov_2 <- stats::cov(mat_2) * (n-1)/n
  cov_12 <- crossprod(mat_1, mat_2)/n
  
  cca_res <- .cca(cov_1, cov_2, cov_12, K = K)
  
  res1 <- t(cca_res$factor_1) %*% cov_12 %*% cca_res$factor_2
  
  score_1 <- mat_1 %*% cca_res$factor_1/sqrt(n)
  score_2 <- mat_2 %*% cca_res$factor_2/sqrt(n)
  
  res2 <- stats::cor(score_1, score_2)
  
  expect_true(abs(res1[1,2]) + abs(res1[2,1]) <= 1e-4)
  expect_true(abs(res2[1,2]) + abs(res2[2,1]) <= 1e-4)
  expect_true(sum(abs(res1 - res2)) <= 1e-4)
  expect_true(sum(abs(diag(res2) - cca_res$obj_vec)) <= 1e-4)
})

test_that(".cca obtains the maximum correlation", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  cov_1 <- stats::cov(mat_1) * (n-1)/n; cov_2 <- stats::cov(mat_2) * (n-1)/n
  cov_12 <- crossprod(mat_1, mat_2)/n
  
  cca_res <- .cca(cov_1, cov_2, cov_12, K = 1)
  
  score_1 <- mat_1 %*% cca_res$factor_1/sqrt(n)
  score_2 <- mat_2 %*% cca_res$factor_2/sqrt(n)
  
  target_cor <- stats::cor(score_1, score_2)
  
  # randomly generate some transformations
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    rand_factor_1 <- stats::runif(p1, min = -1, max = 1)
    rand_factor_2 <-stats::runif(p2, min = -1, max = 1)
    
    rand_score_1 <- mat_1 %*% rand_factor_1
    rand_score_2 <- mat_2 %*% rand_factor_2
    
    rand_cor <- stats::cor(rand_score_1, rand_score_2)
    
    rand_cor < target_cor
  })
  
  expect_true(all(bool_vec))
})

#################################3

## dcca is correct
test_that("dcca works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  res <- dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, rank_12 = K)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_factors", "common_mat_1", "common_mat_2",
                                             "distinct_mat_1", "distinct_mat_2"))))
  expect_true(all(dim(res$common_factors) == c(n, K)))
  expect_true(all(dim(res$common_mat_1) == dim(mat_1)))
  expect_true(all(dim(res$common_mat_2) == dim(mat_2)))
  expect_true(all(dim(res$distinct_mat_1) == dim(mat_1)))
  expect_true(all(dim(res$distinct_mat_2) == dim(mat_2)))
})

test_that("dcca yields uncorrelated distinct matrices", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_res <- dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, rank_12 = K)
  
  res <- crossprod(dcca_res$distinct_mat_1, dcca_res$distinct_mat_2)
  
  expect_true(sum(abs(res)) <= 1e-4)
})


