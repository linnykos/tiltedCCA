context("Test cca")

## .cca is correct

test_that("(Basic) .cca works", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
  svd_1 <- .check_svd(svd_1, dims = 1:min(dim(mat_1)))
  svd_2 <- .check_svd(svd_2, dims = 1:min(dim(mat_2)))
  class(svd_1) <- "svd"; class(svd_2) <- "svd"
  
  res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("loading_1", "loading_2", "obj_vec"))))
  expect_true(all(dim(res$loading_1) == c(p1, length(svd_1$d))))
  expect_true(all(dim(res$loading_2) == c(p2, length(svd_2$d))))
  expect_true(length(res$obj_vec) <= min(c(length(svd_1$d), length(svd_2$d))))
})

test_that("(Coding) .cca preserves colnames", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  colnames(mat_1) <- paste0("a", 1:p1)
  colnames(mat_2) <- paste0("b", 1:p2)
  
  res <- .cca(mat_1, mat_2, dims_1 = 1:2, dims_2 = 1:2, return_scores = F)
  
  expect_true(all(rownames(res$loading_1) == colnames(mat_1)))
  expect_true(all(rownames(res$loading_2) == colnames(mat_2)))
  
  ##
  
  svd_1 <- .svd_safe(mat = mat_1,
                     check_stability = T, 
                     K = min(dim(mat_1)), 
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_2 <- .svd_safe(mat = mat_2,
                     check_stability = T, 
                     K = min(dim(mat_2)), 
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_1 <- .check_svd(svd_1, dims = 1:min(dim(mat_1)))
  svd_2 <- .check_svd(svd_2, dims = 1:min(dim(mat_2)))
  
  res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  
  expect_true(all(rownames(res$loading_1) == colnames(mat_1)))
  expect_true(all(rownames(res$loading_2) == colnames(mat_2)))
})

test_that("(Coding) .cca works when the two matrices have different ranks larger than 1", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  rank_1 <- 2; rank_2 <- 4
  res <- .cca(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 1:rank_2, return_scores = F)
  
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 2)
  
  ####
  
  svd_1 <- .svd_safe(mat = mat_1,
                     check_stability = T, 
                     K = rank_1,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_2 <- .svd_safe(mat = mat_2,
                     check_stability = T, 
                     K = rank_2,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_1 <- .check_svd(svd_1, dims = 1:rank_1)
  svd_2 <- .check_svd(svd_2, dims = 1:rank_2)
  
  res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 2)
})

test_that("(Coding) .cca works when either matrix have rank 1 for matrix inputs", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  rank_1 <- 2; rank_2 <- 4
  
  res <- .cca(mat_1, mat_2, dims_1 = 1, dims_2= 1:rank_2, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 1, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(mat_1, mat_2, dims_1 = 1, dims_2 = 1, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
})

test_that("(Coding) .cca works when either matrix have rank 1 for svd inputs", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  rank_1 <- 2; rank_2 <- 4
  
  svd_1 <- .svd_safe(mat = mat_1,
                     check_stability = T, 
                     K = rank_1,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_2 <- .svd_safe(mat = mat_2,
                     check_stability = T, 
                     K = rank_2,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_1 <- .check_svd(svd_1, dims = 1:rank_1)
  svd_2 <- .check_svd(svd_2, dims = 1:rank_2)
  
  svd_1b <- .svd_safe(mat = mat_1,
                     check_stability = T, 
                     K = 1,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_2b <- .svd_safe(mat = mat_2,
                     check_stability = T, 
                     K = 1,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_1b <- .check_svd(svd_1b, dims = 1)
  svd_2b <- .check_svd(svd_2b, dims = 1)
  
  res <- .cca(svd_1b, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(svd_1, svd_2b, dims_1 = NA, dims_2 = NA, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(svd_1b, svd_2b, dims_1 = NA, dims_2 = NA, return_scores = F)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
})

test_that("(Math) .cca yields empirically uncorrelated canoncical scores", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 150; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    svd_1 <- .check_svd(svd_1, dims = 1:min(dim(mat_1)))
    svd_2 <- .check_svd(svd_2, dims = 1:min(dim(mat_2)))
    class(svd_1) <- "svd"; class(svd_2) <- "svd"
    
    res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    
    bool1 <- sum(abs(t(res$loading_1) %*% (stats::cov(mat_1)*(n-1)/n) %*% res$loading_1 - diag(p1))) <= 1e-4
    bool2 <- sum(abs(t(res$loading_2) %*% (stats::cov(mat_2)*(n-1)/n) %*% res$loading_2 - diag(p2))) <= 1e-4
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .cca obtains the correlation that is the same value as the objective", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 150; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
    
    svd_1 <- .svd_safe(mat = mat_1,
                       check_stability = T, 
                       K = K,
                       mean_vec = NULL, 
                       rescale = F, 
                       scale_max = NULL, 
                       sd_vec = NULL)
    svd_2 <- .svd_safe(mat = mat_2,
                       check_stability = T, 
                       K = K,
                       mean_vec = NULL, 
                       rescale = F, 
                       scale_max = NULL, 
                       sd_vec = NULL)
    
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    
    cov_12 <- crossprod(mat_1, mat_2)/n
    res1 <- t(cca_res$loading_1) %*% cov_12 %*% cca_res$loading_2 
    
    score_1 <- mat_1 %*% cca_res$loading_1
    score_2 <- mat_2 %*% cca_res$loading_2
    
    res2 <- stats::cor(score_1, score_2) 
    
    bool1 <- abs(res1[1,2]) + abs(res1[2,1]) <= 1e-4
    bool2 <- abs(res2[1,2]) + abs(res2[2,1]) <= 1e-4
    bool3 <- sum(abs(res1 - res2)) <= 1e-4
    bool4 <- sum(abs(diag(res2) - cca_res$obj_vec)) <= 1e-4
    
    bool1 & bool2 & bool3 & bool4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .cca obtains the maximum correlation", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- .svd_safe(mat = mat_1,
                     check_stability = T, 
                     K = K,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  svd_2 <- .svd_safe(mat = mat_2,
                     check_stability = T, 
                     K = K,
                     mean_vec = NULL, 
                     rescale = F, 
                     scale_max = NULL, 
                     sd_vec = NULL)
  
  cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  
  score_1 <- mat_1 %*% cca_res$loading_1
  score_2 <- mat_2 %*% cca_res$loading_2
  
  target_cor <- stats::cor(score_1, score_2)[1]
  
  # randomly generate some transformations
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    rand_loading_1 <- stats::runif(p1, min = -1, max = 1)
    rand_loading_2 <- stats::runif(p2, min = -1, max = 1)
    
    rand_score_1 <- mat_1 %*% rand_loading_1
    rand_score_2 <- mat_2 %*% rand_loading_2
    
    rand_cor <- stats::cor(rand_score_1, rand_score_2)
    
    rand_cor < target_cor
  })
  
  expect_true(all(bool_vec))
})
