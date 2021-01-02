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
  
  expect_true(all(dim(mat) == c(nrow(res$u), nrow(res$v))))
  expect_true(length(res$d) == 2)
  expect_true(all(sort(names(res)) == sort(c("d", "d_original", "u", "v"))))
})

################################

## .compute_cca_aggregate_matrix is correct

test_that(".compute_cca_aggregate_matrix works", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  res <- .compute_cca_aggregate_matrix(svd_1, svd_2)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(length(svd_1$d), length(svd_2$d))))
})

test_that(".compute_cca_aggregate_matrix is correct", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 20; p1 <- 50; p2 <- 40
    mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    res <- .compute_cca_aggregate_matrix(svd_1, svd_2)
    
    cov_1 <- stats::cov(mat_1) * (n-1)/n
    cov_2 <- stats::cov(mat_2) * (n-1)/n
    cov_12 <- crossprod(mat_1, mat_2)/n
    
    r1 <- Matrix::rankMatrix(cov_1); r2 <- Matrix::rankMatrix(cov_2)
    eigen_1 <- eigen(cov_1); eigen_2 <- eigen(cov_2)
    cov_1_invhalf <- .mult_mat_vec(eigen_1$vectors[,1:r1], 1/sqrt(eigen_1$values[1:r1]))
    cov_2_invhalf <- .mult_mat_vec(eigen_2$vectors[,1:r2], 1/sqrt(eigen_2$values[1:r2]))
    
    res2 <-  t(cov_1_invhalf) %*% cov_12 %*% cov_2_invhalf
    
    sum(abs(abs(res) - abs(res2))) <= 1e-6
  })
 
  expect_true(all(bool_vec))
})

test_that(".compute_cca_aggregate_matrix is gives the correct svd", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 20; p1 <- 50; p2 <- 40
    mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    res1 <- .compute_cca_aggregate_matrix(svd_1, svd_2)
    svd_res1 <- svd(res1)
    
    cov_1 <- stats::cov(mat_1) * (n-1)/n
    cov_2 <- stats::cov(mat_2) * (n-1)/n
    cov_12 <- crossprod(mat_1, mat_2)/n
    
    r1 <- Matrix::rankMatrix(cov_1); r2 <- Matrix::rankMatrix(cov_2)
    eigen_1 <- eigen(cov_1); eigen_2 <- eigen(cov_2)
    cov_1_invhalf <- .mult_mat_vec(eigen_1$vectors[,1:r1], 1/sqrt(eigen_1$values[1:r1]))
    cov_2_invhalf <- .mult_mat_vec(eigen_2$vectors[,1:r2], 1/sqrt(eigen_2$values[1:r2]))
    
    res2 <-  t(cov_1_invhalf) %*% cov_12 %*% cov_2_invhalf
    svd_res2 <- svd(res2)
    
    bool1 <- sum(abs(tcrossprod(svd_res1$u) - tcrossprod(svd_res2$u))) <= 1e-6
    bool2 <- sum(abs(tcrossprod(svd_res1$v) - tcrossprod(svd_res2$v))) <= 1e-6
    bool3 <- sum(abs(svd_res1$d - svd_res2$d)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

#################################

## .cca is correct

test_that(".cca works", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  res <- .cca(svd_1, svd_2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("loading_1", "loading_2", "obj_vec"))))
  expect_true(all(dim(res$loading_1) == c(p1, length(svd_1$d))))
  expect_true(all(dim(res$loading_2) == c(p2, length(svd_2$d))))
  expect_true(length(res$obj_vec) <= min(c(length(svd_1$d), length(svd_2$d))))
})

test_that(".cca yields empirically uncorrelated canoncical scores", {
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
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    res <- .cca(svd_1, svd_2)
    
    bool1 <- sum(abs(t(res$loading_1) %*% (stats::cov(mat_1)*(n-1)/n) %*% res$loading_1 - diag(min(p1,p2)))) <= 1e-4
    bool2 <- sum(abs(t(res$loading_2) %*% (stats::cov(mat_2)*(n-1)/n) %*% res$loading_2 - diag(min(p1,p2)))) <= 1e-4
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that(".cca obtains the correlation that is the same value as the objective", {
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
    
    svd_1 <- .svd_truncated(mat_1, K); svd_2 <- .svd_truncated(mat_2, K)
    
    cca_res <- .cca(svd_1, svd_2)
    
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

test_that(".cca obtains the maximum correlation", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- .svd_truncated(mat_1, K); svd_2 <- .svd_truncated(mat_2, K)

  cca_res <- .cca(svd_1, svd_2)
  
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

#################################

## .compute_unnormalized_scores is correct

test_that(".compute_unnormalized_scores computes the correct scores", {
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
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    cca_res <- .cca(svd_1, svd_2)
    res <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    rank_val <- Matrix::rankMatrix(.compute_cca_aggregate_matrix(svd_1, svd_2))
    
    bool1 <- sum(crossprod(res$score_1) - n*diag(rank_val)) <= 1e-4
    bool2 <- sum(crossprod(res$score_2) - n*diag(rank_val)) <= 1e-4
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

##################################

## .dcca is correct

test_that(".dcca works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  res <- .dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, 
               apply_shrinkage = T, compute_common_factors = T, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_factors", "score_1", "score_2",
                                             "svd_1", "svd_2", "cca_obj"))))
  
})

#################################3

## dcca_factor is correct

test_that("dcca_factor works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca")
  expect_true(all(sort(names(res)) == sort(c("common_factors", "score_1", "score_2",
                                             "svd_1", "svd_2", "cca_obj"))))
  expect_true(all(dim(res$common_factors) == c(n, K)))
})

test_that("dcca_factor works with meta-cells", {
  set.seed(5)
  n <- 200; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  nc <- 20
  meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                     apply_shrinkage = F, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca_meta")
  expect_true(all(sort(names(res)) == sort(c("common_factors", "score_1", "score_2",
                                             "svd_1", "svd_2", "cca_obj", "meta_score_1",
                                             "meta_score_2", "meta_svd_1", "meta_svd_2"))))
  expect_true(all(dim(res$common_factors) == c(n, K)))
  expect_true(all(dim(res$meta_svd_1$u) == c(nc, K)))
  expect_true(all(dim(res$meta_svd_2$u) == c(nc, K)))
  expect_true(all(dim(res$svd_1$u) == c(n, K)))
  expect_true(all(dim(res$svd_2$u) == c(n, K)))
  expect_true(all(dim(res$meta_svd_1$v) == c(p1, K)))
  expect_true(all(dim(res$meta_svd_2$v) == c(p2, K)))
  expect_true(all(dim(res$svd_1$v) == c(p1, K)))
  expect_true(all(dim(res$svd_2$v) == c(p2, K)))
  expect_true(all(dim(res$meta_score_1) == c(p1, K)))
  expect_true(all(dim(res$meta_score_2) == c(p2, K)))
  expect_true(all(dim(res$score_1) == c(p1, K)))
  expect_true(all(dim(res$score_2) == c(p2, K)))
})

####################################

## dcca_decomposition is correct

test_that("dcca_decomposition yields uncorrelated distinct matrices", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_12 = K, verbose = F)
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("dcca_decomposition yields a low-rank matrix", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_12 = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("dcca_decomposition yields common matrices with the same column space", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_12 = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})


