context("Test D-CCA")

## .spoet is correct

test_that("(Basic) .spoet works", {
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

test_that("(Coding) .spoet preserves rownames and columnnames", {
  set.seed(10)
  p <- 100; n <- 20; K <- 2
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  mat <- abs(MASS::mvrnorm(n = n, mu = c(200,rep(30, p/2-1), rep(5, p/2)), Sigma = cov_mat))
  rownames(mat) <- paste0("a", 1:n)
  colnames(mat) <- paste0("b", 1:p)
  
  res <- .spoet(mat, K = 2)
  
  expect_true(all(rownames(res$u) == rownames(mat)))
  expect_true(all(rownames(res$v) == colnames(mat)))
})

################################

## .compute_cca_aggregate_matrix is correct

test_that("(Basic) .compute_cca_aggregate_matrix works", {
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

test_that("(Coding) .compute_cca_aggregate_matrix works when ranks are not the same", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  svd_1 <- .svd_truncated(mat_1, K = 2); svd_2 <- .svd_truncated(mat_2, K = 3)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  res <- .compute_cca_aggregate_matrix(svd_1, svd_2, augment = F)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(length(svd_1$d), length(svd_2$d))))
})

test_that("(Math) .compute_cca_aggregate_matrix preserves SVD values when augment = T", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 20; p1 <- 50; p2 <- 40
    mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    rank_1 <- sample(1:10,1); rank_2 <- sample(1:10,1); rank_full <- min(c(rank_1, rank_2))
    svd_1 <- .svd_truncated(mat_1, K = rank_1); svd_2 <- .svd_truncated(mat_2, K = rank_2)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    res1 <- .compute_cca_aggregate_matrix(svd_1, svd_2, augment = T)
    svd_res_1 <- svd(res1)
    svd_res_1$u <- svd_res_1$u[1:rank_1, 1:rank_full, drop  = F]
    svd_res_1$v <- svd_res_1$v[1:rank_2, 1:rank_full, drop  = F]
    
    res2 <- .compute_cca_aggregate_matrix(svd_1, svd_2, augment = F)
    svd_res_2 <- svd(res2)
    
    tmp <- svd(crossprod(svd_res_1$u, svd_res_2$u))
    svd_res_1$u <- svd_res_1$u %*% tcrossprod(tmp$u, tmp$v)
    
    tmp <- svd(crossprod(svd_res_1$v, svd_res_2$v))
    svd_res_1$v <- svd_res_1$v %*% tcrossprod(tmp$u, tmp$v)

    bool1 <- sum(abs(svd_res_1$u - svd_res_2$u)) <= 1e-6
    bool2 <- sum(abs(svd_res_1$v - svd_res_2$v)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .compute_cca_aggregate_matrix is correct", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 20; p1 <- 50; p2 <- 40
    mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    Sigma_mat <- matrix(1, p2, p2); diag(Sigma_mat) <- 2
    mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = Sigma_mat), center = T, scale = F)
    
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

test_that("(Math) .compute_cca_aggregate_matrix is gives the correct svd", {
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

test_that("(Basic) .cca works", {
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

test_that("(Coding) .cca preserves colnames", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  colnames(mat_1) <- paste0("a", 1:p1)
  colnames(mat_2) <- paste0("b", 1:p2)

  res <- .cca(mat_1, mat_2, rank_1 = 2, rank_2 = 2)
  
  expect_true(all(rownames(res$loading_1) == colnames(mat_1)))
  expect_true(all(rownames(res$loading_2) == colnames(mat_2)))
  
  ##
  
  svd_1 <- .svd_truncated(mat_1); svd_2 <- .svd_truncated(mat_2)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  res <- .cca(svd_1, svd_2)
  
  expect_true(all(rownames(res$loading_1) == colnames(mat_1)))
  expect_true(all(rownames(res$loading_2) == colnames(mat_2)))
})

test_that("(Coding) .cca works when the two matrices have different ranks larger than 1", {
  set.seed(10)
  n <- 20; p1 <- 50; p2 <- 40
  mat_1 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  rank_1 <- 2; rank_2 <- 4
  res <- .cca(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2)
  
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 2)
  
  ####
  
  svd_1 <- .svd_truncated(mat_1, K = rank_1); svd_2 <- .svd_truncated(mat_2, K = rank_2)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  res <- .cca(svd_1, svd_2)
  
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
  
  res <- .cca(mat_1, mat_2, rank_1 = 1, rank_2 = rank_2)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(mat_1, mat_2, rank_1 = rank_1, rank_2 = 1)
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(mat_1, mat_2, rank_1 = 1, rank_2 = 1)
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
  
  svd_1 <- .svd_truncated(mat_1, K = rank_1); svd_2 <- .svd_truncated(mat_2, K = rank_2)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  svd_1b <- .svd_truncated(mat_1, K = 1); svd_2b <- .svd_truncated(mat_2, K = 1)
  svd_1b <- .check_svd(svd_1b); svd_2b <- .check_svd(svd_2b)
  
  res <- .cca(svd_1b, svd_2)
  expect_true(all(dim(res$loading_1) == c(p1, 1)))
  expect_true(all(dim(res$loading_2) == c(p2, rank_2)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(svd_1, svd_2b)
  expect_true(all(dim(res$loading_1) == c(p1, rank_1)))
  expect_true(all(dim(res$loading_2) == c(p2, 1)))
  expect_true(length(res$obj_vec) == 1)
  
  res <- .cca(svd_1b, svd_2b)
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
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    res <- .cca(svd_1, svd_2)
    
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

test_that("(Math) .cca obtains the maximum correlation", {
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

test_that("(Math) .compute_unnormalized_scores computes the correct scores", {
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
    
    bool1 <- sum(crossprod(res$score_1) - n*diag(p1)) <= 1e-4
    bool2 <- sum(crossprod(res$score_2) - n*diag(p2)) <= 1e-4
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

###################################

## .dcca_common_score is correct

test_that("(Basic) .dcca_common_score works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  cca_res <- .cca(svd_1, svd_2)
  
  res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) .dcca_common_score preserves rownames and colnames", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2)
  
  res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
})

test_that("(Coding) .dcca_common_score works with K=1 for either", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  svd_1 <- .spoet(mat_1, K = 1); svd_2 <- .spoet(mat_2, K = K)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,1)))
  expect_true(all(dim(res$distinct_score_2) == c(n,K)))
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = 1)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,K)))
  expect_true(all(dim(res$distinct_score_2) == c(n,1)))
  
  svd_1 <- .spoet(mat_1, K = 1); svd_2 <- .spoet(mat_2, K = 1)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,1)))
  expect_true(all(dim(res$distinct_score_2) == c(n,1)))
})

test_that("(Math) .dcca_common_score yields uncorrelated residuals", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    
    svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
    
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    cca_res <- .cca(svd_1, svd_2)
    
    res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = F, verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .dcca_common_score yields uncorrelated residuals with meta-cells", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    
    svd_1 <- .svd_truncated(mat_1, K); svd_2 <- .svd_truncated(mat_2, K)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    num_meta <- max(meta_clustering)
    
    mat_1_meta <- t(sapply(1:num_meta, function(x){
      idx <- which(meta_clustering == x)
      apply(mat_1[idx,,drop = F], 2, mean)
    }))
    
    mat_2_meta <- t(sapply(1:num_meta, function(x){
      idx <- which(meta_clustering == x)
      apply(mat_2[idx,,drop = F], 2, mean)
    }))
    
    cca_res <- .cca(mat_1_meta, mat_2_meta, rank_1 = K, rank_2 = K)
    
    res <- .dcca_common_score(svd_1, svd_2, cca_res, check_alignment = T, verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

#################################3

## dcca_factor is correct

test_that("(Basic) dcca_factor works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2",
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) dcca_factor preserves rownames and colnames", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  expect_true(all(rownames(mat_1) == rownames(res$score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$score_2)))
  
  expect_true(all(rownames(mat_1) == rownames(res$svd_1$u)))
  expect_true(all(rownames(mat_1) == rownames(res$svd_2$u)))
  expect_true(all(colnames(mat_1) == rownames(res$svd_1$v)))
  expect_true(all(colnames(mat_2) == rownames(res$svd_2$v)))
})

test_that("(Coding) dcca_factor preserves rownames and colnames with metacells", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  nc <- 20
  meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                     apply_shrinkage = T, verbose = F)
  
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  expect_true(all(rownames(mat_1) == rownames(res$score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$score_2)))
  
  expect_true(all(rownames(mat_1) == rownames(res$svd_1$u)))
  expect_true(all(rownames(mat_1) == rownames(res$svd_2$u)))
  expect_true(all(colnames(mat_1) == rownames(res$svd_1$v)))
  expect_true(all(colnames(mat_2) == rownames(res$svd_2$v)))
})

test_that("(Coding) dcca_factor does not require rank_1 = rank_2", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = 2, rank_2 = 5, verbose = F)
  res2 <- dcca_factor(mat_2, mat_1, rank_1 = 5, rank_2 = 2, verbose = F)

  expect_true(length(res$cca_obj) == 2)
  expect_true(length(res2$cca_obj) == 2)
  expect_true(all(dim(res$distinct_score_1) == c(n, 2)))
  expect_true(all(dim(res$distinct_score_2) == c(n, 5)))
  expect_true(all(dim(res2$distinct_score_1) == c(n, 5)))
  expect_true(all(dim(res2$distinct_score_2) == c(n, 2)))
})

test_that("(Coding) dcca_factor works with meta-cells", {
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
                     apply_shrinkage = T, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

####################################

## dcca_decomposition is correct

test_that("(Basic) dcca_decomposition works", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "cca_obj", 
                                             "distinct_score_1", 
                                             "distinct_score_2",
                                             "common_mat_1", "common_mat_2",
                                             "distinct_mat_1", "distinct_mat_2"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) dcca_decomposition preserves rownames and colnames", {
  set.seed(1)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_2)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_2)))
  expect_true(all(colnames(mat_1) == colnames(res$common_mat_1)))
  expect_true(all(colnames(mat_1) == colnames(res$distinct_mat_1)))
  expect_true(all(colnames(mat_2) == colnames(res$common_mat_2)))
  expect_true(all(colnames(mat_2) == colnames(res$distinct_mat_2)))
})

test_that("(Math) dcca_decomposition yields uncorrelated distinct matrices", {
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
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields uncorrelated distinct matrices with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields a low-rank matrix", {
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
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K+1
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K+1
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields a low-rank matrix with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K+1
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K+1
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields common matrices with the same column space", {
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
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields common matrices with the same column space with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition is a decomposition under no noise", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1; mat_2 <- common_space %*% transform_mat_2 
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, 
                            apply_shrinkage = F, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    # correct common_space
    proj_1 <- common_space %*% solve(crossprod(common_space)) %*% t(common_space)
    Ktmp <- Matrix::rankMatrix(res$common_score)
    proj_2 <- res$common_score[,1:Ktmp] %*% solve(crossprod(res$common_score[,1:Ktmp])) %*% t(res$common_score[,1:Ktmp])
    bool1 <- sum(abs(proj_1 - proj_2)) <= 1e-6
    
    # no distinctive componet
    bool2 <- all(abs(res$distinct_mat_1) <= 1e-6)
    bool3 <- all(abs(res$distinct_mat_2) <= 1e-6)
    
    # exact decomposition
    bool4 <- all(abs(mat_1 - res$common_mat_1) <= 1e-6)
    bool5 <- all(abs(mat_2 - res$common_mat_2) <= 1e-6)
    
    bool1 & bool2 & bool3 & bool4 & bool5
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition is a decomposition under no noise with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1; mat_2 <- common_space %*% transform_mat_2 
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                            apply_shrinkage = F, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    # correct common_space
    proj_1 <- common_space %*% solve(crossprod(common_space)) %*% t(common_space)
    Ktmp <- Matrix::rankMatrix(res$common_score)
    proj_2 <- res$common_score[,1:Ktmp] %*% solve(crossprod(res$common_score[,1:Ktmp])) %*% t(res$common_score[,1:Ktmp])
    bool1 <- sum(abs(proj_1 - proj_2)) <= 1e-6
    
    # no distinctive componet
    bool2 <- all(abs(res$distinct_mat_1) <= 1e-6)
    bool3 <- all(abs(res$distinct_mat_2) <= 1e-6)
    
    # exact decomposition
    bool4 <- all(abs(mat_1 - res$common_mat_1) <= 1e-6)
    bool5 <- all(abs(mat_2 - res$common_mat_2) <= 1e-6)
    
    bool1 & bool2 & bool3 & bool4 & bool5
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition can obtain the same result when fed into itself (i.e., stability/identifiability)", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  set.seed(10)
  p_1 <- 20; p_2 <- 40
  coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
  coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)
  
  set.seed(10)
  dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

  dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K,
                          apply_shrinkage = F, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

  dcca_res2 <- dcca_factor(res$common_mat_1 + res$distinct_mat_1,
                           res$common_mat_2 + res$distinct_mat_2,
                           rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
  res2 <- dcca_decomposition(dcca_res2, rank_c = K, verbose = F)

  expect_true(sum(abs(res$common_mat_1 - res2$common_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$common_mat_2 - res2$common_mat_2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_1 - res2$distinct_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_2 - res2$distinct_mat_2)) <= 1e-6)
})

####################################

## .compute_common_score is correct

test_that("(Basic) .compute_common_score works", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  res <- .compute_common_score(score_1, score_2)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(score_1)))
})

test_that("(Coding) .compute_common_score can handle when the two scores have different dimensions", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  score_1b <- score_1[,1:2]
  res <- .compute_common_score(score_1b, score_2)
  expect_true(all(dim(res) == dim(score_1b)))
  
  score_2b <- score_2[,1:2]
  res <- .compute_common_score(score_1, score_2b)
  expect_true(all(dim(res) == dim(score_2b)))
})

test_that("(Coding) .compute_common_score works with some of them being K=1", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  score_1b <- score_1[, 1,drop = F]
  res <- .compute_common_score(score_1b, score_2)
  expect_true(all(dim(res) == dim(score_1b)))
  
  score_2b <- score_2[, 1,drop = F]
  res <- .compute_common_score(score_1, score_2b)
  expect_true(all(dim(res) == dim(score_2b)))
  
  res <- .compute_common_score(score_1b, score_2b)
  expect_true(all(dim(res) == dim(score_1b)))
})

test_that(".compute_common_score preserves rownames and colnames", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  rownames(score_1) <- paste0("a", 1:n)
  rownames(score_2) <- paste0("a", 1:n)
  
  res <- .compute_common_score(score_1, score_2)
  expect_true(all(dim(res) == dim(score_1)))
  expect_true(all(rownames(res) == rownames(score_1)))
  
  # observe that if rownames(score_1) != rownames(score_2), the names would still be rownames(score_1)
})

#####################################

## .compute_distinct_score is correct

test_that("(Basic) .compute_distinct_score works", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  common_score <- .compute_common_score(score_1, score_2)
  
  res <- .compute_distinct_score(score_1, score_2, common_score)
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_score_1", "distinct_score_2"))))
  expect_true(all(dim(res$common_score) == dim(common_score)))
  expect_true(all(dim(res$distinct_score_1) == dim(score_1)))
  expect_true(all(dim(res$distinct_score_2) == dim(score_2)))
})

test_that("(Math) .compute_distinct_score generates orthogonal distinct matrices", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    B_mat <- matrix(c(0.9, 0.4, 0.1,
                      0.4, 0.9, 0.1,
                      0.1, 0.1, 0.3), 3, 3)
    K <- ncol(B_mat); n_clust <- 50; rho <- 0.1
    
    true_membership_vec <- rep(1:4, each = n_clust)
    n <- length(true_membership_vec)
    
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
    score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
    score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    common_score <- .compute_common_score(score_1, score_2)
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    all(abs(crossprod(res$distinct_score_1, res$distinct_score_2)) <= 1e-6)
  })
  
  expect_true(all(bool_vec))
})














