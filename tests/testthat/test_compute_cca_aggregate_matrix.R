context("Test .compute_cca_aggregate_matrix")

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
  
  svd_1 <- .svd_truncated(mat_1, K = 2,
                          symmetric = F, rescale = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = 3,
                          symmetric = F, rescale = F, K_full_rank = F)
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
    svd_1 <- .svd_truncated(mat_1, K = rank_1,
                            symmetric = F, rescale = F, K_full_rank = F)
    svd_2 <- .svd_truncated(mat_2, K = rank_2,
                            symmetric = F, rescale = F, K_full_rank = F)
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

