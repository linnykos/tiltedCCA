context("Test utility")

## .mult_vec_mat is correct

test_that(".mult_vec_mat works", {
  trials <- 100
  n <- 10
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- round(100*runif(n))
    mat <- matrix(100*runif(n^2), n, n)
    res1 <- .mult_vec_mat(vec, mat)
    res2 <- diag(vec) %*% mat
    
    sum(abs(res1 - res2)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

########################

## .mult_mat_vec is correct

test_that(".mult_mat_vec works", {
  trials <- 100
  n <- 10
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- round(100*runif(n))
    mat <- matrix(100*runif(n^2), n, n)
    res1 <- .mult_mat_vec(mat, vec)
    res2 <- mat %*% diag(vec)
    
    sum(abs(res1 - res2)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that(".mult_mat_vec can be used to compute inv_onehalf", {
  trials <- 100
  p <- 6
  cov_mat <- matrix(0, p, p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  n <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- MASS::mvrnorm(n = 100, mu = rep(0, p), Sigma = cov_mat)
    mat <- stats::cov(dat)
    eigen_res <- eigen(mat)
    res <- .mult_mat_vec(eigen_res$vectors, 1/sqrt(eigen_res$values))
    
    sum(abs(crossprod(res, mat) %*% res - diag(p))) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})
