context("Test heteroPCA")

## heteroPCA is correct

test_that("heteroPCA works", {
  set.seed(10)
  p <- 50; n <- 20
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 0.5
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- 1
  latent_mat <- log(abs(MASS::mvrnorm(n = n, mu = rep(5, p), Sigma = cov_mat)))
  obs_mat <- matrix(NA, nrow = n, ncol = p) 
  for(i in 1:n){
    for(j in 1:p){
      obs_mat[i,j] <- stats::rpois(n = 1, lambda = exp(latent_mat[i,j]))
    }
  }
  
  res <- heteroPCA(obs_mat, K = 2)
  
  expect_true(all(dim(res) == c(p,2)))
})

test_that("heteroPCA improves upon the typical covariance estimator", {
  set.seed(10)
  p <- 100; n <- 20; K <- 2
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  
  trials <- 20
  result_mat <- sapply(1:trials, function(x){
    set.seed(x)
    latent_mat <- abs(MASS::mvrnorm(n = n, mu = c(200,rep(30, p/2-1), rep(5, p/2)), Sigma = cov_mat))
    
    obs_mat <- matrix(NA, nrow = n, ncol = p)
    for(i in 1:n){
      for(j in 1:p){
        obs_mat[i,j] <- stats::rpois(n = 1, lambda = latent_mat[i,j])
      }
    }
    
    target_mat <- stats::cov(latent_mat)
    
    heteroPCA_res <- heteroPCA(obs_mat, K = K)
    naive_res <- stats::cov(obs_mat)
    
    c(max(eigen(heteroPCA_res - target_mat)$values), max(eigen(naive_res - target_mat)$values))
  })
  
  expect_true(median(result_mat[1,]-result_mat[2,]) < 0)
})
