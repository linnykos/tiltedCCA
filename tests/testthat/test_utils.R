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

#######################################

## .reparameterize is correct

test_that(".reparameterize works", {
  set.seed(10)
  mat_1 <- matrix(rnorm(200), 20, 5)
  mat_2 <- matrix(rnorm(200), 20, 5)
  
  res <- .reparameterize(mat_1, mat_2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2", "diag_vec"))))
})

test_that(".reparameterize yields diagonal crossproducts", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    mat_1 <- matrix(rnorm(200), 20, 5)
    mat_2 <- matrix(rnorm(200), 20, 5)
    
    res <- .reparameterize(mat_1, mat_2)
    
    prod_mat <- crossprod(res$mat_1, res$mat_2)
    
    abs(sum(abs(prod_mat)) - sum(abs(diag(prod_mat)))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that(".reparameterize yields diagonal covariances", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    mat_1 <- matrix(rnorm(200), 20, 5)
    mat_2 <- matrix(rnorm(200), 20, 5)
    
    res <- .reparameterize(mat_1, mat_2)
    
    cov_1 <- crossprod(res$mat_1)
    cov_2 <- crossprod(res$mat_2)
    
    bool1 <- abs(sum(abs(cov_1)) - sum(abs(diag(cov_1)))) <= 1e-6
    bool2 <- abs(sum(abs(cov_2)) - sum(abs(diag(cov_2)))) <= 1e-6
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that(".reparameterize preserves the column space", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    mat_1 <- matrix(rnorm(200), 20, 5)
    mat_2 <- matrix(rnorm(200), 20, 5)
    
    res <- .reparameterize(mat_1, mat_2)
    weight_1 <- runif(ncol(res$mat_1)); weight_2 <- runif(ncol(res$mat_2))
    
    vec_1 <- res$mat_1 %*% weight_1; vec_2 <- res$mat_2 %*% weight_2
    
    # if vec_1 were in the column space of mat_1, then the vector should be unchanged
    #  after projection
    tmp_1 <- mat_1 %*% solve(crossprod(mat_1)) %*% t(mat_1) %*% vec_1
    tmp_2 <- mat_2 %*% solve(crossprod(mat_2)) %*% t(mat_2) %*% vec_2
    
    sum(abs(vec_1 - tmp_1)) <= 1e-6 & sum(abs(vec_2 - tmp_2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})


test_that(".reparameterize preserves the column space (2nd test)", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    mat_1 <- matrix(rnorm(200), 20, 5)
    mat_2 <- matrix(rnorm(200), 20, 5)
    
    res <- .reparameterize(mat_1, mat_2)
    
    # if the column space were the same, the projection matrices should be the same
    proj_1 <- mat_1 %*% solve(crossprod(mat_1)) %*% t(mat_1)
    proj_2 <- mat_2 %*% solve(crossprod(mat_2)) %*% t(mat_2)
    proj_1b <- res$mat_1 %*% solve(crossprod(res$mat_1)) %*% t(res$mat_1)
    proj_2b <- res$mat_2 %*% solve(crossprod(res$mat_2)) %*% t(res$mat_2)
    
    sum(abs(proj_1 - proj_1b)) <= 1e-6 & sum(abs(proj_2 - proj_2b)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})
