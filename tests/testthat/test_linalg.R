context("Test linear algabra")

## .representation_2d is correct

test_that("(Math) .representation_2d has the correct representation", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(11*x)
    n <- 50
    vec1 <- rnorm(n); vec2 <- rnorm(n)
    res <- .representation_2d(vec1, vec2)
    
    # ensure basis matrix is a valid basis matrix
    bool1 <- sum(abs(crossprod(res$basis_mat) - diag(2))) <= 1e-6
    
    # ensure coefficients are correct
    bool2 <- (sum(abs(vec1 - res$basis_mat %*% res$rep1)) <= 1e-6 &
                sum(abs(vec2 - res$basis_mat %*% res$rep2)) <= 1e-6)
    
    # ensure all linear combinations can be represented in that basis
    coef_vec <- runif(2, min = -5, max = 5)
    vec <- res$basis_mat %*% coef_vec
    bool3 <- sum(abs(vec - tcrossprod(res$basis_mat)%*%vec)) <= 1e-6
    
    
    bool1 & bool2 & bool3
  })
  
  expect_true(all(bool_vec))
})
