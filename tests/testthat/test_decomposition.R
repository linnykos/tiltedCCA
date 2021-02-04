context("Test decomposition")

## .decomposition_2d is correct

test_that("(Math) .decomposition_2d works", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(13*x)
    vec1 <- runif(2, min = 0, max = 2)
    vec2 <- runif(2, min = 0, max = 2)
    
    res <- .decomposition_2d(vec1, vec2, gridsize = 50, plotting = F)
    
    bool1 <- sum(abs(vec1 - res$common_vec1 - res$distinct_vec1)) <= 1e-6
    bool2 <- sum(abs(vec2 - res$common_vec2 - res$distinct_vec2)) <= 1e-6
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .decomposition_2d yields orthogonal distinct vectors", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(13*x)
    vec1 <- runif(2, min = 0, max = 2)
    vec2 <- runif(2, min = 0, max = 2)
    
    res <- .decomposition_2d(vec1, vec2, gridsize = 50, plotting = F)
    abs(.angle_between_vectors(res$distinct_vec1, res$distinct_vec2)-90) <= 1
  })
  
  expect_true(all(bool_vec))
})
