context("Test generate data")

test_that("generate_data works", {
  set.seed(10)
  n <- 100; rank_12 <- 2; rank_1 <- 2; rank_2 <- 2; p_1 <- 20; p_2 <- 40
  common_loading <- matrix(0, n, rank_12)
  diag(common_loading) <- 1
  distinct_loading_1 <- matrix(0, n, rank_1)
  for(i in 1:2){
    distinct_loading_1[i+2,i] <- 1
  }
  distinct_loading_2 <- matrix(0, n, rank_2)
  for(i in 1:2){
    distinct_loading_2[i+4,i] <- 1
  }
  coef_mat_1 <- matrix(runif(rank_1*p_1), rank_1, p_1)
  coef_mat_2 <- matrix(runif(rank_2*p_2), rank_2, p_2)
  noise_func <- function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}
  
  res <- generate_data(common_loading, distinct_loading_1, distinct_loading_2,
                       coef_mat_1, coef_mat_2, noise_func = noise_func)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2"))))
  expect_true(all(sapply(res, is.matrix)))
  expect_true(all(dim(res$mat_1) == c(n, p_1)))
  expect_true(all(dim(res$mat_2) == c(n, p_2)))
})
