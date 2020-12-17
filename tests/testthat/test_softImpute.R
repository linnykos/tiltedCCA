context("Test softImpute")

## softImpute_diagnostic is correct

test_that("softImpute_diagnostic works", {
  set.seed(10)
  dat <- MASS::mvrnorm(n = 100, mu = rep(0, 20), Sigma = diag(20))
  svd_res <- svd(dat)
  dat <- tcrossprod(.mult_mat_vec(svd_res$u[,1:3], svd_res$d[1:3]), svd_res$v[,1:3])
  
  expect_true(Matrix::rankMatrix(dat) == 3)
  
  res <- softImpute_diagnostic(dat, K = 3)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("training_mat", "testing_mat"))))
  expect_true(all(sapply(res, ncol) == 2))
  expect_true(sum(sapply(res, nrow)) == prod(dim(dat)))
  expect_true(all(sort(colnames(res$training_mat)) == sort(c("observed_val", "predicted_val"))))
  expect_true(all(sort(colnames(res$testing_mat)) == sort(c("observed_val", "predicted_val"))))
})

##################

## plot_prediction_against_observed is correct

test_that("plot_prediction_against_observed does not crash", {
  set.seed(10)
  dat <- MASS::mvrnorm(n = 100, mu = rep(0, 20), Sigma = diag(20))
  svd_res <- svd(dat)
  dat <- tcrossprod(.mult_mat_vec(svd_res$u[,1:3], svd_res$d[1:3]), svd_res$v[,1:3])
  
  expect_true(Matrix::rankMatrix(dat) == 3)
  
  diag_res <- softImpute_diagnostic(dat, K = 3)
  
  res <- plot_prediction_against_observed(diag_res$training_mat)
  
  expect_true(TRUE)
})