context("Test meta-cell formation")

## form_metacells is correct
test_that("form_metacells works", {
  set.seed(10)
  n <- 1000
  p <- 10
  mat <- matrix(runif(n*p), ncol = p, nrow = n)
  res <- form_metacells(mat, verbose = F)
  
  expect_true(is.factor(res))
  expect_true(length(res) == n)
})
