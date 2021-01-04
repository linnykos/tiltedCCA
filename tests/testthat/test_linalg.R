context("Test linear algebra")

## .projection is correct

test_that(".projection works", {
  set.seed(10)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection preserves orthogonal vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  res2 <- .projection(res, vec2)
  
  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection removes parallel vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- 2*vec1
  res <- .projection(vec1, vec2)
  
  expect_true(sum(abs(res)) < 1e-6)
})

test_that(".projection reduces the norm", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  
  expect_true(.l2norm(vec1) >= .l2norm(res))
})


################################

## .orthogonal_vector is correct

test_that(".orthogonal_vector works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .orthogonal_vector(vec, mat)
  
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".orthogonal_vector returns a vector that is orthogonal to rows of mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .orthogonal_vector(vec, mat)
  
  for(i in 1:nrow(mat)){
    expect_true(abs(res %*% mat[i,]) < 1e-6)
  }
})

test_that(".orthogonal_vector reduces the norm", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .orthogonal_vector(vec, mat)
  
  expect_true(.l2norm(vec) >= .l2norm(res))
})

test_that(".orthogonal_vector preserves orthogonal vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .orthogonal_vector(vec, mat)
  res2 <- .orthogonal_vector(res, mat)
  
  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".orthogonal_vector removes parallel vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- 2*mat[1,]
  res <- .orthogonal_vector(vec, mat)
  
  expect_true(sum(abs(res)) < 1e-6)
})
