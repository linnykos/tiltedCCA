context("Test linear algabra")

## .cor_vectors is correct

test_that("(Code) .cor_vectors works with len is 0", {
  res <- .cor_vectors(c(1,0),c(0,0))
  expect_true(is.na(res))
})

######################

## .angle_from_vector is correct

test_that(".angle_from_vector works", {
  res <- .angle_from_vector(rep(1/sqrt(2), 2), 45)
  expect_true(sum(abs(res - c(0,1))) <= 1e-6)
  
  res <- .angle_from_vector(c(0,1), 90)
  expect_true(sum(abs(res - c(-1,0))) <= 1e-6)
  
  res <- .angle_from_vector(c(-1,0), 90)
  expect_true(sum(abs(res - c(0,-1))) <= 1e-6)
})

###########################3#

## .rightmost_vector is correct

test_that(".rightmost_vector works", {
  set.seed(10)
  vec1 <- abs(rnorm(2)); vec2 <- abs(rnorm(2))
  res <- .rightmost_vector(vec1, vec2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("vec_right", "vec_left", "len_right", "len_left"))))
})
