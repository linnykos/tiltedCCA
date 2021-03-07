context("Test spoet")

## .spoet is correct

test_that("(Basic) .spoet works", {
  set.seed(10)
  p <- 100; n <- 20; K <- 2
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  mat <- abs(MASS::mvrnorm(n = n, mu = c(200,rep(30, p/2-1), rep(5, p/2)), Sigma = cov_mat))
  
  res <- .spoet(mat, K = 2)
  
  expect_true(all(dim(mat) == c(nrow(res$u), nrow(res$v))))
  expect_true(length(res$d) == 2)
  expect_true(all(sort(names(res)) == sort(c("d", "d_original", "u", "v"))))
})

test_that("(Coding) .spoet preserves rownames and columnnames", {
  set.seed(10)
  p <- 100; n <- 20; K <- 2
  cov_mat <- matrix(0, nrow = p, ncol = p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  mat <- abs(MASS::mvrnorm(n = n, mu = c(200,rep(30, p/2-1), rep(5, p/2)), Sigma = cov_mat))
  rownames(mat) <- paste0("a", 1:n)
  colnames(mat) <- paste0("b", 1:p)
  
  res <- .spoet(mat, K = 2)
  
  expect_true(all(rownames(res$u) == rownames(mat)))
  expect_true(all(rownames(res$v) == colnames(mat)))
})
