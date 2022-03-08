context("Test svd")

## .svd_truncated is correct

test_that(".svd_truncated works", {
  set.seed(10)
  mat <- MASS::mvrnorm(n = 100, mu = rep(0, 10), Sigma = diag(10))
  rownames(mat) <- paste0("c", 1:100)
  colnames(mat) <- paste0("g", 1:10)
  res <- .svd_truncated(mat = mat, K = 5, symmetric = F, 
                        rescale = F, mean_vec = F, sd_vec = F,
                        K_full_rank = F)
  
  expect_true(all(rownames(res$u) == rownames(mat)) & length(rownames(res$u)) == 100)
  expect_true(all(rownames(res$v) == colnames(mat)) & length(rownames(res$v)) == 10)
})

########################

## create_multiSVD is correct

test_that("create_multiSVD works", {
  load(paste0("../assets/test_data", i, ".RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  res <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                         dims_1 = 1:2, dims_2 = 1:2)
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(sort(names(res)) == sort(c("svd_1", "svd_2", "default_assay", "param")))
})
