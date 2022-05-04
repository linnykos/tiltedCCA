context("Test consensus PCA")

## consensus_pca is correct

test_that("consensus_pca works", {
  # load("tests/assets/test_data4.RData")
  load("../assets/test_data4.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  
  res <- consensus_pca(mat_1 = mat_1, mat_2 = mat_2,
                       dims_1 = 1:2, dims_2 = 1:2,
                       dims_consensus = 1:2)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "consensusPCA"))
  expect_true(all(sort(names(res)) == sort(c("dimred_1", "dimred_2", "dimred_consensus", "param"))))
  expect_true(all(dim(res$dimred_1) == c(nrow(mat_1),2)))
  expect_true(all(dim(res$dimred_2) == c(nrow(mat_1),2)))
  expect_true(all(dim(res$dimred_consensus) == c(nrow(mat_1),2)))
  expect_true(length(rownames(res$dimred_consensus)) > 0)
  expect_true(all(rownames(res$dimred_consensus) == rownames(mat_1)))
  
  # plot(res$dimred_1[,1], res$dimred_1[,2], col = rep(1:4, each = 100), pch = 16, asp = T)
  # plot(res$dimred_2[,1], res$dimred_2[,2], col = rep(1:4, each = 100), pch = 16, asp = T)
  # plot(res$consensus_dimred[,1], res$consensus_dimred[,2], col = rep(1:4, each = 100), pch = 16, asp = T)
})


test_that("consensus_pca can take SVDs as inputs", {
  # load("tests/assets/test_data4.RData")
  load("../assets/test_data4.RData")
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  
  res <- consensus_pca(mat_1 = NULL, mat_2 = NULL,
                       dims_1 = NULL, dims_2 = NULL,
                       dims_consensus = 1:2,
                       svd_1 = svd_1, svd_2 = svd_2)
  
  expect_true(inherits(res, "consensusPCA"))
  expect_true(all(sort(names(res)) == sort(c("dimred_1", "dimred_2", "dimred_consensus", "param"))))
})