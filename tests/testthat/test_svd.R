context("Test svd")

## .svd_safe is correct

test_that(".svd_safe works", {
  set.seed(10)
  mat <- MASS::mvrnorm(n = 100, mu = rep(0, 10), Sigma = diag(10))
  rownames(mat) <- paste0("c", 1:100)
  colnames(mat) <- paste0("g", 1:10)
  res <- .svd_safe(mat = mat,
                   check_stability = T, 
                   K = 5, 
                   mean_vec = F, 
                   rescale = F, 
                   scale_max = NULL, 
                   sd_vec = NULL)
  
  expect_true(all(rownames(res$u) == rownames(mat)) & length(rownames(res$u)) == 100)
  expect_true(all(rownames(res$v) == colnames(mat)) & length(rownames(res$v)) == 10)
})

test_that(".svd_safe works with tall and skinny matrices", {
  set.seed(10)
  mat <- MASS::mvrnorm(n = 100, mu = rep(0, 3), Sigma = diag(3))
  rownames(mat) <- paste0("c", 1:100)
  colnames(mat) <- paste0("g", 1:3)
  res <- .svd_safe(mat = mat,
                   check_stability = T, 
                   K = 3, 
                   mean_vec = F, 
                   rescale = F, 
                   scale_max = NULL, 
                   sd_vec = NULL)
  
  mat2 <- tcrossprod(.mult_mat_vec(res$u, res$d), res$v)
  expect_true(sum(abs(mat - mat2)) <= 1e-5)
})

test_that(".svd_safe works with scale_max", {
  set.seed(10)
  mat <- MASS::mvrnorm(n = 100, mu = rep(0, 10), Sigma = diag(10))*10
  rownames(mat) <- paste0("c", 1:100)
  colnames(mat) <- paste0("g", 1:10)
  res1 <- .svd_safe(mat = mat,
                   check_stability = T, 
                   K = 5, 
                   mean_vec = F, 
                   rescale = F, 
                   scale_max = NULL, 
                   sd_vec = NULL)
  res2 <- .svd_safe(mat = mat,
                    check_stability = T, 
                    K = 5, 
                    mean_vec = F, 
                    rescale = F, 
                    scale_max = 10, 
                    sd_vec = NULL)
  res3 <- .svd_safe(mat = mat,
                    check_stability = T, 
                    K = 5, 
                    mean_vec = F, 
                    rescale = F, 
                    scale_max = max(mat)*2, 
                    sd_vec = NULL)
  
  
  expect_true(sum(abs(res1$d - res2$d)) > 0)
  expect_true(sum(abs(res1$d - res3$d)) < 1e-6)
})

########################

## create_multiSVD is correct

test_that("create_multiSVD works", {
  load(paste0("../assets/test_data1.RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  res <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                         dims_1 = 1:2, dims_2 = 1:2)
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(sort(names(res)) == sort(c("svd_1", "svd_2", "default_assay", "param"))))
})

test_that("create_multiSVD respects scale and center", {
  load(paste0("../assets/test_data1.RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2 + 2
  res <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                         dims_1 = 1:2, dims_2 = 1:2,
                         center_1 = T, center_2 = F,
                         scale_1 = T, scale_2 = F)
  
  res <- .set_defaultAssay(res, assay = 1)
  svd_1 <- .get_SVD(res)
  res <- .set_defaultAssay(res, assay = 2)
  svd_2 <- .get_SVD(res)
  
  expect_true(sum(abs(colMeans(svd_1$u))) <= 1e-6)
  expect_true(sum(abs(colMeans(svd_2$u))) >= 1e-6)
  
  mat_1b <- test_data$mat_1*100
  res2 <- create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2,
                         dims_1 = 1:2, dims_2 = 1:2,
                         center_1 = T, center_2 = F,
                         scale_1 = F, scale_2 = F)
  res2 <- .set_defaultAssay(res2, assay = 1)
  svd_1b <- .get_SVD(res2)
  
  expect_true(sum(svd_1b$d) > sum(svd_1$d))
})
