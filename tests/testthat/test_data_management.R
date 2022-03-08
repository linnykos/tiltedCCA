context("Test svd")

## .get_SVD is correct

test_that(".get_SVD works for matrix", {
  mat <- matrix(1:50, 10, 5)
  rownames(mat) <- paste0("c", 1:10)
  colnames(mat) <- paste0("g", 1:5)
  res <- .get_SVD(mat, center = T, dims = 1:2, scale = T)
  
  expect_true(inherits(res, "svd"))
  expect_true(all(sort(names(res)) == sort(c("u", "d", "v"))))
  expect_true(all(rownames(res$u) == rownames(mat)) & length(rownames(res$u)) == 10)
  expect_true(all(rownames(res$v) == colnames(mat)) & length(rownames(res$v)) == 5)
})


test_that(".get_SVD works for dgCMatrix", {
  set.seed(10)
  mat <- matrix(1:50, 10, 5)
  mat[sample(1:prod(dim(mat)), round(0.8*prod(dim(mat))))] <- 0
  mat <- Matrix::Matrix(mat, sparse = T)
  rownames(mat) <- paste0("c", 1:10)
  colnames(mat) <- paste0("g", 1:5)
  res <- .get_SVD(mat, center = T, dims = 1:2, scale = T)
  
  expect_true(inherits(res, "svd"))
  expect_true(all(sort(names(res)) == sort(c("u", "d", "v"))))
  expect_true(all(rownames(res$u) == rownames(mat)) & length(rownames(res$u)) == 10)
  expect_true(all(rownames(res$v) == colnames(mat)) & length(rownames(res$v)) == 5)
})


test_that(".get_SVD works for multiSVD", {
  load(paste0("../assets/test_data", i, ".RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2)
  
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 1)
  res1 <- .get_SVD(multiSVD_obj)
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 2)
  res2 <- .get_SVD(multiSVD_obj)
  
  expect_true(inherits(res1, "svd"))
  expect_true(all(c("u", "d", "v") %in% names(res1)))
  expect_true(all(rownames(res1$u) == rownames(mat_1)) & length(rownames(res1$u)) == nrow(mat_1))
  expect_true(all(rownames(res1$v) == colnames(mat_1)) & length(rownames(res1$v)) == ncol(mat_1))
  
  expect_true(inherits(res2, "svd"))
  expect_true(all(c("u", "d", "v") %in% names(res2)))
  expect_true(all(rownames(res2$u) == rownames(mat_2)) & length(rownames(res2$u)) == nrow(mat_2))
  expect_true(all(rownames(res2$v) == colnames(mat_2)) & length(rownames(res2$v)) == ncol(mat_2))
})

#########################

## .get_Dimred is correct

test_that(".get_Dimred works for matrix", {
  mat <- matrix(1:50, 10, 5)
  rownames(mat) <- paste0("c", 1:10)
  colnames(mat) <- paste0("g", 1:5)
  res <- .get_Dimred(mat, center = T, dims = 1:2, scale = T,
                     normalize_singular_value = T)
  
  expect_true(inherits(res, "matrix"))
  expect_true(all(rownames(res) == rownames(mat)) & length(rownames(res)) == 10)
})

test_that(".get_Dimred works for svd", {
  mat <- matrix(1:50, 10, 5)
  rownames(mat) <- paste0("c", 1:10)
  colnames(mat) <- paste0("g", 1:5)
  svd_obj <- .get_SVD(mat, center = T, dims = 1:2, scale = T)
  
  res <- .get_Dimred(svd_obj, center = T, dims = 1:2, scale = T,
                     normalize_singular_value = T)
  
  expect_true(inherits(res, "matrix"))
  expect_true(all(rownames(res) == rownames(mat)) & length(rownames(res)) == 10)
})

############################

## .get_postDimred is correct

test_that(".get_postDimred works", {
  load(paste0("../assets/test_data", i, ".RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2)
  
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 1)
  res <- .get_postDimred(multiSVD_obj, averaging_mat = NULL)
  
  expect_true(inherits(res, "matrix"))
  expect_true(all(rownames(res) == rownames(mat_1)) & length(rownames(res)) == nrow(mat_1))
})


test_that(".get_postDimred can differentiate between centering and not", {
  load(paste0("../assets/test_data", i, ".RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2; mat_2 <- mat_2 + 2
  multiSVD_obj1 <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = T, center_2 = F,
                                  normalize_row = F,
                                  normalize_singular_value = T,
                                  recenter_1 = F, recenter_2 = T,
                                  rescale_1 = F, rescale_2 = T,
                                  scale_1 = T, scale_2 = F)
  multiSVD_obj2 <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = T, center_2 = F,
                                   normalize_row = F,
                                   normalize_singular_value = T,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = T, scale_2 = F)
  
  multiSVD_obj1 <- .set_defaultAssay(multiSVD_obj1, assay = 2)
  res1 <- .get_postDimred(multiSVD_obj1, averaging_mat = NULL)
  multiSVD_obj2 <- .set_defaultAssay(multiSVD_obj2, assay = 2)
  res2 <- .get_postDimred(multiSVD_obj2, averaging_mat = NULL)
  
  colmeans_1 <- colMeans(res1)
  colmeans_2 <- colMeans(res2)
  expect_true(sum(abs(colmeans_1)) <= 1e-6)
  expect_true(sum(abs(colmeans_2)) >= 1e-6)
  
  colsd_1 <- matrixStats::colSds(res1)
  colsd_2 <- matrixStats::colSds(res2)
  expect_true(sum(abs(colsd_1 - 1)) <= 1e-6)
  expect_true(sum(abs(colsd_2 - 1)) >= 1e-6)
})


test_that(".get_postDimred respects normalize_row", {
  load(paste0("../assets/test_data", i, ".RData"))
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2; mat_2 <- mat_2 + 2
  multiSVD_obj1 <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = T, center_2 = T,
                                   normalize_row = T,
                                   normalize_singular_value = T,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = T, scale_2 = T)
  multiSVD_obj2 <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = T, center_2 = T,
                                   normalize_row = F,
                                   normalize_singular_value = T,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = T, scale_2 = T)
  
  multiSVD_obj1 <- .set_defaultAssay(multiSVD_obj1, assay = 2)
  res1 <- .get_postDimred(multiSVD_obj1, averaging_mat = NULL)
  multiSVD_obj2 <- .set_defaultAssay(multiSVD_obj2, assay = 2)
  res2 <- .get_postDimred(multiSVD_obj2, averaging_mat = NULL)
  
  expect_true(sum(abs(apply(res1, 1, .l2norm) - 1)) <= 1e-6)
  expect_true(sum(abs(apply(res2, 1, .l2norm) - 1)) >= 1e-6)
})

