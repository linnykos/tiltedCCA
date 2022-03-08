context("Test metacell")

## .intersect_clustering is correct

test_that(".intersect_clustering works", {
  set.seed(10)
  vec <- rep(1:5, each = 50)
  n <- length(vec)
  large_clustering_1 <- factor(vec)
  large_clustering_2 <- factor(sample(vec))
  res <- .intersect_clustering(large_clustering_1 = large_clustering_1,
                               large_clustering_2 = large_clustering_2)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("large_clustering_list", 
                                             "clustering_hierarchy_1",
                                             "clustering_hierarchy_2"))))
  vec2 <- .convert_list2factor(res$large_clustering_list, n = n)
  level_vec <- levels(vec2)
  bool_vec <- sapply(1:length(level_vec), function(i){
    idx <- which(vec2 == level_vec[i])
    uniq_vec <- unique(vec[idx])
    length(uniq_vec) == 1
  })
  expect_true(all(bool_vec))
  expect_true(max(table(unlist(res$large_clustering_list))) == 1)
  expect_true(max(table(unlist(res$clustering_hierarchy_1))) == 1)
  expect_true(max(table(unlist(res$clustering_hierarchy_2))) == 1)
})

###############################

## .compute_metacell_splits is correct

test_that(".compute_metacell_splits works", {
  load(paste0("../assets/test_data3.RData"))
  svd_1 <- test_data$svd_1; svd_2 <- test_data$svd_2
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  n <- nrow(dimred_1)
  set.seed(10)
  large_clustering_1 <- factor(stats::kmeans(dimred_1, centers = 2)$cluster)
  large_clustering_2 <- factor(stats::kmeans(dimred_2, centers = 2)$cluster)
  tmp <- .intersect_clustering(large_clustering_1 = large_clustering_1,
                               large_clustering_2 = large_clustering_2)
  res <- .compute_metacell_splits(tmp$large_clustering_list,
                                  n = n,
                                  num_metacells = 100)
  
  expect_true(is.data.frame(res))
  expect_true(all(sort(colnames(res)) == sort(c("total_size", "size", "num"))))
  expect_true(sum(res$total_size) == n)
  expect_true(sum(res$num) == 100)
})

###############################

## .compute_metacells is correct

test_that(".compute_metacells works", {
  load(paste0("../assets/test_data3.RData"))
  svd_1 <- test_data$svd_1; svd_2 <- test_data$svd_2
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  n <- nrow(dimred_1)
  set.seed(10)
  large_clustering_1 <- factor(stats::kmeans(dimred_1, centers = 2)$cluster)
  large_clustering_2 <- factor(stats::kmeans(dimred_2, centers = 2)$cluster)
  tmp <- .intersect_clustering(large_clustering_1 = large_clustering_1,
                               large_clustering_2 = large_clustering_2)
  idx <- tmp$large_clustering_list[[1]]
  
  res <- .compute_metacells(dimred_1 = dimred_1[idx,],
                            dimred_2 = dimred_2[idx,],
                            k = 10,
                            row_indices = idx)
  
  expect_true(all(sort(unlist(res)) == sort(idx)))
  expect_true(max(table(unlist(res))) == 1)
})

###############################

## .form_metacells is correct

test_that(".form_metacells works", {
  load(paste0("../assets/test_data3.RData"))
  svd_1 <- test_data$svd_1; svd_2 <- test_data$svd_2
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  set.seed(10)
  large_clustering_1 <- factor(stats::kmeans(dimred_1, centers = 2)$cluster)
  large_clustering_2 <- factor(stats::kmeans(dimred_2, centers = 2)$cluster)
  tmp <- .intersect_clustering(large_clustering_1 = large_clustering_1,
                               large_clustering_2 = large_clustering_2)
  res <- .form_metacells(dimred_1 = dimred_1, 
                         dimred_2 = dimred_2,
                         large_clustering_list = tmp$large_clustering_list,
                         num_metacells = 100)
  
  expect_true(length(res) == 100)
  expect_true(all(sort(unlist(res)) == 1:nrow(dimred_1)))
})

##############################

## .generate_averaging_matrix is correct

test_that(".generate_averaging_matrix works", {
  set.seed(10)
  vec <- factor(rep(c("a","b","c","d","e"), each = 5))
  vec[sample(length(vec), 10)] <- NA
  n <- length(vec)
  lis <- .convert_factor2list(vec)
  res <- .generate_averaging_matrix(metacell_clustering_list = lis,
                                    n = n)
  rowsum_vec <- sparseMatrixStats::rowSums2(res)
  idx_list <- sapply(1:ncol(res), function(j){
    .nonzero_col(mat = res, col_idx = j, bool_value = F)
  })
  
  expect_true(sum(abs(rowsum_vec - 1)) <= 1e-6)
  expect_true(all(sapply(idx_list, length) <= 1))
})

#############################

## form_metacells is correct

test_that("form_metacells works", {
  load(paste0("../assets/test_data3.RData"))
  mat_1 <- test_data$mat_1; mat_2 <- test_data$mat_2
  svd_1 <- test_data$svd_1; svd_2 <- test_data$svd_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2)
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  n <- nrow(dimred_1)
  set.seed(10)
  large_clustering_1 <- factor(stats::kmeans(dimred_1, centers = 2)$cluster)
  large_clustering_2 <- factor(stats::kmeans(dimred_2, centers = 2)$cluster)
  
  res <- form_metacells(input_obj = multiSVD_obj,
                        large_clustering_1 = large_clustering_1, 
                        large_clustering_2 = large_clustering_2,
                        num_metacells = 100)
  
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all("metacell_obj" %in% names(res)))
  expect_true(inherits(res$metacell_obj, "metacell"))
  expect_true(all(sort(names(res$metacell_obj)) == sort(c("large_clustering_1",
                                                          "large_clustering_2",
                                                          "metacell_clustering_list"))))
  expect_true(all(res$metacell_obj$large_clustering_1 == large_clustering_1))
  expect_true(all(res$metacell_obj$large_clustering_2 == large_clustering_2))
  expect_true(max(table(unlist(res$metacell_obj$metacell_clustering_list))) == 1)
  expect_true(all(sort(unlist(res$metacell_obj$metacell_clustering_list)) == 1:n))
})