context("Test tiltedCCA_decomposition")

## tiltedCCA_decomposition is correct

test_that("(Basic) tiltedCCA_decomposition works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  res <- tiltedCCA_decomposition(multiSVD_obj)
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("common_mat_1", "common_mat_2", "distinct_mat_1", "distinct_mat_2") %in% names(res)))
  expect_true(all(dim(res$common_mat_1) == dim(mat_1)))
  expect_true(all(rownames(res$common_mat_1) == rownames(mat_1)) & length(rownames(res$common_mat_1)) > 0)
  expect_true(all(rownames(res$common_mat_2) == rownames(mat_2)) & length(rownames(res$common_mat_2)) > 0)
  expect_true(all(colnames(res$common_mat_1) == colnames(mat_1)) & length(colnames(res$common_mat_1)) > 0)
  expect_true(all(colnames(res$common_mat_2) == colnames(mat_2)) & length(colnames(res$common_mat_2)) > 0)
  
  expect_true(all(rownames(res$distinct_mat_1) == rownames(mat_1)) & length(rownames(res$distinct_mat_1)) > 0)
  expect_true(all(rownames(res$distinct_mat_2) == rownames(mat_2)) & length(rownames(res$distinct_mat_2)) > 0)
  expect_true(all(colnames(res$distinct_mat_1) == colnames(mat_1)) & length(colnames(res$distinct_mat_1)) > 0)
  expect_true(all(colnames(res$distinct_mat_2) == colnames(mat_2)) & length(colnames(res$distinct_mat_2)) > 0)
})

test_that("(Math) tiltedCCA_decomposition yields uncorrelated distinct matrices", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    n <- nrow(mat_1)
    large_clustering_1 <- test_data$clustering_1
    large_clustering_2 <- test_data$clustering_2
    multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                    dims_1 = 1:2, dims_2 = 1:2,
                                    center_1 = F, center_2 = F,
                                    normalize_row = T,
                                    normalize_singular_value = F,
                                    recenter_1 = F, recenter_2 = F,
                                    rescale_1 = F, rescale_2 = F,
                                    scale_1 = F, scale_2 = F)
    multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                   large_clustering_1 = large_clustering_1, 
                                   large_clustering_2 = large_clustering_2,
                                   num_metacells = NULL)
    multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                                 latent_k = 2,
                                 num_neigh = 10,
                                 bool_cosine = T,
                                 bool_intersect = T,
                                 min_deg = 1)
    multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                              verbose = F)
    res <- tiltedCCA_decomposition(multiSVD_obj)
    
    res <- .set_defaultAssay(res, assay = 1)
    distinct_mat_1 <- .get_tCCAobj(res, apply_postDimred = F, what = "distinct_mat")
    res <- .set_defaultAssay(res, assay = 2)
    distinct_mat_2 <- .get_tCCAobj(res, apply_postDimred = F, what = "distinct_mat")
    
    tmp <- crossprod(distinct_mat_1, distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition yields a low-rank matrix", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    n <- nrow(mat_1)
    large_clustering_1 <- test_data$clustering_1
    large_clustering_2 <- test_data$clustering_2
    multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                    dims_1 = 1:2, dims_2 = 1:2,
                                    center_1 = F, center_2 = F,
                                    normalize_row = T,
                                    normalize_singular_value = F,
                                    recenter_1 = F, recenter_2 = F,
                                    rescale_1 = F, rescale_2 = F,
                                    scale_1 = F, scale_2 = F)
    multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                   large_clustering_1 = large_clustering_1, 
                                   large_clustering_2 = large_clustering_2,
                                   num_metacells = NULL)
    multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                                 latent_k = 2,
                                 num_neigh = 10,
                                 bool_cosine = T,
                                 bool_intersect = T,
                                 min_deg = 1)
    multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                              verbose = F)
    res <- tiltedCCA_decomposition(multiSVD_obj)
    
    K <- 2
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition yields common matrices with the same column space", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    n <- nrow(mat_1)
    large_clustering_1 <- test_data$clustering_1
    large_clustering_2 <- test_data$clustering_2
    multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                    dims_1 = 1:2, dims_2 = 1:2,
                                    center_1 = F, center_2 = F,
                                    normalize_row = T,
                                    normalize_singular_value = F,
                                    recenter_1 = F, recenter_2 = F,
                                    rescale_1 = F, rescale_2 = F,
                                    scale_1 = F, scale_2 = F)
    multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                   large_clustering_1 = large_clustering_1, 
                                   large_clustering_2 = large_clustering_2,
                                   num_metacells = NULL)
    multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                                 latent_k = 2,
                                 num_neigh = 10,
                                 bool_cosine = T,
                                 bool_intersect = T,
                                 min_deg = 1)
    multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                              verbose = F)
    res <- tiltedCCA_decomposition(multiSVD_obj)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition can obtain the same result when fed into itself (i.e., stability/identifiability)", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  res <- tiltedCCA_decomposition(multiSVD_obj)
  
  multiSVD_obj2 <- create_multiSVD(mat_1 = res$common_mat_1 + res$distinct_mat_1, 
                                   mat_2 = res$common_mat_2 + res$distinct_mat_2, 
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = F, center_2 = F,
                                   normalize_row = T,
                                   normalize_singular_value = F,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = F, scale_2 = F)
  multiSVD_obj2 <- form_metacells(input_obj = multiSVD_obj2,
                                  large_clustering_1 = large_clustering_1, 
                                  large_clustering_2 = large_clustering_2,
                                  num_metacells = NULL)
  multiSVD_obj2 <- compute_snns(input_obj = multiSVD_obj2,
                                latent_k = 2,
                                num_neigh = 10,
                                bool_cosine = T,
                                bool_intersect = T,
                                min_deg = 1)
  multiSVD_obj2 <- tiltedCCA(input_obj = multiSVD_obj2,
                             verbose = F)
  res2 <- tiltedCCA_decomposition(multiSVD_obj2)
  
  expect_true(sum(abs(res$common_mat_1 - res2$common_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$common_mat_2 - res2$common_mat_2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_1 - res2$distinct_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_2 - res2$distinct_mat_2)) <= 1e-6)
})
