context("Test SNN")

## .form_snn_mat is correct
test_that(".form_snn_mat works", {
  load(paste0("../assets/test_data1.RData"))
  svd_1 <- test_data$svd_1
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  res <- .form_snn_mat(dimred_1, 
                       num_neigh = 10,
                       bool_cosine = T, 
                       bool_intersect = T, 
                       min_deg = 5)
  
  expect_true(inherits(res, "dgCMatrix"))
  expect_true(sum(abs(res - Matrix::t(res))) <= 1e-6)
  expect_true(all(rownames(res) == rownames(dimred_1)) && length(rownames(res)) == nrow(res))
  expect_true(all(colnames(res) == rownames(dimred_1)) && length(colnames(res)) == nrow(res))
  
  expect_true(all(Matrix::rowSums(res) >= 5))
})

#############################

## .compute_laplacian_basis is correct

test_that(".compute_laplacian_basis works", {
  load(paste0("../assets/test_data1.RData"))
  svd_1 <- test_data$svd_1
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  snn <- .form_snn_mat(dimred_1, 
                       num_neigh = 10,
                       bool_cosine = T, 
                       bool_intersect = T, 
                       min_deg = 5)
  res <- .compute_laplacian_basis(latent_k = 2,
                                  sparse_mat = snn)
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(300, 2)))
})

##########################

## compute_snns is correct

test_that("compute_snns works", {
  load(paste0("../assets/test_data3.RData"))
  mat_1 <- test_data$mat_1; mat_2 <- test_data$mat_2
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  n <- nrow(mat_1)
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  
  res <- compute_snns(input_obj = multiSVD_obj,
                      latent_k = 2,
                      num_neigh = 10,
                      bool_cosine = T,
                      bool_intersect = T,
                      min_deg = 5)
  
  expect_true(length(grep("^snn*", names(.get_param(res)))) == 5)
  expect_true(inherits(res, "snn"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("snn_list", "laplacian_list") %in% names(res)))
  expect_true(all(sort(names(res$snn_list)) == sort(c("snn_1", "snn_2", "common_snn"))))
  expect_true(all(sort(names(res$laplacian_list)) == sort(c("laplacian_1", "laplacian_2", "common_laplacian"))))
  
  for(i in 1:3){
    snn <- res$snn_list[[i]]
    expect_true(inherits(snn, "dgCMatrix"))
    expect_true(sum(abs(snn - Matrix::t(snn))) <= 1e-6)
    expect_true(all(dim(snn) == n))
    expect_true(all(rownames(snn) == rownames(mat_1)) && length(rownames(snn)) == nrow(snn))
    expect_true(all(colnames(snn) == rownames(mat_1)) && length(colnames(snn)) == nrow(snn))
    
    mat <- res$laplacian_list[[i]]
    expect_true(is.matrix(mat))
    expect_true(all(dim(mat) == c(n, 2)))
  }
})


test_that("compute_snns works with metacells", {
  load(paste0("../assets/test_data3.RData"))
  mat_1 <- test_data$mat_1; mat_2 <- test_data$mat_2
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  n <- nrow(mat_1)
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = 100)
  metacell_names <- .get_metacell(input_obj = multiSVD_obj,
                                  resolution = "metacell", 
                                  type = "factor", 
                                  what = "metacell_clustering")
  
  res <- compute_snns(input_obj = multiSVD_obj,
                      latent_k = 2,
                      num_neigh = 10,
                      bool_cosine = T,
                      bool_intersect = T,
                      min_deg = 5)
  
  expect_true(length(grep("^snn*", names(.get_param(res)))) == 5)
  expect_true(inherits(res, "snn"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("snn_list", "laplacian_list") %in% names(res)))
  expect_true(all(sort(names(res$snn_list)) == sort(c("snn_1", "snn_2", "common_snn"))))
  expect_true(all(sort(names(res$laplacian_list)) == sort(c("laplacian_1", "laplacian_2", "common_laplacian"))))
  
  for(i in 1:3){
    snn <- res$snn_list[[i]]
    expect_true(inherits(snn, "dgCMatrix"))
    expect_true(sum(abs(snn - Matrix::t(snn))) <= 1e-6)
    expect_true(all(dim(snn) == 100))
    expect_true(all(rownames(snn) == as.character(metacell_names)) && length(rownames(snn)) == nrow(snn))
    expect_true(all(colnames(snn) == as.character(metacell_names)) && length(colnames(snn)) == nrow(snn))
    
    mat <- res$laplacian_list[[i]]
    expect_true(is.matrix(mat))
    expect_true(all(dim(mat) == c(100, 2)))
  }
})