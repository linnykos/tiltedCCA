context("Test common SNN")

## .compute_common_snn_hardclustering is correct

test_that(".compute_common_snn_hardclustering works", {
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
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 1)
  dimred_1 <- .get_postDimred(multiSVD_obj, averaging_mat = NULL)
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 2)
  dimred_2 <- .get_postDimred(multiSVD_obj, averaging_mat = NULL)
  
  num_neigh <- 10; bool_cosine <- T; bool_intersect <- T; min_deg <- 5
  snn_1 <- .form_snn_mat(mat = dimred_1,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg)
  snn_2 <- .form_snn_mat(mat = dimred_2,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg)
  
  res <- .compute_common_snn_hardclustering(snn_1 = snn_1,
                                            snn_2 = snn_2,
                                            clustering_1 = large_clustering_1,
                                            clustering_2 = large_clustering_2,
                                            num_neigh = num_neigh)
  
  expect_true(inherits(res, "dgCMatrix"))
  expect_true(sum(abs(res - Matrix::t(res))) <= 1e-6)
  expect_true(all(rownames(res) == rownames(dimred_1)) && length(rownames(res)) == nrow(res))
  expect_true(all(colnames(res) == rownames(dimred_1)) && length(colnames(res)) == nrow(res))
})

############################

## .l2_selection_qp is correct

test_that(".l2_selection_qp works", {
  num_neigh <- 30
  obs_tab <- as.table(matrix(c(10,10,5, 5,5,10), nrow = 3, ncol = 2))
  prior_1 <- c(1/3,1/3,1/3)
  prior_2 <- c(1/2, 1/2)
  res <- .l2_selection_qp(num_neigh = num_neigh,
                          obs_tab = obs_tab,
                          prior_1 = prior_1,
                          prior_2 = prior_2)
  expect_true(is.table(res))
  expect_true(sum(res) == num_neigh)
  expect_true(all(dim(res) == dim(obs_tab)))
})

test_that(".l2_selection_qp works when the requested number is too small", {
  obs_tab <- as.table(matrix(c(10,10,5, 5,5,10), nrow = 3, ncol = 2))
  num_neigh <- 2*sum(obs_tab)
  prior_1 <- c(1/3,1/3,1/3)
  prior_2 <- c(1/2, 1/2)
  res <- .l2_selection_qp(num_neigh = num_neigh,
                          obs_tab = obs_tab,
                          prior_1 = prior_1,
                          prior_2 = prior_2)
  expect_true(is.table(res))
  expect_true(all(res == obs_tab))
})

############################

## .compute_common_snn_softclustering is correct

test_that(".compute_common_snn_softclustering works",{
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
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 1)
  dimred_1 <- .get_postDimred(multiSVD_obj, averaging_mat = NULL)
  multiSVD_obj <- .set_defaultAssay(multiSVD_obj, assay = 2)
  dimred_2 <- .get_postDimred(multiSVD_obj, averaging_mat = NULL)
  
  num_neigh <- 10; bool_cosine <- T; bool_intersect <- T; min_deg <- 5
  snn_1 <- .form_snn_mat(mat = dimred_1,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg)
  snn_2 <- .form_snn_mat(mat = dimred_2,
                         num_neigh = num_neigh,
                         bool_cosine = bool_cosine, 
                         bool_intersect = bool_intersect, 
                         min_deg = min_deg)
  
  res <- .compute_common_snn_softclustering(snn_1 = snn_1,
                                            snn_2 = snn_2,
                                            num_neigh = num_neigh)
  
  expect_true(inherits(res, "dgCMatrix"))
  expect_true(sum(abs(res - Matrix::t(res))) <= 1e-6)
  expect_true(all(rownames(res) == rownames(dimred_1)) && length(rownames(res)) == nrow(res))
  expect_true(all(colnames(res) == rownames(dimred_1)) && length(colnames(res)) == nrow(res))
})
