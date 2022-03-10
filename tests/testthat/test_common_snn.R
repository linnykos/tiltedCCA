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
