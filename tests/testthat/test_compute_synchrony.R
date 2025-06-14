context("Test compute synchrony")

## tiltedCCA_decomposition is correct

test_that("(Basic) compute_synchrony works", {
  # load("tests/assets/test_data4.RData")
  load("../assets/test_data4.RData")
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
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  res <- compute_synchrony(
    multiSVD_obj,
    anchor_modality = 1
  )
  
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(nrow(multiSVD_obj$common_mat_1), 2))
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
  expect_true(length(rownames(res)) == nrow(multiSVD_obj$common_mat_1))
  expect_equal(rownames(res), rownames(multiSVD_obj$common_mat_1))
})