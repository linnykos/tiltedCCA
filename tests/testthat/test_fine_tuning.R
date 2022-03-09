context("Test fine tuning")

## fine_tuning is correct

test_that("(Basic) fine_tuning works", {
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
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
  
  res <- fine_tuning(input_obj = multiSVD_obj)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("tcca_obj", "cca_obj") %in% names(res)))
  expect_true(inherits(res$cca_obj, "cca"))
  expect_true(inherits(res$tcca_obj, "tcca"))
  expect_true(all(sort(names(res$cca_obj)) == sort(c("score_1", "score_2", "cca_obj"))))
  expect_true(all(sort(names(res$tcca_obj)) == sort(c("common_basis",
                                                      "common_score",
                                                      "distinct_score_1",
                                                      "distinct_score_2",
                                                      "df_percentage",
                                                      "tilt_perc"))))
  expect_true(length(grep("^tcca*", names(.get_param(res)))) == 3)
  expect_true(all(dim(res$tcca_obj$common_score) == c(n,2)))
  expect_true(length(res$tcca_obj$tilt_perc) == 2)
  expect_true(all(is.na(res$tcca_obj$df_percentage)))
})