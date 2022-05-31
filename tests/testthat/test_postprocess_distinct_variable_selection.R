context("Test postprocess distinct variable selection")

## postprocess_distinct_variable_selection is correct

test_that("postprocess_distinct_variable_selection works", {
  # load("tests/assets/test_data2.RData")
  load("../assets/test_data2.RData")
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
                               min_deg = 10)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  logpval_vec <- seq(0, 1, length.out = ncol(mat_2))
  names(logpval_vec) <- colnames(mat_2)
  res <- postprocess_distinct_variable_selection(input_obj = multiSVD_obj,
                                                 input_assay = 2,
                                                 logpval_vec = logpval_vec,
                                                 num_variables = 2,
                                                 verbose = 0)
  
  expect_true(class(res) == "varSelect")
  expect_true(all(sort(names(res)) == sort(c("selected_variables", "candidate_list", "logpval_vec",
                                             "cor_threshold", "cor_vec_intial"))))
  expect_true(all(res$selected_variables %in% colnames(mat_2)))
  for(i in 1:length(res$candidate_list)){
    expect_true(all(res$candidate_list[[i]] %in% colnames(mat_2)))
  }
})
