context("Test clisi information")

## .clisi is correct

test_that(".clisi works", {
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
  graph_list <- .construct_snn_from_tcca(multiSVD_obj)
  
  res1 <- .clisi(cell_subidx = 1:nrow(mat_1), 
                 g = graph_list$snn_common, 
                 membership_vec = test_data$true_membership_vec)
  res2 <- .clisi(cell_subidx = 1:nrow(mat_1), 
                 g = graph_list$snn_distinct_1, 
                 membership_vec = test_data$true_membership_vec)
  res3 <- .clisi(cell_subidx = 1:nrow(mat_1), 
                 g = graph_list$snn_distinct_2, 
                 membership_vec = test_data$true_membership_vec)
  
  expect_true(length(res1) == 3)
  expect_true(all(sort(names(res1)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
  
  expect_true(all(sort(colnames(res1$df)) == sort(c("celltype", "value", "sd"))))
  expect_true(all(dim(res1$df) == c(3,3)))
  expect_true(all(dim(res1$clisi_mat) == c(3,3)))
  expect_true(all(dim(res1$clisi_cell_mat) == c(3,nrow(mat_1))))
  expect_true(all(res1$clisi_cell_mat >= 0))
  expect_true(all(res1$clisi_cell_mat <= 1))
  
  expect_true(all(res2$df[,"value"] <= res1$df[,"value"]))
  expect_true(all(res3$df[,"value"] <= res1$df[,"value"]))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
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
  graph_list <- .construct_snn_from_tcca(multiSVD_obj)
  
  res1 <- clisi_information(multiSVD_obj,
                            membership_vec = test_data$true_membership_vec,
                            verbose = F)
  
  expect_true(class(res1) == "clisi")
  expect_true(is.list(res1))
  expect_true(all(sort(names(res1)) == sort(c("clisi_common", "clisi_distinct_1", "clisi_distinct_2"))))
  expect_true(all(sort(names(res1$common_clisi)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
})

#########################

## .clisi_cell is correct

test_that(".clisi_cell works", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
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
  graph_list <- .construct_snn_from_tcca(multiSVD_obj)
  
  res <- .clisi_cell(g = graph_list$snn_common,
                     membership_vec = test_data$true_membership_vec,
                     position = 1)
  expect_true(length(res) == 3)
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
})
