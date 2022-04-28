context("Test clisi information")

# ## .clisi is correct
# 
# test_that(".clisi works", {
#   # load("tests/assets/test_data1.RData")
#   load("../assets/test_data1.RData")
#   mat_1 <- test_data$mat_1
#   mat_2 <- test_data$mat_2
#   target_dimred <- test_data$target_dimred
#   true_membership_vec <- test_data$true_membership_vec
#   K <- 2
#   
#   tilted_res <- tiltedCCA(mat_1, mat_2, 
#                           dims_1 = 1:K, dims_2 = 1:K, 
#                           target_dimred = target_dimred,
#                           snn_k = 2,
#                           snn_min_deg = 1,
#                           snn_num_neigh = 10,
#                           verbose = F)
#   
#   list_g <- construct_frnn(tilted_res, num_neigh = 10, 
#                            verbose = F, bool_matrix = T, 
#                            data_1 = T, data_2 = F)
#   c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
#   res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
#   d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
#   res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
#   
#   expect_true(length(res1) == 3)
#   expect_true(all(sort(names(res1)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
#   
#   expect_true(all(sort(colnames(res1$df)) == sort(c("celltype", "value", "sd"))))
#   expect_true(all(dim(res1$df) == c(3,3)))
#   expect_true(all(dim(res1$clisi_mat) == c(3,3)))
#   expect_true(all(dim(res1$clisi_cell_mat) == c(3,nrow(mat_1))))
#   expect_true(all(res1$clisi_cell_mat >= 0))
#   expect_true(all(res1$clisi_cell_mat <= 1))
#   
#   expect_true(all(res1$df[,"value"] <= 0.1))
#   expect_true(all(res2$df[,"value"] >= 0.7))
# })
# 
# test_that(".clisi works for a different setting", {
#   # load("tests/assets/test_data4.RData")
#   load("../assets/test_data4.RData")
#   mat_1 <- test_data$mat_1
#   mat_2 <- test_data$mat_2
#   target_dimred <- test_data$target_dimred
#   true_membership_vec <- test_data$true_membership_vec
#   K <- 2
#   
#   tilted_res <- tiltedCCA(mat_1, mat_2, 
#                           dims_1 = 1:K, dims_2 = 1:K, 
#                           target_dimred = target_dimred,
#                           snn_k = 2,
#                           snn_min_deg = 1,
#                           snn_num_neigh = 10,
#                           verbose = F)
#   list_g <- construct_frnn(tilted_res, num_neigh = 10, 
#                            verbose = F, bool_matrix = T, 
#                            data_1 = T, data_2 = F)
#   c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
#   res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
#   d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
#   res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
#   
#   expect_true(all(res1$df[,"value"] <= 0.1))
#   expect_true(all(res2$df[,"value"] >= 0.1))
# })
# 
# #########################
# 
# ## clisi_information is correct
# 
# test_that("clisi_information works", {
#   # load("tests/assets/test_data1.RData")
#   load("../assets/test_data1.RData")
#   mat_1 <- test_data$mat_1
#   mat_2 <- test_data$mat_2
#   target_dimred <- test_data$target_dimred
#   true_membership_vec <- test_data$true_membership_vec
#   K <- 2
#   
#   tilted_res <- tiltedCCA(mat_1, mat_2, 
#                           dims_1 = 1:K, dims_2 = 1:K, 
#                           target_dimred = target_dimred,
#                           snn_k = 2,
#                           snn_min_deg = 1,
#                           snn_num_neigh = 10,
#                           verbose = F)
#   list_g <- construct_frnn(tilted_res, num_neigh = 10, 
#                            verbose = F, bool_matrix = T, 
#                            data_1 = T, data_2 = F)
#   res1 <- clisi_information(list_g$c_g, list_g$d_g, 
#                             membership_vec = true_membership_vec,
#                             verbose = F)
#   
#   expect_true(class(res1) == "clisi")
#   expect_true(is.list(res1))
#   expect_true(all(sort(names(res1)) == sort(c("common_clisi", "distinct_clisi"))))
#   expect_true(all(sort(names(res1$common_clisi)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
# })
# 
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
