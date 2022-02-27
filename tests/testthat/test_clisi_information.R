context("Test clisi information")

## .clisi is correct

test_that(".clisi works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  
  list_g <- construct_frnn(tilted_res, num_neigh = 10, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
  res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  
  expect_true(length(res1) == 3)
  expect_true(all(sort(names(res1)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
  
  expect_true(all(sort(colnames(res1$df)) == sort(c("celltype", "value", "sd"))))
  expect_true(all(dim(res1$df) == c(3,3)))
  expect_true(all(dim(res1$clisi_mat) == c(3,3)))
  expect_true(all(dim(res1$clisi_cell_mat) == c(3,nrow(mat_1))))
  expect_true(all(res1$clisi_cell_mat >= 0))
  expect_true(all(res1$clisi_cell_mat <= 1))
  
  expect_true(all(res1$df[,"value"] <= 0.1))
  expect_true(all(res2$df[,"value"] >= 0.7))
})

test_that(".clisi works for a different setting", {
  # load("tests/assets/test_data4.RData")
  load("../assets/test_data4.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  list_g <- construct_frnn(tilted_res, num_neigh = 10, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res1 <- .clisi(c_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  d_g <- .symmetrize_sparse(list_g[[2]], set_ones = T)
  res2 <- .clisi(d_g, true_membership_vec, 1:nrow(mat_1), verbose = F)
  
  expect_true(all(res1$df[,"value"] <= 0.1))
  expect_true(all(res2$df[,"value"] >= 0.1))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  list_g <- construct_frnn(tilted_res, num_neigh = 10, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  res1 <- clisi_information(list_g$c_g, list_g$d_g, 
                            membership_vec = true_membership_vec,
                            verbose = F)
  
  expect_true(class(res1) == "clisi")
  expect_true(is.list(res1))
  expect_true(all(sort(names(res1)) == sort(c("common_clisi", "distinct_clisi"))))
  expect_true(all(sort(names(res1$common_clisi)) == sort(c("df", "clisi_mat", "clisi_cell_mat"))))
})

#########################

## .clisi_cell is correct

test_that(".clisi_cell works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  list_g <- construct_frnn(tilted_res, num_neigh = 10, 
                           verbose = F, bool_matrix = T, 
                           data_1 = T, data_2 = F)
  g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res <- .clisi_cell(g, 
                     membership_vec = true_membership_vec, 
                     position = 1)
  expect_true(length(res) == 3)
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
})
