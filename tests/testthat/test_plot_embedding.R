context("Test embedding functions")

## .extract_matrix_helper is correct

.check_extract_matrix_helper <- function(res, n){
  expect_true(is.matrix(res))
  expect_true(nrow(res) == n)
  expect_true(nrow(res) == n)
  invisible()
}

test_that(".extract_matrix_helper works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  n <- nrow(mat_1)
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  
  res <- .extract_matrix_helper(tilted_res$common_score, tilted_res$distinct_score_1,
                                tilted_res$svd_1, common_bool = T, distinct_bool = T,
                                center = T, renormalize = T)
  .check_extract_matrix_helper(res, n)
})

test_that(".extract_matrix_helper is already centered if input centered", {
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
  
  res <- .extract_matrix_helper(tilted_res$common_score, 
                                tilted_res$distinct_score_1,
                                tilted_res$svd_1, 
                                common_bool = T, 
                                distinct_bool = T,
                                center = F, 
                                renormalize = F)

  expect_true(all(abs(colMeans(res)) <= 1e-6))
})

test_that(".extract_matrix_helper has everything embedding have uncorrelated coordinates", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  true_membership_vec <- test_data$true_membership_vec
  K <- 2
  n <- nrow(mat_1)
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  
  res <- .extract_matrix_helper(tilted_res$common_score, 
                                tilted_res$distinct_score_1,
                                tilted_res$svd_1, 
                                common_bool = T, 
                                distinct_bool = T, 
                                center = F, 
                                renormalize = F)
  
  expect_true(abs(sum(crossprod(res)) - sum(diag(crossprod(res)))) <= 1e-4)
})

test_that(".extract_matrix_helper has correct spectral norm", {
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
  
  res <- .extract_matrix_helper(tilted_res$common_score,
                                tilted_res$distinct_score_1,
                                tilted_res$svd_1, 
                                common_bool = T, 
                                distinct_bool = T,
                                center = F, 
                                renormalize = F)
  
  n <- nrow(mat_1)
  expect_true(abs(max(svd(res)$d) - sqrt(n)) <= 1e-4)
})

