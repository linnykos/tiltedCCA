context("Test .tiltedCCA_common_score")

## .tiltedCCA_common_score is correct

test_that("(Basic) .tiltedCCA_common_score works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  cca_res <- test_data$cca_res
  averaging_mat <- test_data$averaging_mat
  target_dimred <- test_data$target_dimred
  
  res <- .tiltedCCA_common_score(averaging_mat = averaging_mat,
                                 cca_res = cca_res, 
                                 discretization_gridsize = 9,
                                 enforce_boundary = T,
                                 fix_tilt_perc = F, 
                                 snn_bool_cosine = T,
                                 snn_bool_intersect = T,
                                 snn_k = 2,
                                 snn_min_deg = 1,
                                 snn_num_neigh = 10,
                                 svd_1 = svd_1, 
                                 svd_2 = svd_2, 
                                 target_dimred = target_dimred,
                                 verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "df_percentage",
                                             "common_basis", "target_dimred",
                                             "tilt_perc"))))
  expect_true(all(dim(res$common_score) == c(nrow(mat_1), 2)))
})

test_that("(Coding) .tiltedCCA_common_score preserves rownames and colnames", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  cca_res <- test_data$cca_res
  averaging_mat <- test_data$averaging_mat
  target_dimred <- test_data$target_dimred
  
  res <- .tiltedCCA_common_score(averaging_mat = averaging_mat,
                                 cca_res = cca_res, 
                                 discretization_gridsize = 9,
                                 enforce_boundary = T,
                                 fix_tilt_perc = F, 
                                 snn_bool_cosine = T,
                                 snn_bool_intersect = T,
                                 snn_k = 2,
                                 snn_min_deg = 1,
                                 snn_num_neigh = 10,
                                 svd_1 = svd_1, 
                                 svd_2 = svd_2, 
                                 target_dimred = target_dimred,
                                 verbose = F)
  
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
})

test_that("(Math) .tiltedCCA_common_score yields uncorrelated residuals", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    svd_1 <- test_data$svd_1
    svd_2 <- test_data$svd_2
    cca_res <- test_data$cca_res
    averaging_mat <- test_data$averaging_mat
    target_dimred <- test_data$target_dimred
    
    
    res <- .tiltedCCA_common_score(averaging_mat = averaging_mat,
                                   cca_res = cca_res, 
                                   discretization_gridsize = 9,
                                   enforce_boundary = T,
                                   fix_tilt_perc = F, 
                                   snn_bool_cosine = T,
                                   snn_bool_intersect = T,
                                   snn_k = 2,
                                   snn_min_deg = 1,
                                   snn_num_neigh = 10,
                                   svd_1 = svd_1, 
                                   svd_2 = svd_2, 
                                   target_dimred = target_dimred,
                                   verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})
