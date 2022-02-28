context("Test .search_tilt_perc")

test_that(".search_tilt_perc works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  averaging_mat <- test_data$averaging_mat
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  target_dimred <- test_data$target_dimred
  basis_list <- test_data$basis_list
  circle_list <- test_data$circle_list
  discretization_gridsize <- 9
  
  res <- .search_tilt_perc(averaging_mat = averaging_mat,
                           basis_list = basis_list,
                           circle_list = circle_list,
                           discretization_gridsize = discretization_gridsize,
                           enforce_boundary = T,
                           score_1 = score_1,
                           score_2 = score_2,
                           snn_bool_intersect = T,
                           snn_k = 2,
                           snn_min_deg = 1,
                           snn_num_neigh = 10,
                           svd_1 = svd_1,
                           svd_2 = svd_2,
                           target_dimred = target_dimred,
                           tol = 1e-3)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("df", "percentage"))))
  expect_true(all(colnames(res$df) == c("percentage", "ratio_val")))
  expect_true(nrow(res$df) == discretization_gridsize)
  expect_true(length(res$percentage) == 1)
  expect_true(is.numeric(res$percentage))
})

test_that(".search_tilt_perc gives reasonable values across different settings", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(c(1,3,4), function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    score_1 <- test_data$score_1
    score_2 <- test_data$score_2
    svd_1 <- test_data$svd_1
    svd_2 <- test_data$svd_2
    basis_list <- test_data$basis_list
    circle_list <- test_data$circle_list
    averaging_mat <- test_data$averaging_mat
    target_dimred <- test_data$target_dimred
    
    res <- .search_tilt_perc(averaging_mat = averaging_mat,
                             basis_list = basis_list,
                             circle_list = circle_list,
                             discretization_gridsize = 9,
                             enforce_boundary = T,
                             score_1 = score_1,
                             score_2 = score_2,
                             snn_bool_intersect = T,
                             snn_k = 2,
                             snn_min_deg = 1,
                             snn_num_neigh = 10,
                             svd_1 = svd_1,
                             svd_2 = svd_2,
                             target_dimred = target_dimred,
                             tol = 1e-3)
    
    if(i == 1){
      expect_true(res$percentage <= 0.25)
    } else if(i == 3){
      expect_true(res$percentage < 0.5)
    } else if(i == 4){
      expect_true(res$percentage >= 0.25)
      expect_true(res$percentage <= 0.75)
    }
  })
})

