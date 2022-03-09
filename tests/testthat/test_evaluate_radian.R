context("Test .evaluate_radian")

## .evaluate_radian is correct

test_that(".evaluate_radian works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  basis_list <- test_data$basis_list
  circle_list <- test_data$circle_list
  averaging_mat <- test_data$averaging_mat
  target_dimred <- test_data$target_dimred
  
  res <- .evaluate_radian(averaging_mat = averaging_mat,
                          basis_list = basis_list, 
                          circle_list = circle_list,
                          enforce_boundary = T,
                          percentage = T,
                          return_common_score_basis = F,
                          score_1 = score_1,
                          score_2 = score_2,
                          snn_bool_cosine = T,
                          snn_bool_intersect = T,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          svd_1 = svd_1, 
                          svd_2 = svd_2,
                          target_dimred = target_dimred)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  
  res <- .evaluate_radian(averaging_mat = averaging_mat,
                          basis_list = basis_list, 
                          circle_list = circle_list,
                          enforce_boundary = T,
                          percentage = T,
                          return_common_score_basis = T,
                          score_1 = score_1,
                          score_2 = score_2,
                          snn_bool_cosine = T,
                          snn_bool_intersect = T,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          svd_1 = svd_1, 
                          svd_2 = svd_2,
                          target_dimred = target_dimred)
  expect_true(is.matrix(res$common_score))
  expect_true(all(dim(res$common_score) == c(nrow(score_1), 2)))
  expect_true(is.matrix(res$common_basis))
  expect_true(all(dim(res$common_basis) == c(nrow(score_1), 2)))
})

test_that(".evaluate_radian is maximized at appropriate values", {
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
    
    res <- sapply(seq(0,1,length.out=9), function(percentage){
      .evaluate_radian(averaging_mat = averaging_mat,
                       basis_list = basis_list, 
                       circle_list = circle_list,
                       enforce_boundary = T,
                       percentage = percentage,
                       return_common_score_basis = F,
                       score_1 = score_1,
                       score_2 = score_2,
                       snn_bool_cosine = T,
                       snn_bool_intersect = T,
                       snn_k = 2,
                       snn_min_deg = 1,
                       snn_num_neigh = 10,
                       svd_1 = svd_1, 
                       svd_2 = svd_2,
                       target_dimred = target_dimred)
    })
    expect_true(length(res) == 9)
    
    if(i == 1){
      expect_true(res[1] <= max(res[5:9]))
    } else if(i == 3){
      expect_true(max(res[1:4]) <= max(res[5:9]))
    } else if(i == 4){
      expect_true(max(res[4:6]) <= max(res[1:3]))
      expect_true(max(res[4:6]) <= max(res[7:9]))
    }
  })
})

