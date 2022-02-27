context("Test .convert_common_score_to_mat")

## .convert_common_score_to_mat is correct
test_that(".convert_common_score_to_mat returns reasonable values", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  common_score <- test_data$common_score
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  n <- nrow(score_1)
  
  res <- .convert_common_score_to_mat(common_score,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  expect_true(all(dim(res) == c(nrow(score_1), 2*ncol(score_1))))
  
  res <- .convert_common_score_to_mat(score_1,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  rescaling_factor <- max(c(svd_1$d, svd_2$d))
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_1 <- dimred_1/svd_1$d[1]*rescaling_factor
  target <- tcrossprod(score_1) %*% dimred_1/n
  expect_true(sum(abs(res[,1:2] - target)) <= 1e-6)
  expect_true(abs(max(svd(dimred_1)$d) - max(svd(target)$d)) <= 1e-6)
  
  res <- .convert_common_score_to_mat(score_2,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  rescaling_factor <- max(c(svd_1$d, svd_2$d))
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  dimred_2 <- dimred_2/svd_2$d[1]*rescaling_factor
  target <- tcrossprod(score_2) %*% dimred_2/n
  expect_true(sum(abs(res[,3:4] - target)) <= 1e-6)
  expect_true(abs(max(svd(dimred_2)$d) - max(svd(target)$d)) <= 1e-6)
})
