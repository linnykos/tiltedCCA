context("Test .common_decomposition")

test_compute_common_score <- function(score_1, score_2, obj_vec = NA){
  stopifnot(nrow(score_1) == nrow(score_2), nrow(score_1) >= ncol(score_1),
            nrow(score_2) >= ncol(score_2), is.matrix(score_1), is.matrix(score_2))
  
  p <- min(c(ncol(score_1), ncol(score_2)))
  if(ncol(score_1) > p){ score_1 <- score_1[,1:p,drop = F]}
  if(ncol(score_2) > p){ score_2 <- score_2[,1:p,drop = F]}
  
  n <- nrow(score_1)
  obj_vec <- diag(crossprod(score_1, score_2))/n
  
  R_vec <- sapply(obj_vec, function(x){x <- min(1,max(x,0)); 1-sqrt((1-x)/(1+x))})
  common_score <- .mult_mat_vec((score_1+score_2)/2, R_vec)
  
  common_score
}

## .common_decomposition is correct

test_that("(Basic) .common_decomposition works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  averaging_mat <- test_data$averaging_mat
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  target_dimred <- test_data$target_dimred
  
  res <- .common_decomposition(averaging_mat = averaging_mat,
                               discretization_gridsize = 9,
                               enforce_boundary = T,
                               fix_tilt_perc = F,
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
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "common_basis", "tilt_perc", "df_percentage"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(all(dim(res$common_basis) == c(nrow(score_1), 2)))
  expect_true(length(res$tilt_perc) == 1)
})

test_that("(Math) .common_decomposition is correct when fix_tilt_perc = T", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    averaging_mat <- test_data$averaging_mat
    cca_res_obj <- test_data$cca_res$obj
    score_1 <- test_data$score_1
    score_2 <- test_data$score_2
    svd_1 <- test_data$svd_1
    svd_2 <- test_data$svd_2
    target_dimred <- test_data$target_dimred
    
    res1 <- .common_decomposition(averaging_mat = averaging_mat,
                                  discretization_gridsize = 9,
                                  enforce_boundary = T,
                                  fix_tilt_perc = 0.5,
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
    res2 <- test_compute_common_score(score_1, score_2, 
                                      obj_vec = cca_res_obj)
    sum(abs(res1$common_score - res2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Coding) .common_decomposition preserves rownames and colnames", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  averaging_mat <- test_data$averaging_mat
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  svd_1 <- test_data$svd_1
  svd_2 <- test_data$svd_2
  target_dimred <- test_data$target_dimred
  
  res <- .common_decomposition(averaging_mat = averaging_mat,
                               discretization_gridsize = 9,
                               enforce_boundary = T,
                               fix_tilt_perc = F,
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
  expect_true(all(dim(res$common_score) == dim(score_1)))
  expect_true(length(rownames(res$common_score)) > 1)
  expect_true(all(rownames(res$common_score) == rownames(score_1)))
  
  # [[note to self: observe that if rownames(score_1) != rownames(score_2), the names would still be rownames(score_1)]]
})

test_that("(Math) .common_decomposition gives sensible numbers in asymmetric information settings and strong clusters", {
  # load("tests/assets/test_data1.RData")
  bool_vec <- sapply(c(1,3), function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    averaging_mat <- test_data$averaging_mat
    cca_res_obj <- test_data$cca_res$obj
    score_1 <- test_data$score_1
    score_2 <- test_data$score_2
    svd_1 <- test_data$svd_1
    svd_2 <- test_data$svd_2
    target_dimred <- test_data$target_dimred
    
    res <- .common_decomposition(averaging_mat = averaging_mat,
                                 discretization_gridsize = 9,
                                 enforce_boundary = T,
                                 fix_tilt_perc = 0.5,
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
    res$tilt_perc <= 0.5
  })
  
  expect_true(all(bool_vec))
})

########################

## .compute_radian is correct

test_that(".compute_radian works", {
  vec1 <- c(1, 0)
  vec2 <- c(1/sqrt(2), 1/sqrt(2))
  circle <- .construct_circle(vec1, vec2)
  
  res <- .compute_radian(percentage_val = 0,
                         enforce_boundary = F,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  expect_true(sum(abs(vec2 - c(x_coord, y_coord))) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 1,
                         enforce_boundary = F,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  expect_true(sum(abs(vec1 - c(x_coord, y_coord))) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 0.5,
                         enforce_boundary = F,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  vec <- c(x_coord, y_coord)
  expect_true(abs(.l2norm(vec1 - vec) - .l2norm(vec2 - vec)) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 0.3,
                         enforce_boundary = F,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  vec <- c(x_coord, y_coord)
  expect_true(.l2norm(vec1 - vec) > .l2norm(vec2 - vec))
})

test_that(".compute_radian always the same value for enforce_boundary = T if vec1 and vec2 are orthogonal", {
  vec1 <- c(1, 0)
  vec2 <- c(0,1)
  circle <- .construct_circle(vec1, vec2)
  
  percentage_vec <- seq(0, 1, length.out = 5)
  vec <- sapply(percentage_vec, function(percentage){
    .compute_radian(percentage_val = percentage,
                    enforce_boundary = T,
                    vec1 = vec1, vec2 = vec2,
                    circle = circle)
  })
  
  expect_true(abs(diff(range(vec))) <= 1e-6)
})

test_that(".compute_radian returns the correct radian when enforce_boundary = T", {
  vec1 <- c(1, 0)
  vec2 <- c(1/sqrt(2), 1/sqrt(2))
  circle <- .construct_circle(vec1, vec2)
  
  res <- .compute_radian(percentage_val = 1,
                         enforce_boundary = T,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  vec <- .position_from_circle(circle, res)
  expect_true(abs(t(vec) %*% vec1 / (.l2norm(vec) * .l2norm(vec1)) - 1) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 0,
                         enforce_boundary = T,
                         vec1 = vec1, vec2 = vec2,
                         circle = circle)
  vec <- .position_from_circle(circle, res)
  expect_true(abs(t(vec) %*% vec2 / (.l2norm(vec) * .l2norm(vec2)) - 1) <= 1e-6)
})

#####################

## .position_from_circle is correct

test_that(".position_from_circle works", {
  circle <- list(center = c(0,0), radius = 2)
  res <- .position_from_circle(circle, radian = 0)
  expect_true(sum(abs(res - c(2,0))) <= 1e-6)
  
  res <- .position_from_circle(circle, radian = pi/2)
  expect_true(sum(abs(res - c(0,2))) <= 1e-6)
  
  res <- .position_from_circle(circle, radian = pi/4)
  expect_true(sum(abs(res - rep(sqrt(2),2))) <= 1e-6)
})

##################

## .select_minimum is correct

test_that(".select_minimum works", {
  x_val <- c(0, 0.25, 0.5, 0.75, 1)
  y_val <- 1:5
  res <- .select_minimum(minimum = T, x_val, y_val)
  
  expect_true(res == 1)
})

test_that(".select_minimum works when there are multiple similar values", {
  x_val <- c(0, 0.25, 0.5, 0.75, 1)
  y_val <- c(3,3,1,1,1)
  res <- .select_minimum(minimum = T, x_val, y_val)
  expect_true(res == 3)
  
  y_val <- c(3,1,1,1,1)
  res <- .select_minimum(minimum = T, x_val, y_val)
  expect_true(res == 3)
  
  y_val <- c(3,2,1,1,2)
  res <- .select_minimum(minimum = T, x_val, y_val)
  expect_true(res == 3)
})
