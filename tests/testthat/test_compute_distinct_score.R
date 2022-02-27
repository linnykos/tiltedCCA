context("Test .compute_distinct_score")

## .compute_distinct_score is correct

test_that("(Basic) .compute_distinct_score works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  score_1 <- test_data$score_1
  score_2 <- test_data$score_2
  common_score <- test_data$common_score
  
  res <- .compute_distinct_score(score_1, score_2, common_score)
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("distinct_score_1", "distinct_score_2"))))
  expect_true(all(dim(res$distinct_score_1) == dim(score_1)))
  expect_true(all(dim(res$distinct_score_2) == dim(score_2)))
})

test_that("(Math) .compute_distinct_score generates orthogonal distinct matrices", {
  # load("tests/assets/test_data1.RData")
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    score_1 <- test_data$score_1
    score_2 <- test_data$score_2
    common_score <- test_data$common_score
    
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    all(abs(crossprod(res$distinct_score_1, res$distinct_score_2)) <= 1e-6)
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .compute_distinct_score generates equal-lengthed matrices when fix_tilt_perc = T", {
  # load("tests/assets/test_data1.RData")
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    score_1 <- test_data$score_1
    score_2 <- test_data$score_2
    common_score <- test_data$common_score
    
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    diff_vec1 <- sapply(1:ncol(common_score), function(k){
      .l2norm(common_score[,k] - res$distinct_score_1[,k])
    })
    diff_vec2 <- sapply(1:ncol(common_score), function(k){
      .l2norm(common_score[,k] - res$distinct_score_2[,k])
    })
    
    sum(abs(diff_vec1 - diff_vec2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

