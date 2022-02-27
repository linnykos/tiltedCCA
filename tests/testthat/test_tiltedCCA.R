context("Test tiltedCCA")

## tiltedCCA is correct

test_that("(Basic) tiltedCCA works", {
  # load("tests/assets/test_data1.RData")
  load("assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- 2
  
  res <- tiltedCCA(mat_1, mat_2, 
                   dims_1 = 1:K, dims_2 = 1:K, 
                   target_dimred = target_dimred,
                   verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "tiltedCCA")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2",
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "tilt_perc",
                                             "target_dimred",
                                             "averaging_mat",
                                             "param_list",
                                             "df_percentage"))))
  expect_true(all(dim(res$common_score) == c(nrow(mat_1), 2)))
})


test_that("(Basic) tiltedCCA works with variable dimensions", {
  # load("tests/assets/test_data1.RData")
  load("assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  
  set.seed(10)
  mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1))), nrow = nrow(mat_1), ncol = ncol(mat_1))
  mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2))), nrow = nrow(mat_2), ncol = ncol(mat_2))
  
  res <- tiltedCCA(mat_1, mat_2, 
                   dims_1 = 1:4, dims_2 = 2:3, 
                   target_dimred = target_dimred,
                   verbose = F)
  
  n <- nrow(mat_1)
  expect_true(is.list(res))
  expect_true(class(res) == "tiltedCCA")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2",
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "tilt_perc",
                                             "target_dimred",
                                             "averaging_mat",
                                             "param_list",
                                             "df_percentage"))))
  expect_true(all(dim(res$common_score) == c(n,2)))
  expect_true(all(dim(res$distinct_score_1) == c(n,4)))
  expect_true(all(dim(res$distinct_score_2) == c(n,2)))
})

test_that("(Basic) tiltedCCA works with a sparse matrix", {
  # load("tests/assets/test_data1.RData")
  load("assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  
  n <- nrow(mat_1)
  mat_1[sample(1:prod(dim(mat_1)), round(prod(dim(mat_1))/2))] <- 0
  mat_1 <- Matrix::Matrix(mat_1, sparse = T)
  mat_2[sample(1:prod(dim(mat_2)), round(prod(dim(mat_1))/2))] <- 0
  mat_2 <- Matrix::Matrix(mat_2, sparse = T)
  res <- tiltedCCA(mat_1, mat_2, 
                   dims_1 = 1:2, dims_2 = 1:2, 
                   target_dimred = target_dimred,
                   verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "tiltedCCA")
  expect_true(all(dim(res$common_score) == c(n,2)))
  expect_true(all(dim(res$distinct_score_1) == c(n,2)))
  expect_true(all(dim(res$distinct_score_2) == c(n,2)))
})

test_that("(Coding) tiltedCCA preserves rownames and colnames", {
  # load("tests/assets/test_data1.RData")
  load("assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  
  n <- nrow(mat_1)
  p1 <- ncol(mat_1); p2 <- ncol(mat_2)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  K <- 2
  
  res <- tiltedCCA(mat_1, mat_2, 
                   dims_1 = 1:K, dims_2 = 1:K, 
                   target_dimred = target_dimred,
                   verbose = F)
  
  expect_true(length(res$common_score) > 1)
  expect_true(length(res$distinct_score_1) > 1)
  expect_true(length(res$distinct_score_2) > 1)
  expect_true(length(res$score_1) > 1)
  expect_true(length(res$score_2) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  expect_true(all(rownames(mat_1) == rownames(res$score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$score_2)))
  
  expect_true(length(rownames(res$svd_1$u)) > 1)
  expect_true(length(rownames(res$svd_2$u)) > 1)
  expect_true(length(rownames(res$svd_1$v)) > 1)
  expect_true(length(rownames(res$svd_2$v)) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$svd_1$u)))
  expect_true(all(rownames(mat_1) == rownames(res$svd_2$u)))
  expect_true(all(colnames(mat_1) == rownames(res$svd_1$v)))
  expect_true(all(colnames(mat_2) == rownames(res$svd_2$v)))
})

test_that("(Math) tiltedCCA is symmetric if the arguments are flipped", {
  # load("tests/assets/test_data1.RData")
  load("assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- 2
  
  set.seed(10)
  res <- tiltedCCA(mat_1, mat_2, 
                   dims_1 = 1:K, dims_2 = 1:K, 
                   target_dimred = target_dimred,
                   verbose = F)
  set.seed(10)
  res2 <- tiltedCCA(mat_2, mat_1, 
                    dims_1 = 1:K, dims_2 = 1:K,
                    target_dimred = target_dimred, 
                    verbose = F)
  
  expect_true(abs(res$tilt_perc - (1-res2$tilt_perc)) <= 1e-6)
  
  tmp1 <- res$common_score
  tmp2 <- res2$common_score
  for(j in 1:ncol(tmp1)){
    if(sign(sum(tmp1[1:10,j])) != sign(sum(tmp2[1:10,j]))) {
      tmp1[,j] <- -tmp1[,j]
      res$distinct_score_1[,j] <- -res$distinct_score_1[,j]
      res$distinct_score_2[,j] <- -res$distinct_score_2[,j]
    }
  }
  
  expect_true(sum(abs(tmp1 - tmp2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_score_1 - res2$distinct_score_2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_score_2 - res2$distinct_score_1)) <= 1e-6)
})
