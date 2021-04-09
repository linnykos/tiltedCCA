
## dcca_factor is correct

test_that("(Basic) dcca_factor works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2",
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "distinct_perc_1"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) dcca_factor preserves rownames and colnames", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  
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

test_that("(Coding) dcca_factor preserves rownames and colnames with metacells", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  nc <- 20
  meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                     apply_shrinkage = T, verbose = F)
  
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

test_that("(Coding) dcca_factor does not require rank_1 = rank_2", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = 2, rank_2 = 5, verbose = F)
  res2 <- dcca_factor(mat_2, mat_1, rank_1 = 5, rank_2 = 2, verbose = F)
  
  expect_true(length(res$cca_obj) == 2)
  expect_true(length(res2$cca_obj) == 2)
  expect_true(all(dim(res$distinct_score_1) == c(n, 2)))
  expect_true(all(dim(res$distinct_score_2) == c(n, 5)))
  expect_true(all(dim(res2$distinct_score_1) == c(n, 5)))
  expect_true(all(dim(res2$distinct_score_2) == c(n, 2)))
})

test_that("(Coding) dcca_factor works with meta-cells", {
  set.seed(5)
  n <- 200; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  nc <- 20
  meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
  
  res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                     apply_shrinkage = T, verbose = F)
  
  expect_true(is.list(res))
  expect_true(class(res) == "dcca")
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "distinct_perc_1"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})
