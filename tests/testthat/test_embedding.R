context("Test embedding functions")

## .prepare_umap_embedding is correct

test_that(".prepare_umap_embedding works", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  dcca_res <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  dcca_obj <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  res <- .prepare_umap_embedding(dcca_obj)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_score_1", "distinct_score_2", "svd_list"))))
  expect_true(all(dim(res$common_score) == c(nrow(mat_1), K)))
  expect_true(all(dim(res$distinct_score_1) == c(nrow(mat_1), K)))
  expect_true(all(dim(res$distinct_score_2) == c(nrow(mat_1), K)))
  expect_true(all(sort(names(res$svd_list)) == sort(c("e1", "e2"))))
  expect_true(all(sort(names(res$svd_list$e1)) == sort(c("u", "d", "v"))))
  expect_true(all(sort(names(res$svd_list$e2)) == sort(c("u", "d", "v"))))
})

#############################

## .extract_matrix_helper is correct

test_that(".extract_matrix_helper works", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  dcca_res <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  dcca_obj <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  prep_obj <- .prepare_umap_embedding(dcca_obj)
  
  res <- .extract_matrix_helper(prep_obj$common_score, prep_obj$distinct_score_1,
                                prep_obj$svd_list$e1, common_bool = T, distinct_bool = T, add_noise = T,
                                center = T, renormalize = T)
  
  expect_true(is.matrix(res))
  expect_true(nrow(res) == nrow(prep_obj$common_score))
  expect_true(nrow(res) == nrow(prep_obj$distinct_score_1))
})

##############################

## .extract_svd_embedding is correct
test_that(".extract_svd_embedding works for dcca_decomposition", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  dcca_res <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  dcca_obj <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  res <- .extract_svd_embedding(dcca_obj)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("c1", "c2", "d1", "d2", "e1", "e2"))))
})
