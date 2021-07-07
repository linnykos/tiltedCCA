context("Test clisi information")

## .clisi is correct

test_that(".clisi works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space1 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))
  common_space2 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))+20
  common_space <- rbind(common_space1, common_space2)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_obj <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  
  membership_vec1 <- as.factor(rep(c("a","b"), each = n/2))
  list_g <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec1,
                           verbose = F, bool_matrix = T)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res1 <- .clisi(c_g, membership_vec1, 1:n, verbose = F)
  
  membership_vec2 <- as.factor(rep(c("a","b"), times = n/2))
  list_g <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec2,
                           verbose = F, bool_matrix = T)
  c_g <- .symmetrize_sparse(list_g[[1]], set_ones = T)
  res2 <- .clisi(c_g, membership_vec2, 1:n, verbose = F)
  
  expect_true(length(res1) == 2)
  expect_true(all(sort(names(res1)) == sort(c("cell_info", "membership_info"))))
  expect_true(all(sort(names(res1$cell_info)) == sort(c("idx", "celltype", "len", "in_ratio", "clisi_score"))))
  expect_true(all(sapply(1:ncol(res1$cell_info), function(j){class(res1$cell_info[,j])}) == c("integer", "factor", "numeric", "numeric", "numeric")))
  expect_true(all(sapply(1:ncol(res1$membership_info), function(j){class(res1$membership_info[,j])}) == c("character", rep("numeric", 6))))
  expect_true(all(sort(names(res1$membership_info)) == sort(c("celltype", "mean_len", "mean_ratio", "mean_clisi", 
                                                              "sd_len", "sd_ratio", "sd_clisi"))))
  expect_true(all(dim(res1$cell_info) == c(n, 5)))
  expect_true(all(res1$membership_info$mean_ratio >= res2$membership_info$mean_ratio))
  expect_true(all(res1$membership_info$mean_clisi >= res2$membership_info$mean_clisi))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space1 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))
  common_space2 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))+20
  common_space <- rbind(common_space1, common_space2)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_obj <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)
  
  membership_vec1 <- as.factor(rep(c("a","b"), each = n/2))
  list_g <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec1,
                           verbose = F, bool_matrix = T)
  res1 <- clisi_information(list_g$c_g, list_g$d_g, list_g$e_g, membership_vec1,
                            verbose = F)
  
  membership_vec2 <- as.factor(rep(c("a","b"), times = n/2))
  list_g <- construct_frnn(dcca_obj, nn = 25, membership_vec = membership_vec2,
                           verbose = F, bool_matrix = T)
  res2 <- clisi_information(list_g$c_g, list_g$d_g, list_g$e_g, membership_vec2,
                            verbose = F)
  
  expect_true(class(res1) == "clisi")
  expect_true(is.list(res1))
  expect_true(all(sort(names(res1)) == sort(c("common_clisi", "distinct_clisi", "everything_clisi"))))
  expect_true(all(dim(res1$common_clisi$cell_info) == dim(res1$distinct_clisi$cell_info)))
  expect_true(all(dim(res1$common_clisi$cell_info) == dim(res1$everything_clisi$cell_info)))
  
  expect_true(all(dim(res1$common_clisi$membership_info) == dim(res1$distinct_clisi$membership_info)))
  expect_true(all(dim(res1$common_clisi$membership_info) == dim(res1$everything_clisi$membership_info)))
  
  expect_true(all(res1$common_clisi$membership_info$mean_clisi >= 
                    res2$common_clisi$membership_info$mean_clisi))
  expect_true(all(res1$distinct_clisi$membership_info$mean_clisi >= 
                    res2$distinct_clisi$membership_info$mean_clisi))
  expect_true(all(res1$everything_clisi$membership_info$mean_clisi >= 
                    res2$everything_clisi$membership_info$mean_clisi))
})
