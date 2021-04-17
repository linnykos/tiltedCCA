context("Test clisi information")

## .compute_radius is correct

test_that(".compute_radius works", {
  set.seed(10)
  mat <- matrix(rnorm(100), 20, 5)
  res <- .compute_radius(mat, 2, radius_quantile = 0.95)
  
  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
})

test_that(".compute_radius increases as neighbor increases", {
  set.seed(10)
  n <- 100; p <- 5
  mat <- matrix(rnorm(n*p), n, p)
  k_vec <- 5:50
  rad_vec <- sapply(k_vec, function(k){
    .compute_radius(mat, k, radius_quantile = 0.95)
  })
  
  expect_true(all(diff(rad_vec) >= 0))
})

############################

## .construct_frnn is correct

test_that(".construct_frnn works", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  rad <- .compute_radius(mat, 5, radius_quantile = 0.95)
  res <- .construct_frnn(mat, radius = rad, frnn_approx = 0)
  
  expect_true(is.list(res))
  expect_true(length(res) == n)
  for(i in 1:n){
    if(length(res[[i]] > 0)){
      expect_true(all(res[[i]] %% 1 == 0))
      expect_true(!i %in% res[[i]])
    }
  }
})

######################

## .clisi is correct

test_that(".clisi works", {
  set.seed(10)
  n <- 100; p <- 2
  mat1 <- matrix(rnorm(n*p), n, p)
  mat2 <- matrix(rnorm(n*p), n, p)+50
  mat <- rbind(mat1, mat2)
  rad <- .compute_radius(mat, 10, radius_quantile = 0.95)
  g <- .construct_frnn(mat, rad, frnn_approx = 0)
  membership_vec1 <- as.factor(rep(c("a","b"), each = n))
  membership_vec2 <- as.factor(rep(c("a","b"), times = n))
  cell_subidx1 <- .construct_celltype_subsample(membership_vec1, 50)
  cell_subidx2 <- .construct_celltype_subsample(membership_vec2, 50)
  
  res1 <- .clisi(g, membership_vec1, cell_subidx1)
  res2 <- .clisi(g, membership_vec2, cell_subidx2)
  
  expect_true(class(res1) == "clisi")
  expect_true(length(res1) == 2)
  expect_true(all(sort(names(res1)) == sort(c("cell_info", "membership_info"))))
  expect_true(all(sort(names(res1$cell_info)) == sort(c("idx", "celltype", "len", "in_ratio", "clisi_score"))))
  expect_true(all(sapply(1:ncol(res1$cell_info), function(j){class(res1$cell_info[,j])}) == c("integer", "factor", "numeric", "numeric", "numeric")))
  expect_true(all(sapply(1:ncol(res1$membership_info), function(j){class(res1$membership_info[,j])}) == c("character", rep("numeric", 6))))
  expect_true(all(sort(names(res1$membership_info)) == sort(c("celltype", "mean_len", "mean_ratio", "mean_clisi", 
                                                              "sd_len", "sd_ratio", "sd_clisi"))))
  expect_true(all(dim(res1$cell_info) == c(length(cell_subidx1), 5)))
  expect_true(all(res1$membership_info$mean_ratio >= res2$membership_info$mean_ratio))
  expect_true(all(res1$membership_info$mean_clisi >= res2$membership_info$mean_clisi))
})

#########################

## clisi_information is correct

test_that("clisi_information works", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  p1 <- 10; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  res <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                           membership_vec = as.factor(sample(c(1,2), size = n, replace = T)),
                           rank_c = p1, rank_d = p1, nn = 10, frnn_approx = 0, 
                           verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_clisi", "distinct_clisi", "everything_clisi"))))
})