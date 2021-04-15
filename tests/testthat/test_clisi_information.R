context("Test clisi information")

## .compute_radius is correct

test_that(".compute_radius works", {
  set.seed(10)
  mat <- matrix(rnorm(100), 20, 5)
  res <- .compute_radius(mat, 2)
  
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
    .compute_radius(mat, k)
  })
  
  expect_true(all(diff(rad_vec) >= 0))
})

############################

## .construct_frnn is correct

test_that(".construct_frnn works with debugging = T", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  rad <- .compute_radius(mat, 5)
  res <- .construct_frnn(mat, radius = rad, frnn_approx = 0,
                         subsampling_rate = 1, min_subsample = 1,
                         debug = T)
  
  expect_true(is.list(res))
  expect_true(length(res) == n)
  for(i in 1:n){
    if(length(res[[i]] > 0)){
      expect_true(all(res[[i]] %% 1 == 0))
      expect_true(!i %in% res[[i]])
    }
  }
})

test_that(".construct_frnn always outputs a graph with all nodes", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  mat <- rbind(mat, c(100,100))
  rad <- .compute_radius(mat, 50)
  res <- .construct_frnn(mat, radius = rad, frnn_approx = 0,
                         subsampling_rate = 1, min_subsample = 1,
                         debug = F)
  
  expect_true(igraph::vcount(res) == n+1)
})

########################

## .convert_frnn2igraph is correct

test_that(".convert_frnn2igraph is correct", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  rad <- .compute_radius(mat, 50)
  frnn_obj <- .construct_frnn(mat, radius = rad, frnn_approx = 0,
                         subsampling_rate = 1, min_subsample = 1,
                         debug = T)
  res <- .convert_frnn2igraph(frnn_obj)
  
  expect_true(class(res) == "igraph")
  expect_true(igraph::vcount(res) == n)
})

########################

## .connect_graph is correct

test_that(".connect_graph works", {
  set.seed(10)
  n <- 100; p <- 2
  mat1 <- matrix(rnorm(n*p), n, p)
  mat2 <- matrix(rnorm(n*p), n, p)+50
  mat <- rbind(mat1, mat2)
  rad <- .compute_radius(mat, 50)
  frnn_obj <- dbscan::frNN(mat, eps = rad, sort = F, approx = 0)$id
  g <- .convert_frnn2igraph(frnn_obj)
  expect_true(igraph::components(g)$no == 2)
  
  res <- .connect_graph(g, mat, subsampling_rate = 1, min_subsample = 1)
  expect_true(igraph::components(res)$no == 1)
})

######################

## .clisi is correct

test_that(".clisi works", {
  set.seed(10)
  n <- 100; p <- 2
  mat1 <- matrix(rnorm(n*p), n, p)
  mat2 <- matrix(rnorm(n*p), n, p)+50
  mat <- rbind(mat1, mat2)
  rad <- .compute_radius(mat, 10)
  g <- .construct_frnn(mat, rad, frnn_approx = 0, subsampling_rate = 1,
                       min_subsample = 1)
  
  res1 <- .clisi(g, rep(c(1,2), each = n))
  res2 <- .clisi(g, rep(c(1,2), times = n))
  
  expect_true(length(res1) == 2)
  expect_true(all(sort(names(res1)) == sort(c("cell_info", "membership_info"))))
  expect_true(all(sort(names(res1$cell_info)) == sort(c("len", "in_ratio", "clisi_score"))))
  expect_true(all(sort(names(res1$membership_info)) == sort(c("mean_len", "mean_ratio", "mean_clisi", 
                                                              "sd_len", "sd_ratio", "sd_clisi"))))
  expect_true(all(dim(res1$cell_info) == c(2*n, 3)))
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
                           dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1, 
                           membership_vec = sample(c(1,2), size = n, replace = T),
                           p = 10, nn = 10, frnn_approx = 0, 
                           subsampling_rate = 0.1, min_subsample = 10)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_clisi", "distinct_clisi", "everything_clisi"))))
})