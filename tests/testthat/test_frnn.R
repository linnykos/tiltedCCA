context("Test frNN")

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
  rad <- .compute_radius(mat, 5, radius_quantile = 0.5)
  res <- .construct_frnn(mat, 
                         num_neigh = 5, 
                         frnn_approx = 0,
                         radius = rad, 
                         resolve_isolated_nodes = T,
                         radius_quantile = NA)
  
  expect_true(is.list(res))
  expect_true(length(res$id) == n)
  expect_true(length(res$dist) == n)
  for(i in 1:n){
    expect_true(length(res$id[[i]]) > 0)
    expect_true(all(res$id[[i]] %% 1 == 0))
    expect_true(!i %in% res$id[[i]])
  }
})

test_that(".construct_frnn works in the presence of outliers", {
  set.seed(10)
  n <- 100; p <- 2
  mat <- matrix(rnorm(n*p), n, p)
  mat <- rbind(mat, matrix(c(10,-10,10,-10), 2, 2))
  n <- nrow(mat)
  rad <- .compute_radius(mat, 5, radius_quantile = 0.5)
  res <- .construct_frnn(mat, 
                         num_neigh = 5, 
                         frnn_approx = 0,
                         radius = rad, 
                         resolve_isolated_nodes = T,
                         radius_quantile = NA)
  
  expect_true(is.list(res))
  expect_true(length(res$id) == n)
  expect_true(length(res$dist) == n)
  for(i in 1:n){
    expect_true(length(res$id[[i]]) > 0)
    expect_true(all(res$id[[i]] %% 1 == 0))
    expect_true(!i %in% res$id[[i]])
  }
})

#####################

## construct_frnn is correct

test_that("construct_frnn works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- test_data$K
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  
  res <- construct_frnn(tilted_res, 
                        num_neigh = 10, 
                        verbose = F, bool_matrix = T)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("c_g", "d_g", "membership_vec", "original_radius"))))
  expect_true(inherits(res$c_g, "dgCMatrix"))
  expect_true(inherits(res$d_g, "dgCMatrix"))
  
  res <- construct_frnn(tilted_res, 
                        num_neigh = 10, 
                        bool_matrix = F,
                        verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("c_g", "d_g", "membership_vec", "original_radius"))))
  expect_true(inherits(res$c_g, "frNN"))
  expect_true(inherits(res$d_g, "frNN"))
})

################################

## combine_frnn is correct

test_that("combine_frnn works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- test_data$K
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  
  set.seed(10)
  list_g_1 <- construct_frnn(tilted_res, 
                             num_neigh = 5, 
                             membership_vec = NA, 
                             data_1 = T, data_2 = F,
                             verbose = F, bool_matrix = T)
  set.seed(10)
  list_g_2 <- construct_frnn(tilted_res, num_neigh = 5, 
                             membership_vec = NA, 
                             data_1 = F, data_2 = T,
                             verbose = F, bool_matrix = T)
  
  set.seed(10)
  res <- combine_frnn(tilted_res, list_g_1$c_g, list_g_2$c_g, num_neigh = 5)
  
  expect_true(all(dim(res) == nrow(mat_1)))
  expect_true(inherits(res, "dgCMatrix"))
})
