context("Test generate data")

test_that("generate_data works", {
  set.seed(10)
  n <- 100; rank_12 <- 2; rank_1 <- 2; rank_2 <- 2; p_1 <- 20; p_2 <- 40
  common_loading <- matrix(0, n, rank_12)
  diag(common_loading) <- 1
  distinct_loading_1 <- matrix(0, n, rank_1)
  for(i in 1:2){
    distinct_loading_1[i+2,i] <- 1
  }
  distinct_loading_2 <- matrix(0, n, rank_2)
  for(i in 1:2){
    distinct_loading_2[i+4,i] <- 1
  }
  coef_mat_1 <- matrix(runif(rank_1*p_1), rank_1, p_1)
  coef_mat_2 <- matrix(runif(rank_2*p_2), rank_2, p_2)
  noise_func <- function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}
  
  res <- generate_data(common_loading, distinct_loading_1, distinct_loading_2,
                       coef_mat_1, coef_mat_2, noise_func = noise_func)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2"))))
  expect_true(all(sapply(res, is.matrix)))
  expect_true(all(dim(res$mat_1) == c(n, p_1)))
  expect_true(all(dim(res$mat_2) == c(n, p_2)))
})

test_that("generate_data can handle different ranks", {
  set.seed(10)
  n <- 100; rank_12 <- 2; rank_1 <- 4; rank_2 <- 5; p_1 <- 20; p_2 <- 40
  common_loading <- matrix(0, n, rank_12)
  diag(common_loading) <- 1
  distinct_loading_1 <- matrix(0, n, rank_1)
  for(i in 1:rank_1){
    distinct_loading_1[i+rank_12,i] <- 1
  }
  distinct_loading_2 <- matrix(0, n, rank_2)
  for(i in 1:rank_2){
    distinct_loading_2[i+rank_12+rank_1,i] <- 1
  }
  coef_mat_1 <- matrix(runif(rank_1*p_1), rank_1, p_1)
  coef_mat_2 <- matrix(runif(rank_2*p_2), rank_2, p_2)
  noise_func <- function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}
  
  res <- generate_data(common_loading, distinct_loading_1, distinct_loading_2,
                       coef_mat_1, coef_mat_2, noise_func = noise_func)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2"))))
  expect_true(all(sapply(res, is.matrix)))
  expect_true(all(dim(res$mat_1) == c(n, p_1)))
  expect_true(all(dim(res$mat_2) == c(n, p_2)))
})

#################

# generate_sbm_orthogonal is correct

test_that("generate_sbm_orthogonal works", {
  set.seed(10)
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  rho <- 0.5
  membership_vec <- c(rep(1,20), rep(2,20), rep(3,20))
  
  res <- generate_sbm_orthogonal(rho*B, membership_vec)
  
  expect_true(all(dim(res) == c(60, 3)))
  tmp <- crossprod(res)
  expect_true(abs(sum(abs(tmp)) - sum(abs(diag(tmp)))) <= 1e-6)
})

##################


## .generate_membership_matrix is correct

test_that(".generate_membership_matrix works", {
  res <- .generate_membership_matrix(sample(1:5, size = 100, replace = T))
  expect_true(all(dim(res) == c(100, 5)))
  expect_true(all(apply(res, 1, function(x){sum(x) == 1 & length(which(x==1)) == 1})))
})

## .compute_prob_mat is correct

test_that(".compute_prob_mat works", {
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  rho <- 0.5
  
  n <- 500
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
  
  prob_mat <- .compute_prob_mat(rho*B, membership_vec)
  
  expect_true(nrow(prob_mat) == ncol(prob_mat))
  expect_true(sum(abs(prob_mat - t(prob_mat))) <= 1e-4)
})

## .generate_adjaceny_mat is correct

test_that(".generate_adjaceny_mat works", {
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  rho <- 0.5
  
  n <- 500
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
  
  prob_mat <- .compute_prob_mat(rho*B, membership_vec)
  adj_mat <- .generate_adjaceny_mat(prob_mat)
  
  expect_true(nrow(adj_mat) == ncol(adj_mat))
  expect_true(sum(abs(adj_mat - t(adj_mat))) <= 1e-4)
  expect_true(all(adj_mat %in% c(0,1)))
})

