context("Test generate data")

test_that("generate_data works", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
  
  set.seed(10)
  p_1 <- 20; p_2 <- 40
  coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
  coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)
  
  set.seed(10)
  res <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("mat_1", "mat_2", "common_score", 
                                             "distinct_score_1", "distinct_score_2",
                                             "rank_c"))))
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

