context("Test .compute_unnormalized_scores")

## .compute_unnormalized_scores is correct

test_that("(Math) .compute_unnormalized_scores computes the correct scores", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 150; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
    res <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    rank_val <- Matrix::rankMatrix(.compute_cca_aggregate_matrix(svd_1, svd_2, augment = T))
    
    bool1 <- sum(crossprod(res$score_1) - n*diag(p1)) <= 1e-4
    bool2 <- sum(crossprod(res$score_2) - n*diag(p2)) <= 1e-4
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})