context("Test .dcca_common_score")

## .dcca_common_score is correct

test_that("(Basic) .dcca_common_score works", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  
  res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                            num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                            fix_distinct_perc = F,
                            check_alignment = F, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "svd_1", "svd_2", 
                                             "score_1", "score_2", 
                                             "cca_obj", "distinct_score_1", 
                                             "distinct_score_2", "distinct_perc_1"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) .dcca_common_score preserves rownames and colnames", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  
  res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                            num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                            fix_distinct_perc = F,
                            check_alignment = F, verbose = F)
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
})

test_that("(Coding) .dcca_common_score works with K=1 for either", {
  set.seed(5)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  svd_1 <- .spoet(mat_1, K = 1); svd_2 <- .spoet(mat_2, K = K)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                            num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                            fix_distinct_perc = F,
                            check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,1)))
  expect_true(all(dim(res$distinct_score_2) == c(n,K)))
  
  svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = 1)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                            num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                            fix_distinct_perc = F,
                            check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,K)))
  expect_true(all(dim(res$distinct_score_2) == c(n,1)))
  
  svd_1 <- .spoet(mat_1, K = 1); svd_2 <- .spoet(mat_2, K = 1)
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                            num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                            fix_distinct_perc = F,
                            check_alignment = F, verbose = F)
  expect_true(all(dim(res$common_score) == c(n,1)))
  expect_true(all(dim(res$distinct_score_1) == c(n,1)))
  expect_true(all(dim(res$distinct_score_2) == c(n,1)))
})

test_that("(Math) .dcca_common_score yields uncorrelated residuals", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    
    svd_1 <- .spoet(mat_1, K = K); svd_2 <- .spoet(mat_2, K = K)
    
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
    
    res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                              num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                              fix_distinct_perc = F,
                              check_alignment = F, verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .dcca_common_score yields uncorrelated residuals with meta-cells", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    
    svd_1 <- .svd_truncated(mat_1, K,
                            symmetric = F, rescale = F, K_full_rank = F)
    svd_2 <- .svd_truncated(mat_2, K,
                            symmetric = F, rescale = F, K_full_rank = F)
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
    
    num_meta <- max(meta_clustering)
    
    mat_1_meta <- t(sapply(1:num_meta, function(x){
      idx <- which(meta_clustering == x)
      apply(mat_1[idx,,drop = F], 2, mean)
    }))
    
    mat_2_meta <- t(sapply(1:num_meta, function(x){
      idx <- which(meta_clustering == x)
      apply(mat_2[idx,,drop = F], 2, mean)
    }))
    
    cca_res <- .cca(mat_1_meta, mat_2_meta, rank_1 = K, rank_2 = K, return_scores = F)
    
    res <- .dcca_common_score(svd_1, svd_2, cca_res, 
                              num_neigh = max(round(nrow(svd_1$u)/20), 40), 
                              fix_distinct_perc = F,
                              check_alignment = T, verbose = F)
    
    prod_mat <- t(res$distinct_score_1) %*% res$distinct_score_2
    
    sum(abs(prod_mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})