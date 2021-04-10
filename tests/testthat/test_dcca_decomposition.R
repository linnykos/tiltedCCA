context("Test dcca_decomposition")

## dcca_decomposition is correct

test_that("(Basic) dcca_decomposition works", {
  set.seed(1)
  n <- 100; K <- 3
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "cca_obj", 
                                             "distinct_score_1", 
                                             "distinct_score_2",
                                             "common_mat_1", "common_mat_2",
                                             "distinct_mat_1", "distinct_mat_2", "distinct_perc_2",
                                             "svd_1", "svd_2", "score_1", "score_2"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) dcca_decomposition preserves rownames and colnames", {
  set.seed(1)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  rownames(mat_1) <- paste0("a", 1:n); rownames(mat_2) <- paste0("a", 1:n)
  colnames(mat_1) <- paste0("b", 1:p1)
  colnames(mat_2) <- paste0("c", 1:p2)
  
  dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  expect_true(length(rownames(res$common_score)) > 1)
  expect_true(length(rownames(res$distinct_score_1)) > 1)
  expect_true(length(rownames(res$distinct_score_2)) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  
  expect_true(length(rownames(res$common_mat_1)) > 1)
  expect_true(length(rownames(res$common_mat_2)) > 1)
  expect_true(length(rownames(res$distinct_mat_1)) > 1)
  expect_true(length(rownames(res$distinct_mat_2)) > 1)
  expect_true(length(colnames(res$common_mat_1)) > 1)
  expect_true(length(colnames(res$common_mat_2)) > 1)
  expect_true(length(colnames(res$distinct_mat_1)) > 1)
  expect_true(length(colnames(res$distinct_mat_2)) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_2)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_2)))
  expect_true(all(colnames(mat_1) == colnames(res$common_mat_1)))
  expect_true(all(colnames(mat_1) == colnames(res$distinct_mat_1)))
  expect_true(all(colnames(mat_2) == colnames(res$common_mat_2)))
  expect_true(all(colnames(mat_2) == colnames(res$distinct_mat_2)))
})

test_that("(Math) dcca_decomposition yields uncorrelated distinct matrices", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields uncorrelated distinct matrices with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields a low-rank matrix", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K+1
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K+1
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields a low-rank matrix with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K+1
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K+1
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields common matrices with the same column space", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition yields common matrices with the same column space with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
    mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K+1, meta_clustering = meta_clustering,
                            verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition is a decomposition under no noise", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 100; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1; mat_2 <- common_space %*% transform_mat_2 
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, 
                            apply_shrinkage = F, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    # correct common_space
    proj_1 <- common_space %*% solve(crossprod(common_space)) %*% t(common_space)
    Ktmp <- Matrix::rankMatrix(res$common_score)
    proj_2 <- res$common_score[,1:Ktmp] %*% solve(crossprod(res$common_score[,1:Ktmp])) %*% t(res$common_score[,1:Ktmp])
    bool1 <- sum(abs(proj_1 - proj_2)) <= 1e-6
    
    # no distinctive componet
    bool2 <- all(abs(res$distinct_mat_1) <= 1e-6)
    bool3 <- all(abs(res$distinct_mat_2) <= 1e-6)
    
    # exact decomposition
    bool4 <- all(abs(mat_1 - res$common_mat_1) <= 1e-6)
    bool5 <- all(abs(mat_2 - res$common_mat_2) <= 1e-6)
    
    bool1 & bool2 & bool3 & bool4 & bool5
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition is a decomposition under no noise with meta-cells", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 200; K <- 2
    common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
    
    p1 <- 5; p2 <- 10
    transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
    transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
    
    mat_1 <- common_space %*% transform_mat_1; mat_2 <- common_space %*% transform_mat_2 
    nc <- 20
    meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster
    
    dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, meta_clustering = meta_clustering,
                            apply_shrinkage = F, verbose = F)
    res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
    
    # correct common_space
    proj_1 <- common_space %*% solve(crossprod(common_space)) %*% t(common_space)
    Ktmp <- Matrix::rankMatrix(res$common_score)
    proj_2 <- res$common_score[,1:Ktmp] %*% solve(crossprod(res$common_score[,1:Ktmp])) %*% t(res$common_score[,1:Ktmp])
    bool1 <- sum(abs(proj_1 - proj_2)) <= 1e-6
    
    # no distinctive componet
    bool2 <- all(abs(res$distinct_mat_1) <= 1e-6)
    bool3 <- all(abs(res$distinct_mat_2) <= 1e-6)
    
    # exact decomposition
    bool4 <- all(abs(mat_1 - res$common_mat_1) <= 1e-6)
    bool5 <- all(abs(mat_2 - res$common_mat_2) <= 1e-6)
    
    bool1 & bool2 & bool3 & bool4 & bool5
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) dcca_decomposition can obtain the same result when fed into itself (i.e., stability/identifiability)", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3)
  K <- ncol(B_mat)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec)
  rho <- 1
  svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  
  set.seed(10)
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- generate_random_orthogonal(p_2, K-1)
  
  set.seed(10)
  dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  
  dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K,
                          apply_shrinkage = F, verbose = F)
  res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
  
  dcca_res2 <- dcca_factor(res$common_mat_1 + res$distinct_mat_1,
                           res$common_mat_2 + res$distinct_mat_2,
                           rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
  res2 <- dcca_decomposition(dcca_res2, rank_c = K, verbose = F)
  
  expect_true(sum(abs(res$common_mat_1 - res2$common_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$common_mat_2 - res2$common_mat_2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_1 - res2$distinct_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_2 - res2$distinct_mat_2)) <= 1e-6)
})
