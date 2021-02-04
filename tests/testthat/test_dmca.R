context("Test D-MCA")

## .mca_coef_mat is correct

test_that("(Math) .mca_coef_mat actually recovers the original matrix", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    low_rank <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat <- tcrossprod(low_rank)
    svd_res <- .svd_truncated(mat, K = 2)
    
    low_rank2 <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat2 <- tcrossprod(low_rank2)
    svd_res2 <- .svd_truncated(mat2, K = 2)
    
    res <- .mca_coef_mat(svd_res$v, svd_res2$u)
    
    proj_mat <- mat %*% svd_res2$u
    
    sum(abs(proj_mat %*% res - mat)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

########################

## .mca is correct

test_that("(Math) .mca yields correct answer", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    low_rank1 <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat1 <- tcrossprod(low_rank1)
    svd_1 <- .svd_truncated(mat1, K = 2)
    mat1 <- .recover_mat_from_svd(svd_1)
    
    low_rank2 <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat2 <- tcrossprod(low_rank2)
    svd_2 <- .svd_truncated(mat2, K = 2)
    mat2 <- .recover_mat_from_svd(svd_2)
    
    res <- .mca(svd_1, svd_2, rank = 2)
    res2 <- .svd_truncated(t(mat1) %*% mat2, K = 2)
    
    # rotation for u
    svd_u <- svd(t(res$u) %*% res2$u)
    bool1 <- sum(abs(res$u %*% svd_u$u %*% t(svd_u$v) - res2$u)) <= 1e-6
    
    # rotation for u
    svd_v <- svd(t(res$v) %*% res2$v)
    bool2 <- sum(abs(res$v %*% svd_v$u %*% t(svd_v$v) - res2$v)) <= 1e-6
    
    bool3 <- sum(abs(res$d - res2$d)) <= 1e-6
    
    bool1 & bool2 & bool3
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .mca yields uncorrelated scores", {
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
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    
    res <- .mca(svd_1, svd_2, rank = 2)
    
    tmp1 <- mat_1 %*% res$u
    tmp2 <- mat_2 %*% res$v
    prod_mat <- t(tmp1) %*% tmp2
    
    abs(sum(abs(prod_mat)) - sum(diag(abs(prod_mat)))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

##################################

## .mca_common_score is correct

test_that(".mca_common_score works", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    low_rank1 <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat1 <- tcrossprod(low_rank1)
    svd_1 <- .svd_truncated(mat1, K = 2)
    mat1 <- .recover_mat_from_svd(svd_1)
    
    low_rank2 <- matrix(rnorm(50), nrow = 25, ncol = 2)
    mat2 <- tcrossprod(low_rank2)
    svd_2 <- .svd_truncated(mat2, K = 2)
    mat2 <- .recover_mat_from_svd(svd_2)
    
    mca_res <- .mca(svd_1, svd_2, rank = 2)
    res <- .mca_common_score(svd_1, svd_2, mca_res, verbose = F)
    
    bool1 <- all(dim(res$score_1) == dim(res$common_mat_1)) & all(dim(res$score_2) == dim(res$common_mat_2))
    bool2 <- nrow(res$score_1) == nrow(svd_1$u) & nrow(res$score_2) == nrow(svd_2$u)
    
    bool1 & bool2
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .mca_common_score yields uncorrelated scores", {
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
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    
    mca_res <- .mca(svd_1, svd_2, rank = 2)
    res <- .mca_common_score(svd_1, svd_2, mca_res, verbose = F)
    
    prod_mat <- t(res$score_1) %*% res$score_2
    
    abs(sum(abs(prod_mat)) - sum(diag(abs(prod_mat)))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .mca_common_score yields orthogonal distinct vectors", {
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
    svd_1 <- svd(mat_1); svd_2 <- svd(mat_2)
    
    mca_res <- .mca(svd_1, svd_2, rank = 2)
    res <- .mca_common_score(svd_1, svd_2, mca_res, verbose = F)
    
    distinct_mat1 <- res$score_1 - res$common_mat_1
    distinct_mat2 <- res$score_2 - res$common_mat_2
    
    ang_vec <- sapply(1:ncol(distinct_mat1), function(j){
      abs(.angle_between_vectors(distinct_mat1[,j], distinct_mat2[,j]) - 90)
    })
    
    all(ang_vec <= 1)
  })
  
  expect_true(all(bool_vec))
})

###############################

## dmca_factor is correct

test_that("(Basic) dmca_factor works", {
  set.seed(1)
  low_rank1 <- matrix(rnorm(50), nrow = 25, ncol = 2)
  mat1 <- tcrossprod(low_rank1) 
  mat1 <- mat1 + matrix(rnorm(prod(dim(mat1)), sd = 0.1), nrow = nrow(mat1), ncol = ncol(mat1))
  low_rank2 <- matrix(rnorm(50), nrow = 25, ncol = 2)
  mat2 <- cbind(tcrossprod(low_rank2), matrix(rnorm(nrow(low_rank2)*5), nrow(low_rank2), 5))
  mat2 <- mat2 + matrix(rnorm(prod(dim(mat2)), sd = 0.1), nrow = nrow(mat2), ncol = ncol(mat2))
  n <- nrow(mat1); p1 <- ncol(mat1); p2 <- ncol(mat2)
  rownames(mat1) <- paste0("n", 1:n)
  rownames(mat2) <- paste0("n", 1:n)
  colnames(mat1) <- paste0("p", 1:p1)
  colnames(mat2) <- paste0("d", 1:p2)
  
  res <- dmca_factor(mat1, mat2, rank_1 = 2, rank_2 = 2, verbose = F, apply_shrinkage = F)
  
  expect_true(class(res) == "dmca")
  expect_true(all(c("u","d","v") %in% names(res$svd_1)))
  expect_true(all(c("u","d","v") %in% names(res$svd_2)))
  expect_true(all(dim(res$svd_1$d) == 2))
  expect_true(all(dim(res$svd_2$d) == 2))
  expect_true(all(dim(res$svd_1$u) == c(n, 2)))
  expect_true(all(dim(res$svd_2$u) == c(n, 2)))
  expect_true(all(dim(res$svd_1$v) == c(p1, 2)))
  expect_true(all(dim(res$svd_2$v) == c(p2, 2)))
  expect_true(all(dim(res$score_1) == c(n, 2)))
  expect_true(all(dim(res$score_2) == c(n, 2)))
  expect_true(all(dim(res$common_mat_1) == c(n, 2)))
  expect_true(all(dim(res$common_mat_2) == c(n, 2)))
  
  expect_true(all(rownames(res$svd_1$u) == rownames(mat1)))
  expect_true(all(rownames(res$svd_2$u) == rownames(mat2)))
  expect_true(all(rownames(res$svd_1$v) == colnames(mat1)))
  expect_true(all(rownames(res$svd_2$v) == colnames(mat2)))
  expect_true(all(rownames(res$score_1) == rownames(mat1)))
  expect_true(all(rownames(res$score_2) == rownames(mat2)))
  expect_true(all(rownames(res$common_mat_1) == rownames(mat1)))
  expect_true(all(rownames(res$common_mat_2) == rownames(mat2)))
})

##############################

## dmca_decomposition is correct

test_that("dmca_decomposition works", {
  set.seed(1)
  n <- 100; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)
  
  rownames(mat_1) <- paste0("n", 1:n)
  rownames(mat_2) <- paste0("n", 1:n)
  colnames(mat_1) <- paste0("p", 1:p1)
  colnames(mat_2) <- paste0("d", 1:p2)
  
  dmca_res <- dmca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = F, apply_shrinkage = F)
  res <- dmca_decomposition(dmca_res, verbose = F)
  
  expect_true(class(res) == "dmca_decomp")
  expect_true(all(sort(names(res)) == sort(c("common_mat_1", "common_mat_2", "distinct_mat_1", "distinct_mat_2", "mca_obj"))))
  expect_true(all(dim(res$common_mat_1) == dim(mat_1)))
  expect_true(all(dim(res$common_mat_2) == dim(mat_2)))
  expect_true(all(dim(res$distinct_mat_1) == dim(mat_1)))
  expect_true(all(dim(res$distinct_mat_2) == dim(mat_2)))
  
  expect_true(all(rownames(res$common_mat_1) == rownames(mat_1)))
  expect_true(all(rownames(res$common_mat_2) == rownames(mat_2)))
  expect_true(all(rownames(res$distinct_mat_1) == rownames(mat_1)))
  expect_true(all(rownames(res$distinct_mat_2) == rownames(mat_2)))
  
  expect_true(all(colnames(res$common_mat_1) == colnames(mat_1)))
  expect_true(all(colnames(res$common_mat_2) == colnames(mat_2)))
  expect_true(all(colnames(res$distinct_mat_1) == colnames(mat_1)))
  expect_true(all(colnames(res$distinct_mat_2) == colnames(mat_2)))
})
