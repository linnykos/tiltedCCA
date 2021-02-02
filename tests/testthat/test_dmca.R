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
