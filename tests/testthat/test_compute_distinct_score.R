context("Test .compute_distinct_score")

## .compute_distinct_score is correct

test_that("(Basic) .compute_distinct_score works", {
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
  
  tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  nn_1 <- RANN::nn2(score_1, k = 50)$nn.idx
  nn_2 <- RANN::nn2(score_1, k = 50)$nn.idx
  
  tmp <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)
  common_score <- tmp$common_score; distinct_perc_1 <- tmp$distinct_perc_1
  
  res <- .compute_distinct_score(score_1, score_2, common_score)
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("distinct_score_1", "distinct_score_2"))))
  expect_true(all(dim(res$distinct_score_1) == dim(score_1)))
  expect_true(all(dim(res$distinct_score_2) == dim(score_2)))
})

test_that("(Math) .compute_distinct_score generates orthogonal distinct matrices", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    B_mat <- matrix(c(0.9, 0.4, 0.1,
                      0.4, 0.9, 0.1,
                      0.1, 0.1, 0.3), 3, 3)
    K <- ncol(B_mat); n_clust <- 50; rho <- 0.1
    
    true_membership_vec <- rep(1:4, each = n_clust)
    n <- length(true_membership_vec)
    
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
    score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
    score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    nn_1 <- RANN::nn2(score_1, k = 50)$nn.idx
    nn_2 <- RANN::nn2(score_1, k = 50)$nn.idx
    
    tmp <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)
    common_score <- tmp$common_score; distinct_perc_1 <- tmp$distinct_perc_1
    
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    all(abs(crossprod(res$distinct_score_1, res$distinct_score_2)) <= 1e-6)
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .compute_distinct_score generates equal-lengthed matrices when fix_distinct_perc = T", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    B_mat <- matrix(c(0.9, 0.4, 0.1,
                      0.4, 0.9, 0.1,
                      0.1, 0.1, 0.3), 3, 3)
    K <- ncol(B_mat); n_clust <- 50; rho <- 0.1
    
    true_membership_vec <- rep(1:4, each = n_clust)
    n <- length(true_membership_vec)
    
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
    score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
    score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
    
    tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    tmp <- .common_decomposition(score_1, score_2, nn_1 = NA, nn_2 = NA, fix_distinct_perc = T)
    common_score <- tmp$common_score; distinct_perc_1 <- tmp$distinct_perc_1
    
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    diff_vec1 <- sapply(1:ncol(common_score), function(k){
      .l2norm(common_score[,k] - res$distinct_score_1[,k])
    })
    diff_vec2 <- sapply(1:ncol(common_score), function(k){
      .l2norm(common_score[,k] - res$distinct_score_2[,k])
    })
    
    sum(abs(diff_vec1 - diff_vec2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})














