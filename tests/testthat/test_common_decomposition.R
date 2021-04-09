context("Test .common_decomposition")

test_compute_common_score <- function(score_1, score_2, obj_vec = NA){
  stopifnot(nrow(score_1) == nrow(score_2), nrow(score_1) >= ncol(score_1),
            nrow(score_2) >= ncol(score_2), is.matrix(score_1), is.matrix(score_2))
  
  p <- min(c(ncol(score_1), ncol(score_2)))
  if(ncol(score_1) > p){ score_1 <- score_1[,1:p,drop = F]}
  if(ncol(score_2) > p){ score_2 <- score_2[,1:p,drop = F]}
  
  n <- nrow(score_1)
  obj_vec <- diag(crossprod(score_1, score_2))/n
  
  R_vec <- sapply(obj_vec, function(x){x <- min(1,max(x,0)); 1-sqrt((1-x)/(1+x))})
  common_score <- .mult_mat_vec((score_1+score_2)/2, R_vec)
  
  common_score
}


test_that("(Basic) .common_decomposition works", {
  set.seed(5)
  n <- 200; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  svd_1 <- .svd_truncated(mat_1, p1,
                          symmetric = F, rescale = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, p2,
                          symmetric = F, rescale = F, K_full_rank = F) 
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  nn_1 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_1$v), k = 50)$nn.idx
  nn_2 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v), k = 50)$nn.idx
  
  res <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_perc_1"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(length(res$distinct_perc_1) == ncol(res$common_score))
  
  res <- .common_decomposition(score_1, score_2, nn_1 = NA, nn_2 = NA, fix_distinct_perc = T)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_perc_1"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(length(res$distinct_perc_1) == ncol(res$common_score))
})

test_that("(Test) .common_decomposition is correct when fix_distinct_perc = T", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
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
    res <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
    
    svd_1 <- .svd_truncated(res$mat_1, p_1,
                            symmetric = F, rescale = F, K_full_rank = F)
    svd_2 <- .svd_truncated(res$mat_2, p_2,
                            symmetric = F, rescale = F, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    res1 <- .common_decomposition(score_1, score_2, nn_1 = NA, nn_2 = NA, fix_distinct_perc = T)
    res2 <- test_compute_common_score(score_1, score_2, obj_vec = cca_res$obj_vec)
    
    sum(abs(res1$common_score - res2)) <= 1e-6
  })
 
  expect_true(all(bool_vec))
})

test_that("(Coding) .compute_common_score preserves rownames and colnames", {
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
  res <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  
  svd_1 <- .svd_truncated(res$mat_1, p_1,
                          symmetric = F, rescale = F, K_full_rank = F)
  svd_2 <- .svd_truncated(res$mat_2, p_2,
                          symmetric = F, rescale = F, K_full_rank = F) 
  cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  rownames(score_1) <- paste0("a", 1:n)
  rownames(score_2) <- paste0("a", 1:n)
  
  nn_1 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_1$v), k = 50)$nn.idx
  nn_2 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v), k = 50)$nn.idx
  
  res <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)
  expect_true(all(dim(res$common_score) == dim(score_1)))
  expect_true(length(rownames(res$common_score)) > 1)
  expect_true(all(rownames(res$common_score) == rownames(score_1)))
  
  # [[note to self: observe that if rownames(score_1) != rownames(score_2), the names would still be rownames(score_1)]]
})

test_that("(Math) .common_decomposition gives sensible numbers in asymmetric information settings", {
  trials <- 25
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n_clust <- 100
    B_mat1 <- matrix(c(0.9, 0, 0, 
                       0, 0.9, 0,
                       0, 0, 0.9), 3, 3, byrow = T)
    B_mat2 <- matrix(c(0.9, 0.85, 0, 
                       0.85, 0.9, 0,
                       0, 0, 1), 3, 3, byrow = T)
    K <- ncol(B_mat1)
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
    n <- length(membership_vec); true_membership_vec <- membership_vec
    svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
    svd_u_2 <- generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)
    
    set.seed(10)
    p_1 <- 20; p_2 <- 40
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- generate_random_orthogonal(p_1, K-1)
    svd_v_2 <- generate_random_orthogonal(p_2, K-1)
    
    set.seed(10)
    dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, 
                         svd_v_1, svd_v_2, noise_val = 0.1)
    
    svd_1 <- .svd_truncated(dat$mat_1, 2,
                            symmetric = F, rescale = F, K_full_rank = F)
    svd_2 <- .svd_truncated(dat$mat_2, 2,
                            symmetric = F, rescale = F, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, rank_1 = NA, rank_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    num_neigh <- 40
    nn_1 <- RANN::nn2(.mult_mat_vec(svd_1$u, svd_1$d), k = num_neigh)$nn.idx
    nn_2 <- RANN::nn2(.mult_mat_vec(svd_2$u, svd_2$d), k = num_neigh)$nn.idx
    
    res <- .common_decomposition(score_1, score_2, nn_1 = nn_1, nn_2 = nn_2, 
                                 fix_distinct_perc = F)
    
    all(res$distinct_perc_1 <= 0.5)
  })
  
  expect_true(all(bool_vec))
})

###############################

## .sigmoid_ratio is correct

test_that("(Math) .sigmoid_ratio gives an appropriate value w.r.t 0.5 ", {
  trials <- 50
  
  a_vec <- seq(0, 10, length.out = 11)
  bool_vec <- sapply(a_vec, function(a){
    all(sapply(1:trials, function(x){
      b <- runif(1, min = 0, max = 2)*abs(a+rnorm(1, sd = 0.01))
      val <- .sigmoid_ratio(a, b)
      if(a < b) val < 0.5 else val > 0.5
    }))
  })
  
  expect_true(all(bool_vec))
  
  b_vec <- seq(0, 10, length.out = 11)
  bool_vec <- sapply(b_vec, function(b){
    all(sapply(1:trials, function(x){
      a <- runif(1, min = 0, max = 2)*abs(b+rnorm(1, sd = 0.01))
      val <- .sigmoid_ratio(a, b)
      if(a < b) val < 0.5 else val > 0.5
    }))
  })
  
  expect_true(all(bool_vec))
})
