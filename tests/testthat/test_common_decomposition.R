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

## .common_decomposition is correct

test_that("(Basic) .common_decomposition works", {
  set.seed(5)
  n <- 200; K <- 2
  common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)
  
  p1 <- 5; p2 <- 10
  transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
  transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)
  
  mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
  mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
  
  svd_1 <- .svd_truncated(mat_1, p1, symmetric = F, rescale = F,
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, p2, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
  cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  nn_1 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_1$v), k = 50)$nn.idx
  nn_2 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v), k = 50)$nn.idx
  
  res <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_perc_2"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(length(res$distinct_perc_2) == ncol(res$common_score))
  
  res <- .common_decomposition(score_1, score_2, nn_1 = NA, nn_2 = NA, fix_distinct_perc = T)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "distinct_perc_2"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(length(res$distinct_perc_2) == ncol(res$common_score))
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
    
    svd_1 <- .svd_truncated(res$mat_1, p_1, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    svd_2 <- .svd_truncated(res$mat_2, p_2, symmetric = F, rescale = F,
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    res1 <- .common_decomposition(score_1, score_2, nn_1 = NA, nn_2 = NA, fix_distinct_perc = T)
    res2 <- test_compute_common_score(score_1, score_2, obj_vec = cca_res$obj_vec)
    
    sum(abs(res1$common_score - res2)) <= 1e-6
  })
 
  expect_true(all(bool_vec))
})

test_that("(Coding) .common_decomposition preserves rownames and colnames", {
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
  
  svd_1 <- .svd_truncated(res$mat_1, p_1, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
  svd_2 <- .svd_truncated(res$mat_2, p_2, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
  cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
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
    
    svd_1 <- .svd_truncated(dat$mat_1, 2, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    svd_2 <- .svd_truncated(dat$mat_2, 2, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    num_neigh <- 40
    nn_1 <- RANN::nn2(.mult_mat_vec(svd_1$u, svd_1$d), k = num_neigh)$nn.idx
    nn_2 <- RANN::nn2(.mult_mat_vec(svd_2$u, svd_2$d), k = num_neigh)$nn.idx
    
    res <- .common_decomposition(score_1, score_2, nn_1 = nn_1, nn_2 = nn_2, 
                                 fix_distinct_perc = F)
    
    bool1 <- all(res$distinct_perc_2 <= 0.5)
    
    distinct_1 <- score_1 - res$common_score
    distinct_2 <- score_2 - res$common_score
    
    bool2_vec <- rep(NA, 2)
    for(i in 1:2){
      d1_norm <- .l2norm(distinct_1[,i]); d2_norm <- .l2norm(distinct_2[,i])
      bool2_vec[i] <- abs(d2_norm/(d1_norm + d2_norm) - res$distinct_perc_2[i]) <= 1e-2
    }
    
    bool1 & all(bool2_vec)
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

#####################

## .sigmoid_ratio is correct

test_that(".sigmoid_ratio lies between 0 and 1", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    a <- abs(rnorm(1)); b <- abs(rnorm(1))
    res <- .sigmoid_ratio(a,b)
    res >= 0 & res <= 1
  })
  
  expect_true(all(bool_vec))
})

test_that(".sigmoid_ratio is correctly above/below 0.5", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    a <- abs(rnorm(1)); b <- abs(rnorm(1))
    res <- .sigmoid_ratio(a,b)
    
    if(a < b) return(res < 0.5)
    if(a > b) return(res > 0.5)
  })
  
  expect_true(all(bool_vec))
})

test_that(".sigmoid_ratio is symmetric", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    a <- abs(rnorm(1)); b <- abs(rnorm(1))
    res1 <- .sigmoid_ratio(a,b)
    res2 <- .sigmoid_ratio(b,a)
    
    abs(res1 - (1-res2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that(".sigmoid_ratio appropriately scales", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    
    a <- abs(rnorm(1)); b <- abs(rnorm(1))
    res1 <- .sigmoid_ratio(a,b)
    if(a < b){
      res2 <- .sigmoid_ratio(a/2,b)
      res3 <- .sigmoid_ratio(a,2*b)
      
      bool2 <- res2 <= res1
      bool3 <- res3 <= res1
    } else {
      res2 <- .sigmoid_ratio(a,b/2)
      res3 <- .sigmoid_ratio(2*a,b)
      
      bool2 <- res2 >= res1
      bool3 <- res3 >= res1
    }
    
    res4 <- .sigmoid_ratio(2*a,2*b)
    bool1 <- abs(res1 - res4) <= 1e-6
  
    bool1 & bool2 & bool3
  })
  
  expect_true(all(bool_vec))
})

#######################

## .latent_distinct_perc_2 is correct

test_that(".latent_distinct_perc_2 works", {
  set.seed(10)
  n <- 100; k <- 5
  score_vec_1 <- rnorm(n); score_vec_1 <- score_vec_1/.l2norm(score_vec_1)
  score_vec_2 <- rnorm(n); score_vec_2 <- score_vec_2/.l2norm(score_vec_2)
  nn_1 <- matrix(sample(1:n, size = n*k, replace = T), n, k)
  nn_2 <- matrix(sample(1:n, size = n*k, replace = T), n, k)
  
  res <- .latent_distinct_perc_2(score_vec_1, score_vec_2, nn_1, nn_2)
  
  expect_true(is.numeric(res))
})

test_that(".latent_distinct_perc_2 roughly has the correct magnitude", {
  trials <- 25
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    n <- 150; k <- 5
    score_vec_1 <- rnorm(n)
    score_vec_1 <- score_vec_1-mean(score_vec_1); score_vec_1 <- score_vec_1/.l2norm(score_vec_1)
    score_vec_2 <- rnorm(n)
    score_vec_2 <- score_vec_2-mean(score_vec_2); score_vec_2 <- score_vec_2/.l2norm(score_vec_2)
    nn_1 <- matrix(sample(1:n, size = n*k, replace = T), n, k)
    nn_2 <- matrix(sample(1:n, size = n*k, replace = T), n, k)
    res1 <- .latent_distinct_perc_2(score_vec_1, score_vec_2, nn_1, nn_2)
    
    score_vec_3 <- c(rnorm(50, mean = -5), rnorm(50, mean = 0), rnorm(50, mean = 5))
    score_vec_3 <- score_vec_3-mean(score_vec_3); score_vec_3 <- score_vec_3/.l2norm(score_vec_3)
    nn_3 <- do.call(rbind, lapply(1:3, function(x){
      matrix(sample(1:50, size = 50*k, replace = T)+(x-1)*50, 50, k)
    }))
    res2 <- .latent_distinct_perc_2(score_vec_1, score_vec_3, nn_1, nn_3)
    res3 <- .latent_distinct_perc_2(score_vec_3, score_vec_1, nn_3, nn_1)
    
    res1 < res2 & res1 > res3
  })
  
  expect_true(all(bool_vec))
})

################################

## .binary_search_radian is correct

test_that(".binary_search_radian works", {
  circle <- .construct_circle(c(0,1), c(1,0))
  res <- .binary_search_radian(circle, left_radian = 0, right_radian = pi/2, distinct_perc_2 = 0.5)
  expect_true(is.numeric(res))
})

test_that(".binary_search_radian gives the correct value", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- abs(rnorm(2)); vec1 <- vec1/.l2norm(vec1)
    vec2 <- abs(rnorm(2)); vec2 <- vec2/.l2norm(vec2)
    circle <- .construct_circle(vec1, vec2)
    tmp <- .rightmost_vector(vec1, vec2)
    vec1 <- tmp$vec_right; vec2 <- tmp$vec_left
    right_radian <- .find_radian(circle, tmp$vec_right)
    left_radian <- .find_radian(circle, tmp$vec_left)
    stopifnot(right_radian < left_radian)
    
    distinct_perc_2 <- runif(1)
    right_radian <- right_radian + 2*pi
    res <- .binary_search_radian(circle, left_radian, right_radian, 
                                 distinct_perc_2, max_iter = 50, tol = 1e-6)
    common_vec <- .position_from_circle(circle, res)
    
    distinct1 <- vec1 - common_vec
    distinct2 <- vec2 - common_vec
    ratio <- .l2norm(distinct2)/(.l2norm(distinct1) + .l2norm(distinct2))
    
    abs(distinct_perc_2 - ratio) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that(".binary_search_radian gives the correct value after basis transformation", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- abs(rnorm(2)); vec1 <- vec1/.l2norm(vec1)
    vec2 <- abs(rnorm(2)); vec2 <- vec2/.l2norm(vec2)
    n <- 100
    basis_mat <- svd(matrix(rnorm(2*n), n, 2))$u[,1:2]
    vec1 <- as.numeric(basis_mat%*%vec1); vec2 <- as.numeric(basis_mat%*%vec2)
    basis_res <- .representation_2d(vec1, vec2)
    
    circle <- .construct_circle(basis_res$rep1, basis_res$rep2)
    tmp <- .rightmost_vector(basis_res$rep1, basis_res$rep2)
    right_radian <- .find_radian(circle, tmp$vec_right)
    left_radian <- .find_radian(circle, tmp$vec_left)
    stopifnot(right_radian < left_radian)
    
    distinct_perc_2 <- runif(1)
    right_radian <- right_radian + 2*pi
    res <- .binary_search_radian(circle, left_radian, right_radian, 
                                 distinct_perc_2, max_iter = 50, tol = 1e-6)
    common_vec <- basis_res$basis_mat %*% .position_from_circle(circle, res)
    
    distinct1 <- vec1 - common_vec
    distinct2 <- vec2 - common_vec
    ratio <- .l2norm(distinct2)/(.l2norm(distinct1) + .l2norm(distinct2))
    
    abs(distinct_perc_2 - ratio) <= 1e-3
  })
  
  expect_true(all(bool_vec))
})


