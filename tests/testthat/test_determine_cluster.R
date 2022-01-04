context("Test determine cluster")

# .grassmann_distance is correct

test_that(".grassmann_distance works", {
  set.seed(10)
  mat_1 <- matrix(rnorm(75), ncol = 3, nrow = 25)
  mat_2 <- matrix(rnorm(75), ncol = 3, nrow = 25)
  orthonormal_1 <- svd(mat_1)$u[,1:3]
  orthonormal_2 <- svd(mat_2)$u[,1:3]
  
  res <- .grassmann_distance(orthonormal_1 = orthonormal_1,
                             orthonormal_2 = orthonormal_2)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".grassmann_distance gives reasonable distances", {
  trials <- 25
  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                       0.1, 0.9, 0.1,
                       0.1, 0.1, 0.9), ncol = 3, nrow = 3)
    B_mat2 <- matrix(runif(9), ncol = 3, nrow = 3)
    B_mat2 <- (B_mat2 + t(B_mat2))/2
    membership_vec <- rep(1:3, each = 25)
    orthonormal_1a <- generate_sbm_orthogonal(B_mat = B_mat1,
                                              membership_vec = membership_vec)
    orthonormal_1b <- generate_sbm_orthogonal(B_mat = B_mat1,
                                              membership_vec = membership_vec)
    orthonormal_2 <- generate_sbm_orthogonal(B_mat = B_mat2,
                                             membership_vec = membership_vec)
    
    res1 <- .grassmann_distance(orthonormal_1 = orthonormal_1a,
                                orthonormal_2 = orthonormal_1b)
    res2 <- .grassmann_distance(orthonormal_1 = orthonormal_1a,
                                orthonormal_2 = orthonormal_2)
    
    res1 <= res2
  })
  
  expect_true(all(bool_vec))
})

#####################

# .determine_cluster is correct

test_that(".determine_cluster works", {
  
})

##################

# .kl_divergence is correct

test_that(".kl_divergence works", {
  query_dist <- c(.5, .4, .1)
  reference_dist <- c(.3, .3, .4)
  val1 <- .kl_divergence(query_dist, reference_dist)
  
  query_dist2 <- c(.35, .25, .4)
  val2 <- .kl_divergence(query_dist2, reference_dist)
  
  expect_true(is.numeric(val1))
  expect_true(is.numeric(val2))
  expect_true(val1 >= val2)
})

test_that(".kl_divergence works with values of 0", {
  query_dist <- c(.25, .25, 0, .25, .1, 0, .15)
  reference_dist <- c(.15, .15, .15, .15, .4, 0, 0)
  val1 <- .kl_divergence(query_dist, reference_dist)
  
  expect_true(is.numeric(val1))
})

############################

## .compute_factor_diversity is correct

test_that(".compute_factor_diversity works", {
  n <- 100
  set.seed(10)
  dat1 <- c(stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3))
  dat2 <- stats::rnorm(2*n, mean = 3)
  dat3 <- c(stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3))
  dat4 <- stats::rnorm(2*n, mean = 3)
  
  factor_anchor <- factor(c(rep("A", n), rep("B", n)))
  factor_other <- factor(rep("AB", 2*n))
  
  
  list_nn <- RANN::nn2(dat3, k = 11)$nn.idx
  list_nn <- lapply(1:nrow(list_nn), function(i){
    list_nn[i,-1]
  })
  val1 <- .compute_factor_diversity(factor_anchor = factor_anchor,
                                    factor_other = factor_other,
                                    list_nn = list_nn)
  expect_true(abs(val1) <= 1e-3)
  val2 <- .compute_factor_diversity(factor_anchor = factor_other,
                                    factor_other = factor_anchor,
                                    list_nn = list_nn)
  expect_true(is.numeric(val2))
  
  
  list_nn <- RANN::nn2(dat4, k = 11)$nn.idx
  list_nn <- lapply(1:nrow(list_nn), function(i){
    list_nn[i,-1]
  })
  val3 <- .compute_factor_diversity(factor_anchor = factor_anchor,
                                    factor_other = factor_other,
                                    list_nn = list_nn)
  expect_true(abs(val3) <= 1e-3)
  val4 <- .compute_factor_diversity(factor_anchor = factor_other,
                                    factor_other = factor_anchor,
                                    list_nn = list_nn)
  expect_true(is.numeric(val4))
  
  expect_true(val2 < val4)
})

test_that(".compute_factor_diversity works with NAs in the factors", {
  n <- 100
  set.seed(10)
  dat1 <- c(stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3))
  dat2 <- stats::rnorm(2*n, mean = 3)
  dat3 <- c(stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3))
  
  factor_anchor <- factor(c(rep("A", n), rep("B", n)))
  factor_anchor[sample(2*n, 40)] <- NA
  factor_other <- factor(rep("AB", 2*n))
  factor_other[sample(2*n, 40)] <- NA
  
  list_nn <- RANN::nn2(dat3, k = 11)$nn.idx
  list_nn <- lapply(1:nrow(list_nn), function(i){
    list_nn[i,-1]
  })
  val1 <- .compute_factor_diversity(factor_anchor = factor_anchor,
                                    factor_other = factor_other,
                                    list_nn = list_nn)
  expect_true(abs(val1) <= 1e-3)
  val2 <- .compute_factor_diversity(factor_anchor = factor_other,
                                    factor_other = factor_anchor,
                                    list_nn = list_nn)
  expect_true(val2 < 0)
})

test_that(".compute_factor_diversity works for a more interesting example", {
  n <- 100
  set.seed(10)
  dat1 <- c(stats::rnorm(2*n, mean = 3), stats::rnorm(2*n, mean = -3))
  dat2 <- c(stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3),
            stats::rnorm(n, mean = 3), stats::rnorm(n, mean = -3))
  
  factor_anchor <- factor(c(rep("A", 2*n), rep("B", 2*n)))
  factor_other <- factor(c(rep("C", n), rep("D", n), rep("C", n), rep("D", n)))
  
  candidate_list <- list(c(stats::rnorm(2*n, mean = 3), stats::rnorm(2*n, mean = -3)),
                         stats::rnorm(4*n),
                         c(stats::rnorm(n, mean = -9), stats::rnorm(n, mean = -3),
                           stats::rnorm(n, mean = 3), stats::rnorm(n, mean = 9)))
  
  val_vec <- sapply(1:length(candidate_list), function(j){
    list_nn <- RANN::nn2(candidate_list[[j]], k = 11)$nn.idx
    list_nn <- lapply(1:nrow(list_nn), function(i){
      list_nn[i,-1]
    })
    min(.compute_factor_diversity(factor_anchor = factor_anchor,
                                  factor_other = factor_other,
                                  list_nn = list_nn),
        .compute_factor_diversity(factor_anchor = factor_other,
                                  factor_other = factor_anchor,
                                  list_nn = list_nn))
  })
  
  expect_true(all(val_vec <= 0))
  expect_true(which.max(val_vec) == 2)
})
