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

########################

## .grab_previous_values is correct

test_that(".grab_previous_values works", {
  percentage_grid <- c(0, 1/4, 1/2)
  percentage_grid_all <- c(0, 1/2, 1)
  value_vec_all <- c(1, 2, 3)
  
  res <- .grab_previous_values(percentage_grid,
                               percentage_grid_all,
                               value_vec_all)
  
  expect_true(all(res[c(1,3)] == c(1, 2)))
  expect_true(is.na(res[2]))
})

###############

## .update_values is correct

test_that(".update_values works", {
  percentage_grid <- c(0, 1/4, 1/2)
  percentage_grid_all <- c(0, 1/2, 1)
  value_vec <- c(1, 5, 2)
  value_vec_all <- c(1, 2, 3)
  
  res <- .update_values(percentage_grid, percentage_grid_all,
                        value_vec, value_vec_all)
  
  expect_true(all(sort(names(res)) == sort(c("percentage_grid_all",
                                             "value_vec_all"))))
  expect_true(all(res$percentage_grid_all == c(0, 1/4, 1/2, 1)))            
  expect_true(all(res$value_vec_all == c(1, 5, 2, 3)))            
})

########################

## .compute_radian is correct

test_that(".compute_radian works", {
  vec1 <- c(1, 0)
  vec2 <- c(1/sqrt(2), 1/sqrt(2))
  circle <- .construct_circle(vec1, vec2)
  
  res <- .compute_radian(percentage_val = 0,
                         vec1, vec2,
                         circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  expect_true(sum(abs(vec2 - c(x_coord, y_coord))) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 1,
                         vec1, vec2,
                         circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  expect_true(sum(abs(vec1 - c(x_coord, y_coord))) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 0.5,
                         vec1, vec2,
                         circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  vec <- c(x_coord, y_coord)
  expect_true(abs(.l2norm(vec1 - vec) - .l2norm(vec2 - vec)) <= 1e-6)
  
  res <- .compute_radian(percentage_val = 0.3,
                         vec1, vec2,
                         circle)
  x_coord <- circle$center[1] + circle$radius*cos(res)
  y_coord <- circle$center[2] + circle$radius*sin(res)
  vec <- c(x_coord, y_coord)
  expect_true(.l2norm(vec1 - vec) > .l2norm(vec2 - vec))
})

#####################

## .position_from_circle is correct

test_that(".position_from_circle works", {
  circle <- list(center = c(0,0), radius = 2)
  res <- .position_from_circle(circle, radian = 0)
  expect_true(sum(abs(res - c(2,0))) <= 1e-6)
  
  res <- .position_from_circle(circle, radian = pi/2)
  expect_true(sum(abs(res - c(0,2))) <= 1e-6)
  
  res <- .position_from_circle(circle, radian = pi/4)
  expect_true(sum(abs(res - rep(sqrt(2),2))) <= 1e-6)
})
