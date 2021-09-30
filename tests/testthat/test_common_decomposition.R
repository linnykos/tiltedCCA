context("Test .common_decomposition")

compute_common_decomposition_ingredients <- function(setting = 1){
  # setting 1 has modality 2 having no distinct information
  if(setting == 1){
    n_clust <- 100
    high <- 0.9; low <- 0.05
    B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                       0.1, 0.9, 0.1,
                       0.1, 0.1, 0.9), 3, 3, byrow = T)
    K <- ncol(B_mat1)
    membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
    n <- length(membership_vec); true_membership_vec <- membership_vec
    svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
    svd_u_2 <- multiomicCCA::generate_random_orthogonal(n, 2, centered = T)
    
    p_1 <- 20; p_2 <- 40
    svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
  } else if(setting == 2){
    # setting 2 is where both modalities are the same
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:3, function(i){
      if(i == 1){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else if(i == 2){
        MASS::mvrnorm(n = n_each, mu = c(0,12), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:3, function(i){
      if(i == 1){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else if(i == 2){
        MASS::mvrnorm(n = n_each, mu = c(0,12), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else if(setting == 3){
    # setting 3 is where modality 2 has information, but not as much
    n_each <- 100
    true_membership_vec <- rep(1:3, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:3, function(i){
      if(i %in% c(1,2)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(9,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:3, function(i){
      if(i %in% c(1,3)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(3,0), Sigma = diag(2)) 
      }
    }))
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  } else {
    # setting 4 is two modalities with high distinct information
    n_each <- 100
    true_membership_vec <- rep(1:4, each = n_each)
    mat_1 <- do.call(rbind, lapply(1:4, function(i){
      if(i %in% c(1,2)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_2 <- do.call(rbind, lapply(1:4, function(i){
      if(i %in% c(1,3)){
        MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
      } else {
        MASS::mvrnorm(n = n_each, mu = c(12,0), Sigma = diag(2)) 
      }
    }))
    
    mat_1 <- scale(mat_1, center = T, scale = F)
    mat_2 <- scale(mat_2, center = T, scale = F)
    svd_1 <- svd(mat_1)
    svd_2 <- svd(mat_2)
    
    p_1 <- 40; p_2 <- 40
    svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
    svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
    
    mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
    mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  }
  
  svd_1 <- .svd_truncated(mat_1, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = c(1:2))
  svd_2 <- .check_svd(svd_2, dims = c(1:2))
  
  cca_res <- .cca(svd_1, svd_2, 
                  dims_1 = NA, dims_2 = NA, 
                  return_scores = F)
  
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  list(score_1 = score_1,
       score_2 = score_2,
       svd_1 = svd_1, 
       svd_2 = svd_2)
}

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
  
  res <- suppressWarnings(.common_decomposition(discretization_gridsize = 9,
                                                fix_tilt_perc = F,
                                                score_1 = score_1,
                                                score_2 = score_2,
                                                svd_1 = svd_1, 
                                                svd_2 = svd_2,
                                                trials = 100))
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", "tilt_perc", "df_percentage"))))
  expect_true(nrow(res$common_score) == nrow(score_1))
  expect_true(ncol(res$common_score) == min(ncol(score_1), ncol(score_2)))
  expect_true(length(res$tilt_perc) == 1)
})

test_that("(Test) .common_decomposition is correct when fix_tilt_perc = T", {
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
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1) 
    mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1))), ncol = ncol(mat_1), nrow = nrow(mat_1))
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
    mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2))), ncol = ncol(mat_2), nrow = nrow(mat_2))
    
    svd_1 <- .svd_truncated(mat_1, p_1, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    svd_2 <- .svd_truncated(mat_2, p_2, symmetric = F, rescale = F,
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    res1 <- .common_decomposition(discretization_gridsize = NA,
                                  fix_tilt_perc = T,
                                  score_1 = score_1,
                                  score_2 = score_2,
                                  svd_1 = svd_1, 
                                  svd_2 = svd_2)
    res2 <- test_compute_common_score(score_1, score_2, obj_vec = cca_res$obj_vec)
    
    sum(abs(res1$common_score - res2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Coding) .common_decomposition preserves rownames and colnames", {
  set.seed(10)
  tmp <- compute_common_decomposition_ingredients(setting = 1)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  svd_1 <- tmp$svd_1; svd_2 <- tmp$svd_2
  rownames(score_1) <- paste0("a", 1:n)
  rownames(score_2) <- paste0("a", 1:n)
  
  res <- .common_decomposition(discretization_gridsize = 9,
                               fix_tilt_perc = F,
                               score_1 = score_1,
                               score_2 = score_2,
                               svd_1 = svd_1, 
                               svd_2 = svd_2,
                               trials = 100)
  expect_true(all(dim(res$common_score) == dim(score_1)))
  expect_true(length(rownames(res$common_score)) > 1)
  expect_true(all(rownames(res$common_score) == rownames(score_1)))
  
  # [[note to self: observe that if rownames(score_1) != rownames(score_2), the names would still be rownames(score_1)]]
})

test_that("(Math) .common_decomposition gives sensible numbers in asymmetric information settings", {
  trials <- 10
  
  bool_vec <- sapply(1:trials, function(x){
    print(x)
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
    mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1) 
    mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1)), sd = 0.1), ncol = ncol(mat_1), nrow = nrow(mat_1))
    mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
    mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2)), sd = 0.1), ncol = ncol(mat_2), nrow = nrow(mat_2))
    
    svd_1 <- .svd_truncated(mat_1, 2, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    svd_2 <- .svd_truncated(mat_2, 2, symmetric = F, rescale = F, 
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F) 
    cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)
    tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    
    res <- .common_decomposition(discretization_gridsize = 9,
                                 fix_tilt_perc = F,
                                 score_1 = score_1,
                                 score_2 = score_2,
                                 svd_1 = svd_1, 
                                 svd_2 = svd_2,
                                 trials = 50)
    
    res$tilt_perc <= 0.5
  })
  
  expect_true(all(bool_vec))
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

##################

## .select_minimum is correct

test_that(".select_minimum works", {
  x_val <- c(0, 0.25, 0.5, 0.75, 1)
  y_val <- 1:5
  res <- .select_minimum(x_val, y_val)
  
  expect_true(res == 1)
})

test_that(".select_minimum works when there are multiple similar values", {
  x_val <- c(0, 0.25, 0.5, 0.75, 1)
  y_val <- c(3,3,1,1,1)
  res <- .select_minimum(x_val, y_val)
  expect_true(res == 5)
  
  y_val <- c(3,1,1,1,1)
  res <- .select_minimum(x_val, y_val)
  expect_true(res == 5)
  
  y_val <- c(3,2,1,1,2)
  res <- .select_minimum(x_val, y_val)
  expect_true(res == 4)
  
  expect_warning(.select_minimum(x_val, c(3,1,1,1,3)))
})
