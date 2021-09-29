context("Test .convert_common_score_to_mat")

compute_common_score_ingredients <- function(){
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
  
  n <- nrow(mat_1)
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  svd_1 <- svd(mat_1)
  svd_2 <- svd(mat_2)
  
  p_1 <- 40; p_2 <- 40
  svd_v_1 <- generate_random_orthogonal(p_1, 2)
  svd_v_2 <- generate_random_orthogonal(p_2, 2)
  
  mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
  mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
  
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
  
  r <- min(ncol(score_1), ncol(score_2))
  basis_list <- lapply(1:r, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  circle_list <- lapply(1:r, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  radian_vec <- sapply(1:r, function(k){
    .compute_radian(percentage_val = 0.5, 
                    vec1 = basis_list[[k]]$rep1,
                    vec2 = basis_list[[k]]$rep2,
                    circle = circle_list[[k]])
  })
  
  
  common_representation <- sapply(1:r, function(k){
    .position_from_circle(circle_list[[k]], radian_vec[k])
  })
  
  common_score <- sapply(1:r, function(k){
    basis_list[[k]]$basis_mat %*% common_representation[,k]
  })
  
  list(common_score = common_score,
       score_1 = score_1,
       score_2 = score_2,
       svd_1 = svd_1, 
       svd_2 = svd_2)
}

## .convert_common_score_to_mat is correct
test_that(".convert_common_score_to_mat returns reasonable values", {
  set.seed(10)
  tmp <- compute_common_score_ingredients()
  common_score <- common_score
  score_1 <- score_1
  score_2 <- score_2
  svd_1 <- svd_1
  svd_2 <- svd_2
  n <- nrow(score_1)
  
  res <- .convert_common_score_to_mat(common_score,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  expect_true(all(dim(res) == c(nrow(score_1), 2*ncol(score_1))))
  
  res <- .convert_common_score_to_mat(score_1,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  rescaling_factor <- max(c(svd_1$d, svd_2$d))
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d)
  dimred_1 <- dimred_1/svd_1$d[1]*rescaling_factor
  target <- tcrossprod(score_1) %*% dimred_1/n
  expect_true(sum(abs(res[,1:2] - target)) <= 1e-6)
  expect_true(abs(max(svd(dimred_1)$d) - max(svd(target)$d)) <= 1e-6)
  
  res <- .convert_common_score_to_mat(score_2,
                                      score_1,
                                      score_2,
                                      svd_1, 
                                      svd_2)
  rescaling_factor <- max(c(svd_1$d, svd_2$d))
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d)
  dimred_2 <- dimred_2/svd_2$d[1]*rescaling_factor
  target <- tcrossprod(score_2) %*% dimred_2/n
  expect_true(sum(abs(res[,3:4] - target)) <= 1e-6)
  expect_true(abs(max(svd(dimred_2)$d) - max(svd(target)$d)) <= 1e-6)
})

#########################

## .computer_overlap is correct
test_that(".computer_overlap works", {
  set.seed(10)
  n <- 1000
  mat <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(2))
  snn_target <- .form_snn(mat, num_neigh = 30)
  
  mat <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(2))
  snn_query <- .form_snn(mat, num_neigh = 30)
  
  res <- .computer_overlap(snn_target, snn_query)
  expect_true(res >= 0)
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})




