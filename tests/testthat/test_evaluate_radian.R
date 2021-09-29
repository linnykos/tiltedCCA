context("Test .evaluate_radian")

compute_evaluate_radian_ingredients <- function(){
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
  
  svd_1 <- .svd_truncated(mat_1, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  svd_2 <- .svd_truncated(mat_2, K = 2, symmetric = F, rescale = F, 
                          mean_vec = T, sd_vec = F, K_full_rank = F)
  
  svd_1 <- .check_svd(svd_1, dims = c(1:2))
  svd_2 <- .check_svd(svd_2, dims = c(1:2))
  
  num_neigh <- 30
  frnn_union <- .compute_frnn_union(svd_1,
                                    svd_2,
                                    num_neigh,
                                    radius_quantile = 0.5)
  
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
  
  list(basis_list = basis_list, 
       frnn_union = frnn_union,
       num_neigh = num_neigh,
       score_1 = score_1,
       score_2 = score_2,
       svd_1 = svd_1, 
       svd_2 = svd_2,
       circle_list = circle_list)
}

## .evaluate_radian is correct

test_that(".evaluate_radian works", {
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients()
  basis_list <- tmp$basis_list; frnn_union <- tmp$frnn_union
  num_neigh <- tmp$num_neigh; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  
  res <- .evaluate_radian(percentage_val = 0.5,
                          basis_list = basis_list, 
                          frnn_union = frnn_union,
                          num_neigh = num_neigh,
                          score_1 = score_1,
                          score_2 = score_2,
                          svd_1 = svd_1, 
                          svd_2 = svd_2,
                          circle_list = circle_list,
                          return_common_score = F,
                          return_snn = F)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})

test_that(".evaluate_radian is maximized reasonably at 0", {
  set.seed(10)
  tmp <- compute_evaluate_radian_ingredients()
  basis_list <- tmp$basis_list; frnn_union <- tmp$frnn_union
  num_neigh <- tmp$num_neigh; score_1 <- tmp$score_1
  score_2 <- tmp$score_2; svd_1 <- tmp$svd_1
  svd_2 <- tmp$svd_2; circle_list <- tmp$circle_list
  
  res <- sapply(seq(0, 1, length.out = 9), function(percentage_val){
    .evaluate_radian(percentage_val = percentage_val,
                     basis_list = basis_list, 
                     frnn_union = frnn_union,
                     num_neigh = num_neigh,
                     score_1 = score_1,
                     score_2 = score_2,
                     svd_1 = svd_1, 
                     svd_2 = svd_2,
                     radius_quantile = 0.5,
                     circle_list = circle_list,
                     return_common_score = F,
                     return_frnn = F)
  })
  
  expect_true(length(res) == 9)
  expect_true(all(res >= 0))
})