context("Test .compute_distinct_score")

compute_compute_distinct_score_ingredients <- function(setting = 1){
  # setting 1 has modality 2 having no distinct information
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
  
  cca_res <- .cca(svd_1, svd_2, 
                  dims_1 = NA, dims_2 = NA, 
                  return_scores = F)
  
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  tmp <- .form_snns(num_neigh = 30, svd_1 = svd_1, svd_2 = svd_2)
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
  
  res <- .common_decomposition(discretization_gridsize = 9,
                               enforce_boundary = T,
                               fix_tilt_perc = T,
                               metacell_clustering_1 = metacell_clustering_1,
                               metacell_clustering_2 = metacell_clustering_2,
                               n_idx = 1:nrow(score_1),
                               num_neigh = 30,
                               score_1 = score_1,
                               score_2 = score_2,
                               svd_1 = svd_1, 
                               svd_2 = svd_2)
  
  list(score_1 = score_1,
       score_2 = score_2,
       common_score = res$common_score)
}

## .compute_distinct_score is correct

test_that("(Basic) .compute_distinct_score works", {
  set.seed(10)
  tmp <- compute_compute_distinct_score_ingredients()
  score_1 <- tmp$score_1
  score_2 <- tmp$score_2
  common_score <- tmp$common_score
  
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
    tmp <- compute_compute_distinct_score_ingredients()
    score_1 <- tmp$score_1
    score_2 <- tmp$score_2
    common_score <- tmp$common_score
    
    res <- .compute_distinct_score(score_1, score_2, common_score)
    
    all(abs(crossprod(res$distinct_score_1, res$distinct_score_2)) <= 1e-6)
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) .compute_distinct_score generates equal-lengthed matrices when fix_distinct_perc = T", {
  trials <- 20
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    tmp <- compute_compute_distinct_score_ingredients()
    score_1 <- tmp$score_1
    score_2 <- tmp$score_2
    common_score <- tmp$common_score
    
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

