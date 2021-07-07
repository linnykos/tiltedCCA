context("Test plot embedding 2")

test_that("plot_embeddings2 works", {
  set.seed(10)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.2, 0.1, 
                    0.2, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3, byrow = T)
  K <- ncol(B_mat); rho <- 1
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  
  set.seed(10)
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- generate_random_orthogonal(p_2, K-1)
  
  set.seed(10)
  dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  dcca_obj <- dcca_factor(dat$mat_1, dat$mat_2, dims_1 = 1:(K-1), dims_2 = 1:(K-1), 
                          verbose = F)
  membership_vec <- as.factor(membership_vec)
  list_g <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
                        verbose = F, bool_matrix = T)
  
  set.seed(10)
  res <- plot_embeddings2(list_g[[1]], list_g[[2]], list_g[[3]],
                          only_embedding = T)
  
})