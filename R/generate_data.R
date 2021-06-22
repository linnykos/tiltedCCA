# WARNING: currently requires the latent dimension to be same, ie: ncol(score_1) == ncol(score_2)
#' Generate data
#'
#' @param svd_u_1 \code{n} by \code{rank_1} orthogonal matrix
#' @param svd_u_2 \code{n} by \code{rank_2} orthogonal matrix
#' @param svd_d_1 vector of length \code{rank_1}
#' @param svd_d_2 vector of length \code{rank_2}
#' @param svd_v_1 \code{p1} by \code{rank_1} orthogonal matrix
#' @param svd_v_2 \code{p2} by \code{rank_2} orthogonal matrix
#' @param num_neigh number of neighbors to consider to computed the common percentage 
#' @param noise_val numeric for scalar of centered Gaussian noise
#'
#' @return list of outputs
#' @export
generate_data <- function(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                          num_neigh = max(round(nrow(svd_u_1)/20), 40), 
                          noise_val = 1){
  stopifnot(nrow(svd_u_1) == nrow(svd_u_2), ncol(svd_u_1) == ncol(svd_u_2), # simplification for now
            ncol(svd_u_1) == length(svd_d_1), ncol(svd_u_2) == length(svd_d_2), 
            ncol(svd_v_1) == length(svd_d_1), ncol(svd_v_2) == length(svd_d_2), 
            num_neigh <= min(nrow(svd_u_1), nrow(svd_u_2)))

  n <- nrow(svd_u_1); p1 <- nrow(svd_v_1); p2 <- nrow(svd_v_2)
  agg_mat <- crossprod(svd_u_1, svd_u_2)
  cca_svd <- svd(agg_mat)
  score_1 <- svd_u_1 %*% cca_svd$u
  score_2 <- svd_u_2 %*% cca_svd$v
 
  mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
  mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
  nn_1 <- RANN::nn2(.mult_mat_vec(svd_u_1, svd_d_1), k = num_neigh)$nn.idx
  nn_2 <- RANN::nn2(.mult_mat_vec(svd_u_2, svd_d_2), k = num_neigh)$nn.idx
  
  res <- .common_decomposition(score_1, score_2, nn_1, nn_2, fix_distinct_perc = F)
  common_score <- res$common_score; distinct_perc_2 <- res$distinct_perc_2
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  tmp1 <- crossprod(common_score); tmp2 <- crossprod(distinct_score_1); tmp3 <- crossprod(distinct_score_2)
  
  stopifnot(abs(sum(abs(tmp1)) - sum(abs(diag(tmp1)))) <= 1e-6,  
            abs(sum(abs(tmp2)) - sum(abs(diag(tmp2)))) <= 1e-6,
            abs(sum(abs(tmp3)) - sum(abs(diag(tmp3)))) <= 1e-6)
  
  common_mat_1 <- common_score %*% crossprod(score_1, mat_1)
  distinct_mat_1 <- distinct_score_1 %*% crossprod(score_1, mat_1)
  common_mat_2 <- common_score %*% crossprod(score_2, mat_2)
  distinct_mat_2 <- distinct_score_2 %*% crossprod(score_2, mat_2)

  noise_1 <- matrix(stats::rnorm(n*p1, mean = 0, sd = noise_val), nrow = n, ncol = p1)
  noise_1 <- scale(noise_1, center = T, scale = F)
  noise_2 <- matrix(stats::rnorm(n*p2, mean = 0, sd = noise_val), nrow = n, ncol = p2)
  noise_2 <- scale(noise_2, center = T, scale = F)
  mat_1 <- mat_1 + noise_1
  mat_2 <- mat_2 + noise_2
  
  n <- nrow(mat_1)
  rownames(mat_1) <- paste0("n", 1:n); rownames(mat_2) <- paste0("n", 1:n)
  rownames(common_mat_1) <- paste0("n", 1:n); rownames(common_mat_2) <- paste0("n", 1:n)
  rownames(distinct_mat_1) <- paste0("n", 1:n); rownames(distinct_mat_2) <- paste0("n", 1:n)
  rownames(common_score) <- paste0("n", 1:n)
  rownames(distinct_score_1) <- paste0("n", 1:n); rownames(distinct_score_2) <- paste0("n", 1:n)
  
  p_1 <- ncol(mat_1); p_2 <- ncol(mat_2)
  colnames(mat_1) <- paste0("g", 1:p_1)
  colnames(common_mat_1) <- paste0("g", 1:p_1); colnames(distinct_mat_1) <- paste0("g", 1:p_1)
  colnames(mat_2) <- paste0("p", 1:p_2)
  colnames(common_mat_2) <- paste0("p", 1:p_2); colnames(distinct_mat_2) <- paste0("p", 1:p_2)
  
  structure(list(mat_1 = mat_1, mat_2 = mat_2, 
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2,
       distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
       common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2,
       distinct_perc_2 = distinct_perc_2), class = "dcca_data")
}

form_seurat_obj <- function(mat_1, mat_2){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1); p_1 <- ncol(mat_1); p_2 <- ncol(mat_2)

  obj <- Seurat::CreateSeuratObject(counts = t(mat_1), assay = "mode1")
  obj[["mode2"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
  
  obj
}

#' Generate orthogonal matrices correspond to an SBM
#'
#' @param B_mat symmetric matrix
#' @param membership_vec vector of positive integers
#' @param centered boolean, where if \code{TRUE}, one less dimension
#' is provided but the singular vectors are centered
#'
#' @return matrix
#' @export
generate_sbm_orthogonal <- function(B_mat, membership_vec, centered = T){
  stopifnot(ncol(B_mat) == nrow(B_mat), all(B_mat >= 0), all(B_mat <= 1),
            sum(abs(B_mat - t(B_mat))) <= 1e-6,
            all(membership_vec > 0), all(membership_vec %% 1 == 0), 
            max(membership_vec) == length(unique(membership_vec)))
  
  K <- max(membership_vec)
  
  prob_mat <- .compute_prob_mat(B_mat, membership_vec)
  adj_mat <- .generate_adjaceny_mat(prob_mat)
  .svd_truncated(adj_mat, K = ifelse(centered, K-1, K), 
                 symmetric = F, rescale = F, mean_vec = centered, sd_vec = NULL, K_full_rank = F)$u
}

#' Generate orthogonal matrices via Gaussian noise
#'
#' @param n integer
#' @param K integer
#' @param centered boolean, where if \code{TRUE}, one less dimension
#' is provided but the singular vectors are centered
#'
#' @return matrix
#' @export
generate_random_orthogonal <- function(n, K, centered = F){
  stopifnot(K+1 <= n)
  mat <- matrix(stats::rnorm(n^2), n, n)
  mat <- mat + t(mat)
  .svd_truncated(mat, K = K, symmetric = F, rescale = F, mean_vec = centered, 
                 sd_vec = NULL, K_full_rank = F)$u
}

###########################

.generate_membership_matrix <- function(membership_vector){
  K <- max(membership_vector)
  
  sapply(1:K, function(x){
    as.numeric(membership_vector == x)
  })
}

#' Compute probability matrix
#'
#' @param B_mat symmetric connectivity matrix
#' @param membership_vec vector containing values \code{1} through \code{ncol(B_mat)}
#'
#' @return symmetric matrix of dimension \code{length(membership_vec)}
#' @export
.compute_prob_mat <- function(B_mat, membership_vec){
  membership_mat <- .generate_membership_matrix(membership_vec)
  prob_mat <- membership_mat %*% B_mat %*% t(membership_mat)
  
  prob_mat
}

#' Simulate adjacency matrix
#' 
#' The matrix is not symmetric.
#'
#' @param prob_mat probability matrix
#'
#' @return adjacency matrix
#' @export
.generate_adjaceny_mat <- function(prob_mat){
  n <- nrow(prob_mat)
  val <- stats::rbinom(n^2, 1, prob_mat)
  matrix(val, ncol = n, nrow = n)
}
