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