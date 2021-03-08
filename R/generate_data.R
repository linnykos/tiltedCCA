# WARNING: currently requires the latent dimension to be same, ie: ncol(score_1) == ncol(score_2)
#' Generate data
#'
#' @param score_1 \code{n} by \code{rank_1} orthogonal matrix
#' @param score_2 \code{n} by \code{rank_2} orthogonal matrix
#' @param coef_mat_1 \code{rank_1} by \code{p_1} matrix
#' @param coef_mat_2 \code{rank_2} by \code{p_2} matrix
#' @param num_neigh number of neighbors to consider to computed the common percentage 
#' @param noise_func function that takes in a matrix and outputs a matrix of the same dimension
#'
#' @return list of two matrices, one of dimension \code{n} by \code{p_1} and another of dimension \code{n} by \code{p_2}
#' @export
generate_data <- function(score_1, score_2, coef_mat_1, coef_mat_2, 
                               num_neigh = max(round(nrow(score_1)/20), 40),
                               noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}){
  stopifnot(nrow(score_1) == nrow(score_2), nrow(score_1) > ncol(score_1),
            nrow(score_2) > ncol(score_2), ncol(score_1) == nrow(coef_mat_1),
            ncol(score_2) == nrow(coef_mat_2), num_neigh <= min(nrow(score_1), nrow(score_2)))

  
  tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  
  mat_1 <- score_1 %*% coef_mat_1; mat_2 <- score_2 %*% coef_mat_2
  nn_1 <- RANN::nn2(mat_1, k = num_neigh)$nn.idx
  nn_2 <- RANN::nn2(mat_2, k = num_neigh)$nn.idx
  
  full_rank <- length(tmp$obj_vec)
  common_score <- .common_decomposition(score_1, score_2, nn_1, nn_2)$common_score
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  tmp1 <- crossprod(common_score); tmp2 <- crossprod(distinct_score_1); tmp3 <- crossprod(distinct_score_2)
  
  stopifnot(abs(sum(abs(tmp1)) - sum(abs(diag(tmp1)))) <= 1e-6,  
            abs(sum(abs(tmp2)) - sum(abs(diag(tmp2)))) <= 1e-6,
            abs(sum(abs(tmp3)) - sum(abs(diag(tmp3)))) <= 1e-6)
  
  rank_c <- ncol(common_score)
  common_mat_1 <- common_score %*% coef_mat_1[1:rank_c,,drop = F] 
  distinct_mat_1 <- distinct_score_1 %*% coef_mat_1
  common_mat_2 <- common_score %*% coef_mat_2[1:rank_c,,drop = F] 
  distinct_mat_2 <- distinct_score_2 %*% coef_mat_2

  mat_1 <- noise_func(mat_1); mat_2 <- noise_func(mat_2)

  list(mat_1 = mat_1, mat_2 = mat_2, 
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2,
       distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
       common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2)
}

form_seurat_obj <- function(mat_1, mat_2){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  n <- nrow(mat_1); p_1 <- ncol(mat_1); p_2 <- ncol(mat_2)
  rownames(mat_1) <- paste0("n", 1:n); rownames(mat_2) <- paste0("n", 1:n)
  colnames(mat_1) <- paste0("g", 1:p_1); colnames(mat_2) <- paste0("p", 1:p_2)
  
  obj <- Seurat::CreateSeuratObject(counts = t(mat_1), assay = "mode1")
  obj[["mode2"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
  
  obj
}

generate_sbm_orthogonal <- function(B_mat, membership_vec){
  stopifnot(ncol(B_mat) == nrow(B_mat), all(B_mat >= 0), all(B_mat <= 1),
            sum(abs(B_mat - t(B_mat))) <= 1e-6,
            all(membership_vec > 0), all(membership_vec %% 1 == 0), 
            max(membership_vec) == length(unique(membership_vec)))
  
  K <- max(membership_vec)
  
  prob_mat <- .compute_prob_mat(B_mat, membership_vec)
  adj_mat <- .generate_adjaceny_mat(prob_mat)
  .svd_truncated(adj_mat, K = K)$u
}

generate_random_orthogonal <- function(n, K){
  stopifnot(K <= n)
  mat <- matrix(stats::rnorm(n^2), n, n)
  mat <- mat + t(mat)
  eigen(mat)$vectors[, 1:K, drop = F]
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
#' @param prob_mat symmetric probability matrix
#'
#' @return symmetric adjacency matrix
#' @export
.generate_adjaceny_mat <- function(prob_mat){
  upper_tri_idx <- upper.tri(prob_mat)
  prob_upper <- prob_mat[upper_tri_idx]
  
  n <- nrow(prob_mat)
  adj_upper <- stats::rbinom(n*(n-1)/2, 1, prob_upper)
  adj_mat <- matrix(0, ncol = n, nrow = n)
  adj_mat[upper_tri_idx] <- adj_upper
  adj_mat <- adj_mat + t(adj_mat)
  
  adj_mat
}
