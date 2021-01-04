#' Generate data
#'
#' @param common_loading \code{n} by \code{rank_12} orthogonal matrix
#' @param distinct_loading_1 \code{n} by \code{rank_1} orthogonal matrix where \code{rank_1} is larger than \code{rank_12}
#' @param distinct_loading_2 \code{n} by \code{rank_2} orthogonal matrix where \code{rank_2} is larger than \code{rank_12}
#' @param coef_mat_1 \code{rank_1} by \code{p_1} matrix
#' @param coef_mat_2 \code{rank_2} by \code{p_2} matrix
#' @param noise_func function that takes in a matrix and outputs a matrix of the same dimension
#'
#' @return list of two matrices, one of dimension \code{n} by \code{p_1} and another of dimension \code{n} by \code{p_2}
#' @export
generate_data <- function(common_loading, distinct_loading_1, distinct_loading_2,
                          coef_mat_1, coef_mat_2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat), nrow(mat), ncol(mat))}){
  stopifnot(nrow(common_loading) == nrow(distinct_loading_1), nrow(common_loading) == nrow(distinct_loading_2),
            ncol(common_loading) <= ncol(distinct_loading_1), ncol(common_loading) <= ncol(distinct_loading_2),
            ncol(distinct_loading_1) == nrow(coef_mat_1), ncol(distinct_loading_2) == nrow(coef_mat_2))
  tmp1 <- crossprod(common_loading); tmp2 <- crossprod(distinct_loading_1); tmp3 <- crossprod(distinct_loading_2)
  stopifnot(abs(sum(abs(tmp1)) - sum(abs(diag(tmp1)))) <= 1e-6,  
            abs(sum(abs(tmp2)) - sum(abs(diag(tmp2)))) <= 1e-6,
            abs(sum(abs(tmp3)) - sum(abs(diag(tmp3)))) <= 1e-6)
  
  rank_12 <- ncol(common_loading)
  mat_1 <- common_loading %*% coef_mat_1[1:rank_12,,drop = F] + distinct_loading_1 %*% coef_mat_1
  mat_2 <- common_loading %*% coef_mat_2[1:rank_12,,drop = F] + distinct_loading_2 %*% coef_mat_2
  
  mat_1 <- noise_func(mat_1); mat_2 <- noise_func(mat_2)
  
  list(mat_1 = mat_1, mat_2 = mat_2)
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
  svd_res <- .svd_truncated(adj_mat, K = K)
  .mult_mat_vec(svd_res$u, svd_res$d)
}

generate_random_orthogonal <- function(n, K){
  stopifnot(K <= n)
  mat <- matrix(stats::rnorm(n^2), n, n)
  mat <- mat + t(mat)
  eigen(mat)$vectors[, 1:K, drop = F]
}

equalize_norm <- function(mat_1, mat_2){
  vec_1 <- apply(mat_1, 2, .l2norm)
  vec_2 <- apply(mat_2, 2, .l2norm)
  
  mat_1 <- .mult_mat_vec(mat_1, vec_2/vec_1)
  mat_1
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

  
  
  
  