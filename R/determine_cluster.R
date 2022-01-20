.grassmann_distance <- function(orthonormal_1,
                                orthonormal_2,
                                tol = 1e-3){
  stopifnot(all(nrow(orthonormal_1) == nrow(orthonormal_2)),
            ncol(orthonormal_1) <= nrow(orthonormal_1),
            ncol(orthonormal_2) <= nrow(orthonormal_2))
  
  l2_vec <- apply(orthonormal_1, 2, .l2norm)
  orthonormal_1 <- .mult_mat_vec(orthonormal_1, 1/l2_vec)
  l2_vec <- apply(orthonormal_2, 2, .l2norm)
  orthonormal_2 <- .mult_mat_vec(orthonormal_2, 1/l2_vec)
  
  stopifnot(sum(abs(crossprod(orthonormal_1) - diag(ncol(orthonormal_1)))) <= tol,
            sum(abs(crossprod(orthonormal_2) - diag(ncol(orthonormal_2)))) <= tol)
  
  k <- ncol(orthonormal_1)
  crossprod_mat <- crossprod(orthonormal_1, orthonormal_2)
  svd_res <- svd(crossprod_mat)
  sing_vec <- svd_res$d
  theta_vec <- acos(sing_vec)
  .l2norm(theta_vec)
}