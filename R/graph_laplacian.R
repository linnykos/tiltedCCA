#' Compute the eigenbasis of a Laplacian based on a sparse matrix
#'
#' @param g_mat sparse matrix of class \code{dgCMatrix}
#' @param k_max positive integer
#' @param normalize boolean
#' @param rowname_vec vector of rownames
#' @param colname_vec vector of colnames
#'
#' @return matrix of dimension \code{nrow(g_mat)} by \code{k_max}
#' @export
compute_laplacian <- function(g_mat, k_max = 100, normalize = T,
                        rowname_vec, colname_vec, verbose = T){
  stopifnot(inherits(g_mat, "dgCMatrix"))
  
  n <- nrow(g_mat)
  g_mat <- .symmetrize_sparse(g_mat, set_ones = T)
  
  if(normalize){
    deg_vec <- sparseMatrixStats::rowSums2(g_mat)
    invdeg_mat <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec)
    g_mat <- invdeg_mat %*% g_mat %*% invdeg_mat
  }
  
  deg_vec2 <- sparseMatrixStats::rowSums2(g_mat)
  invdeg_mat2 <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec2)
  lap_mat <- invdeg_mat2 %*% g_mat
  
  eigenbasis <- .laplacian_eigenvectors(lap_mat, k_max = k_max, reorient = T,
                                        print_approximation = verbose)
  rownames(eigenbasis) <- rowname_vec
  colnames(eigenbasis) <- colname_vec
  
  eigenbasis
}

#' Compute smoothed vector based on eigenbases
#'
#' @param vec vector
#' @param eigenbasis matrix with rows equal to \code{length(vec)}
#'
#' @return a list
#' @export
compute_smooth_signal <- function(vec, eigenbasis){
  lm_res <- stats::lm(vec ~ eigenbasis) 
  list(smoothed_vec = lm_res$fitted.values, 
       variance = stats::var(lm_res$fitted.values),
       r_squared = summary(lm_res)$adj.r.squared)
}

####################3

.laplacian_eigenvectors <- function(mat, k_max, reorient, print_approximation){
  res <- RSpectra::eigs(mat, k = k_max) # it's important to use eigs since this matrix isn't symmetric
  res$values <- Re(res$values)
  res$vectors <- Re(res$vectors)
  eigenbasis <- .mult_mat_vec(res$vectors[,-1], res$values[-1])
  
  if(print_approximation){
    inv_vectors <- MASS::ginv(res$vectors)
    approx_mat <- .mult_mat_vec(res$vectors, res$values) %*% inv_vectors
    print(paste0("Approximation quality: ", round(sqrt(sum((mat - approx_mat)^2))/sqrt(sum(mat^2)),2)))
  }
 
  if(reorient){
    for(j in 1:ncol(eigenbasis)){
      max_pos <- max(pmax(eigenbasis[,j], 0))
      max_neg <- abs(min(pmin(eigenbasis[,j],0)))
      if(max_neg > max_pos) eigenbasis[,j] <- -eigenbasis[,j]
    }
  }
  
  eigenbasis
}

