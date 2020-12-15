.l2norm <- function(x){sqrt(sum(x^2))}

##############

.svd_truncated <- function(mat, K){
  stopifnot(min(dim(mat)) >= K + 2, nrow(mat) == ncol(mat))
  
  if(nrow(mat) > K+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = K + 2)
    }, error = function(e){
      svd(mat)
    })
  } else {
    res <- svd(mat)
  }
  
  res
}

#' Do an SVD projection
#'
#' Uses \code{RSpectra::svds} to compute the \code{k} leading singular vectors, but
#' sometimes there are numerical instability issues. In case of crashes, the code
#' then uses the default \code{svd} function.
#'
#' @param mat numeric matrix with \code{n} rows and \code{n} columns
#' @param K positive integer less than \code{n}
#' @param weighted boolean
#'
#' @return numeric matrix
.svd_projection <- function(mat, K, weighted = F){
  res <- .svd_truncated(mat, K)
  
  if(weighted){
    diag_mat <- .diag_matrix(sqrt(abs(res$d[1:K])))
    res$u[,1:K, drop = F] %*% diag_mat
  } else {
    res$u[,1:K,drop = F]
  }
  
  res
}

.diag_matrix <- function(vec){
  K <- length(vec)
  if(K == 1) {
    matrix(vec, 1, 1)
  } else {
    diag(vec)
  }
}

#######################

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == ncol(mat))
  mat * rep(vec, rep(nrow(mat), length(vec)))
}
