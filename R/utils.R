.l2norm <- function(x){sqrt(sum(x^2))}

##############

.svd_truncated <- function(mat, K = min(dim(mat))){
  stopifnot(min(dim(mat)) >= K)
  
  if(min(dim(mat)) > K+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = K + 2)
    }, error = function(e){
      svd(mat)
    })
  } else {
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  res
}

.check_svd <- function(svd_res, tol = 1e-6){
  idx <- which(svd_res$d > tol)
  if(length(idx) == length(svd_res$d)) return(svd_res)
  
  svd_res$u <- svd_res$u[, idx, drop = F]
  svd_res$v <- svd_res$v[, idx, drop = F]
  svd_res$d <- svd_res$d[idx]
  
  svd_res
}

# reparameterize according to CCA: so t(mat_1) %*% mat_2 is diagonal with descending values,
# and t(mat_1)%*%mat_1 and t(mat_2)%*%mat_2 are the identity
# this essentially runs CCA
# we should move/integrate this into a more appropriate function later
.reparameterize <- function(mat_1, mat_2){
  stopifnot(all(dim(mat_1) == dim(mat_2)), Matrix::rankMatrix(mat_1) == ncol(mat_1),
            Matrix::rankMatrix(mat_2) == ncol(mat_2))
  n <- nrow(mat_1); p <- nrow(mat_2)
  
  res <- .cca(mat_1, mat_2, rank_1 = ncol(mat_1), rank_2 = ncol(mat_2))
  mat_1b <- mat_1 %*% res$loading_1; mat_2b <- mat_2 %*% res$loading_2
  
  list(mat_1 = mat_1b, mat_2 = mat_2b, diag_vec = res$obj_vec)
}

.equalize_norm <- function(mat_1, mat_2){
  vec_1 <- apply(mat_1, 2, .l2norm)
  vec_2 <- apply(mat_2, 2, .l2norm)
  
  mat_1 <- .mult_mat_vec(mat_1, vec_2/vec_1)
  mat_1
}

###############################################

#' Projection of vector onto another vector
#'
#' Returns the component of \code{vec1} that is orthogonal to \code{vec2}
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param tol small positive number
#'
#' @return vector
.projection <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))
  
  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  as.numeric((diag(d) - vec2%*%t(vec2))%*%vec1)
}

#' Projection of vector onto rows of a matrix
#'
#' Returns the component of \code{vec} that is orthogonal to all the rows
#' of \code{mat}
#'
#' @param vec vector
#' @param mat matrix
#'
#' @return vector
.orthogonal_vector <- function(vec, mat){
  stopifnot(length(vec) == nrow(mat), ncol(mat) <= nrow(mat))
  n <- length(vec)
  proj_mat <- diag(n) - mat %*% tcrossprod(solve(crossprod(mat)), mat)
  
  as.numeric(proj_mat %*% vec)
}

.orthogonal_matrix <- function(mat1, mat2){
  stopifnot(nrow(mat1) == nrow(mat2), ncol(mat1) <= nrow(mat1), ncol(mat2) <= nrow(mat2))
  n <- nrow(mat1)
  proj_mat <- diag(n) - mat2 %*% tcrossprod(solve(crossprod(mat2)), mat2)
  
  proj_mat %*% mat1
}

.orthogonalize <- function(mat){
  if(ncol(mat) == 1) return(mat)
  
  mat2 <- mat
  for(j in 2:ncol(mat)){
    mat2[,j] <- .orthogonal_vector(mat[,j], mat[,1:(j-1),drop = F])
  }
  
  mat2
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
