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