## given 2 vectors, find the hyperplane 2D representation
.representation_2d <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  
  unit1 <- vec1/.l2norm(vec1)
  tmp <- .orthogonal_vec2vec(vec2, unit1)
  unit2 <- tmp/.l2norm(tmp)
  
  basis_mat <- cbind(unit1, unit2)
  rep1 <- c(.l2norm(vec1), 0)
  rep2 <- as.numeric(crossprod(basis_mat, vec2))
  names(rep2) <- NULL
  
  list(basis_mat = basis_mat, rep1 = rep1, rep2 = rep2)
}

#############

.cor_vectors <- function(vec1, vec2){
  (vec1 %*%vec2)/(.l2norm(vec1)*.l2norm(vec2))
}

.angle_between_vectors <- function(vec1, vec2){
  acos(.cor_vectors(vec1, vec2))*180/pi
}

.angle_from_vector <- function(vec, angle){
  ang <- .angle_between_vectors(vec, c(1,0))
  ang <- ang + angle
  
  rad <- ang * pi/180
  c(cos(rad), sin(rad))
}

.rightmost_vector <- function(vec1, vec2){
  stopifnot(length(vec1) == 2, length(vec2) == 2, all(c(vec1,vec2) >= 0))
  
  ang1 <- .angle_between_vectors(vec1, c(1,0))
  ang2 <- .angle_between_vectors(vec2, c(1,0))
  
  if(ang1 < ang2) {
    list(vec_right = vec1, vec_left = vec2, len_right = .l2norm(vec1), len_left = .l2norm(vec2))
  } else {
    list(vec_right = vec2, vec_left = vec1, len_right = .l2norm(vec2), len_left = .l2norm(vec1))
  }
}

#' Projection of vector onto another vector
#'
#' Returns the component of \code{vec1} that is orthogonal to \code{vec2}
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param tol small positive number
#'
#' @return vector
.orthogonal_vec2vec <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))
  
  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  as.numeric((diag(d) - vec2%*%t(vec2))%*%vec1)
}

