## given 2 vectors, find the hyperplane 2D representation
.representation_2d <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2), .l2norm(vec1) > tol, .l2norm(vec2) > tol)
  
  unit1 <- vec1/.l2norm(vec1)
  tmp <- .project_vec2vec(vec2, unit1, orthogonal = T)
  if(.l2norm(tmp) < tol){
    return(list(basis_mat = cbind(unit1, 0), rep1 = c(.l2norm(vec1), 0), 
                rep2 = c(.l2norm(vec2), 0)))
  }
  unit2 <- tmp/.l2norm(tmp)
  
  basis_mat <- cbind(unit1, unit2)
  rep1 <- c(.l2norm(vec1), 0)
  rep2 <- as.numeric(crossprod(basis_mat, vec2))
  names(rep2) <- NULL
  
  list(basis_mat = basis_mat, rep1 = rep1, rep2 = rep2)
}

#############

.cor_vectors <- function(vec1, vec2, tol = 1e-3){
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  if(len1 <= tol | len2 <= tol) return(NA)
  (vec1 %*%vec2)/(len1*len2)
}

.angle_between_vectors <- function(vec1, vec2){
  cor_val <- .cor_vectors(vec1, vec2)
  if(is.na(cor_val)) return(NA)
  ang <- acos(min(max(cor_val,-1),1))*180/pi

  ang
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
#' if \code{orthogonal} is \code{TRUE}, and the
#' component of \code{vec1} that is parallel to \code{vec2} if 
#' \code{orthogonal} is \code{FALSE}.
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param orthogonal boolean
#' @param tol small positive number
#'
#' @return vector
.project_vec2vec <- function(vec1, vec2, orthogonal, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))
  
  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  if(orthogonal){
    vec1 - vec2 %*% crossprod(vec2, vec1)
  } else {
    vec2 %*% crossprod(vec2, vec1)
  }
}

############################
# both points are lie on opposite ends of the circle
.construct_circle <- function(endpoint1, endpoint2){
  stopifnot(length(endpoint1) == 2, length(endpoint2) == 2)
  
  center <- c(endpoint1 + endpoint2)/2
  radius <- .l2norm(endpoint1 - center)
  
  list(center = center, radius = radius)
}

.find_radian <- function(circle, point, tol = 1e-6){
  stopifnot(length(point) == 2, length(circle$center) == 2,
            .l2norm(circle$center - point) <= circle$radius+tol)
  
  atan2(point[2]-circle$center[2], point[1]-circle$center[1])
}

.position_from_circle <- function(circle, radian){
  circle$radius *c(cos(radian),  sin(radian)) + circle$center
}
