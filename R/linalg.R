## given 2 vectors, find the hyperplane 2D representation

.angle_between_vectors <- function(vec1, vec2){
  theta <- (vec1 %*%vec2)/(.l2norm(vec1)*.l2norm(vec2))
  acos(theta)*180/pi
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