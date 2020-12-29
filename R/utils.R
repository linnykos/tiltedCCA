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
