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
  stopifnot(all(dim(mat_1) == dim(mat_2)))
  n <- nrow(mat_1); p <- nrow(mat_2)
  
  res <- .cca(mat_1, mat_2, rank_1 = ncol(mat_1), rank_2 = ncol(mat_2))
  mat_1b <- mat_1 %*% res$loading_1; mat_2b <- mat_2 %*% res$loading_2
  
  list(mat_1 = mat_1b, mat_2 = mat_2b, diag_vec = res$obj_vec)
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
