.l2norm <- function(x){sqrt(sum(x^2))}

##############

.svd_truncated <- function(mat, K = min(dim(mat))){
  stopifnot(min(dim(mat)) >= K)
  
  if(min(dim(mat)) > K+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = K+2)
    }, error = function(e){
      irlba::irlba(mat, nv = K+2)
    })
  } else {
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  
  # pass row-names and column-names
  if(length(rownames(mat)) != 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) != 0) rownames(res$v) <- colnames(mat)
  
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

.equalize_norm <- function(mat_1, mat_2){
  vec_1 <- apply(mat_1, 2, .l2norm)
  vec_2 <- apply(mat_2, 2, .l2norm)
  
  mat_1 <- .mult_mat_vec(mat_1, vec_2/vec_1)
  mat_1
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

###########################

# return idx such that vec1[idx] == vec2
.matching_idx <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  
  ord_1 <- order(vec1, decreasing = F)
  rank_2 <- rank(vec2)
  
  ord_1[rank_2]
}
