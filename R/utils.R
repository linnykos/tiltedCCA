.l2norm <- function(x){sqrt(sum(x^2))}

##############

.svd_truncated <- function(mat, K, symmetric, rescale,
                           mean_vec, sd_vec,
                           K_full_rank){
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  if(K == min(dim(mat))) K_full_rank <- T
  
  ## [[note to self: probably factor this out and add sd_vec]]
  if(length(mean_vec) == 1 && !is.null(mean_vec)){
    if(mean_vec){
      if(inherits(x = mat, what = c('dgCMatrix', 'dgTMatrix'))){
        mean_vec <- sparseMatrixStats::colMeans2(mat)
      } else {
        mean_vec <- matrixStats::colMeans2(mat)
      }
    } else{
      mean_vec <- NULL
    }
  }
  
  if(min(dim(mat)) > 2*(K+2)){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      if(symmetric){
        tmp <- irlba::partial_eigen(mat, n = ifelse(K_full_rank, K, K+2),
                                    center = mean_vec, scale = sd_vec)
        list(u = tmp$vectors, d = tmp$values, v = tmp$vectors)
      } else {
        irlba::irlba(mat, nv = ifelse(K_full_rank, K, K+2),
                     center = mean_vec, scale = sd_vec)
      }
    }, warning = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    }, error = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    })
  } else {
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  
  # pass row-names and column-names
  if(length(rownames(mat)) != 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) != 0) rownames(res$v) <- colnames(mat)
  
  # useful only if your application requires only the singular vectors
  # if the number of rows or columns is too large, the singular vectors themselves
  # are often a bit too small numerically
  if(rescale){
    n <- nrow(mat); p <- ncol(mat)
    res$u <- res$u * sqrt(n)
    res$v <- res$v * sqrt(p)
    res$d <- res$d / (sqrt(n)*sqrt(p))
  }
  
  res
}

.check_svd <- function(svd_res, dims, tol = 1e-6){
  idx <- sort(intersect(which(svd_res$d > tol), dims), decreasing = F)
  if(length(idx) == length(svd_res$d)) return(svd_res)
  
  svd_res$u <- svd_res$u[, idx, drop = F]
  svd_res$v <- svd_res$v[, idx, drop = F]
  svd_res$d <- svd_res$d[idx]
  
  svd_res
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

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= nrow(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  mat@i[(val1+1):val2]+1
}

## from https://stackoverflow.com/questions/11832170/element-wise-max-operation-on-sparse-matrices-in-r
## seems to be an error from https://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package
## this means we need to include "data.table" into the "Depends"
.pmax_sparse <- function(..., na.rm = FALSE) {
  
  # check that all matrices have conforming sizes
  num_rows <- unique(sapply(list(...), nrow))
  num_cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num_rows) == 1, length(num_cols) == 1)
  
  cat_summary <- data.table::rbindlist(lapply(list(...), Matrix::summary)) 
  out_summary <- cat_summary[, list(x = max(x)), by = c("i", "j")]
  
  Matrix::sparseMatrix(i = out_summary[,i],
                       j = out_summary[,j],
                       x = out_summary[,x],
                       dims = c(num_rows, num_cols))
}
