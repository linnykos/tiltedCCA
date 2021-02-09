#' Weight based on D-CCA
#'
#' @param mat Desired observed matrix
#' @param common_mat Desired common matrix
#' @param tol numeric
#' @param iter_max numeric
#' @param verbose boolean
#'
#' @return three vectors as a list
#' @export
dcca_information_weight <- function(mat, common_mat, tol = 0.05, iter_max = 100,
                                    verbose = T){
  if(verbose) print("D-CCA: Initializing weights")
  ratio_mat <- common_mat/mat
  ratio_mat <- pmin(pmax(ratio_mat, tol), 1)
  svd_res <- .svd_truncated(ratio_mat, K = 1)
  if(sum(svd_res$u<0) + sum(svd_res$v<0) > sum(svd_res$u>0) + sum(svd_res$v>0)){
    svd_res$u <- -svd_res$u; svd_res$v <- -svd_res$v
  }
  
  alpha_vec <- as.numeric(svd_res$u*sqrt(svd_res$d))
  alpha_vec <- pmin(pmax(alpha_vec, tol), 1)
  beta_vec <- as.numeric(svd_res$v*sqrt(svd_res$d))
  beta_vec <- pmin(pmax(beta_vec, tol), 1)
  
  # now alternating minimization
  if(verbose) print("D-CCA: Starting alternating minimization")
  n <- length(alpha_vec); p <- length(beta_vec)
  obj_vec <- rep(NA, iter_max)
  for(iter in 1:iter_max){
    if(verbose && iter %% floor(iter_max/10) == 0) cat('*')
    
    obj_vec[iter] <- .l2norm(common_mat - .mult_vec_mat(alpha_vec, .mult_mat_vec(mat, beta_vec)))^2/(n*p)
    if(iter > 1 && (abs(obj_vec[iter] - obj_vec[iter-1])/ obj_vec[iter-1]) <= 1e-6) break()
    
    tmp <- .mult_mat_vec(mat, beta_vec)
    num_vec <- sapply(1:n, function(i){sum(common_mat[i,]*mat[i,])})
    denom_vec <- sapply(1:n, function(i){sum(mat[i,]^2)})
    alpha_vec <- num_vec/denom_vec
    alpha_vec <- pmin(pmax(alpha_vec, tol), 1)
    
    tmp <- .mult_vec_mat(alpha_vec, mat)
    num_vec <- sapply(1:p, function(j){sum(common_mat[,j]*mat[,j])})
    denom_vec <- sapply(1:p, function(j){sum(mat[,j]^2)})
    beta_vec <- num_vec/denom_vec
    beta_vec <- pmin(pmax(beta_vec, tol), 1)
  }

  list(alpha_vec = alpha_vec, beta_vec = beta_vec, obj_vec = obj_vec)
}