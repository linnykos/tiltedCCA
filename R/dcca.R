# from Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.

#' D-CCA Factorization
#'
#' @param mat_1 data matrix 1
#' @param mat_2 data matrix 2
#' @param rank_1 desired rank of data matrix 1
#' @param rank_2 desired rank of data matrix 2
#' @param meta_clustering optional clustering
#' @param apply_shrinkage boolean 
#' @param reorthogonalize boolean
#' @param verbose boolean
#'
#' @return list of class \code{dcca}
#' @export
dcca_factor <- function(mat_1, mat_2, rank_1, rank_2, meta_clustering = NA,
                        apply_shrinkage = T, reorthogonalize = F, verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2), 
            rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
  n <- nrow(mat_1)
  if(verbose) print("D-CCA: Rescaling matrices")
  mat_1 <- scale(mat_1, center = T, scale = F)
  mat_2 <- scale(mat_2, center = T, scale = F)
  
  if(verbose) print(paste0("D-CCA: Starting matrix shrinkage"))
  if(apply_shrinkage) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
  if(apply_shrinkage) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)
  
  svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  
  if(all(!is.na(meta_clustering))){
    # apply D-CCA to meta-cells
    msg <- " (meta-cells)"
    
    stopifnot(all(meta_clustering > 0), all(meta_clustering %% 1 == 0),
              max(meta_clustering) == length(unique(meta_clustering)),
              length(meta_clustering) == nrow(mat_1))
    num_meta <- max(meta_clustering)
    
    if(verbose) print(paste0("D-CCA", msg, ": Constructing meta-cells for matrix 1"))
    mat_1_meta <- t(sapply(1:num_meta, function(x){
      if(verbose & x %% floor(num_meta/10) == 0) cat('*')
      idx <- which(meta_clustering == x)
      apply(mat_1[idx,,drop = F], 2, mean)
    }))
    
    if(verbose) print(paste0("D-CCA", msg, ": Constructing meta-cells for matrix 2"))
    mat_2_meta <- t(sapply(1:num_meta, function(x){
      if(verbose & x %% floor(num_meta/10) == 0) cat('*')
      idx <- which(meta_clustering == x)
      apply(mat_2[idx,,drop = F], 2, mean)
    }))
    
    if(verbose) print(paste0("D-CCA", msg, ": Computing CCA"))
    # note, since we're already doing averaging, we don't further shrink the spectrum
    cca_res <- .cca(mat_1_meta, mat_2_meta, rank_1 = rank_1, rank_2 = rank_2)
  } else {
    
    # alternatively, apply D-CCA to all cells
    msg <- " (all cells)"
    if(verbose) print(paste0("D-CCA", msg, ": Computing CCA"))
    cca_res <- .cca(svd_1, svd_2)
  }
  
  res <- .dcca_common_score(svd_1, svd_2, cca_res,
                            check_alignment = all(!is.na(meta_clustering)),
                            reorthogonalize = reorthogonalize,
                            verbose = verbose, msg = msg)

  class(res) <- "dcca"
  res
}

#' D-CCA Decomposition
#'
#' @param dcca_res output from \code{dcca_factor}
#' @param rank_c desired rank of cross-correlation matrix between \code{mat_1} and \code{mat_2} when running \code{dcca_factor}
#' @param verbose boolean
#'
#' @return list of class \code{dcca_decomp}
#' @export
dcca_decomposition <- function(dcca_res, rank_c, verbose = T){
  stopifnot(class(dcca_res) == "dcca")
  n <- nrow(dcca_res$svd_1$u)
  full_rank <- length(dcca_res$cca_obj)
  
  if(verbose) print("D-CCA: Form denoised observation matrices")
  mat_1 <- tcrossprod(.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d) , dcca_res$svd_1$v)
  mat_2 <- tcrossprod(.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d) , dcca_res$svd_2$v)
  
  if(verbose) print("D-CCA: Computing common matrices")
  coef_mat_1 <- crossprod(dcca_res$score_1, mat_1)/n
  coef_mat_2 <- crossprod(dcca_res$score_2, mat_2)/n
  
  common_mat_1 <- dcca_res$common_score[,1:rank_c, drop = F] %*% coef_mat_1[1:rank_c,,drop = F]
  common_mat_2 <- dcca_res$common_score[,1:rank_c, drop = F] %*% coef_mat_2[1:rank_c,,drop = F]
   
  if(verbose) print("D-CCA: Computing distinctive matrices")
  distinct_mat_1 <- dcca_res$distinct_score_1 %*% coef_mat_1
  distinct_mat_2 <- dcca_res$distinct_score_2 %*% coef_mat_2
  
  if(verbose) print("D-CCA: Done")
  structure(list(common_score = dcca_res$common_score[,1:rank_c, drop = F],
       distinct_score_1 = dcca_res$distinct_score_1,
       distinct_score_2 = dcca_res$distinct_score_2,
       common_mat_1 = common_mat_1, common_mat_2 = common_mat_2, 
       distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2,
       cca_obj = dcca_res$cca_obj), class = "dcca_decomp")
}


#################################

.dcca_common_score <- function(svd_1, svd_2, cca_res, check_alignment = T,
                               reorthogonalize = F, verbose = T, msg = ""){
  full_rank <- length(cca_res$obj_vec)
  n <- nrow(svd_1$u)
  
  if(verbose) print(paste0("D-CCA", msg, ": Computing unnormalized scores"))
  tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
  score_1 <- tmp$score_1; score_2 <- tmp$score_2
  stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
            nrow(score_1) == nrow(score_2))

  if(verbose) print(paste0("D-CCA", msg, ": Computing common factors"))
  if(check_alignment){
    # reparameterize the scores
    tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
    score_1 <- tmp$score_1; score_2 <- tmp$score_2
    stopifnot(is.matrix(score_1), is.matrix(score_2))
    obj_vec <- diag(crossprod(score_1, score_2))/n
  } else {
    obj_vec <- cca_res$obj_vec
  }
  
  # threshold x for numerical stability
  common_score <- .compute_common_score(score_1, score_2, obj_vec = obj_vec)
  stopifnot(all(dim(common_score) == c(nrow(score_1), full_rank)))
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score, reorthogonalize = reorthogonalize)
  common_score <- tmp$common_score
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
 
  if(verbose) print(paste0("D-CCA", msg, ": Done"))
  list(common_score = common_score, 
       distinct_score_1 = distinct_score_1,
       distinct_score_2 = distinct_score_2,
       score_1 = score_1, score_2 = score_2, 
       svd_1 = svd_1, svd_2 = svd_2, cca_obj = obj_vec)
}

# score is a matrix with n rows. loadings is a mat with p rows
.compute_common_score <- function(score_1, score_2, obj_vec = NA){
  stopifnot(nrow(score_1) == nrow(score_2), nrow(score_1) >= ncol(score_1),
            nrow(score_2) >= ncol(score_2), is.matrix(score_1), is.matrix(score_2))
  
  p <- min(c(ncol(score_1), ncol(score_2)))
  if(ncol(score_1) > p){ score_1 <- score_1[,1:p,drop = F]}
  if(ncol(score_2) > p){ score_2 <- score_2[,1:p,drop = F]}
  
  n <- nrow(score_1)
  obj_vec <- diag(crossprod(score_1, score_2))/n
  
  R_vec <- sapply(obj_vec, function(x){x <- min(1,max(x,0)); 1-sqrt((1-x)/(1+x))})
  common_score <- .mult_mat_vec((score_1+score_2)/2, R_vec)
  
  common_score
}

.compute_distinct_score <- function(score_1, score_2, common_score, reorthogonalize = F){
  full_rank <- ncol(common_score)
  distinct_score_1 <- score_1[,1:full_rank, drop = F] - common_score
  distinct_score_2 <- score_2[,1:full_rank, drop = F] - common_score
  
  if(reorthogonalize){
    tmp <- .reorthgonalize_scores(common_score, distinct_score_1, distinct_score_2)
    common_score <- tmp$common_score
    distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  }
  
  # handle different ranks
  distinct_score_1 <- cbind(distinct_score_1, score_1[,-(1:full_rank), drop = F])
  distinct_score_2 <- cbind(distinct_score_2, score_2[,-(1:full_rank), drop = F])
  
  list(common_score = common_score, distinct_score_1 = distinct_score_1, 
       distinct_score_2 = distinct_score_2)
}

.reorthgonalize_scores <- function(common_score, distinct_score_1, distinct_score_2, check = F){
  tmp <- solve(crossprod(common_score))
  tmp2 <- tcrossprod(tmp, common_score)
  
  common_comp1 <- common_score %*% (tmp2 %*% distinct_score_1)
  
  if(check){
    common_comp2 <- common_score %*% (tmp2 %*% distinct_score_2)
    stopifnot(all(abs(common_comp1 - common_comp2) <= 1e-6))
  }
 
  common_score <- common_score + common_comp1
  distinct_score_1 <- distinct_score_1 - common_comp1
  distinct_score_2 <- distinct_score_2 - common_comp1
  
  list(common_score = common_score, distinct_score_1 = distinct_score_1, 
       distinct_score_2 = distinct_score_2)
}

##############################################

.spoet <- function(mat, K){
  n <- nrow(mat); p <- ncol(mat); m <- min(n, p)
  target_full_dim <- min(c(nrow(mat), ncol(mat), K*10))
  svd_res <- .svd_truncated(mat, target_full_dim)
  tau <- sum((svd_res$d[(K+1):length(svd_res$d)])^2)/(n*p - n*K - p*K)
  sing_vec <- sapply(svd_res$d[1:K], function(x){
    sqrt(max(c(x^2-tau*p, 0)))
  })
  
  svd_res$d_original <- svd_res$d[1:K]
  svd_res$d <- sing_vec
  svd_res$u <- svd_res$u[,1:K,drop = F]
  svd_res$v <- svd_res$v[,1:K,drop = F]
  
  # pass row-names and column-names
  if(length(rownames(mat)) != 0) rownames(svd_res$u) <- rownames(mat)
  if(length(colnames(mat)) != 0) rownames(svd_res$v) <- colnames(mat)
  
  svd_res
}

# takes either two matrices or two SVDs
.cca <- function(input_1, input_2, rank_1 = NA, rank_2 = NA, 
                 return_scores = F, tol = 1e-6){
  if(is.matrix(input_1) & is.matrix(input_2)){
    stopifnot(nrow(input_1) == nrow(input_2), 
              rank_1 <= ncol(input_1), rank_2 <= ncol(input_2))
    svd_1 <- .svd_truncated(input_1, rank_1)
    svd_2 <- .svd_truncated(input_2, rank_2)
    
    svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
  } else {
    stopifnot(is.list(input_1), is.list(input_2), 
              all(c("u","d","v") %in% names(input_1)),
              all(c("u","d","v") %in% names(input_2)),
              all(input_1$d >= 0), all(input_2$d >= 0), 
              all(diff(input_1$d) <= 1e-6), all(diff(input_2$d) <= 1e-6),
              nrow(input_1$u) == nrow(input_2$u), is.na(rank_1), is.na(rank_2))
    
    rank_1 <- length(which(input_1$d >= 1e-6)); rank_2 <- length(which(input_2$d >= 1e-6))
    svd_1 <- input_1; svd_2 <- input_2
  }
  
  n <- nrow(svd_1$u)
  
  # perform CCA
  cov_1_invhalf <- .mult_mat_vec(svd_1$v, sqrt(n)/svd_1$d)
  cov_2_invhalf <- .mult_mat_vec(svd_2$v, sqrt(n)/svd_2$d)
                                    
  agg_mat <-  .compute_cca_aggregate_matrix(svd_1, svd_2, augment = T)
  svd_res <- svd(agg_mat)
  full_rank <- min(c(rank_1, rank_2))
  
  loading_1 <- cov_1_invhalf %*% svd_res$u[1:ncol(cov_1_invhalf), 1:rank_1, drop = F]
  loading_2 <- cov_2_invhalf %*% svd_res$v[1:ncol(cov_2_invhalf), 1:rank_2, drop = F]
  
  # return
  if(!return_scores){
    list(loading_1 = loading_1, loading_2 = loading_2, obj_vec = svd_res$d[1:full_rank])
  } else {
    if(!is.matrix(input_1) | !is.matrix(input_2)){
      input_1 <- tcrossprod(.mult_mat_vec(input_1$u, input_1$d), input_1$v) 
      input_2 <- tcrossprod(.mult_mat_vec(input_2$u, input_2$d), input_2$v) 
    }
    
    list(score_1 = input_1 %*% loading_1, score_2 = input_2 %*% loading_2, obj_vec = svd_res$d[1:full_rank])
  }
  
}

.compute_cca_aggregate_matrix <- function(svd_1, svd_2, augment = T){
  res <- crossprod(svd_1$u, svd_2$u)
  
  if(augment){
    if(nrow(res) < ncol(res)) res <- rbind(res, matrix(0, ncol(res)-nrow(res), ncol(res)))
    if(ncol(res) < nrow(res)) res <- cbind(res, matrix(0, nrow(res), nrow(res)-ncol(res)))
  }

  res
}

# unnoramlized means that while score_1 is orthogonal, it is not orthonormal
.compute_unnormalized_scores <- function(svd_1, svd_2, cca_res){
  score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, cca_res$loading_1) 
  score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, cca_res$loading_2)
  
  list(score_1 = score_1, score_2 = score_2)
}
