#' Construct fixed-radius NN graph
#'
#' @param common_score output from \code{dcca_decomposition}
#' @param distinct_score output from \code{dcca_decomposition}
#' @param svd_e list containing the SVD of the full matrix 
#' @param cell_subidx vector of integers between 1 and \code{nrow(common_score)}
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param frnn_approx small non-negative number
#' @param radius_quantile value between 0 and 1
#' @param bool_matrix boolean. If \code{TRUE}, output the graphs as a sparse matrix.
#' If \code{FALSE}, output the graphs as a list where each element of the list
#' corresponds with the element's neighbors
#' @param include_diag boolean on whether or not the diagonal is included in the 
#' graph (currently only impactful if \code{bool_matrix=TRUE}).
#' @param verbose boolean
#'
#' @return list, depends on \code{bool_matrix}
#' @export
construct_frnn <- function(common_score, distinct_score, svd_e,
                           membership_vec,
                           max_subsample_frnn = nrow(common_score),
                           nn, frnn_approx = 0, radius_quantile = 0.9,
                           bool_matrix = F, include_diag = T, verbose = T){
  
  stopifnot(nrow(common_score) == nrow(distinct_score))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1)
  
  # compute the radius
  # [[note to self: the code below is quite messy -- fix it when i finalized the method]]
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
  if(length(cell_subidx) < nrow(common_score)) {
    membership_vec <- membership_vec[cell_subidx]
  }
  
  # extract embeddings
  c_embedding <- .extract_matrix_helper(common_score, distinct_score, svd_e,
                                        common_bool = T, distinct_bool = F, add_noise = F,
                                        center = F, renormalize = F)
  d_embedding <- .extract_matrix_helper(common_score, distinct_score, svd_e,
                                        common_bool = F, distinct_bool = T, add_noise = F,
                                        center = F, renormalize = F)
  e_embedding <- .extract_matrix_helper(common_score, distinct_score, svd_e,
                                        common_bool = T, distinct_bool = T, add_noise = F,
                                        center = F, renormalize = F)
  
  c_embedding <- c_embedding[cell_subidx,,drop = F]
  d_embedding <- d_embedding[cell_subidx,,drop = F]
  e_embedding <- e_embedding[cell_subidx,,drop = F]
  n <- nrow(c_embedding)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- common"))
  c_rad <- .compute_radius(c_embedding, nn, radius_quantile, 1:n)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- distinct"))
  d_rad <- .compute_radius(d_embedding, nn, radius_quantile, 1:n)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- everything"))
  e_rad <- .compute_radius(e_embedding, nn, radius_quantile, 1:n)
  sub_rad <- max(c_rad, d_rad)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- common"))
  c_g <- .construct_frnn(c_embedding, radius = sub_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- distinct"))
  d_g <- .construct_frnn(d_embedding, radius = sub_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- everything"))
  e_g <- .construct_frnn(e_embedding, radius = e_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  
  if(bool_matrix){
    c_g <- .nnlist_to_matrix(c_g, include_diag)
    d_g <- .nnlist_to_matrix(d_g, include_diag)
    e_g <- .nnlist_to_matrix(e_g, include_diag)
  }
  
  return(list(c_g = c_g, d_g = d_g, e_g = e_g, 
              membership_vec = membership_vec))
}

########################

.construct_celltype_subsample <- function(membership_vec, max_subsample_cell){
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    stopifnot(length(idx) > 2)
    
    if(length(idx) <= max_subsample_cell) return(idx)
    
    sample(idx, max_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}

.compute_radius <- function(mat, nn, radius_quantile, cell_subidx){
  res <- RANN::nn2(mat[cell_subidx,,drop = F], k = nn)
  as.numeric(stats::quantile(res$nn.dists[,nn], probs = radius_quantile))[1]
}

.construct_frnn <- function(mat, radius, nn, frnn_approx, verbose = F){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing frNN"))
  res <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)
  
  idx <- which(sapply(res$id, length) < nn)
  if(verbose) print(paste0(Sys.time(),": cLISI: ", length(idx), " cells with too few neighbors"))
  if(length(idx) > 0){
    res2 <- RANN::nn2(mat, query = mat[idx,,drop = F], k = nn+1, eps = frnn_approx)
    
    if(verbose) print(paste0(Sys.time(),": cLISI: Plugging in kNN neighbors"))
    for(i in 1:length(idx)){
      tmp_id <- c(res$id[[idx[i]]], res2$nn.idx[i,-1])
      tmp_dist <- c(res$dist[[idx[i]]], res2$nn.dists[i,-1])
      duplicates <- !duplicated(tmp_id)
      
      res$id[[idx[i]]] <- tmp_id[duplicates]
      res$dist[[idx[i]]] <- tmp_dist[duplictes]
    }
  }
  
  res
}

.nnlist_to_matrix <- function(rann_obj, include_diag){
  if(!include_diag){
    for(i in 1:length(rann_obj$id)){
      idx <- rann_obj$id[[i]]
      rann_obj$id[[i]] <- rann_obj$id[[i]][idx != i]
      rann_obj$dist[[i]] <- rann_obj$dist[[i]][idx != i]
    }
  }
  
  j_vec <- unlist(rann_obj$id)
  i_vec <- unlist(lapply(1:length(rann_obj$id), function(i){
    rep(i, length(rann_obj$id[[i]]))
  }))
  x_vec <- unlist(rann_obj$dist)
  Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)
}
