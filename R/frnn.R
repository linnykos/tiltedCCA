#' Construct fixed-radius NN graph
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param membership_vec factor vector
#' @param data_1 boolean
#' @param data_2 boolean
#' @param max_subsample_frnn positive integer, used for determining number of cells to 
#' compute cLISI for
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
construct_frnn <- function(obj, nn, membership_vec, data_1 = T, data_2 = F,
                           max_subsample_frnn = nrow(obj$common_score),
                           frnn_approx = 0, radius_quantile = 0.9,
                           bool_matrix = T, include_diag = T, verbose = T){
  stopifnot(frnn_approx >= 0, frnn_approx <= 1,
            length(membership_vec) == nrow(obj$common_score),
            is.factor(membership_vec))
  
  embedding <- .prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                   add_noise = F)
  n <- nrow(embedding[[1]])
  
  # construct subsamples
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
  if(length(cell_subidx) < n) {
    membership_vec <- membership_vec[cell_subidx]
  }
  for(i in 1:3){
    embedding[[i]] <- embedding[[i]][cell_subidx,,drop = F]
  }
  n <- nrow(embedding[[1]])
  
  # compute the radius
  vec_print <- c("common", "distinct", "everything")
  vec_rad <- sapply(1:3, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- ", vec_print[i]))
    .compute_radius(embedding[[i]], nn, radius_quantile)
  })
  vec_rad_org <- vec_rad
  vec_rad[1:2] <- max(vec_rad[1:2])
  
  # construct frnn
  list_g <- lapply(1:3, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- ", vec_print[i]))
    .construct_frnn(embedding[[i]], radius = vec_rad[i], nn = nn, 
                    frnn_approx = frnn_approx, verbose = verbose)
  }) 
 
  # convert to matrix if needed
  if(bool_matrix){
    for(i in 1:3){
      list_g[[i]] <- .nnlist_to_matrix(list_g[[i]], include_diag)
      
      if(length(rownames(obj$common_score)) != 0){
        rownames(list_g[[i]]) <- rownames(obj$common_score)
        colnames(list_g[[i]]) <- rownames(obj$common_score)
      }
    }
  }
  
  return(list(c_g = list_g[[1]], d_g = list_g[[2]], e_g = list_g[[3]], 
              membership_vec = membership_vec,
              original_radius = vec_rad_org))
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

.compute_radius <- function(mat, nn, radius_quantile){
  res <- RANN::nn2(mat, k = nn)
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
      res$dist[[idx[i]]] <- tmp_dist[duplicates]
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
