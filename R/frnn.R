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
#' @param verbose boolean
#'
#' @return list, depends on \code{bool_matrix}
#' @export
construct_frnn <- function(obj, nn, membership_vec, data_1 = T, data_2 = F,
                           max_subsample_frnn = nrow(obj$common_score),
                           frnn_approx = 0, radius_quantile = 0.5,
                           bool_matrix = T, verbose = T){
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
  names(vec_rad_org) <- c("common", "distinct", "everything")
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
      list_g[[i]] <- .nnlist_to_matrix(list_g[[i]])
      
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

combine_frnn <- function(dcca_obj, g_1, g_2, nn, 
                         common_1 = T, common_2 = T, keep_n = T,
                         verbose = T){
  stopifnot(all(dim(g_1) == dim(g_2)))
  
  # extract the relevant embeddings from dcca_obj
  embedding_1 <- .prepare_embeddings(dcca_obj, data_1 = T, data_2 = F, 
                                    add_noise = F)
  if(common_1){
    embedding_1 <- embedding_1$common
  } else {
    embedding_1 <- embedding_1$distinct
  }
  
  embedding_2 <- .prepare_embeddings(dcca_obj, data_1 = F, data_2 = T, 
                                     add_noise = F)
  if(common_2){
    embedding_2 <- embedding_2$common 
  } else {
    embedding_2 <- embedding_2$distinct
  }
  
  # symmetrize g_1 and g_2
  g_1 <- .symmetrize_sparse(g_1, set_ones = F)
  g_2 <- .symmetrize_sparse(g_2, set_ones = F)
  
  # prepare
  n <- nrow(g_1)
  nn_idx_1 <- lapply(1:n, function(j){.nonzero_col(g_1, j, bool_value = F)})
  nn_dist_1 <- lapply(1:n, function(j){.nonzero_col(g_1, j, bool_value = T)})
  nn_idx_2 <- lapply(1:n, function(j){.nonzero_col(g_2, j, bool_value = F)})
  nn_dist_2 <- lapply(1:n, function(j){.nonzero_col(g_2, j, bool_value = T)})
  
  # apply the following procedure for each cell n
  nn_idx_all <- vector("list", n); nn_dist_all <- vector("list", n)
  for(i in 1:n){
    # intersect
    idx_all <- unique(c(nn_idx_1[[i]], nn_idx_2[[i]]))
    idx_intersect <- intersect(nn_idx_1[[i]], nn_idx_2[[i]])
    if(length(idx_intersect) < nn){
      idx_subset <- setdiff(idx_all, idx_intersect)
      idx_intersect <- c(idx_intersect, sample(idx_subset, size = nn - length(idx_intersect)))
    }
    
    # start tabulating nn_dist
    nn_idx_all[[i]] <- idx_intersect
    nn_dist_all[[i]] <- .compute_distance_from_idx(embedding_1, embedding_2, 
                                                   nn_idx_1[[i]], nn_idx_2[[i]],
                                                   nn_dist_1[[i]], nn_dist_2[[i]],
                                                   start_idx = i, end_idx_vec = idx_intersect)
    
    # find the nn's
    if(length(nn_idx_1[[i]]) < nn){
      order_1 <- nn_idx_1[[i]]
    } else {
      order_1 <- order(nn_dist_1[[i]], decreasing = F)[1:nn]
    }
    if(length(nn_idx_2[[i]]) < nn){
      order_2 <- nn_idx_2[[i]]
    } else {
      order_2 <- order(nn_dist_2[[i]], decreasing = F)[1:nn]
    }
   
    tmp_idx <- unique(c(nn_idx_1[[i]][order_1], nn_idx_2[[i]][order_2]))
    tmp_dist <- .compute_distance_from_idx(embedding_1, embedding_2, 
                                           nn_idx_1[[i]], nn_idx_2[[i]],
                                           nn_dist_1[[i]], nn_dist_2[[i]],
                                           start_idx = i, end_idx_vec = tmp_idx)
    zz <- intersect(order(tmp_dist, decreasing = F)[1:nn], which(!tmp_idx %in% nn_idx_all[[i]]))
    nn_idx_all[[i]] <- c(nn_idx_all[[i]], tmp_idx[zz])
    nn_dist_all[[i]] <- c(nn_dist_all[[i]], tmp_dist[zz])
  }
  
  tmp_list <- list(id = nn_idx_all, dist = nn_dist_all)
  res <- .nnlist_to_matrix(tmp_list)
  
  if(length(rownames(dcca_obj$common_score)) != 0){
    rownames(res) <- rownames(dcca_obj$common_score)
    colnames(res) <- rownames(dcca_obj$common_score)
  }
  
  res
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

.nnlist_to_matrix <- function(rann_obj){
  stopifnot(all(c("id", "dist") %in% names(rann_obj)))
  
  for(i in 1:length(rann_obj$id)){
    idx <- rann_obj$id[[i]]
    rann_obj$id[[i]] <- rann_obj$id[[i]][idx != i]
    rann_obj$dist[[i]] <- rann_obj$dist[[i]][idx != i]
  }
  
  j_vec <- unlist(rann_obj$id)
  i_vec <- unlist(lapply(1:length(rann_obj$id), function(i){
    rep(i, length(rann_obj$id[[i]]))
  }))
  x_vec <- unlist(rann_obj$dist)
  Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)
}

.compute_distance_from_idx <- function(embedding_1, embedding_2, 
                                       nn_idx_1_vec, nn_idx_2_vec,
                                       nn_dist_1_vec, nn_dist_2_vec,
                                       start_idx, end_idx_vec){
  stopifnot(!start_idx %in% end_idx_vec)
  
  res <- sapply(end_idx_vec, function(j){
    tmp <- 0
    z1 <- which(nn_idx_1_vec == j)
    
    if(length(z1) == 1) {
      tmp <- tmp + nn_dist_1_vec[z1]
    } else {
      tmp <- tmp + .l2norm(embedding_1[start_idx,] - embedding_1[j,])
    }
    
    z2 <- which(nn_idx_2_vec == j)
    if(length(z2) == 1) {
      tmp <- tmp + nn_dist_2_vec[z2]
    } else {
      tmp <- tmp + .l2norm(embedding_2[start_idx,] - embedding_2[j,])
    }
  })
  
  res
}
