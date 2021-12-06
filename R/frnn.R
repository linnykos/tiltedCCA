#' Construct fixed-radius NN graph
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec factor vector
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param data_1 boolean
#' @param data_2 boolean
#' @param normalization_type character among \code{"cosine_itself"},
#' \code{"cosine_everything"}, \code{"signac_itself"}, \code{"signac_everything"},
#' and \code{"none"}
#' @param max_subsample_frnn positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param frnn_approx small non-negative number
#' @param radius_quantile value between 0 and 1
#' @param bool_matrix boolean. If \code{TRUE}, output the graphs as a sparse matrix.
#' If \code{FALSE}, output the graphs as a list where each element of the list
#' corresponds with the element's neighbors
#' @param center boolean
#' @param renormalize boolean
#' @param symmetrize boolean
#' @param verbose boolean
#'
#' @return list, depends on \code{bool_matrix}
#' @export
construct_frnn <- function(obj, 
                           membership_vec = NA, 
                           nn = 30, 
                           data_1 = T, 
                           data_2 = F,
                           normalization_type = "none",
                           max_subsample_frnn = nrow(obj$common_score),
                           frnn_approx = 0, 
                           radius_quantile = 0.5,
                           bool_matrix = T, 
                           center = T,
                           renormalize = F, 
                           symmetrize = F,
                           verbose = T){
  stopifnot(frnn_approx >= 0, frnn_approx <= 1)
  if(!all(is.na(membership_vec))){
    stopifnot(length(membership_vec) == nrow(obj$common_score),
              is.factor(membership_vec))
  }
  
  embedding <- .prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                   center = center, 
                                   renormalize = renormalize)
  embedding <- .normalize_embeddings(embedding, normalization_type)
  n <- nrow(embedding[[1]])
  
  # construct subsamples
  if(!all(is.na(membership_vec))){
    cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
    if(length(cell_subidx) < n) {
      membership_vec <- membership_vec[cell_subidx]
    }
    for(i in 1:3){
      embedding[[i]] <- embedding[[i]][cell_subidx,,drop = F]
    }
  }
  n <- nrow(embedding[[1]])
  
  # compute the radius
  vec_print <- c("common", "distinct")
  vec_rad <- sapply(1:2, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- ", vec_print[i]))
    .compute_radius(embedding[[i]], nn, radius_quantile)
  })
  vec_rad_org <- vec_rad
  names(vec_rad_org) <- c("common", "distinct")
  vec_rad[1:2] <- max(vec_rad[1:2])
  
  # construct frnn
  list_g <- lapply(1:2, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- ", vec_print[i]))
    .construct_frnn(embedding[[i]], radius = vec_rad[i], nn = nn, 
                    frnn_approx = frnn_approx, 
                    resolve_isolated_nodes = T,
                    radius_quantile = NA,
                    verbose = verbose)
  }) 
  
  for(i in 1:2){
    if(verbose) print(paste0(Sys.time(),": cLISI: Converting into matrix -- ", vec_print[i]))
    list_g[[i]] <- .nnlist_to_matrix(list_g[[i]], set_to_one = F)
    if(symmetrize){
      list_g[[i]] <- .symmetrize_sparse(list_g[[i]], set_ones = F)
    }
    
    # convert back to list form if needed
    if(bool_matrix){
      if(length(rownames(obj$common_score)) != 0){
        rownames(list_g[[i]]) <- rownames(obj$common_score)
        colnames(list_g[[i]]) <- rownames(obj$common_score)
      }
    } else {
      # [[note to self: add a test to make sure this conversion is bijective]]
      list_g[[i]] <- .matrix_to_nnlist(list_g[[i]])
    }
  }
  
  # [[note to self: I'm not sure if we need to output e_g]]
  return(list(c_g = list_g[[1]], d_g = list_g[[2]],
              membership_vec = membership_vec,
              original_radius = vec_rad_org))
}

#' Combining two frNN graphs
#'
#' @param dcca_obj output of \code{dcca_decomposition}
#' @param g_1 graph object of class \code{dgCMatrix} where the non-zero entries represent distances
#' between nearest cells
#' @param g_2 graph object of class \code{dgCMatrix} where the non-zero entries represent distances
#' between nearest cells
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param common_1 boolean
#' @param common_2 boolean
#' @param keep_nn boolean
#' @param sampling_type string
#' @param center boolean
#' @param renormalize boolean
#' @param verbose numeric from 0 to 2
#'
#' @return object of class \code{dgCMatrix}
#' @export
combine_frnn <- function(dcca_obj, 
                         g_1, 
                         g_2, 
                         nn, 
                         common_1 = T, 
                         common_2 = T, 
                         keep_nn = T,
                         sampling_type = "adaptive_gaussian",
                         center = T, 
                         renormalize = F,
                         verbose = 0){
  stopifnot(all(dim(g_1) == dim(g_2)))
  
  # extract the relevant embeddings from dcca_obj
  if(verbose > 0) print(paste0(Sys.time(),": Preparing first embedding"))
  embedding_1 <- .prepare_embeddings(dcca_obj, data_1 = T, data_2 = F, 
                                     center = center, 
                                     renormalize = renormalize)
  if(common_1){
    embedding_1 <- embedding_1$common
  } else {
    embedding_1 <- embedding_1$distinct
  }
  
  if(verbose > 0) print(paste0(Sys.time(),": Preparing second embedding"))
  embedding_2 <- .prepare_embeddings(dcca_obj, data_1 = F, data_2 = T, 
                                     center = center, 
                                     renormalize = renormalize)
  if(common_2){
    embedding_2 <- embedding_2$common 
  } else {
    embedding_2 <- embedding_2$distinct
  }
  
  # symmetrize g_1 and g_2
  if(verbose > 0) print(paste0(Sys.time(),": Symmetrizing"))
  g_1 <- .symmetrize_sparse(g_1, set_ones = F)
  g_2 <- .symmetrize_sparse(g_2, set_ones = F)
  
  # prepare
  if(verbose > 0) print(paste0(Sys.time(),": Converting matrices to list"))
  n <- nrow(g_1)
  nn_idx_1 <- lapply(1:n, function(j){.nonzero_col(g_1, j, bool_value = F)})
  nn_dist_1 <- lapply(1:n, function(j){.nonzero_col(g_1, j, bool_value = T)})
  nn_idx_2 <- lapply(1:n, function(j){.nonzero_col(g_2, j, bool_value = F)})
  nn_dist_2 <- lapply(1:n, function(j){.nonzero_col(g_2, j, bool_value = T)})
  
  # apply the following procedure for each cell n
  nn_idx_all <- vector("list", n); nn_dist_all <- vector("list", n)
  for(i in 1:n){
    if(verbose == 2) {
      print(i)
    } else if(verbose == 1 && n > 10 && i %% floor(n/10) == 0) {
      cat('*')
    }
    
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
    # sample these entries
    tmp <- .embedding_resampling(nn_idx_all, nn_dist_all, nn = nn, 
                                 sampling_type = sampling_type, keep_nn = F)
    nn_idx_all <- tmp$id; nn_dist_all <- tmp$dist
    
    # find the nn's
    if(keep_nn){
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
  }
  
  tmp_list <- list(id = nn_idx_all, dist = nn_dist_all)
  res <- .nnlist_to_matrix(tmp_list, set_to_one = F)
  res <- .symmetrize_sparse(res, set_ones = F)
  
  if(length(rownames(dcca_obj$common_score)) != 0){
    rownames(res) <- rownames(dcca_obj$common_score)
    colnames(res) <- rownames(dcca_obj$common_score)
  }
  
  res
}

########################

.normalize_embeddings <- function(embedding, normalization_type){
  stopifnot(normalization_type %in% c("cosine_itself", "cosine_everything",
                                      "signac_itself", "signac_everything",
                                      "none"),
            sort(names(embedding)) == sort(c("common", "distinct", "everything")))
  
  len <- length(embedding)
  
  if(normalization_type == "cosine_everything"){
    l2_vec <- apply(embedding[["everything"]], 1, .l2norm)
    
    for(i in 1:len){
      embedding[[i]] <-.mult_vec_mat(1/l2_vec, embedding[[i]])
    }
    
  } else if(normalization_type == "cosine_itself"){
    for(i in 1:len){
      l2_vec <- apply(embedding[[i]], 1, .l2norm)
      embedding[[i]] <- .mult_vec_mat(1/l2_vec, embedding[[i]])
    }
    
  } else if(normalization_type == "signac_everything"){
    center_vec <- matrixStats::colMeans2(embedding[["everything"]])
    sd_vec <- matrixStats::colSds(embedding[["everything"]])
    embedding[["everything"]] <- sweep(embedding[["everything"]], 
                                       MARGIN = 2, 
                                       STATS = center_vec, 
                                       FUN = "-")
    embedding[["everything"]] <- sweep(embedding[["everything"]], 
                                       MARGIN = 2, 
                                       STATS = sd_vec, 
                                       FUN = "/")
    l2_vec <- apply(embedding[["everything"]], 1, .l2norm)
    embedding[["everything"]] <- .mult_vec_mat(1/l2_vec, embedding[["everything"]])
    
    for(i in c("common", "distinct")){
      embedding[[i]] <- sweep(embedding[[i]], 
                              MARGIN = 2, 
                              STATS = center_vec, 
                              FUN = "-")
      embedding[[i]] <- sweep(embedding[[i]], 
                              MARGIN = 2, 
                              STATS = sd_vec, 
                              FUN = "/")
      embedding[[i]] <- .mult_vec_mat(1/l2_vec, embedding[[i]])
    }
    
  } else if(normalization_type == "signac_itself"){
    for(i in 1:len){
      center_vec <- matrixStats::colMeans2(embedding[[i]])
      sd_vec <- matrixStats::colSds(embedding[[i]])
      embedding[[i]] <- sweep(embedding[[i]], 
                              MARGIN = 2, 
                              STATS = center_vec, 
                              FUN = "-")
      embedding[[i]] <- sweep(embedding[[i]], 
                              MARGIN = 2, 
                              STATS = sd_vec, 
                              FUN = "/")
      l2_vec <- apply(embedding[[i]], 1, .l2norm)
      embedding[[i]] <- .mult_vec_mat(1/l2_vec, embedding[[i]])
    }
    
  } 
  
  embedding
}

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

.construct_frnn <- function(mat, 
                            radius, 
                            nn, 
                            frnn_approx, 
                            resolve_isolated_nodes,
                            radius_quantile, 
                            verbose = F){
  if(is.na(radius)){
    radius <- .compute_radius(mat, 
                              nn = nn, 
                              radius_quantile = radius_quantile)
  }
  
  if(verbose) print(paste0(Sys.time(),": Computing frNN"))
  res <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)
  
  if(resolve_isolated_nodes){
    idx <- which(sapply(res$id, length) < nn)
    if(verbose) print(paste0(Sys.time(),": ", length(idx), " cells with too few neighbors"))
    if(length(idx) > 0){
      res2 <- RANN::nn2(mat, query = mat[idx,,drop = F], k = nn+1, eps = frnn_approx)
      
      if(verbose) print(paste0(Sys.time(),": Plugging in kNN neighbors"))
      for(i in 1:length(idx)){
        tmp_id <- c(res$id[[idx[i]]], res2$nn.idx[i,-1])
        tmp_dist <- c(res$dist[[idx[i]]], res2$nn.dists[i,-1])
        duplicates <- !duplicated(tmp_id)
        
        res$id[[idx[i]]] <- tmp_id[duplicates]
        res$dist[[idx[i]]] <- tmp_dist[duplicates]
      }
    }
  }
  
  res
}

.nnlist_to_matrix <- function(rann_obj, set_to_one){
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
  if(set_to_one){
    x_vec <- rep(1, length = length(i_vec))
  } else {
    x_vec <- unlist(rann_obj$dist)
  }
  
  Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)
}

.matrix_to_nnlist <- function(g_mat){
  n <- nrow(g_mat)
  nn_idx <- lapply(1:n, function(j){.nonzero_col(g_mat, j, bool_value = F)})
  nn_dist <- lapply(1:n, function(j){.nonzero_col(g_mat, j, bool_value = T)})
  
  structure(list(id = nn_idx, dist = nn_dist), class = "frNN")
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
      tmp <- tmp + (nn_dist_1_vec[z1])^2
    } else {
      tmp <- tmp + .l2norm(embedding_1[start_idx,] - embedding_1[j,])^2
    }
    
    z2 <- which(nn_idx_2_vec == j)
    if(length(z2) == 1) {
      tmp <- tmp + (nn_dist_2_vec[z2])^2
    } else {
      tmp <- tmp + .l2norm(embedding_2[start_idx,] - embedding_2[j,])^2
    }
    
    sqrt(tmp)
  })
  
  res
}

.symmetrize_sparse <- function(g_mat, set_ones){
  stopifnot(inherits(g_mat, "dgCMatrix"))
  
  tmp <- Matrix::t(g_mat)
  g_mat <- g_mat + tmp - sqrt(g_mat * tmp)
  if(set_ones) g_mat@x <- rep(1, length(g_mat@x))
  
  g_mat
}
