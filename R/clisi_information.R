#' Compute cLISI information
#'
#' @param common_mat output from \code{dcca_decomposition}
#' @param distinct_mat output from \code{dcca_decomposition}
#' @param membership_vec factor
#' @param rank_c rank of \code{common_mat}
#' @param rank_d rank of \code{distinct_mat}
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param frnn_approx small non-negative number
#' @param subsampling_rate value between 0 and 1, used for ensuring frNN graph is connected
#' @param min_subsample positive integer, used for ensuring frNN graph is connected
#' @param subsampling_rate_cell value between 0 and 1, used for determining number of cells to 
#' compute cLISI for
#' @param min_subsample_cell positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param verbose boolean
#'
#' @return two lists
#' @export
clisi_information <- function(common_mat, distinct_mat,
                              membership_vec, rank_c, rank_d, nn, 
                              frnn_approx = 0, 
                              subsampling_rate = 0.1, min_subsample = 50, 
                              subsampling_rate_cell = 0.1, 
                              min_subsample_cell = 50,
                              verbose = T){
  
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(common_mat))
  # stopifnot(all(membership_vec %% 1 == 0), all(membership_vec >= 1),
  #           length(unique(membership_vec)) == max(membership_vec),
  #           max(membership_vec) <= nrow(common_mat))
  stopifnot(all(dim(common_mat) == dim(distinct_mat)),
            rank_c <= ncol(common_mat), rank_d <= ncol(distinct_mat))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1,
            subsampling_rate >= 0, subsampling_rate <= 1)
  
  everything_mat <- common_mat + distinct_mat
  
  # compute the 3 matrices
  if(verbose) print("cLISI: Computing SVD -- common")
  tmp <- .svd_truncated(common_mat, K = rank_c, 
                        symmetric = F, rescale = F, K_full_rank = T)
  c_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  if(verbose) print("cLISI: Computing SVD -- distinct")
  tmp <- .svd_truncated(distinct_mat, K = rank_d, 
                        symmetric = F, rescale = F, K_full_rank = T)
  d_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  if(verbose) print("cLISI: Computing SVD -- everything")
  tmp <- .svd_truncated(everything_mat, K = rank_d, 
                        symmetric = F, rescale = F, K_full_rank = T)
  e_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  
  # compute the radius
  if(verbose) print("cLISI: Computing radius -- common")
  c_rad <- .compute_radius(c_embedding, nn)
  if(verbose) print("cLISI: Computing radius -- distinct")
  d_rad <- .compute_radius(d_embedding, nn)
  if(verbose) print("cLISI: Computing radius -- everything")
  e_rad <- .compute_radius(e_embedding, nn)
  sub_rad <- max(c_rad, d_rad)
  
  if(verbose) print("cLISI: Construct graph -- common")
  c_g <- .construct_frnn(c_embedding, radius = sub_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample, verbose = verbose)
  if(verbose) print("cLISI: Construct graph -- distinct")
  d_g <- .construct_frnn(d_embedding, radius = sub_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample, verbose = verbose)
  if(verbose) print("cLISI: Construct graph -- everything")
  e_g <- .construct_frnn(e_embedding, radius = e_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample, verbose = verbose)
  
  cell_subidx <- .construct_celltype_subsample(membership_vec, subsampling_rate_cell, min_subsample_cell)
  if(verbose) print("cLISI: Compute cLISI -- common")
  c_score <- .clisi(c_g, membership_vec, cell_subidx)
  if(verbose) print("cLISI: Compute cLISI -- distinct")
  d_score <- .clisi(d_g, membership_vec, cell_subidx)
  if(verbose) print("cLISI: Compute cLISI -- everything")
  e_score <- .clisi(e_g, membership_vec, cell_subidx)
  
  list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score)
}

#############

.compute_radius <- function(mat, nn){
  res <- RANN::nn2(mat, k = nn)
  stats::median(res$nn.dists[,nn])
}

.construct_frnn <- function(mat, radius, frnn_approx, subsampling_rate, min_subsample,
                            debug = F, verbose = F){
  if(verbose) print("cLISI: Computing frNN graph")
  frnn_obj <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)$id
  if(debug) return(frnn_obj)
  
  for(i in 1:length(frnn_obj)){
    if(length(frnn_obj[[i]]) == 0) frnn_obj[[i]] <- i
  }
  
  if(verbose) print("cLISI: Converting frNN into igraph")
  g <- .convert_frnn2igraph(frnn_obj)
  
  if(verbose) print("cLISI: Connecting igraph")
  .connect_graph(g, mat, subsampling_rate, min_subsample)
}

.convert_frnn2igraph <- function(frnn_obj){
  n <- length(frnn_obj)
  g <- igraph::graph.empty(n = n, directed = F)
  for(i in 1:n){
    stopifnot(length(frnn_obj[[i]]) > 0)
    
    edge_mat <- rbind(i, frnn_obj[[i]])
    g <- igraph::add_edges(g, edges = edge_mat)
  }
  g <- igraph::simplify(g)
  
  g
}

.connect_graph <- function(g, mat, subsampling_rate, min_subsample){
  res <- igraph::components(g)
  
  ## [[note to self: this could definitely be improved]]
  while(res$no > 1){
    anchor_id <- which.max(res$csize)
    k <- res$no
    id_list <- lapply(1:k, function(x){
      tmp <- which(res$membership == x)
      len <- length(tmp)
      if(len >= min_subsample){
        tmp <- tmp[sample(1:len, size = max(ceiling(len*subsampling_rate), min_subsample), replace = F)]
      }
      tmp
    })
    
    target_id <- id_list[[anchor_id]]
    source_id <- unlist(id_list[-anchor_id])
    
    nn_res <- RANN::nn2(mat[target_id,,drop = F], query = mat[source_id,,drop = F], k = 1)
    tmp <- which.min(nn_res$nn.dists[,1])
    idx_from <- source_id[tmp]
    idx_to <- target_id[nn_res$nn.idx[tmp,1]]
    
    g <- igraph::add_edges(g, c(idx_from, idx_to))
    res <- igraph::components(g)
  }
  
  igraph::simplify(g)
}

.construct_celltype_subsample <- function(membership_vec, subsampling_rate_cell, 
                                          min_subsample_cell){
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    if(length(idx) <= min_subsample_cell) return(idx)
    
    sample(idx, max(ceiling(subsampling_rate_cell*length(idx)), min_subsample_cell), replace = F)
  })
  
  sort(unlist(res))
}

.clisi <- function(g, membership_vec, cell_subidx){
  stopifnot(all(table(membership_vec[cell_subidx]) > 0))
  n <- igraph::vcount(g)
  bg_prop <- as.numeric(table(membership_vec))/n
  
  clisi_info <- sapply(cell_subidx, function(i){
    neigh <- igraph::neighbors(g, v = i)
    len <- length(neigh)
    mem_vec <- membership_vec[neigh]
    target_mem <- membership_vec[i]
    in_len <- length(which(mem_vec == target_mem))
    
    target_bg <- bg_prop[target_mem]
    clisi_score <- max((in_len/len - target_bg)/(1-target_bg), 0)
    
    c(len = len, in_ratio = in_len/len, clisi_score = clisi_score)
  })
  
  # convert into df
  clisi_info <- data.frame(len = clisi_info["len",], 
                           in_ratio = clisi_info["in_ratio",], 
                           clisi_score = clisi_info["clisi_score",])
  
  res <- sapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    mean_vec <- colMeans(clisi_info[idx,])
    sd_vec <- apply(clisi_info[idx,], 2, stats::sd)
    
    c(celltype = x,
      mean_len = as.numeric(mean_vec["len"]), 
      mean_ratio = as.numeric(mean_vec["in_ratio"]),
      mean_clisi = as.numeric(mean_vec["clisi_score"]),
      sd_len = as.numeric(sd_vec["len"]), 
      sd_ratio = as.numeric(sd_vec["in_ratio"]),
      sd_clisi = as.numeric(sd_vec["clisi_score"]))
  })
  
  res <- as.data.frame(t(res))
  
  list(cell_info = clisi_info, membership_info = res)
}