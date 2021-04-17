#' Compute cLISI information
#'
#' @param common_mat output from \code{dcca_decomposition}
#' @param distinct_mat output from \code{dcca_decomposition}
#' @param membership_vec factor
#' @param rank_c rank of \code{common_mat}
#' @param rank_d rank of \code{distinct_mat}
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param ensure_connected boolean
#' @param frnn_approx small non-negative number
#' @param radius_quantile value between 0 and 1
#' @param subsampling_rate value between 0 and 1, used for ensuring frNN graph is connected
#' @param min_subsample positive integer, used for ensuring frNN graph is connected
#' @param min_subsample_cell positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param verbose boolean
#'
#' @return two lists
#' @export
clisi_information <- function(common_mat, distinct_mat,
                              membership_vec, rank_c, rank_d, nn, 
                              frnn_approx = 0, radius_quantile = 0.95,
                              ensure_connected = F, 
                              subsampling_rate = 0.1, min_subsample = 50, 
                              min_subsample_cell = 50,
                              verbose = T){
  
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(common_mat))
  stopifnot(all(dim(common_mat) == dim(distinct_mat)),
            rank_c <= ncol(common_mat), rank_d <= ncol(distinct_mat))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1,
            subsampling_rate >= 0, subsampling_rate <= 1)
  
  everything_mat <- common_mat + distinct_mat
  
  # compute the 3 matrices
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing SVD -- common"))
  tmp <- .svd_truncated(common_mat, K = rank_c, 
                        symmetric = F, rescale = F, K_full_rank = T)
  c_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing SVD -- distinct"))
  tmp <- .svd_truncated(distinct_mat, K = rank_d, 
                        symmetric = F, rescale = F, K_full_rank = T)
  d_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing SVD -- everything"))
  tmp <- .svd_truncated(everything_mat, K = rank_d, 
                        symmetric = F, rescale = F, K_full_rank = T)
  e_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  
  # compute the radius
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- common"))
  c_rad <- .compute_radius(c_embedding, nn, radius_quantile)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- distinct"))
  d_rad <- .compute_radius(d_embedding, nn, radius_quantile)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- everything"))
  e_rad <- .compute_radius(e_embedding, nn, radius_quantile)
  sub_rad <- max(c_rad, d_rad)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- common"))
  c_g <- .construct_frnn(c_embedding, radius = sub_rad,
                         frnn_approx = frnn_approx,
                         ensure_connected = ensure_connected,
                         subsampling_rate = subsampling_rate,
                         min_subsample = min_subsample, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- distinct"))
  d_g <- .construct_frnn(d_embedding, radius = sub_rad,
                         frnn_approx = frnn_approx,
                         ensure_connected = ensure_connected,
                         subsampling_rate = subsampling_rate,
                         min_subsample = min_subsample, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- everything"))
  e_g <- .construct_frnn(e_embedding, radius = e_rad,
                         frnn_approx = frnn_approx,
                         ensure_connected = ensure_connected,
                         subsampling_rate = subsampling_rate,
                         min_subsample = min_subsample, verbose = verbose)
  
  cell_subidx <- .construct_celltype_subsample(membership_vec, min_subsample_cell)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- common"))
  c_score <- .clisi(c_g, membership_vec, cell_subidx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- distinct"))
  d_score <- .clisi(d_g, membership_vec, cell_subidx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- everything"))
  e_score <- .clisi(e_g, membership_vec, cell_subidx, verbose = verbose)
  
  list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score)
}

#############

.compute_radius <- function(mat, nn, radius_quantile){
  res <- RANN::nn2(mat, k = nn)
  as.numeric(stats::quantile(res$nn.dists[,nn], probs = radius_quantile))[1]
}

.construct_frnn <- function(mat, radius, frnn_approx, ensure_connected, subsampling_rate, min_subsample,
                            debug = F, verbose = F){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing frNN graph"))
  frnn_obj <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)$id
  # if(debug) return(frnn_obj)
  # 
  # for(i in 1:length(frnn_obj)){
  #   if(length(frnn_obj[[i]]) == 0) frnn_obj[[i]] <- i
  # }
  # 
  # if(verbose) print(paste0(Sys.time(),": cLISI: Converting frNN into igraph"))
  # g <- .convert_frnn2igraph(frnn_obj, verbose = verbose)
  # 
  # if(ensure_connected){
  #   if(verbose) print(paste0(Sys.time(),": cLISI: Connecting igraph"))
  #   .connect_graph(g, mat, subsampling_rate, min_subsample, verbose = verbose)
  # } else {
  #   if(verbose) {
  #     print(paste0(Sys.time(),": cLISI: Computing number of connected components"))
  #     res <- igraph::components(g)
  #     print(paste0(Sys.time(),": cLISI: Number of connected components: ", res$no))
  #   }
  #   
  #   g
  # }
  
}

.convert_frnn2igraph <- function(frnn_obj, verbose = F){
  n <- length(frnn_obj)
  g <- igraph::graph.empty(n = n, directed = F)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Constructing massive list of edges"))
  edge_mat <- do.call(cbind, lapply(1:n, function(i){
    rbind(i, frnn_obj[[i]])
  }))
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Adding edges into igraph"))
  g <- igraph::add_edges(g, edges = edge_mat)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Simplifying graph"))
  g <- igraph::simplify(g)
  
  g
}

.connect_graph <- function(g, mat, subsampling_rate, min_subsample,
                           verbose = F){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing number of components"))
  res <- igraph::components(g)
  
  ## [[note to self: this could definitely be improved]]
  while(res$no > 1){
    anchor_id <- which.max(res$csize)
    if(verbose) print(paste0(Sys.time(),": Number of components: ", res$no))
    k <- res$no
    id_list <- lapply(1:k, function(x){
      tmp <- which(res$membership == x)
      len <- length(tmp)
      if(len >= min_subsample){
        tmp <- tmp[sample(1:len, size = max(ceiling(len*subsampling_rate), min_subsample), 
                          replace = F)]
      }
      tmp
    })
    
    target_id <- id_list[[anchor_id]]
    source_id <- unlist(id_list[-anchor_id])
    
    if(verbose) print(paste0(Sys.time(),": cLISI: Computing nn2"))
    nn_res <- RANN::nn2(mat[target_id,,drop = F], query = mat[source_id,,drop = F], k = 1)
    tmp <- which.min(nn_res$nn.dists[,1])
    idx_from <- source_id[tmp]
    idx_to <- target_id[nn_res$nn.idx[tmp,1]]
    
    if(verbose) print(paste0(Sys.time(),": cLISI: Recomputing components"))
    g <- igraph::add_edges(g, c(idx_from, idx_to))
    res <- igraph::components(g)
  }
  
  igraph::simplify(g)
}

# .reduce_components <- function(res, idx1, idx2){
#   min_idx <- min(idx1, idx2)
#   max_idx <- max(idx1, idx2)
#   res$no <- res$no - 1
#   
# }

.construct_celltype_subsample <- function(membership_vec, min_subsample_cell){
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    if(length(idx) <= min_subsample_cell) return(idx)
    
    sample(idx, min_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}

.clisi <- function(g, membership_vec, cell_subidx, verbose = F){
  # stopifnot(all(table(membership_vec[cell_subidx]) > 0))
  # n <- igraph::vcount(g)
  # bg_prop <- as.numeric(table(membership_vec))/n
  # 
  # if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-wise cLISI"))
  # clisi_info <- sapply(1:length(cell_subidx), function(i){
  #   if(verbose && length(cell_subidx) > 10 && i %% floor(length(cell_subidx)/10) == 0) cat('*')
  #   
  #   neigh <- igraph::neighbors(g, v = cell_subidx[i])
  #   len <- length(neigh)
  #   mem_vec <- membership_vec[neigh]
  #   target_mem <- membership_vec[cell_subidx[i]]
  #   in_len <- length(which(mem_vec == target_mem))
  #   
  #   target_bg <- bg_prop[target_mem]
  #   clisi_score <- max((in_len/len - target_bg)/(1-target_bg), 0)
  #   
  #   c(len = len, in_ratio = in_len/len, clisi_score = clisi_score)
  # })
  
  n <- length(g)
  bg_prop <- as.numeric(table(membership_vec))/n
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-wise cLISI"))
  
  clisi_info <- sapply(1:length(cell_subidx), function(i){
    if(verbose && length(cell_subidx) > 10 && i %% floor(length(cell_subidx)/10) == 0) cat('*')
    neigh <- g[[cell_subidx[i]]]
    len <- length(neigh)
    if(len == 0){
      return(c(len = 0, in_ratio = 0, clisi_score = 0))
    }
    target_mem <- membership_vec[cell_subidx[i]]
    mem_vec <- membership_vec[neigh]
    in_len <- length(which(mem_vec == target_mem))
    
    target_bg <- bg_prop[target_mem]
    clisi_score <- max((in_len/len - target_bg)/(1-target_bg), 0)
    
    c(len = len, in_ratio = in_len/len, clisi_score = clisi_score)
  })
  
  # convert into df
  clisi_info <- data.frame(len = clisi_info["len",], 
                           in_ratio = clisi_info["in_ratio",], 
                           clisi_score = clisi_info["clisi_score",])
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-type cLISI"))
  res <- sapply(levels(membership_vec), function(x){
    idx <- which(membership_vec[cell_subidx] == x)
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