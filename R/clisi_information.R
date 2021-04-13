clisi_information <- function(common_mat, distinct_mat, everything_mat, 
                              membership_vec, p, nn, 
                              frnn_approx = 0, 
                              subsampling_rate = 0.1, min_subsample = 50){
  
  stopifnot(all(membership_vec %% 1 == 0), all(membership_vec >= 1),
            length(unique(membership_vec)) == max(membership_vec),
            max(membership_vec) <= nrow(common_mat))
  stopifnot(all(dim(common_mat) == dim(distinct_mat)),
            all(dim(common_mat) == dim(everything_mat)),
            p <= ncol(common_mat))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1,
            subsampling_rate >= 0, subsampling_rate <= 1)
  
  # compute the 3 matrices
  tmp <- .svd_truncated(common_mat, K = p, 
                        symmetric = F, rescale = F, K_full_rank = F)
  c_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  tmp <- .svd_truncated(distinct_mat, K = p, 
                        symmetric = F, rescale = F, K_full_rank = F)
  d_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  tmp <- .svd_truncated(everything_mat, K = p, 
                        symmetric = F, rescale = F, K_full_rank = F)
  e_embedding <- .mult_mat_vec(tmp$u, tmp$d)
  
  # compute the radius
  c_rad <- .compute_radius(c_embedding, nn)
  d_rad <- .compute_radius(d_embedding, nn)
  e_rad <- .compute_radius(e_embedding, nn)
  sub_rad <- max(c_rad, d_rad)
  
  c_g <- .construct_frnn(c_embedding, radius = sub_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample)
  d_g <- .construct_frnn(d_embedding, radius = sub_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample)
  e_g <- .construct_frnn(e_embedding, radius = e_rad,
                              frnn_approx = frnn_approx,
                              subsampling_rate = subsampling_rate,
                              min_subsample = min_subsample)
  
  c_score <- .clisi(c_g, membership_vec)
  d_score <- .clisi(d_g, membership_vec)
  e_score <- .clisi(e_g, membership_vec)
  
  list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score)
}

#############

.compute_radius <- function(mat, nn){
  res <- RANN::nn2(mat, k = nn)
  stats::median(res$nn.dists[,nn])
}

.construct_frnn <- function(mat, radius, frnn_approx, subsampling_rate, min_subsample,
                            debug = F){
  frnn_obj <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)$id
  if(debug) return(frnn_obj)
  
  for(i in 1:length(frnn_obj)){
    if(length(frnn_obj[[i]]) == 0) frnn_obj[[i]] <- i
  }
  
  g <- .convert_frnn2igraph(frnn_obj)
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
        tmp <- tmp[sample(1:len, size = ceiling(len*subsampling_rate), replace = F)]
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

.clisi <- function(g, membership_vec){
  n <- igraph::vcount(g)
  clisi_score <- sapply(1:n, function(i){
    neigh <- igraph::neighbors(g, v = i)
    len <- length(neigh)
    mem_vec <- membership_vec[neigh]
    target_mem <- membership_vec[i]
    in_val <- length(which(mem_vec == target_mem))
    
    in_val/len
  })
  
  k <- max(membership_vec)
  sapply(1:k, function(x){
    idx <- which(membership_vec == x)
    mean(clisi_score[idx])
  })
}