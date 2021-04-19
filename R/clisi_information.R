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
#' @param radius_quantile value between 0 and 1
#' @param max_subsample_frnn positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param max_subsample_clisi positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param verbose boolean
#'
#' @return two lists
#' @export
clisi_information <- function(common_mat, distinct_mat,
                              membership_vec, rank_c, rank_d, nn, 
                              frnn_approx = 0, radius_quantile = 0.9,
                              max_subsample_frnn = nrow(common_mat),
                              max_subsample_clisi = 50,
                              verbose = T){
  
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(common_mat))
  stopifnot(all(dim(common_mat) == dim(distinct_mat)),
            rank_c <= ncol(common_mat), rank_d <= ncol(distinct_mat))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1, max_subsample_clisi <= max_subsample_frnn)
  
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
  # [[note to self: the code below is quite messy -- fix it when i finalized the method]]
  cell_subidx1 <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
  n <- length(cell_subidx1)
  if(n < nrow(c_embedding)){
    c_embedding <- c_embedding[cell_subidx1,,drop = F]
    d_embedding <- d_embedding[cell_subidx1,,drop = F]
    e_embedding <- e_embedding[cell_subidx1,,drop = F]
    membership_vec <- membership_vec[cell_subidx1]
  }
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
  
  cell_subidx2 <- .construct_celltype_subsample(membership_vec, max_subsample_clisi)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- common"))
  c_score <- .clisi(c_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- distinct"))
  d_score <- .clisi(d_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- everything"))
  e_score <- .clisi(e_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  
  structure(list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score), class = "clisi")
}

#############

.compute_radius <- function(mat, nn, radius_quantile, cell_subidx){
  res <- RANN::nn2(mat[cell_subidx,,drop = F], k = nn)
  as.numeric(stats::quantile(res$nn.dists[,nn], probs = radius_quantile))[1]
}

.construct_frnn <- function(mat, radius, nn, frnn_approx, verbose = F){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing frNN"))
  res <- dbscan::frNN(mat, eps = radius, sort = F, approx = frnn_approx)$id
  
  idx <- which(sapply(res, length) < nn)
  if(verbose) print(paste0(Sys.time(),": cLISI: ", length(idx), " cells with too few neighbors"))
  if(length(idx) > 0){
    res2 <- RANN::nn2(mat, query = mat[idx,,drop = F], k = nn+1, eps = frnn_approx)
    
    if(verbose) print(paste0(Sys.time(),": cLISI: Plugging in kNN neighbors"))
    for(i in 1:length(idx)){
      res[[idx[i]]] <- unique(c(res[[idx[i]]], res2$nn.idx[i,-1]))
    }
  }
  
  res
}

.construct_celltype_subsample <- function(membership_vec, min_subsample_cell){
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    stopifnot(length(idx) > 2)
    
    if(length(idx) <= min_subsample_cell) return(idx)
    
    sample(idx, min_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}

.clisi <- function(g, membership_vec, cell_subidx, full_idx, verbose = F){
  stopifnot(is.list(g), is.factor(membership_vec), length(g) == length(full_idx))
  stopifnot(all(full_idx %% 1 == 0), all(full_idx > 0), length(full_idx) == length(unique(full_idx)))
  stopifnot(all(cell_subidx %% 1 == 0), all(cell_subidx > 0), length(cell_subidx) == length(unique(cell_subidx)),
            max(cell_subidx) <= length(full_idx))
  
  
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
  
  clisi_info <- as.data.frame(t(clisi_info))
  clisi_info <- cbind(idx = full_idx[cell_subidx], celltype = membership_vec[cell_subidx], 
                      clisi_info)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-type cLISI"))
  res <- sapply(levels(membership_vec), function(x){
    idx <- which(clisi_info$celltype == x)
    mean_vec <- apply(clisi_info[idx,c("len", "in_ratio", "clisi_score")], 2, stats::median)
    sd_vec <- apply(clisi_info[idx,c("len", "in_ratio", "clisi_score")], 2, function(x){
      diff(stats::quantile(x, probs = c(0.25, 0.75)))
    })
    
    c(mean_len = as.numeric(mean_vec["len"]), 
      mean_ratio = as.numeric(mean_vec["in_ratio"]),
      mean_clisi = as.numeric(mean_vec["clisi_score"]),
      sd_len = as.numeric(sd_vec["len"]), 
      sd_ratio = as.numeric(sd_vec["in_ratio"]),
      sd_clisi = as.numeric(sd_vec["clisi_score"]))
  })
  
  res <- as.data.frame(t(res))
  res <- cbind(celltype = levels(membership_vec), res)
  
  list(cell_info = clisi_info, membership_info = res)
}