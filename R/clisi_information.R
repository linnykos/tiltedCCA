#' Compute cLISI information
#'
#' @param common_score output from \code{dcca_decomposition}
#' @param distinct_score output from \code{dcca_decomposition}
#' @param membership_vec factor
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
clisi_information <- function(common_score, distinct_score,
                              membership_vec, nn, 
                              frnn_approx = 0, radius_quantile = 0.9,
                              max_subsample_frnn = nrow(common_score),
                              max_subsample_clisi = min(500, nrow(common_score)),
                              verbose = T){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(common_score),
            max_subsample_clisi <= max_subsample_frnn)
  
  # compute the radius
  # [[note to self: the code below is quite messy -- fix it when i finalized the method]]
  cell_subidx1 <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
  if(length(cell_subidx1) < nrow(common_score)) membership_vec <- membership_vec[cell_subidx1]
  list_g <- construct_frnn(common_score[cell_subidx1,,drop = F], 
                           distinct_score[cell_subidx1,,drop = F],
                           nn, frnn_approx = frnn_approx, radius_quantile = radius_quantile,
                           bool_matrix = F, verbose = verbose)
  
  cell_subidx2 <- .construct_celltype_subsample(membership_vec, max_subsample_clisi)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- common"))
  c_score <- .clisi(list_g$c_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- distinct"))
  d_score <- .clisi(list_g$d_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- everything"))
  e_score <- .clisi(list_g$e_g, membership_vec, cell_subidx2, full_idx = cell_subidx1, verbose = verbose)
  
  structure(list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score), class = "clisi")
}

#' Construct fixed-radius NN graph
#'
#' @param common_score output from \code{dcca_decomposition}
#' @param distinct_score output from \code{dcca_decomposition}
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#'  @param frnn_approx small non-negative number
#' @param radius_quantile value between 0 and 1
#' @param bool_matrix boolean. If \code{TRUE}, output the graphs as a sparse matrix.
#' If \code{FALSE}, output the graphs as a list where each element of the list
#' corresponds with the element's neighbors
#' @param verbose boolean
#'
#' @return list, depends on \code{bool_matrix}
#' @export
construct_frnn <- function(common_score, distinct_score,
                           nn, frnn_approx = 0, radius_quantile = 0.9,
                           bool_matrix = F, verbose = T){
  
  stopifnot(nrow(common_score) == nrow(distinct_score))
  stopifnot(frnn_approx >= 0, frnn_approx <= 1)
  
  n <- nrow(common_score)
  everything_score <- .add_two_matrices(common_score, distinct_score)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- common"))
  c_rad <- .compute_radius(common_score, nn, radius_quantile, 1:n)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- distinct"))
  d_rad <- .compute_radius(distinct_score, nn, radius_quantile, 1:n)
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- everything"))
  e_rad <- .compute_radius(everything_score, nn, radius_quantile, 1:n)
  sub_rad <- max(c_rad, d_rad)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- common"))
  c_g <- .construct_frnn(common_score, radius = sub_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- distinct"))
  d_g <- .construct_frnn(distinct_score, radius = sub_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- everything"))
  e_g <- .construct_frnn(everything_score, radius = e_rad, nn = nn, 
                         frnn_approx = frnn_approx, verbose = verbose)
  
  if(bool_matrix){
    c_g <- .nnlist_to_matrix(c_g)
    d_g <- .nnlist_to_matrix(d_g)
    e_g <- .nnlist_to_matrix(e_g)
  }
  
  return(list(c_g = c_g, d_g = d_g, e_g = e_g))
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

.nnlist_to_matrix <- function(nn_list){
  j_vec <- unlist(nn_list)
  i_vec <- unlist(lapply(1:length(nn_list), function(i){rep(i, length(nn_list[[i]]))}))
  mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = 1)
  mat <- mat + Matrix::t(mat)
  mat@x <- rep(1, length(mat@i))
  mat
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