#' Compute cLISI information
#'
#' @param c_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the common embedding, where the non-zero entries represent distances
#' @param d_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding, where the non-zero entries represent distances
#' @param e_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the everything embedding, where the non-zero entries represent distances
#' @param membership_vec factor vector
#' @param max_subsample_clisi positive integer, used for determining number of cells to 
#' compute cLISI for
#' @param verbose boolean
#'
#' @return three lists
#' @export
clisi_information <- function(c_g, d_g, e_g, membership_vec, 
                              max_subsample_clisi = min(1000, nrow(c_g)),
                              verbose = T){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(c_g),
            all(dim(c_g) == dim(d_g)), all(dim(c_g) == dim(e_g)))
  n <- length(membership_vec)
  
  # remove all distance information and symmetrize
  if(verbose) print(paste0(Sys.time(),": cLISI: Symmetrizing matrices"))
  c_g <- .symmetrize_sparse(c_g, set_ones = T)
  d_g <- .symmetrize_sparse(d_g, set_ones = T)
  e_g <- .symmetrize_sparse(e_g, set_ones = T)
  
  # compute clisi scores
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_clisi)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- common"))
  c_score <- .clisi(c_g, membership_vec, cell_subidx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- distinct"))
  d_score <- .clisi(d_g, membership_vec, cell_subidx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- everything"))
  e_score <- .clisi(e_g, membership_vec, cell_subidx, verbose = verbose)
  
  structure(list(common_clisi = c_score, distinct_clisi = d_score,
       everything_clisi = e_score), class = "clisi")
}

############

.clisi <- function(g, membership_vec, cell_subidx, tol = 1e-3, verbose = F){
  stopifnot(is.factor(membership_vec))
  stopifnot(all(cell_subidx %% 1 == 0), all(cell_subidx > 0), length(cell_subidx) == length(unique(cell_subidx)))
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-wise cLISI"))
  
  clisi_info <- sapply(1:length(cell_subidx), function(i){
    if(verbose && length(cell_subidx) > 10 && i %% floor(length(cell_subidx)/10) == 0) cat('*')
    
    target_mem <- membership_vec[cell_subidx[i]]
    idx <- which(membership_vec == target_mem)
    .clisi_cell(g, idx, position = which(idx == i), tol = tol)
  })
  
  clisi_info <- as.data.frame(t(clisi_info))
  clisi_info <- cbind(idx = cell_subidx, celltype = membership_vec[cell_subidx], 
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

.clisi_cell <- function(g, idx, position, tol = 1e-3){
  stopifnot(position > 0, position <= length(idx), position %% 1 == 0)
  
  n <- nrow(g)
  target_bg <- length(idx)/n
  neigh <- .nonzero_col(g, idx[position], bool_value = F)
  len <- length(neigh)
  if(len == 0){
    return(c(len = 0, in_ratio = 0, clisi_score = 0))
  }
  in_len <- length(which(neigh %in% idx))
  
  clisi_score <- max((in_len/len - target_bg + tol)/(1-target_bg+tol), 0)
  c(len = len, in_ratio = in_len/len, clisi_score = clisi_score)
}