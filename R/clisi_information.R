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
clisi_information <- function(c_g, d_g, membership_vec, 
                              max_subsample_clisi = min(1000, nrow(c_g)),
                              verbose = T){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(c_g),
            all(dim(c_g) == dim(d_g)), all(dim(c_g) == dim(e_g)))
  n <- length(membership_vec)
  
  # remove all distance information and symmetrize
  if(verbose) print(paste0(Sys.time(),": cLISI: Symmetrizing matrices"))
  c_g <- .symmetrize_sparse(c_g, set_ones = T)
  d_g <- .symmetrize_sparse(d_g, set_ones = T)
  
  # compute clisi scores
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_clisi)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- common"))
  c_score <- .clisi(c_g, membership_vec, cell_subidx, verbose = verbose)
  if(verbose) print(paste0(Sys.time(),": cLISI: Compute cLISI -- distinct"))
  d_score <- .clisi(d_g, membership_vec, cell_subidx, verbose = verbose)
  
  structure(list(common_clisi = c_score, distinct_clisi = d_score), class = "clisi")
}

############

.clisi <- function(g, membership_vec, cell_subidx, tol = 1e-3, verbose = F){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(g))
  stopifnot(all(cell_subidx %% 1 == 0), all(cell_subidx > 0), 
            length(cell_subidx) == length(unique(cell_subidx)),
            length(cell_subidx) <= length(membership_vec))
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-wise cLISI"))
  
  clisi_mat <- sapply(1:length(cell_subidx), function(i){
    if(verbose && length(cell_subidx) > 10 && i %% floor(length(cell_subidx)/10) == 0) cat('*')
    
    .clisi_cell(g, membership_vec, position = i, tol = tol)
  })
  rownames(clisi_mat) <- levels(membership_vec)
  
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing cell-type cLISI"))
  res <- lapply(levels(membership_vec), function(celltype){
    idx <- which(membership_vec[cell_subidx] == celltype)
    tmp_mat <- clisi_mat[,idx,drop = F]
    mean_vec <- matrixStats::rowMedians(tmp_mat)
    names(mean_vec) <- rownames(tmp_mat)
    
    mean_val <- mean_vec[which(levels(membership_vec) == celltype)]
      
    c(vec = mean_vec, 
      value = mean_val)
  })
  names(res) <- levels(membership_vec)
  
  # reorganize values
  df <- data.frame(celltype = levels(membership_vec), value = sapply(res, function(x){x$value}))
  mat <- sapply(res, function(x){x$vec})
  colnames(mat) <- paste0("from_", levels(membership_vec))
  rownames(mat) <- paste0("to_", levels(membership_vec))
  
  list(df = df, clisi_mat = mat)
}

.clisi_cell <- function(g, membership_vec, position, tol = 1e-3){
  stopifnot(position > 0, position <= nrow(g), position %% 1 == 0,
            ncol(g) == nrow(g), is.factor(membership_vec), 
            length(membership_vec) == nrow(g))
  
  n <- nrow(g)
  celltype_vec <- levels(membership_vec)
  k <- length(celltype_vec)
  target_bg <- table(membership_vec)/n
  neigh <- .nonzero_col(g, position, bool_value = F)
  len <- length(neigh)
  if(len == 0){
    return(rep(0, k))
  }
  in_len <- sapply(celltype_vec, function(celltype){
    length(which(membership_vec[neigh] == celltype))
  })
    
  clisi_vec <- pmax((in_len/len-target_bg+tol) / (1-target_bg+tol), 0)
  names(clisi_vec) <- celltype_vec
    
  clisi_vec
}