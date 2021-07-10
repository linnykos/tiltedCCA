#' Plot embeddings via frNN
#'
#' @param c_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the common embedding, where the non-zero entries represent distances
#' @param d_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding, where the non-zero entries represent distances
#' @param nn integer of number of nearest neighbors to determine the appropriate radius
#' for the frNN graph
#' @param membership_vec factor vector
#' @param col_vec vector of colors
#' @param only_embedding boolean
#' @param main_addition additional string to append to main of each plot
#' @param verbose boolean
#' @param ... additional parameters for \code{Seurat:::RunUMAP.Graph}
#'
#' @return depends on \code{only_embedding}. 
#' If \code{TRUE}, returns three matrices as a list.
#' If \code{FALSE}, shows a plot but returns nothing
#' @export
plot_embeddings2 <- function(dcca_res, nn, data_1 = T, data_2 = F, c_g = NA, d_g = NA, 
                             membership_vec = NA,
                             col_vec = scales::hue_pal()(length(levels(membership_vec))),
                             only_embedding = F, main_addition = "",
                             verbose = T, ...){
  stopifnot(!data_1 | !data_2)
  
  if(all(is.na(c_g)) || all(is.na(d_g))){
    rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                             data_1 = data_1, data_2 = data_2,
                                             bool_matrix = T, include_diag = F, verbose = verbose)
    c_g <- rna_frnn$c_g; d_g <- rna_frnn$d_g
  }
  
  if(data_1){
    everything_embedding <- .mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d)
  } else if(data_2) {
    everything_embedding <- .mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d)
  }
 
  list_g <- list(c_g = c_g, d_g = d_g)
  n <- nrow(c_g)
  list_output <- vector("list", 3)

  # Use Seurat::RunUMAP.Graph for the common and distinct embeddings
  for(i in 1:length(list_g)){
    mat <- .symmetrize_sparse(list_g[[i]], set_ones = F)
    nn_idx <- lapply(1:n, function(j){.nonzero_col(mat, j, bool_value = F)})
    nn_dist <- lapply(1:n, function(j){.nonzero_col(mat, j, bool_value = T)})
    
    # remove edges randomly
    for(j in 1:n){
      if(length(nn_idx[[j]]) <= nn) next()
      idx <- sample(1:length(nn_idx[[j]]), size = nn)
      nn_idx[[j]] <- nn_idx[[j]][idx]
      nn_dist[[j]] <- nn_dist[[j]][idx]
    }
    
    rann_obj <- list(id = nn_idx, dist = nn_dist)
    mat <- .nnlist_to_matrix(rann_obj, include_diag = F)
    
    # symmetrize
    mat <- .symmetrize_sparse(mat, set_ones = F)
    if(length(rownames(list_g[[i]])) > 0){
      rownames(mat) <- rownames(list_g[[i]])
      colnames(mat) <- rownames(list_g[[i]])
    }
    graph_obj <- SeuratObject::as.Graph(mat)
    
    # use Seurat
    list_output[[i]]  <- Seurat::RunUMAP(graph_obj, verbose = verbose, assay = "RNA", ...)@cell.embeddings
  }
  
  # run Seurat::RunUMAP.Default on the everything
  list_output[[3]] <- Seurat::RunUMAP(everything_embedding, metric = "euclidean", 
                                      verbose = verbose, assay = "RNA", ...)@cell.embeddings
  
  if(!only_embedding) {
    if(all(is.na(membership_vec))){
      col_cells <- rep(col_vec[1], n)
    } else{
      col_cells <- col_vec[as.numeric(membership_vec)]
    }
    
    main_vec <- c("Common view", "Distinct view", "Everything view")
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
    xlim <- range(sapply(list_output, function(x){x[,1]}))
    ylim <- range(sapply(list_output, function(x){x[,2]}))
    n_idx <- sample(1:nrow(list_output))
    
    graphics::par(mfrow = c(1,length(list_g)))
    for(i in 1:length(list_g)){
      graphics::plot(list_output[[i]][n_idx,1], list_output[[i]][n_idx,2],
                     asp = T, pch = 16, col = col_vec[as.numeric(membership_vec)][n_idx], 
                     main = paste0(main_vec[i], main_addition),
                     xlab = label1, ylab = label2, xlim = xlim, ylim = ylim)
    }
    invisible()
  } else {
    return(list_output)
  }
}