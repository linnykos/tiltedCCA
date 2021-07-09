#' Plot embeddings via frNN
#'
#' @param c_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the common embedding, where the non-zero entries represent distances
#' @param d_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding, where the non-zero entries represent distances
#' @param e_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the everything embedding, where the non-zero entries represent distances
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
plot_embeddings2 <- function(c_g, d_g, e_g, nn, membership_vec = NA,
                             col_vec = scales::hue_pal()(length(levels(membership_vec))),
                             only_embedding = F, main_addition = "",
                             verbose = T, ...){
  list_g <- list(c_g = c_g, d_g = d_g, e_g = e_g)
  n <- nrow(c_g)
  list_output <- vector("list", 3)

  for(i in 1:3){
    nn_idx <- lapply(1:n, function(j){.nonzero_col(list_g[[i]], j, bool_value = F)})
    nn_dist <- lapply(1:n, function(j){.nonzero_col(list_g[[i]], j, bool_value = T)})
    
    # remove edges randomly
    for(j in 1:length(n)){
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
    
    graphics::par(mfrow = c(1,3))
    for(i in 1:3){
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