#' Plot embeddings via frNN
#'
#' @param c_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the common embedding, where the non-zero entries represent distances
#' @param d_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding, where the non-zero entries represent distances
#' @param e_g sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the everything embedding, where the non-zero entries represent distances
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
plot_embeddings2 <- function(c_g, d_g, e_g, membership_vec = NA,
                             col_vec = scales::hue_pal()(length(levels(membership_vec))),
                             only_embedding = F, main_addition = "",
                             verbose = T, ...){
  list_g <- list(c_g = c_g, d_g = d_g, e_g = e_g)
  n <- nrow(c_g)
  list_output <- vector("list", 3)

  for(i in 1:3){
    # symmetrize
    list_g[[i]] <- .symmetrize_sparse(list_g[[i]], set_ones = F)
  
    # convert into kernels
    list_g[[i]] <- .distance_to_kernel(list_g[[i]])
    list_g[[i]] <- SeuratObject::as.Graph(list_g[[i]])
    
    # use Seurat
    tmp <- Seurat::RunUMAP(list_g[[i]], verbose = verbose, ...)
    list_output[[i]] <- tmp@cell.embeddings
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

######################

.distance_to_kernel <- function(mat_g){
  stopifnot(inherits(mat_g, "dgCMatrix"))
  
  n <- nrow(mat_g)
  x_val <- unlist(lapply(1:n, function(col_idx){
    val1 <- mat_g@p[col_idx]+1
    val2 <- mat_g@p[col_idx+1]+1
    
    if(val1 == val2) return(numeric(0))
    vec <- mat_g@x[(val1+1):val2]
    min_val <- min(vec)
    max_val <- max(vec)
    exp(-(vec - min_val)/max_val)
  }))
  
  mat_g@x <- x_val
  mat_g
}