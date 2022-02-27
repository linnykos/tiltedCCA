#' Side-by-side UMAPs of the common, distinct and everything matrix
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec factor vector
#' @param data_1 boolean
#' @param data_2 boolean
#' @param col_vec vector of colors
#' @param pca boolean. If \code{TRUE}, plot the PCA embedding with the leading 2 components. 
#' If \code{FALSE}, plot the UMAP embedding.
#' @param only_embedding boolean
#' @param metric character for the distance metric used in \code{Seurat::RunUMAP},
#' typically \code{"cosine"} or \code{"euclidean"}
#' @param main_addition additional string to append to main of each plot
#' @param verbose boolean
#'
#' @return depends on \code{only_embedding}. 
#' If \code{TRUE}, returns three matrices as a list.
#' If \code{FALSE}, shows a plot but returns nothing
#' @export
plot_embeddings <- function(obj, 
                            membership_vec = NA, 
                            data_1 = T, 
                            data_2 = F, 
                            col_vec = scales::hue_pal()(length(levels(membership_vec))),
                            pca = F, 
                            only_embedding = F,
                            metric = "cosine",
                            main_addition = "", 
                            verbose = F){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(obj$common_score),
            class(obj) %in% c("dcca", "dcca_decomp"))
  stopifnot(data_1 | data_2)
  
  n <- nrow(obj$common_score)
  embedding <- .prepare_embeddings(obj, 
                                   data_1 = data_1, 
                                   data_2 = data_2, 
                                   center = F, 
                                   renormalize = T)
  if(pca) {
    embedding <- .extract_pca_embedding(embedding)
  } else {
    embedding <- .extract_umap_embedding(embedding, only_embedding = T, metric = metric)
  }
  
  if(!only_embedding) {
    if(all(is.na(membership_vec))){
      col_cells <- rep(col_vec[1], n)
    } else{
      col_cells <- col_vec[as.numeric(membership_vec)]
    }
    
    main_vec <- c("Common view", "Distinct view", "Everything view")
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
    xlim <- range(sapply(embedding, function(x){x[,1]}))
    ylim <- range(sapply(embedding, function(x){x[,2]}))
    n_idx <- sample(1:nrow(obj$common_score))
    
    graphics::par(mfrow = c(1,3))
    for(i in 1:3){
      graphics::plot(embedding[[i]][n_idx,1], embedding[[i]][n_idx,2],
                     asp = T, pch = 16, col = col_cells[n_idx], 
                     main = paste0(main_vec[i], main_addition),
                     xlab = label1, ylab = label2, xlim = xlim, ylim = ylim)
    }
    invisible()
  } else {
    return(embedding)
  }
}

#############################

#' Extract UMAP embedding
#'
#' @param embedding return object from \code{.prepare_embeddings}
#' @param only_embedding boolean
#' @param metric character for the distance metric used in \code{Seurat::RunUMAP},
#' typically \code{"cosine"} or \code{"euclidean"}
#' @param reduction_key string for \code{Seurat::RunUMAP}
#'
#' @return list of three 2-column matrix or \code{Seurat} object, one for the
#' common, distinct and everything matrix
.extract_umap_embedding <- function(embedding, only_embedding, metric, reduction_key = "UMAP"){
  stopifnot(length(embedding) == 3)
  for(i in 1:3){
    if(only_embedding){
      embedding[[i]] <- Seurat::RunUMAP(embedding[[i]], metric = metric, 
                                        verbose = F)@cell.embeddings
    } else {
      embedding[[i]] <- Seurat::RunUMAP(embedding[[i]], metric = metric, 
                                        reduction.key = reduction_key, verbose = F)
    }
  }
  
  embedding
}

.extract_pca_embedding <- function(embedding){
  stopifnot(length(embedding) == 3)
  for(i in 1:3){
    tmp <- .svd_truncated(embedding[[i]], K = 2, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    embedding[[i]] <- .mult_mat_vec(tmp$u, tmp$d)
  }
  
  embedding
}

###################################
