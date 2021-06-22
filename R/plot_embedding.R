#' Side-by-side UMAPs of the common, distinct and everything matrix
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec factor vector
#' @param data_1 boolean
#' @param data_2 boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct "on the same scale" by adding appropriately-scaled Gaussian noise
#' @param col_vec vector of colors
#' @param pca boolean. If \code{TRUE}, plot the PCA embedding with the leading 2 components. 
#' If \code{FALSE}, plot the UMAP embedding.
#' @param only_embedding boolean
#' @param main_addition additional string to append to main of each plot
#' @param verbose boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_embeddings <- function(obj, membership_vec, data_1 = T, data_2 = T, 
                            add_noise = T, 
                            col_vec = scales::hue_pal()(length(levels(membership_vec))),
                            pca = F, only_embedding = F,
                            main_addition = "", verbose = F){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(obj$common_score),
            class(obj) %in% c("dcca", "dcca_decomp"))
  
  stopifnot(data_1 | data_2)
  n_idx <- sample(1:nrow(obj$common_score))
  
  if(pca){
    embedding <- .plot_embeddings_pca(obj, data_1 = data_1, data_2 = data_2, col_vec = col_vec, 
                         main_addition = main_addition)
  } else {
    embedding <- .plot_embeddings_umap(obj, data_1 = data_1, data_2 = data_2, 
                                       add_noise = add_noise)
  }
  
  if(!only_embedding) {
    main_vec <- c("Common view", "Distinct view", "Everything view")
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
    xlim <- range(sapply(embedding, function(x){x[,1]}))
    ylim <- range(sapply(embedding, function(x){x[,2]}))
    
    graphics::par(mfrow = c(1,3))
    for(i in 1:3){
      graphics::plot(embedding[[i]][n_idx,1], embedding[[i]][n_idx,2],
                     asp = T, pch = 16, col = col_vec[as.numeric(membership_vec)][n_idx], 
                     main = paste0(main_vec[i], main_addition),
                     xlab = label1, ylab = label2, xlim = xlim, ylim = ylim)
    }
    invisible()
  } else {
    return(embedding)
  }
}

#################

.plot_embeddings_pca <- function(obj, data_1, data_2){
  stopifnot(!data_1 | !data_2) # [[note to self: restriction for now]]
  label1 <- "PCA 1"; label2 <- "PCA 2"
  embedding <- vector("list", 3)
  names(embedding) <- c("common", "distinct", "everything")
  if(data_1){
    distinct_score <- obj$distinct_score_1
    svd_res <- obj$svd_1
  } else {
    distinct_score <- obj$distinct_score_2
    svd_res <- obj$svd_2
  }

  embedding[[1]] <- .extract_matrix_helper(obj$common_score, distinct_score,
                                           svd_res, common_bool = T, distinct_bool = F,
                                           center = F, renormalize = F, add_noise = F)
  embedding[[2]] <- .extract_matrix_helper(obj$common_score, distinct_score,
                                           svd_res, common_bool = F, distinct_bool = T,
                                           center = F, renormalize = F, add_noise = F)
  embedding[[3]] <- .extract_matrix_helper(obj$common_score, distinct_score,
                                           svd_res, common_bool = T, distinct_bool = T,
                                           center = F, renormalize = F, add_noise = F)
  
  for(i in 1:3){
    tmp <- .svd_truncated(embedding[[i]], K = 2, symmetric = F, rescale = F, 
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
    embedding[[i]] <- .mult_mat_vec(tmp$u, tmp$d)
  }
  
  embedding
}

.plot_embeddings_umap <- function(obj, data_1, data_2, add_noise){
  if(verbose) print(paste0(Sys.time(),": Plotting: Preparing objects"))
  prep_obj <- .prepare_umap_embedding(obj)
  embedding <- vector("list", 3)
  names(embedding) <- c("common", "distinct", "everything")

  set.seed(10)
  if(verbose) print(paste0(Sys.time(),": Plotting: UMAP for common matrix"))
  embedding[[1]] <- .extract_umap_embedding(prep_obj, common_1 = data_1, common_2 = data_2, 
                                            distinct_1 = F, distinct_2 = F,
                                            add_noise = add_noise, only_embedding = T)
  set.seed(10)
  if(verbose) print(paste0(Sys.time(),": Plotting: UMAP for distinct matrix"))
  embedding[[2]] <- .extract_umap_embedding(prep_obj, common_1 = F, common_2 = F, 
                                            distinct_1 = data_1, distinct_2 = data_2,
                                            add_noise = add_noise, only_embedding = T)
  
  set.seed(10)
  if(verbose) print(paste0(Sys.time(),": Plotting: UMAP for entire matrix"))
  embedding[[3]] <- .extract_umap_embedding(prep_obj, common_1 = data_1, common_2 = data_2, 
                                            distinct_1 = data_1, distinct_2 = data_2,
                                            add_noise = add_noise, only_embedding = T)
  
  embedding
}