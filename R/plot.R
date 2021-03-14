#' Heatmap of the data
#' 
#' If \code{reserve_zero = T}, then reserve zero for white.
#' If \code{reserve_zero = F}, then all the greens are negative and all the reds are positive
#'
#' @param dat matrix
#' @param luminosity boolean
#' @param asp numeric
#' @param reserve_zero boolean
#' @param ... additional graphical parameters
#'
#' @return shows a plot but returns nothing
#' @export
plot_heatmat <- function(dat, luminosity = T, asp = nrow(dat)/ncol(dat), 
                         reserve_zero = T, ...){
  if(reserve_zero){
    col_vec <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19,
                                 luminosity = luminosity)
    
    col_vec <- c("white", col_vec)
    
    tmp <- as.numeric(dat)
    tmp <- tmp[tmp!=0]
    
    break_vec <- stats::quantile(tmp, probs = seq(0, 1, length.out = 20))
    break_vec <- c(-5, break_vec)
  } else {
    # green is negative
    col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), 5,
                                 luminosity = luminosity)
    tmp <- as.numeric(dat[dat < 0])
    break_vec_neg <- rev(-stats::quantile(abs(tmp), probs = 1-(1-seq(0, 1, length.out = 6))^1.5))
    
    # red is positive
    col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), 5,
                                     luminosity = luminosity)
    tmp <- as.numeric(dat[dat > 0])
    break_vec_pos <- stats::quantile(tmp, probs = 1-(1-seq(0, 1, length.out = 6))^1.5)
    
    # combine the two
    break_vec <- c(break_vec_neg, break_vec_pos)
    col_vec <- c(col_vec_neg, "white", col_vec_pos)
  }
 
  graphics::image(.rotate(dat), breaks = break_vec, 
                  col = col_vec, asp = asp, ...)
  
  invisible()
}

#' Side-by-side plot of the canonical scores, colored by membership
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec integer vector
#' @param decomposition boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_scores <- function(obj, membership_vec, decomposition = F){
  if(decomposition){
    score_list <- vector("list", 3)
    
    score_list[[1]] <- obj$common_score
    score_list[[2]] <- obj$distinct_score_1
    score_list[[3]] <- obj$distinct_score_2
    
  } else {
    score_list <- vector("list", 2)
    
    score_list[[1]] <- obj$common_score
    if(ncol(score_list[[1]]) < ncol(obj$distinct_score_1)){
      score_list[[1]] <- cbind(score_list[[1]], matrix(0, nrow = nrow(score_list[[1]]), ncol = ncol(obj$distinct_score_1) - ncol(score_list[[1]])))
    }
    score_list[[1]] <- score_list[[1]]+obj$distinct_score_1
    
    score_list[[2]] <- obj$common_score
    if(ncol(score_list[[2]]) < ncol(obj$distinct_score_2)){
      score_list[[2]] <- cbind(score_list[[2]], matrix(0, nrow = nrow(score_list[[2]]), ncol = ncol(obj$distinct_score_2) - ncol(score_list[[2]])))
    }
    score_list[[2]] <- score_list[[2]]+obj$distinct_score_2
  }
 
  n <- nrow(score_list[[1]])
  max_col <- max(sapply(score_list, ncol)); xlim <- range(do.call(rbind, score_list))
  
  if(decomposition){
    main_vec <- c("Common", "Distinct 1", "Distinct 2")
  } else {
    main_vec <- c("Dataset 1", "Dataset 2")
  }
  
  graphics::par(mfrow = c(1, length(score_list)))
  for(k in 1:length(score_list)){
    graphics::plot(NA, xlim = xlim, ylim = c(0.5, max_col+.5), main = main_vec[k], xlab = "Value", ylab = "Dimension")
    graphics::lines(rep(0, 2), max_col*c(-10,10), col = "red", lty = 2)
    
    for(i in 1:ncol(score_list[[k]])){
      graphics::points(x = score_list[[k]][,i], y = stats::runif(n, min = i-.2, max = i+.2), col = membership_vec, pch = 16)
    }
  }
  
  invisible()
}


#' Side-by-side plot of the canonical scores as heatmaps
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec integer vector
#' @param num_col positive integers for number of distinct colors
#' @param luminosity boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_scores_heatmap <- function(obj, membership_vec = NA, num_col = 10, luminosity = F){
  n <- nrow(obj$common_score)
  zlim <- range(c(obj$common_score, obj$distinct_score_1, obj$distinct_score_2))
  
  # construct colors. green is negative
  max_val <- max(abs(zlim))
  col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), num_col,
                                   luminosity = luminosity)
  break_vec_neg <- seq(-max_val, 0, length.out = num_col+2)
  break_vec_neg <- break_vec_neg[-length(break_vec_neg)]
  
  # red is positive
  col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), num_col,
                                   luminosity = luminosity)
  break_vec_pos <- seq(0, max_val, length.out = num_col+2)
  break_vec_pos <- break_vec_pos[-1]
  
  # combine the two
  break_vec <- c(break_vec_neg, break_vec_pos)
  col_vec <- c(col_vec_neg, "white", col_vec_pos)
  
  if(!all(is.na(membership_vec))){
    stopifnot(all(membership_vec %% 1 == 0), all(membership_vec > 0),
              max(membership_vec) == length(unique(membership_vec)))
    
    idx <- order(membership_vec, decreasing = F)
    breakpoints <- 1-which(abs(diff(sort(membership_vec, decreasing = F))) >= 1e-6)/n
  } else {
    idx <- 1:n
  }
  
  line_func <- function(){
    if(!all(is.na(membership_vec))){
      for(i in 1:length(breakpoints)){
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.1, col = "white")
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
      }
    }
  }
  
  graphics::par(mfrow = c(1,3))
  graphics::image(.rotate(obj$common_score[idx,,drop = F]), main = "Common score",
                  col = col_vec, breaks = break_vec)
  line_func()
  
  graphics::image(.rotate(obj$distinct_score_1[idx,,drop = F]), main = "Distinct score 1",
                  col = col_vec, breaks = break_vec)
  line_func()
  
  graphics::image(.rotate(obj$distinct_score_2[idx,,drop = F]), main = "Distinct score 2",
                  col = col_vec, breaks = break_vec)
  line_func()
  
  invisible()
}

#' Side-by-side UMAPs of the common, distinct and everything matrix
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec integer vector
#' @param data_1 boolean
#' @param data_2 boolean
#' @param add_noise boolean, intended (if \code{TRUE}) to put the common and 
#' distinct "on the same scale" by adding appropriately-scaled Gaussian noise
#' @param pca boolean. If \code{TRUE}, plot the PCA embedding with the leading 2 components. 
#' If \code{FALSE}, plot the UMAP embedding.
#'
#' @return shows a plot but returns nothing
#' @export
plot_embeddings <- function(obj, membership_vec, data_1 = T, data_2 = T, add_noise = T, pca = F){
  stopifnot(data_1 | data_2)
  if(pca){
    label1 <- "PCA 1"; label2 <- "PCA 2"
  } else {
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
  }
  
  if(pca){
    stopifnot(!data_1 | !data_2) # restriction for now
    embedding <- vector("list", 3)
    if(data_1){
      embedding[[1]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_1,
                                               obj$svd_1, common_bool = T, distinct_bool = F,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
      embedding[[2]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_1,
                                               obj$svd_1, common_bool = F, distinct_bool = T,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
      embedding[[3]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_1,
                                               obj$svd_1, common_bool = T, distinct_bool = T,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
    } else {
      embedding[[1]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_2,
                                               obj$svd_2, common_bool = T, distinct_bool = F,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
      embedding[[2]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_2,
                                               obj$svd_2, common_bool = F, distinct_bool = T,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
      embedding[[3]] <- .extract_matrix_helper(obj$common_score, obj$distinct_score_2,
                                               obj$svd_2, common_bool = T, distinct_bool = T,
                                               center = F, renormalize = F, add_noise = F, premultiply = T)
    }
    
    for(i in 1:3){
      tmp <- .svd_truncated(embedding[[i]], K = 2)
      embedding[[i]] <- .mult_mat_vec(tmp$u, tmp$d)
    }
    xlim <- range(sapply(embedding, function(x){x[,1]}))
    ylim <- range(sapply(embedding, function(x){x[,2]}))
    main_vec <- c("Common view", "Distinct view", "Entire view")
    
    for(i in 1:3){
      graphics::plot(embedding[[i]][,1], embedding[[i]][,2], asp = T, pch = 16, col = membership_vec, main = main_vec[i],
                     xlab = label1, ylab = label2, xlim = xlim, ylim = ylim)
    }
  } else {
    prep_obj <- .prepare_umap_embedding(obj)
    
    graphics::par(mfrow = c(1,3))
    set.seed(10)
    tmp <- .extract_umap_embedding(prep_obj, common_1 = data_1, common_2 = data_2, distinct_1 = F, distinct_2 = F,
                                   add_noise = add_noise)
    graphics::plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Common view",
                   xlab = label1, ylab = label2)
    
    set.seed(10)
    tmp <- .extract_umap_embedding(prep_obj, common_1 = F, common_2 = F, distinct_1 = data_1, distinct_2 = data_2,
                                   add_noise = add_noise)
    graphics::plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct views",
                   xlab = label1, ylab = label2)
    
    set.seed(10)
    tmp <- .extract_umap_embedding(prep_obj, common_1 = data_1, common_2 = data_2, distinct_1 = data_1, distinct_2 = data_2,
                                   add_noise = add_noise)
    graphics::plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Entire view",
                   xlab = label1, ylab = label2)
  }
  
  invisible()
}

#' Side-by-side UMAPs of the data
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec integer vector
#' @param observed boolean. This should \code{TRUE} if \code{obj} is the output of
#' \code{generate_data}, and you want to plot the "observed" data.
#' Otherwise, this is \code{FALSE} by default, meaning that if \code{obj} is the output of
#' \code{generate_data}, you are plotting the "true" data (i.e., not affected by noise),
#' or if \code{obj} is the output of \code{dcca_decomposition}, you are plotting
#' the estimated denoised observed matrix.
#' @param pca boolean. If \code{TRUE}, plot the PCA embedding with the leading 2 components. 
#' If \code{FALSE}, plot the UMAP embedding.
#'
#' @return shows a plot but returns nothing
#' @export
plot_data <- function(obj, membership_vec, observed = F, pca = F){
  if(pca){
    label1 <- "PCA 1"; label2 <- "PCA 2"
  } else {
    label1 <- "UMAP 1"; label2 <- "UMAP 2"
  }
  
  # plot the noise-affected/"observed" data
  if(observed){
    svd_list <- list(.svd_truncated(obj$mat_1, K = ncol(obj$distinct_score_1)), 
                     .svd_truncated(obj$mat_2, K = ncol(obj$distinct_score_2)))
    
    embedding <- lapply(svd_list, function(svd_res){
      tmp <- .mult_mat_vec(svd_res$u, svd_res$d)
      if(pca){
        stopifnot(ncol(tmp) >= 2)
        tmp[,1:2]
      } else {
        set.seed(10)
        Seurat::RunUMAP(tmp, verbose = F)@cell.embeddings
      }
    })
    
    graphics::par(mfrow = c(1,2))
    graphics::plot(embedding[[1]][,1], embedding[[1]][,2], asp = T, pch = 16, col = membership_vec, main = "Obs. dataset 1",
         xlab = label1, ylab = label2)
    graphics::plot(embedding[[2]][,1], embedding[[2]][,2], asp = T, pch = 16, col = membership_vec, main = "Obs. dataset 2",
         xlab = label1, ylab = label2)
    
  } else {
    # plot the denoised/"true" data
    prep_list <- .prepare_umap_embedding(obj)
    
    graphics::par(mfrow = c(1,2))
    if(pca){
      tmp <- .mult_mat_vec(prep_list$svd_list[[1]]$u, prep_list$svd_list[[1]]$d)[,1:2]
    } else {
      set.seed(10)
      tmp <- .extract_umap_embedding(prep_list, common_1 = T, common_2 = F, distinct_1 = T, distinct_2 = F)
    }
    graphics::plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Dataset 1",
         xlab = label1, ylab = label2)
    
    if(pca){
      tmp <- .mult_mat_vec(prep_list$svd_list[[2]]$u, prep_list$svd_list[[2]]$d)[,1:2]
    } else {
      set.seed(10)
      tmp <- .extract_umap_embedding(prep_list, common_1 = F, common_2 = T, distinct_1 = F, distinct_2 = T)
    }
    graphics::plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Dataset 2",
         xlab = label1, ylab = label2)
  }
  
  invisible()
}

#######################################

.colorRamp_custom <- function(vec1, vec2, length, luminosity = T){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }
  
  if(luminosity){
    luminosity_vec <- apply(mat, 1, function(x){
      0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
    })
    
    target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))
    
    mat <- t(sapply(1:nrow(mat), function(x){
      factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
      mat[x,] * factor
    }))
  }
  
  apply(mat, 1, function(x){
    grDevices::rgb(x[1], x[2], x[3])
  })
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
