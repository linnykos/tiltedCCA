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
plot_heatmat <- function(dat, luminosity = F, asp = nrow(dat)/ncol(dat), 
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

#' Plot summary of D-CCA
#'
#' @param obj output of \code{dcca_decomposition}
#' @param xlab \code{xlab} graphical parameter
#' @param ylab1 \code{ylab} graphical parameter of left axis
#' @param ylab2 \code{ylab} graphical parameter of right axis
#' @param main \code{main} graphical parameter
#' @param pch \code{pch} graphical parameter
#' @param cex \code{cex} graphical parameter
#' @param lwd \code{lwd} graphical parameter
#' @param ... additional graphical parameter
#'
#' @return shows a plot but returns nothing
#' @export
plot_summary <- function(obj, xlab = "Latent dimension", 
                         ylab1 = "CCA objective", ylab2 = "Distinct % for Modality 2", 
                         main = "",
                         pch = 16, cex = 1, lwd = 1, ...){
  stopifnot(inherits(obj, c("dcca_decomp", "dcca")), length(obj$cca_obj) == length(obj$distinct_perc_2))
  
  k <- length(obj$cca_obj)
  graphics::plot(x = 1:k, y = obj$cca_obj, xlim = c(1,k), ylim = c(0,1), 
                 xlab = xlab, ylab = ylab1, pch = 16, cex = cex, col = "black",
                 main = main,
                 ...)
  graphics::lines(x = 1:k, y = obj$cca_obj, lwd = lwd)
  
  graphics::par(new = T)
  graphics::plot(x = 1:k, y = obj$distinct_perc_2, xlim = c(1,k), ylim = c(0,1), 
                 xlab = "", ylab = "", pch = 16, cex = cex, col = "red",
                 axes = F, ...)
  graphics::lines(x = 1:k, y = obj$distinct_perc_2, col = "red", lty = 2, lwd = lwd)
  graphics::axis(side = 4, col = "red", col.axis = "red")
  graphics::mtext(ylab2, side = 4, line = 3, col = "red")
  
  graphics::lines(x = c(-10*k,10*k), y = rep(0.5,2), col = "red", 
                  lty = 3, lwd = 1)

  invisible()
}

#' Side-by-side plot of the canonical scores, colored by membership
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param membership_vec factor vector
#' @param col_vec vector of colors
#' @param xlim custom \code{xlim} graphical argument
#' @param decomposition boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_scores <- function(obj, membership_vec, col_vec = scales::hue_pal()(length(levels(membership_vec))), 
                        xlim = NA, decomposition = F){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(obj$common_score))
  
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
  max_col <- max(sapply(score_list, ncol))
  if(all(is.na(xlim))) xlim <- range(do.call(cbind, score_list))
  
  if(decomposition){
    main_vec <- c("Common", "Distinct 1", "Distinct 2")
  } else {
    main_vec <- c("Dataset 1", "Dataset 2")
  }
  
  n_idx <- sample(1:n)
  graphics::par(mfrow = c(1, length(score_list)))
  for(k in 1:length(score_list)){
    graphics::plot(NA, xlim = xlim, ylim = c(0.5, max_col+.5), main = main_vec[k], xlab = "Value", ylab = "Dimension")
    graphics::lines(rep(0, 2), max_col*c(-10,10), col = "red", lty = 2)
    
    for(i in 1:ncol(score_list[[k]])){
      graphics::points(x = score_list[[k]][n_idx,i], y = stats::runif(n, min = i-.2, max = i+.2), 
                       col = col_vec[as.numeric(membership_vec)][n_idx], pch = 16)
    }
  }
  
  invisible()
}


#' Side-by-side plot of the canonical scores as heatmaps
#'
#' @param obj output of either \code{generate_data} or \code{dcca_decomposition}
#' @param main_vec vector of characters for the title of the plots
#' @param membership_vec factor vector
#' @param num_col positive integers for number of distinct colors
#' @param log_scale boolean
#' @param scaling_power positive numeric
#' @param luminosity boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_scores_heatmap.dcca <- function(obj, main_vec = c("Common score", "Distinct score 1", "Distinct score 2"),
                                     membership_vec = NA, num_col = 10, 
                                     log_scale = F, scaling_power = 1, luminosity = F){
  lis <- list(obj$common_score, obj$distinct_score_1, obj$distinct_score_2)
  
  plot_scores_heatmap.list(lis, 
                           main_vec = main_vec, 
                           membership_vec = membership_vec,
                           num_col = num_col,
                           log_scale = log_scale,
                           scaling_power = scaling_power,
                           luminosity = luminosity)
  
}

#' Side-by-side plot of the canonical scores as heatmaps
#'
#' @param obj \code{list}
#' @param main_vec vector of characters for the title of the plots
#' @param membership_vec factor vector
#' @param num_col positive integers for number of distinct colors
#' @param log_scale boolean
#' @param scaling_power positive numeric
#' @param luminosity boolean
#'
#' @return shows a plot but returns nothing
#' @export
plot_scores_heatmap.list <- function(obj, main_vec = NA, membership_vec = NA, num_col = 10, 
                                     log_scale = F, scaling_power = 1, luminosity = F){
  stopifnot(is.list(obj), all(sapply(obj, is.matrix)), length(unique(sapply(obj, nrow))) == 1)
  if(!all(is.na(main_vec))) stopifnot(length(obj) == length(main_vec))
  
  n <- nrow(obj[[1]])
  if(log_scale) obj <- lapply(obj, function(x){log(abs(x)+1)*sign(x)})
  zlim <- range(unlist(obj))
  
  # construct colors. green is negative
  max_val <- max(abs(zlim)); min_val <- max(min(abs(zlim)), 1e-3)
  col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), num_col,
                                   luminosity = luminosity)
  break_vec_neg <- -max_val*seq(1, 0, length.out = num_col+2)^scaling_power
  break_vec_neg <- break_vec_neg[-length(break_vec_neg)]
  
  # red is positive
  col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), num_col,
                                   luminosity = luminosity)
  break_vec_pos <- max_val*seq(0, 1, length.out = num_col+2)^scaling_power
  break_vec_pos <- break_vec_pos[-1]
  
  # combine the two
  break_vec <- c(break_vec_neg, break_vec_pos)
  col_vec <- c(col_vec_neg, "white", col_vec_pos)
  
  if(!all(is.na(membership_vec))){
    stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(obj[[1]]))
    
    membership_vec <- as.numeric(membership_vec) ## convert to integers
    idx <- order(membership_vec, decreasing = F)
    breakpoints <- 1-which(abs(diff(sort(membership_vec, decreasing = F))) >= 1e-6)/n
  } else {
    idx <- 1:n
  }
  
  
  line_func <- function(membership_vec, breakpoints){
    if(!all(is.na(membership_vec))){
      for(i in 1:length(breakpoints)){
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.1, col = "white")
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
      }
    }
  }
  
  graphics::par(mfrow = c(1,length(obj)))
  for(i in 1:length(obj)){
    graphics::image(.rotate(obj[[i]][idx,,drop = F]), 
                    main = ifelse(!all(is.na(main_vec)), main_vec[i], ""),
                    col = col_vec, breaks = break_vec)
    line_func(membership_vec, breakpoints)
  }
  
  invisible()
}
  

#######################################

.colorRamp_custom <- function(vec1, vec2, length, luminosity){
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