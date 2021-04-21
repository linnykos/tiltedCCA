#' Making the plot for cLISI
#' 
#' See https://www.r-bloggers.com/2012/06/two-tips-adding-title-for-graph-with-multiple-plots-add-significance-asterix-onto-a-boxplot/
#'
#' @param clisi_1 output of \code{clisi_information} on one Modality
#' @param clisi_2 output of \code{clisi_information} on the other Modality
#' @param col_vec vector of colors
#' @param cell_max number of cells of each cell-type to plot
#' @param par_mar \code{mar} parameter for \code{par}
#' @param par_oma \code{oma} parameter for \code{par}
#' @param asp graphical \code{asp} parameter
#' @param pch_main \code{pch} parameter for cell-types
#' @param cex_main \code{cex} parameter for cell-types
#' @param pch_bg \code{pch} parameter for individual cells
#' @param cex_bg \code{cex} parameter for individual cells
#' @param alpha_bg \code{alpha} parameter (color) for individual cells
#' @param l_bg \code{l} parameter (luminosity) for individual cells
#' @param c_bg \code{c} parameter (chroma, i.e., brightness) for individual cells
#' @param xlim \code{xlim} graphical parameter
#' @param ylim \code{ylim} graphical parameter
#' @param gridsize number of grid ticks
#' @param col_grid \code{col} for grid
#' @param lty_grid \code{lty} for grid
#' @param lwd_grid \code{lwd} for grid
#' @param col_diag \code{col} for diagonal line
#' @param lty_diag \code{lty} for diagonal line
#' @param lwd_diag \code{lwd} for diagonal line
#' @param xlab1 \code{xlab} for Modality 1
#' @param xlab2 \code{xlab} for Modality 2
#' @param ylab \code{ylab} for the shared modality
#' @param ylab_dist distance of \code{ylab} from y-axis
#' @param main \code{main} parameter
#' @param cex_text_main \code{cex} for \code{main} parameter
#'
#' @return nothing
#' @export
plot_clisi <- function(clisi_1, clisi_2,
                       col_vec = scales::hue_pal()(nrow(clisi_1$common_clisi$membership_info)),
                       cell_max = 5000,
                       par_mar = c(4,2.5,0.5,0.5), par_oma = c(0,0,2,0),
                       asp = T,
                       pch_main = 16, cex_main = 1.5,
                       pch_bg = 16, cex_bg = 1, alpha_bg = 0.5,
                       l_bg = 95, c_bg = 50,
                       xlim = c(0,1), ylim = c(0,1), 
                       gridsize = 5, col_grid = grDevices::rgb(0.8,0.8,0.8),
                       lty_grid = 3, lwd_grid = 1,
                       col_diag = "firebrick", lty_diag = 2, lwd_diag = 2,
                       xlab1 = "Distinct information 1",
                       xlab2 = "Distinct information 2",
                       ylab = "Common information", ylab_dist = 0.5,
                       main = "cLISI Information", cex_text_main = 1.5){
  stopifnot(class(clisi_1) == "clisi", class(clisi_2) == "clisi",
            all(dim(clisi_1$common_clisi$membership_info) == dim(clisi_2$common_clisi$membership_info)))
  stopifnot(length(col_vec) == nrow(clisi_1$common_clisi$membership_info))
  
  bg_col_vec <- .adjust_colors(col_vec, l_bg, c_bg, alpha_bg)
  graphics::par(mfrow = c(1,2), mar = par_mar, oma = par_oma)
  x_vec <- seq(xlim[1], xlim[2], length.out = gridsize)
  y_vec <- seq(ylim[1], ylim[2], length.out = gridsize)
  
  graphics::plot(NA, xlim = sort(-1*xlim), ylim = ylim, asp = asp, xlab = xlab1,
                 ylab = "", yaxt = 'n', xaxt = 'n')
  graphics::title(ylab = ylab, line = ylab_dist)
  graphics::axis(1, at = -1*x_vec, labels = as.character(round(x_vec, 2)))
  graphics::axis(4, at = y_vec, labels = NA)
  .draw_grid(x_vec, y_vec, xlim, ylim, col_grid, lty_grid, lwd_grid,
             flip = T)
  
  .plot_clisi_cell(clisi_1, bg_col_vec, pch_bg, cex_bg, cell_max, flip = T)
  graphics::lines(c(0,-1),c(0,1), col = col_diag, lty = lty_diag, lwd = lwd_diag)
  .plot_clisi_type(clisi_1, col_vec, pch_main, cex_main, flip = T)
  
  ####
  
  graphics::plot(NA, xlim = xlim, ylim = ylim, asp = asp, xlab = xlab2,
                 ylab = "", yaxt = 'n', xaxt = 'n')
  graphics::axis(1, at = x_vec, labels = as.character(round(x_vec, 2)))
  graphics::axis(2, at = y_vec, labels = as.character(round(y_vec, 2)))
  .draw_grid(x_vec, y_vec, xlim, ylim, col_grid, lty_grid, lwd_grid,
             flip = F)
  
  .plot_clisi_cell(clisi_2, bg_col_vec, pch_bg, cex_bg, cell_max, flip = F)
  graphics::lines(c(0,1),c(0,1), col = col_diag, lty = lty_diag, lwd = lwd_diag)
  .plot_clisi_type(clisi_2, col_vec, pch_main, cex_main, flip = F)
  
  ####
  
  graphics::mtext(main, outer = TRUE, cex = cex_text_main)
  
  invisible()
}

#' Plot the cLISI legend
#'
#' @param clisi_obj output of \code{clisi_information} 
#' @param col_vec vector of colors
#' @param percent_coverage numeric
#' @param pch \code{pch} parameter
#' @param cex_point \code{cex} parameter for points
#' @param cex_text \code{cex} parameter for the text
#' @param text_nudge x-axis offset for the text
#' @param xlim \code{xlim} parameter for the plot
#' @param ... additional graphical parameters
#'
#' @return nothing
#' @export
plot_clisi_legend <- function(clisi_obj, col_vec = scales::hue_pal()(nrow(clisi_obj$common_clisi$membership_info)),
                              percent_coverage = 1, pch = 16, cex_point = 1,
                              cex_text = 1, text_nudge = 0, xlim = c(0,1), ...){
  stopifnot(length(col_vec) == nrow(clisi_obj$common_clisi$membership_info))
  
  graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
  graphics::plot(NA, xlim = xlim, ylim = c(0,1), yaxt = "n", xaxt = "n", bty = "n", 
                 xlab = "", ylab = "", ...)
  
  # plot the colors
  n <- length(col_vec)
  graphics::points(x = rep(0,n), y = seq(1,0,length.out=n), pch = pch, cex = cex_point,
                   col = col_vec)
  graphics::text(x = rep(0+text_nudge,n), y = seq(1,0,length.out=n), 
                 labels = sort(clisi_obj$common_clisi$membership_info$celltype, decreasing = F),
                 pos = 4)
  
  invisible()
}

##################################

.draw_grid <- function(x_vec, y_vec, xlim, ylim, 
                       col_grid, lty_grid, lwd_grid,
                       flip){
  s <- ifelse(flip, -1, 1)
  for(i in 1:length(x_vec)){
    graphics::lines(s*rep(x_vec[i],2), ylim, lty = lty_grid, lwd = lwd_grid, col = col_grid)
  }
  
  for(i in 1:length(y_vec)){
    graphics::lines(s*xlim, rep(y_vec[i],2), lty = lty_grid, lwd = lwd_grid, col = col_grid)
  }
  
  invisible()
}

.plot_clisi_cell <- function(clisi_obj, bg_col_vec, pch_bg, cex_bg, cell_max, flip){
  stopifnot(class(clisi_obj) == "clisi")
  stopifnot(length(bg_col_vec) == nrow(clisi_obj$common_clisi$membership_info))
  
  s <- ifelse(flip, -1, 1)
  n <- nrow(clisi_obj$common_clisi$cell_info)
  n_idx <- sample(1:n, size = min(n, cell_max))
  graphics::points(s*clisi_obj$distinct_clisi$cell_info$clisi_score[n_idx],
                   clisi_obj$common_clisi$cell_info$clisi_score[n_idx],
                   col = bg_col_vec[as.numeric(clisi_obj$common_clisi$cell_info$celltype)][n_idx],
                   pch = pch_bg, cex = cex_bg)
  
  invisible()
}

.plot_clisi_type <- function(clisi_obj, col_vec, pch_main, cex_main, flip){
  stopifnot(class(clisi_obj) == "clisi")
  stopifnot(length(col_vec) == nrow(clisi_obj$common_clisi$membership_info))
  
  s <- ifelse(flip, -1, 1)
  n <- nrow(clisi_obj$common_clisi$membership_info)
  n_idx <- sample(1:n)
  graphics::points(s*clisi_obj$distinct_clisi$membership_info$mean_clisi[n_idx],
                   clisi_obj$common_clisi$membership_info$mean_clisi[n_idx],
                   col = "black",
                   pch = pch_main, cex = 1.25*cex_main)
  graphics::points(s*clisi_obj$distinct_clisi$membership_info$mean_clisi[n_idx],
                   clisi_obj$common_clisi$membership_info$mean_clisi[n_idx],
                   col = "white",
                   pch = pch_main, cex = 1.2*cex_main)
  graphics::points(s*clisi_obj$distinct_clisi$membership_info$mean_clisi[n_idx],
                   clisi_obj$common_clisi$membership_info$mean_clisi[n_idx],
                   col = col_vec[order(clisi_obj$common_clisi$membership_info$celltype, decreasing = F)][n_idx],
                   pch = pch_main, cex = cex_main)
  
  invisible()
}

#############

.adjust_colors <- function(col_vec, l_bg, c_bg, alpha_bg){
  # find if conversion to hex is necessary
  idx <- which(sapply(col_vec, function(x){substring(x,1,1) != "#"}))
  if(length(idx) > 0){
    # convert if needed
    col_vec[idx] <- sapply(col_vec[idx], .name_to_colhex)
  }
  
  # detect non-gray colors
  idx <- which(!sapply(col_vec, .detect_gray))
  # apply alpha
  res <- scales::alpha(col_vec, alpha = alpha_bg*(100-c_bg)/100)
  
  # apply luminosity and chroma
  if(length(idx) > 0){
    col_vec[idx] <- scales::col2hcl(col_vec[idx], l = l_bg, c = c_bg)
    res[idx] <- scales::alpha(col_vec[idx], alpha = alpha_bg)
  }
 
  res
}

.name_to_colhex <- function(val){
  tmp <- as.numeric(grDevices::col2rgb(val))
  grDevices::rgb(tmp[1], tmp[2], tmp[3], maxColorValue=255)
}

.detect_gray <- function(str){
  stopifnot(is.character(str), substring(str,1,1) == "#")
  val1 <- substring(str,2,3)
  val2 <- substring(str,4,5)
  val3 <- substring(str,6,7)
  
  val1 == val2 & val2 == val3
}


