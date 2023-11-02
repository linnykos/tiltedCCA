#' Plot alignment of variables
#'
#' @param rsquare_vec output of \code{tiltedCCA::postprocess_modality_alignment}
#' @param logpval_vec output of \code{tiltedCCA::postprocess_depvalue}
#' @param bool_hide_points boolean (graphical parameter)
#' @param bool_mark_ymedian boolean (graphical parameter)
#' @param bool_polygon_mean boolean (graphical parameter)
#' @param bool_truncate_xaxis boolean (graphical parameter)
#' @param bool_white_bg boolean (graphical parameter)
#' @param bty character (graphical parameter)
#' @param cex_axis positive number (graphical parameter)
#' @param cex_gene_highlight_inner positive number (graphical parameter)
#' @param cex_gene_highlight_outer positive number (graphical parameter)
#' @param cex_lab positive number (graphical parameter)
#' @param cex_points positive number (graphical parameter)
#' @param col_points character for a color (graphical parameter)
#' @param col_gene_highlight character for a color (graphical parameter)
#' @param col_gene_highlight_border character for a color (graphical parameter)
#' @param col_grid character for a color (graphical parameter)
#' @param density positive number, \code{NA}, or \code{NULL} (graphical parameter)
#' @param gene_names vector of charaters for which genes (in \code{names(rsquare_vec)}) to highlight
#' @param lty_grid_major positive number (graphical parameter)
#' @param lty_grid_minor positive number (graphical parameter)
#' @param lty_polygon positive number (graphical parameter)
#' @param lwd_axis positive number (graphical parameter)
#' @param lwd_axis_ticks positive number (graphical parameter)
#' @param lwd_grid_major positive number (graphical parameter)
#' @param lwd_grid_minor positive number (graphical parameter)
#' @param lwd_polygon positive number (graphical parameter)
#' @param lwd_polygon_bold positive number (graphical parameter)
#' @param mark_median_xthres boolean (graphical parameter)
#' @param pch positive integer (graphical parameter)
#' @param xaxt_num_ticks positive integer (graphical parameter)
#' @param xaxt_grid_spacing positive integer (graphical parameter)
#' @param xlab character (graphical parameter)
#' @param xlim vector of 2 numerics, or \code{NULL} (graphical parameter)
#' @param yaxt_by positive number (graphical parameter)
#' @param yaxt_num_ticks positive integer (graphical parameter)
#' @param ylab character (graphical parameter)
#' @param ylim vector of 2 numerics, or \code{NULL} (graphical parameter)
#' @param verbose non-negative integer
#' @param ... 
#'
#' @return makes a plot but does not return anything
#' @export
plot_alignment <- function(rsquare_vec,
                           logpval_vec,
                           bool_hide_points = F,
                           bool_mark_ymedian = F,
                           bool_polygon_mean = T,
                           bool_truncate_xaxis = T,
                           bool_white_bg = T,
                           bty = "n",
                           cex_axis = 1,
                           cex_gene_highlight_inner = 3,
                           cex_gene_highlight_outer = 4,
                           cex_lab = 1,
                           cex_points = 2,
                           col_points = grDevices::rgb(0.5, 0.5, 0.5, 0.1),
                           col_gene_highlight = 2,
                           col_gene_highlight_border = "white",
                           col_grid = "gray",
                           density = 20,
                           gene_names = NULL,
                           lty_grid_major = 2,
                           lty_grid_minor = 3,
                           lty_polygon = 2,
                           lwd_axis = 1,
                           lwd_axis_ticks = 1,
                           lwd_grid_major = 1,
                           lwd_grid_minor = 0.5,
                           lwd_polygon = 1,
                           lwd_polygon_bold = 2,
                           mark_median_xthres = 1,
                           pch = 16,
                           xaxt_num_ticks = 10,
                           xaxt_grid_spacing = 2,
                           xlab = "Separability (Median -Log10(p-value))",
                           xlim = NULL,
                           yaxt_by = 0.1,
                           yaxt_num_ticks = 2,
                           ylab = "Alignment with common space (R^2)",
                           ylim = NULL,
                           verbose = T, ...){
  stopifnot(length(rsquare_vec) == length(logpval_vec),
            all(rsquare_vec >= 0), all(rsquare_vec <= 1), 
            all(logpval_vec >= 0),
            length(names(rsquare_vec)) > 0,
            length(names(logpval_vec)) > 0)
  rsquare_vec <- rsquare_vec[names(logpval_vec)]
  
  logpval_vec_org <- logpval_vec
  logpval_vec[logpval_vec == 0] <- min(logpval_vec[logpval_vec != 0])
  logpval_vec <- log10(logpval_vec)
  if(all(is.null(xlim))){
    xmin <- floor(min(logpval_vec)); xmax <- ceiling(max(logpval_vec))
    if(bool_truncate_xaxis){
      xlim <- range(logpval_vec)
    } else {
      xlim <- c(xmin, xmax)
    }
  } else {
    xmin <- xlim[1]; xmax <- xlim[2]
  }
  if(all(is.null(ylim))) ylim <- c(0, 1)
  
  graphics::plot(NA,
                 bty = bty,
                 cex.lab = cex_lab,
                 xaxt = "n",
                 yaxt = "n",
                 xlab = xlab,
                 xlim = xlim,
                 ylab = ylab,
                 ylim = ylim, ...)
  
  # x ticks
  tmp <- seq(xmin, xmax, by = 1)
  xaxt_vec <- sort(unique(unlist(sapply(1:(length(tmp)-1), function(i){
    lower <- 10^(tmp[i]); upper <- 10^(tmp[i+1])
    seq_vec <- seq(lower, upper, length.out = xaxt_num_ticks)
    log10(seq_vec)
  }))))
  if(bool_truncate_xaxis){
    tmp_subset <- tmp[intersect(which(tmp >= xlim[1]), 
                                which(tmp <= xlim[2]))]
    xaxt_vec_subset <- xaxt_vec[intersect(which(xaxt_vec >= xlim[1]), 
                                          which(xaxt_vec <= xlim[2]))]
  } else {
    tmp_subset <- tmp; xaxt_vec_subset <- xaxt_vec
  }
  graphics::axis(2,
                 cex.axis = cex_axis,
                 lwd = lwd_axis,
                 lwd.ticks = lwd_axis_ticks)
  graphics::axis(1, at = xaxt_vec_subset, 
                 labels = rep("", length(xaxt_vec_subset)),
                 cex.axis = cex_axis,
                 lwd = lwd_axis,
                 lwd.ticks = lwd_axis_ticks)
  graphics::axis(1, at = tmp_subset, 
                 labels = as.character(10^tmp_subset),
                 cex.axis = cex_axis,
                 lwd = lwd_axis,
                 lwd.ticks = lwd_axis_ticks)
  
  # y lines
  y_vec <- seq(0, 1, by = yaxt_by)
  y_vec_major <- y_vec[seq(1, length(y_vec), by = yaxt_num_ticks)]
  y_vec_minor <- setdiff(y_vec, y_vec_major)
  for(y_val in y_vec_minor){
    graphics::lines(c(-2,2)*max(abs(xlim)), 
                    rep(y_val, 2), 
                    col = col_grid, 
                    lty = lty_grid_minor, 
                    lwd = lwd_grid_minor)
  }
  for(y_val in y_vec_major){
    graphics::lines(c(-2,2)*max(abs(xlim)), 
                    rep(y_val, 2), 
                    col = col_grid, 
                    lty = lty_grid_major, 
                    lwd = lwd_grid_major)
  }
  
  # x lines 
  tmp <- pmax(seq(0, length(xaxt_vec), by = xaxt_grid_spacing), 1)
  major <- c(1, tmp[tmp %% xaxt_num_ticks == 0])
  minor <- tmp[!tmp %in% major]
  for(x_val in xaxt_vec[minor]){
    if(x_val >= xlim[1] & x_val <= xlim[2]){
      graphics::lines(x = rep(x_val,2), 
                      y = c(-2,2)*max(abs(ylim)),
                      col = col_grid, 
                      lty = lty_grid_minor, 
                      lwd = lwd_grid_minor)
    }
  }
  for(x_val in xaxt_vec[major]){
    if(x_val >= xlim[1] & x_val <= xlim[2]){
      graphics::lines(x = rep(x_val,2), 
                      y = c(-2,2)*max(abs(ylim)),
                      col = col_grid, 
                      lty = lty_grid_major, 
                      lwd = lwd_grid_major)
    }
  }
  
  
  if(!bool_hide_points & bool_mark_ymedian){
    quantile_vec <- stats::quantile(rsquare_vec[logpval_vec_org >= mark_median_xthres],
                                    probs = c(0.25, 0.75))

    graphics::polygon(cbind(c(-2, 2, 2, -2)*max(abs(logpval_vec)),
                            c(rep(quantile_vec[1], 2), rep(quantile_vec[2], 2))),
                      col = col_gene_highlight_border,
                      density = NA,
                      lwd = lwd_polygon)
  }
  
  if(!bool_hide_points & bool_white_bg){
    graphics::points(x = logpval_vec,
                     y = rsquare_vec,
                     cex = cex_gene_highlight_inner + .5*abs(cex_gene_highlight_inner - cex_gene_highlight_outer),
                     col = "black",
                     pch = pch)
    graphics::points(x = logpval_vec,
                     y = rsquare_vec,
                     cex = cex_gene_highlight_inner,
                     col = "white",
                     pch = pch)
  }
  
  if(!bool_hide_points & !all(is.null(gene_names)) & bool_polygon_mean){
    mean_val <- mean(rsquare_vec[gene_names])
    range_vec <- mean_val + c(-1,1)*stats::sd(rsquare_vec[gene_names])
    
    graphics::lines(c(-2,2)*max(abs(xlim)), rep(mean_val, 2),
                    col = col_gene_highlight,
                    lty = lty_polygon,
                    lwd = lwd_polygon_bold)
    graphics::polygon(cbind(c(-2, 2, 2, -2)*max(abs(logpval_vec)),
                            c(rep(range_vec[1], 2), rep(range_vec[2], 2))),
                      col = col_gene_highlight,
                      density = density,
                      lwd = lwd_polygon)
  }
  
  if(!bool_hide_points){
    graphics::points(x = logpval_vec,
                     y = rsquare_vec,
                     cex = cex_points,
                     col = col_points,
                     pch = pch)
    
    if(!all(is.null(gene_names))){
      graphics::points(x = logpval_vec[gene_names],
                       y = rsquare_vec[gene_names],
                       cex = cex_gene_highlight_outer,
                       col = col_gene_highlight_border,
                       pch = pch)
      graphics::points(x = logpval_vec[gene_names],
                       y = rsquare_vec[gene_names],
                       cex = cex_gene_highlight_inner,
                       col = col_gene_highlight,
                       pch = pch)
    }
    
    if(bool_mark_ymedian){
      med_val <- stats::median(rsquare_vec[logpval_vec_org >= mark_median_xthres])
      if(verbose) print(paste0("Median value of: ", round(med_val, 3)))
      
      graphics::lines(c(-2,2)*max(abs(xlim)), rep(med_val, 2),
                      col = col_gene_highlight,
                      lty = lty_polygon,
                      lwd = lwd_polygon_bold)
    }
  }
  
  invisible()
}