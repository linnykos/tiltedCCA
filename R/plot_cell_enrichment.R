#' Plot cell enrichment
#'
#' @param cell_enrichment_res the result of \code{tiltedCCA::postprocess_cell_enrichment}
#' @param cex_axis positive number (graphical parameter)
#' @param cex_lab positive number (graphical parameter)
#' @param lwd_diag positive number (graphical parameter)
#' @param lwd_grid positive number (graphical parameter)
#' @param mar vector of 4 positive numbers (graphical parameter)
#' @param xlab_1 character (graphical parameter)
#' @param xlab_2 character (graphical parameter)
#' @param ylab character (graphical parameter)
#'
#' @return makes a plot but does not return anything
#' @export
plot_cell_enrichment <- function(cell_enrichment_res,
                                 cex_axis = 1,
                                 cex_lab = 1,
                                 lwd_diag = 3,
                                 lwd_grid = 2,
                                 mar = c(4, 5, 0.5, 0.5),
                                 xlab_1 = "Modality 1 Distinct",
                                 xlab_2 = "Modality 2 Distinct",
                                 ylab = "Common"){
  par(mfrow = c(1,2), mar = mar)
  
  y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
  names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
  x_vec <- -cell_enrichment_res$enrichment_distinct_1$df[,"value"]
  names(x_vec) <- cell_enrichment_res$enrichment_distinct_1$df[,"celltype"]
  
  graphics::plot(NA, xlim = c(-1,0), ylim = c(0,1),
       main = "", xlab = xlab_1, ylab = "",
       xaxt = "n", yaxt = "n", bty = "n", asp = T, cex.lab = cex_lab)
  for(x in seq(-1,0,by=0.1)){
    graphics::lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = lwd_diag)
  }
  for(y in seq(0,1,by=0.1)){
    graphics::lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = lwd_diag)
  }
  graphics::lines(c(-10,10), c(10,-10), col = 2, lty = 2, lwd = lwd_diag)
  
  col_vec <- col_palette[names(x_vec)]
  graphics::points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
  graphics::points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
  graphics::points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)
  
  graphics::axis(1, cex.axis = cex_axis)
  
  y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
  names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
  x_vec <- cell_enrichment_res$enrichment_distinct_2$df[,"value"]
  names(x_vec) <- cell_enrichment_res$enrichment_distinct_2$df[,"celltype"]
  
  graphics::plot(NA, xlim = c(0,1), ylim = c(0,1),
       main = "", xlab = xlab_2, ylab = ylab,
       xaxt = "n", yaxt = "n", bty = "n", asp = T, cex.lab = cex_lab)
  for(x in seq(0,1,by=0.1)){
    graphics::lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = lwd_diag)
  }
  for(y in seq(0,1,by=0.1)){
    graphics::lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = lwd_diag)
  }
  graphics::lines(c(-10,10), c(-10,10), col = 2, lty = 2, lwd = lwd_diag)
  
  col_vec <- col_palette[names(x_vec)]
  graphics::points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
  graphics::points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
  graphics::points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)
  
  graphics::axis(1, cex.axis = cex_axis)
  graphics::axis(2, cex.axis = cex_axis)
  
  invisible()
}