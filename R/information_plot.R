#' Constructing information plot
#' 
#' \code{weight_mat} has one row for each unit (ex: cell)
#'
#' @param weight_mat data frame with 3 columns, where the first two columns
#' denote the row-specific weight for distinct information (values between 0
#' and 1) and the third column is the label
#' @param main character
#' @param reorder_types boolean
#' @param plot_individual boolean
#' @param col_vec vector of colors
#' @param y_offset numeric
#'
#' @return none
#' @export
information_plot <- function(weight_mat, main = "", 
                             reorder_types = T, plot_individual = T,
                             col_vec = scales::hue_pal()(length(unique(weight_mat[,3]))),
                             y_offset = 0.1){
  stopifnot(is.data.frame(weight_mat), ncol(weight_mat) == 3, colnames(weight_mat)[3] == "cell_type",
            all(weight_mat[,1:2] <= 1), all(weight_mat[,1:2] >= 0))
  
  # compute group-wise weight matrix
  uniq_types <- unique(weight_mat$cell_type)
  summary_mat <- t(sapply(uniq_types, function(i){
    idx <- which(weight_mat$cell_type == i)
    colMeans(weight_mat[idx,1:2,drop = F])
  }))
  summary_mat <- data.frame(mode_1 = summary_mat[,1], mode_2 = summary_mat[,2], cell_type = uniq_types)
  colnames(summary_mat) <- colnames(weight_mat)
  
  # order the weight matrix
  diff_val <- summary_mat[,1] - summary_mat[,2]
  if(reorder_types){
    summary_mat <- summary_mat[order(diff_val, decreasing = T),]
  }
  
  # take the colors
  type_col_mat <- data.frame(cell_type = as.character(sort(uniq_types)), col = col_vec)
  
  # draw the plot
  len <- length(uniq_types)
  graphics::plot(NA, xlim = c(0, 1.5*(len-1)+1), ylim = c(-1,1), xlab = "", ylab = "Ratio explained",
       xaxt = "n", main = main)
  graphics::title(xlab = "Cell type", mgp = c(4,1,0))
  
  y_vec <- c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)
  for(i in 1:length(y_vec)){
    graphics::lines(c(-10*len, 10*len), rep(y_vec[i], 2), col = "red", lwd = 0.75, lty = 3)
  }
  
  for(i in 1:len){
    xlim <- c((i-1)*1.5, (i-1)*1.5+1)
   
    col <- type_col_mat$col[which(type_col_mat$cell_type == summary_mat$cell_type[i])]
    graphics::rect(xlim[1], -(1-summary_mat[i,2]), xlim[2], 1-summary_mat[i,1], col = col)
    
    if(plot_individual){
      idx <- which(weight_mat$cell_type == summary_mat$cell_type[i])
      if(length(idx) > 50){ idx <- idx[sample(1:length(idx), 50)] }
      
      offsets <- seq(0, 1, length.out = length(idx)+10)
      offsets <- offsets[-c(1:5, (length(offsets)-4):(length(offsets)))]
      
      for(j in 1:length(idx)){
        graphics::lines(rep(xlim[1]+offsets[j], 2), c(1-weight_mat[idx[j],1], -(1-weight_mat[idx[j],2])),
                        col = "gray", lwd = 0.3)
        graphics::points(rep(xlim[1]+offsets[j], 2), c(1-weight_mat[idx[j],1], -(1-weight_mat[idx[j],2])),
                         col = "gray", cex = 0.2, pch = 16)
      }
    }
  }
  graphics::lines(c(-2*len, 2*len), rep(0,2), col = "black", lwd = 3, lty = 2)
  
  graphics::axis(side = 1, at = (((1:len)-1)*1.5+1.5/2), labels = FALSE)
  graphics::text(x = (((1:len)-1)*1.5+1.5/2),
       y = graphics::par("usr")[3]-y_offset,
       labels = as.character(summary_mat$cell_type),
       xpd = NA, srt = 45, adj = 0.965, cex = 0.8)
  
  invisible()
}