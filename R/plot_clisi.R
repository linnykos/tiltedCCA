#' Making the plot for local enrichment
#' 
#'
#' @param local_1 output of \code{clisi_information} on one Modality
#' @param local_2 output of \code{clisi_information} on the other Modality
#' @param col_vec vector of colors
#' @param l_bg \code{l} parameter (luminosity, i.e., brightness) for individual cells
#' @param c_bg \code{c} parameter (chroma, i.e., color intensity) for individual cells
#' @param alpha_bg \code{alpha} parameter (color) for individual cells
#' @param xlab1 \code{xlab} for Modality 1
#' @param xlab2 \code{xlab} for Modality 2
#' @param ylab \code{ylab} for the shared modality
#' @param main1 Title for plot corresponding to Modality 1
#' @param main2 Title for plot corresponding to Modality 2
#' @param ... extra parameters for \code{ggrepel::geom_text_repel}
#'
#' @return List of two \code{gg} objects
#' @export
plot_clisi <- function(local_1, local_2,
                       col_vec = scales::hue_pal()(nrow(local_1$common_clisi$membership_info)),
                       l_bg = 75, c_bg = 50, alpha_bg = 0.5, 
                       xlab1 = "Distinct enrichment",
                       xlab2 = "Distinct enrichment",
                       ylab = "Common enrichment",
                       main1 = "Modality 1", main2 = "Modality 2", ...){
  stopifnot(class(local_1) == "clisi", class(local_2) == "clisi",
            all(dim(local_1$common_clisi$membership_info) == dim(local_2$common_clisi$membership_info)))
  stopifnot(length(col_vec) == nrow(local_1$common_clisi$membership_info))
  
  # setup
  n <- nrow(local_1$common_clisi$cell_info)
  k <- nrow(local_1$common_clisi$membership_info)
  local_lis <- list(local_1, local_2)
  plot_lis <- vector("list", length = 2)
  
  # construct colors
  bg_col_vec <- .adjust_colors(col_vec, l_bg = l_bg, c_bg = c_bg, alpha_bg = alpha_bg)
  all_col_vec <- c(col_vec, bg_col_vec)
  tmp <- local_lis[[1]]$common_clisi$membership_info$celltype
  names(all_col_vec) <- c(tmp, paste0(tmp, "0"))
  custom_colors <- ggplot2::scale_colour_manual(values = all_col_vec)
  category = celltype = common = distinct = NULL # for appeasing R CHECK
  
  for(i in 1:2){
    df <- data.frame("celltype" = as.factor(c(paste0(as.character(local_lis[[i]]$common_clisi$cell_info$celltype), "0"), 
                                            as.character(local_lis[[i]]$common_clisi$membership_info$celltype))), 
                     "common" = c(local_lis[[i]]$common_clisi$cell_info$clisi_score, local_lis[[i]]$common_clisi$membership_info$mean_clisi),
                     "distinct" = c(local_lis[[i]]$distinct_clisi$cell_info$clisi_score, local_lis[[i]]$distinct_clisi$membership_info$mean_clisi),
                     "category" = as.factor(c(rep(0, n), rep(1, k))))
    
    plot1 <- ggplot2::ggplot(data = subset(df, category == 0), ggplot2::aes(x = distinct, y = common, color = celltype))
    plot1 <- plot1 + ggplot2::geom_point()
    if(i == 1){
      plot1 <- plot1 + ggplot2::xlim(1, 0) + ggplot2::ylim(0, 1)
      plot1 <- plot1 + ggplot2::geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed")
    } else {
      plot1 <- plot1 + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
      plot1 <- plot1 + ggplot2::geom_abline(intercept = 1, slope = -1, color = "red", linetype = "dashed")
    }
    plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                         ggplot2::aes(x = distinct, y = common), 
                                         size = 3, color = "black")
    plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                         ggplot2::aes(x = distinct, y = common), 
                                         size = 2.5, color = "white")
    plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                         ggplot2::aes(x = distinct, y = common,
                                                      color = celltype), 
                                         size = 2)
    plot1 <- plot1 + custom_colors
    plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, category == 1), ggplot2::aes(label = celltype),
                                              color = "black",
                                              segment.color = "grey50",
                                              size = 2, ...)
    plot1 <- plot1 + ggplot2::ylab(ylab)
    if(i == 1){
      plot1 <- plot1 + ggplot2::xlab(xlab1)
      plot1 <- plot1 + ggplot2::ggtitle(main1)
    } else {
      plot1 <- plot1 + ggplot2::xlab(xlab2)
      plot1 <- plot1 + ggplot2::ggtitle(main2)
    }
    plot1 <- plot1 + Seurat::NoLegend()
    plot_lis[[i]] <- plot1
  }
  
  plot_lis
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


