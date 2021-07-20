#' Plot signal with respect to Laplacians
#'
#' @param seurat_obj \code{Seurat} object, which is used to apply the plotting functions to.
#' We require 3 reductions available, \code{"common"}, \code{"distinct"}
#' and \code{"everything"}
#' @param var_name string, used for plotting aesthetics
#' @param prefix string, used for plotting aesthetics
#' @param e_vec vector representing the variable (i.e. column) that's the sum of \code{common_mat}
#' and \code{distinct_mat} from the \code{dcca_decomposition} function
#' @param c_vec vector representing the variable (i.e. column) in \code{common_mat}
#'  from the \code{dcca_decomposition} function
#' @param d_vec vector representing the variable (i.e. column) in \code{distinct_mat}
#'  from the \code{dcca_decomposition} function
#' @param e_res result applying \code{e_vec} to \code{compute_smooth_signal}
#' @param c_res result applying \code{c_vec} to \code{compute_smooth_signal}
#' @param d_res result applying \code{d_vec} to \code{compute_smooth_signal}
#' @param sig_digits integer
#'
#' @return list of 6 \code{gg} objects
#' @export
plot_laplacian <- function(seurat_obj, var_name, prefix = "RNA", 
                           e_vec, c_vec, d_vec, 
                           e_res, c_res, d_res,
                           sig_digits = 2){
  gene_mat <- cbind(e_vec, c_vec, d_vec, 
                    e_res$smoothed_vec, c_res$smoothed_vec, d_res$smoothed_vec)
  rownames(gene_mat) <- colnames(seurat_obj)
  colnames(gene_mat) <- paste0("gene", 1:6)
  seurat_obj[["gene"]] <- Seurat::CreateDimReducObject(embedding = gene_mat, 
                                                       key = "gene", assay = "RNA")
  
  e_range <- range(c(e_vec, e_res$smoothed_vec))
  c_range <- range(c(c_vec, c_res$smoothed_vec))
  d_range <- range(c(d_vec, d_res$smoothed_vec))
  
  plot1 <- Seurat::FeaturePlot(seurat_obj, features = "gene_1", reduction = "everything")
  plot1 <- plot1 + ggplot2::labs(title = paste0(var_name, "\n", prefix, " Denoised, Everything"), x = "Everything 1", y = "Everything 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = e_range)
  plot2 <- Seurat::FeaturePlot(seurat_obj, features = "gene_2", reduction = "common")
  plot2 <- plot2 + ggplot2::labs(title = paste0(var_name, "\n", prefix, " Denoised, Common"), x = "Common 1", y = "Common 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c_range)
  plot3 <- Seurat::FeaturePlot(seurat_obj, features = "gene_3", reduction = "distinct")
  plot3 <- plot3 + ggplot2::labs(title = paste0(var_name, "\n", prefix, " Denoised, Distinct"), x = "Distinct 1", y = "Distinct 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = d_range)
  
  plot4 <- Seurat::FeaturePlot(seurat_obj, features = "gene_4", reduction = "everything")
  plot4 <- plot4 + ggplot2::labs(title = paste0("Smoothed Everything\nVar: ", round(e_res$variance, sig_digits), ", R2: ", round(e_res$r_squared, 2)), 
                                 x = "Everything 1", y = "Everything 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = e_range)
  plot5 <- Seurat::FeaturePlot(seurat_obj, features = "gene_5", reduction = "common")
  plot5 <- plot5 + ggplot2::labs(title = paste0("Smoothed Common\nVar: ", round(c_res$variance, sig_digits), ", R2: ", round(c_res$r_squared, 2)), 
                                 x = "Common 1", y = "Common 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c_range)
  plot6 <- Seurat::FeaturePlot(seurat_obj, features = "gene_6", reduction = "distinct")
  plot6 <- plot6 + ggplot2::labs(title = paste0("Smoothed Distinct\nVar: ", round(d_res$variance, sig_digits), ", R2: ", round(d_res$r_squared, 2)), 
                                 x = "Distinct 1", y = "Distinct 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = d_range)
  
  list(plot1, plot2, plot3, plot4, plot5, plot6)
  
  # p <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)
  # cowplot::save_plot(filename = filename, p, ncol = 3, nrow = 2, base_asp = 1.2, device = "png")
  # 
  # invisible()
}


#' Plot the results of compute_smooth_signal
#'
#' @param val_vec vector of values to plot
#' @param name_vec name of the entries of \code{val_vec} as another vector of strings
#' of the same length
#' @param factor_vec vector of factors (among values \code{"0"}, \code{"1"} and \code{"2"})
#' of the same length as \code{val_vec} where entries with value \code{"1"}
#' represent the entries to be labeled via \code{ggrepel::geom_text_repel},
#'  entries with value \code{"2"} represent other entries to be colored (if any),
#'  and entries with value \code{"0"} all remaining entries
#' @param col_vec vector of 2 or 3 colors, depending on the number of levels in
#' \code{factor_vec}
#' @param xlab string
#' @param ylab string
#' @param main string
#' @param baseline_col color for the horizontal line
#' @param ... extra parameters for \code{ggrepel::geom_text_repel}
#'
#' @return \code{gg} object
#' @export
plot_laplacian_variables <- function(val_vec, name_vec, factor_vec, col_vec,
                                     xlab = "Order of variables", ylab, main,
                                     baseline_col = "orange", ...){
  stopifnot(length(val_vec) == length(name_vec), length(val_vec) == length(factor_vec),
            length(col_vec) %in% c(1,2,3), is.factor(factor_vec),
            length(levels(factor_vec)) == length(col_vec),
            all(as.character(factor_vec) %in% c("0", "1", "2")))
  p <- length(val_vec)
  if(length(names(col_vec)) == 0) names(col_vec) <- as.character(c(1:length(col_vec)))
  custom_colors <- ggplot2::scale_colour_manual(values = col_vec)
  
  # construct df
  size_vec <- rep(1, p); size_vec[which(factor_vec == 1)] <- 2
  df <- data.frame("val" = val_vec, "idx" = rank(val_vec),
                   "factor" = factor_vec,
                   "size" = size_vec, "name" = name_vec)
  val = idx = factor = size = name = NULL # for appeasing R CHECK
  
  # reorder df
  df <- df[c(which(df$factor == 0), which(df$factor == 2), which(df$factor == 1)),]
  
  # plot
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = idx, y = val))
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = factor, size = size))
  plot1 <- plot1 + ggplot2::scale_size_continuous(range = c(1, 2))
  plot1 <- plot1 + custom_colors
  plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, factor == 1), 
                                            ggplot2::aes(label = name, color = factor),
                                            box.padding = ggplot2::unit(0.5, 'lines'),
                                            point.padding = ggplot2::unit(1.6, 'lines'),
                                            size = 2, ...)
  plot1 <- plot1 + ggplot2::geom_hline(yintercept = 0, linetype="dashed",
                                       color = baseline_col)
  plot1 <- plot1 + ggplot2::xlab(xlab)
  plot1 <- plot1 + ggplot2::ylab(ylab)
  plot1 <- plot1 + ggplot2::ggtitle(main)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1
}