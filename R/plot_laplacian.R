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
#' @param filename string, filename to save the plot to
#'
#' @return nothing
#' @export
plot_laplacian <- function(seurat_obj, var_name, prefix = "RNA", 
                           e_vec, c_vec, d_vec, 
                           e_res, c_res, d_res,
                           filename){
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
  plot4 <- plot4 + ggplot2::labs(title = paste0("Smoothed Everything\nVar: ", round(e_res$variance, 2), ", R2: ", round(e_res$r_squared, 2)), 
                                 x = "Everything 1", y = "Everything 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = e_range)
  plot5 <- Seurat::FeaturePlot(seurat_obj, features = "gene_5", reduction = "common")
  plot5 <- plot5 + ggplot2::labs(title = paste0("Smoothed Common\nVar: ", round(c_res$variance, 2), ", R2: ", round(c_res$r_squared, 2)), 
                                 x = "Common 1", y = "Common 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c_range)
  plot6 <- Seurat::FeaturePlot(seurat_obj, features = "gene_6", reduction = "distinct")
  plot6 <- plot6 + ggplot2::labs(title = paste0("Smoothed Distinct\nVar: ", round(d_res$variance, 2), ", R2: ", round(d_res$r_squared, 2)), 
                                 x = "Distinct 1", y = "Distinct 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = d_range)
  
  p <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)
  cowplot::save_plot(filename = filename, p, ncol = 3, nrow = 2, base_asp = 1.2, device = "png")
  
  invisible()
}