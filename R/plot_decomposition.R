plot_decomposition_2d <- function(vec1, vec2, common_vec, 
                                  gridsize = 100, col_levels = 21, 
                                  plot_bg = T,
                                  renormalize = T, 
                                  xlim = range(c(0, 1.1*c(vec1, vec2, common_vec))),
                                  ylim = range(c(0, 1.1*c(vec1, vec2, common_vec))), 
                                  length_arrow = 0.1,
                                  lwd_arrow_black = 1.2,
                                  lwd_arrow_redgreen = 2,
                                  lwd_arrow_white = 3, 
                                  ...){
  if(renormalize){
    scalar <- .l2norm(vec1);
    vec1 <- vec1/scalar; vec2 <- vec2/scalar; common_vec <- common_vec/scalar
  }
 
  basis <- .representation_2d(vec1, vec2)
  rep_vec1 <- basis$rep1; rep_vec2 <- basis$rep2
  rep_common <- as.numeric(crossprod(basis$basis_mat, common_vec))
  
  stopifnot(sum(abs(basis$basis_mat %*% rep_common - common_vec)) <= 1e-6)
  
  .plot_decomposition_2d_inner(rep_vec1, rep_vec2, rep_common,
                               gridsize = gridsize,
                               col_levels = col_levels, 
                               plot_bg = plot_bg, 
                               xlim = xlim,
                               ylim = ylim,
                               length_arrow = length_arrow,
                               lwd_arrow_black = lwd_arrow_black,
                               lwd_arrow_redgreen = lwd_arrow_redgreen,
                               lwd_arrow_white = lwd_arrow_white, 
                               ...)
}

#################3

# vec1, vec2 and common_vec are 2d
.plot_decomposition_2d_inner <- function(vec1, vec2, common_vec,
                                  xlim = range(c(0, 1.1*c(vec1, vec2, common_vec))),
                                  ylim = range(c(0, 1.1*c(vec1, vec2, common_vec))), 
                                  gridsize = 100, col_levels = 21,
                                  plot_bg = T, 
                                  length_arrow = 0.1,
                                  lwd_arrow_black = 1.2,
                                  lwd_arrow_redgreen = 2,
                                  lwd_arrow_white = 3, 
                                  ...){
  stopifnot(length(vec1) == 2, length(vec2) == 2,
            all(c(vec1, vec2) >= 0))
  
  # relabel so vec1 is always shorter
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2); flipped <- F
  if(len1 > len2){ flipped <- T; tmp <- vec2; vec2 <- vec1; vec1 <- tmp }
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  ratio <- len1/len2
  
  x_seq <- seq(xlim[1], xlim[2], length.out = gridsize)
  y_seq <- seq(ylim[1], ylim[2], length.out = gridsize)
  grid_mat <- .construct_grid(x_seq, y_seq)
  values <- sapply(1:nrow(grid_mat), function(i){
    .compute_attempted_decomposition(vec1, vec2, common_vec = grid_mat[i,])
  })
  mat <- .arrange_gridvalues(x_seq, y_seq, values)
  
  side_length <- floor(col_levels/2); max_val <- max(abs(mat), na.rm = T)
  col_vec <- c(grDevices::hcl.colors(side_length, palette = "YlGnBu"),
               "gray75",
               rev(grDevices::hcl.colors(side_length, palette = "OrRd")))
  graphics::plot(NA, xlim = range(x_seq), ylim = range(y_seq), 
                 asp = T, ...)
  if(plot_bg){
    graphics::.filled.contour(x = x_seq, y = y_seq,
                              z = mat, levels = c(rev(seq(0, -max_val, length.out = side_length+2)[-1]), seq(0, max_val, length.out = side_length+2)[-1]),
                              col = col_vec)
  }
  graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], 
                   length = length_arrow, col = "white", lwd = lwd_arrow_white)
  graphics::arrows(x0 = common_vec[1], y0 = common_vec[2], x1 = vec1[1], y1 = vec1[2], 
                   length = length_arrow, col = "white", lwd = lwd_arrow_white)
  graphics::arrows(x0 = common_vec[1], y0 = common_vec[2], x1 = vec2[1], y1 = vec2[2], 
                   length = length_arrow, col = "white", lwd = lwd_arrow_white)
  
  graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], 
                   length = length_arrow, col = 3, lwd = lwd_arrow_redgreen)
  graphics::arrows(x0 = common_vec[1], y0 = common_vec[2], x1 = vec1[1], y1 = vec1[2], 
                   length = length_arrow, col = 2, lwd = lwd_arrow_redgreen)
  graphics::arrows(x0 = common_vec[1], y0 = common_vec[2], x1 = vec2[1], y1 = vec2[2], 
                   length = length_arrow, col = 2, lwd = lwd_arrow_redgreen)
  
  graphics::arrows(x0 = 0, y0 = 0, x1 = vec1[1], y1 = vec1[2], 
                   length = length_arrow, lwd = lwd_arrow_black)
  graphics::arrows(x0 = 0, y0 = 0, x1 = vec2[1], y1 = vec2[2], 
                   length = length_arrow, lwd = lwd_arrow_black)
  
  invisible()
}

.construct_grid <- function(x_seq, y_seq){
  as.matrix(expand.grid(x_seq, y_seq))
}

.compute_attempted_decomposition <- function(vec1, vec2, common_vec){
  stopifnot(all(c(vec1, vec2) >= 0), length(vec1) == 2, length(vec2) == 2,
            length(common_vec) == 2)
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  stopifnot(len1 <= len2)
  ratio <- len1/len2
  
  common_vec2 <- common_vec
  common_vec1 <- common_vec2*ratio
  
  dis_vec1 <- vec1 - common_vec1; dis_vec2 <- vec2 - common_vec2
  
  tmp <- .angle_between_vectors(dis_vec1, dis_vec2)
  if(is.na(tmp)) return(NA)
  tmp - 90
}

.arrange_gridvalues <- function(x_seq, y_seq, values){
  stopifnot(length(values) == length(x_seq)*length(y_seq))
  
  gridsize <- length(x_seq)
  mat <- matrix(values, nrow = gridsize, ncol = gridsize)
  rownames(mat) <- x_seq
  colnames(mat) <- y_seq
  
  mat
}

