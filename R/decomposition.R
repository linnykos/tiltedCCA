## construct the plot given 2 vectors in 2D
## blue is angle close to 0, white is angle close to 90, red is angle close to 180
.decomposition_2d <- function(vec1, vec2, 
                             xlim = c(0, 1.1*max(c(vec1, vec2))),
                             ylim = c(0, 1.1*max(c(vec1, vec2))),
                             plotting = T, gridsize = 100, col_levels = 21){
  stopifnot(length(vec1) == 2, length(vec2) == 2,
            all(c(vec1, vec2) >= 0))
  
  # relabel so vec1 is always shorter
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2); flipped <- F
  if(len1 > len2){ flipped <- T; tmp <- vec2; vec2 <- vec1; vec1 <- tmp }
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  rescaling_factor <- max(c(len1, len2))
  vec1 <- vec1/rescaling_factor; vec2 <- vec2/rescaling_factor
  ratio <- len1/len2

  x_seq <- seq(xlim[1], xlim[2], length.out = gridsize)
  y_seq <- seq(ylim[1], ylim[2], length.out = gridsize)
  grid_mat <- .construct_grid(x_seq, y_seq)
  values <- sapply(1:nrow(grid_mat), function(i){
    .compute_attempted_decomposition(vec1, vec2, common_vec = grid_mat[i,])
  })
  mat <- .arrange_gridvalues(x_seq, y_seq, values)
  
  # find the common vector
  common_vec <- .search_common_vec(vec1, vec2, mat, grid_mat)
  
  # refine common vector using zero-order optimization to ensure distinct vectors are orthogonal
  common_vec <- .optimize_common_vector(vec1, vec2, common_vec)

  # restore original scale
  vec1 <- vec1*rescaling_factor; vec2 <- vec2*rescaling_factor
  common_vec <- common_vec*rescaling_factor
  x_seq <- x_seq*rescaling_factor; y_seq <- y_seq*rescaling_factor
  
  # plotting
  if(!plotting) { 
    common_vec2 <- common_vec; common_vec1 <- ratio*common_vec2
    distinct_vec1 <- vec1 - common_vec1; distinct_vec2 <- vec2 - common_vec2
    
    if(!flipped){
      list(common_vec1 = common_vec1, common_vec2 = common_vec2,
           distinct_vec1 = distinct_vec1, distinct_vec2 = distinct_vec2)
    } else {
      list(common_vec1 = common_vec2, common_vec2 = common_vec1,
           distinct_vec1 = distinct_vec2, distinct_vec2 = distinct_vec1)
    }
    
  } else {
    side_length <- floor(col_levels/2); max_val <- max(abs(mat), na.rm = T)
    col_vec <- c(grDevices::hcl.colors(side_length, palette = "YlGnBu"),
                 "gray75",
                 rev(grDevices::hcl.colors(side_length, palette = "OrRd")))
    graphics::plot(NA, xlim = range(x_seq), ylim = range(y_seq), 
                   xlab = "Dimension 1", ylab = "Dimension 2",
                   asp = T)
    graphics::.filled.contour(x = x_seq, y = y_seq,
                              z = mat, levels = c(rev(seq(0, -max_val, length.out = side_length+2)[-1]), seq(0, max_val, length.out = side_length+2)[-1]),
                              col = col_vec)
    
    graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], length = 0.1, col = 3, lwd = 2)
    graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1]*len1/len2, y1 = common_vec[2]*len1/len2, length = 0.1, col = 3, lwd = 2)
    graphics::arrows(x0 = common_vec[1]*len1/len2, y0 = common_vec[2]*len1/len2, x1 = vec1[1], y1 = vec1[2], length = 0.1, col = 2, lwd = 2)
    graphics::arrows(x0 = common_vec[1], y0 = common_vec[2], x1 = vec2[1], y1 = vec2[2], length = 0.1, col = 2, lwd = 2)
    
    graphics::arrows(x0 = 0, y0 = 0, x1 = vec1[1], y1 = vec1[2], length = 0.1, lwd = 1.2)
    graphics::arrows(x0 = 0, y0 = 0, x1 = vec2[1], y1 = vec2[2], length = 0.1, lwd = 1.2)
    
    invisible()
  }
}

#####################################

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

.search_common_vec <- function(vec1, vec2, mat, grid_mat){
  stopifnot(nrow(grid_mat) == prod(dim(mat)), .l2norm(vec1) <= .l2norm(vec2))
  num_val <- nrow(grid_mat)
  len2 <- .l2norm(vec2)
  
  # fill out the boolean matrix
  bool_array <- array(NA, dim = c(nrow(mat), ncol(mat), 4))
  
  # which vectors are not too long
  bool_array[,,1] <- matrix(sapply(1:num_val, function(i){
    .l2norm(grid_mat[i,]) <= len2
  }), nrow = nrow(mat), ncol = ncol(mat))
  
  # which vectors are within the angle
  ang1 <- .angle_between_vectors(vec1, c(1,0))
  ang2 <- .angle_between_vectors(vec2, c(1,0))
  bool_array[,,2] <- matrix(sapply(1:num_val, function(i){
    ang <- .angle_between_vectors(grid_mat[i,], c(1,0))
    ang <= max(c(ang1, ang2)) & ang >= min(c(ang1, ang2))
  }), nrow = nrow(mat), ncol = ncol(mat))
  
  # which vectors have close-to-orthogonal values (try 5 different values to determine the range)
  ang_vec <- seq(min(c(ang1, ang2)), max(c(ang1, ang2)), length.out = 7)[3:5]
  val_vec <- sapply(ang_vec, function(ang){
    .search_min_along_angle(mat, ang)
  })
  bool_array[,,3] <- (abs(mat) <= max(val_vec))
  
  bool_array[,,4] <- !is.na(mat)
  bool_mat <- apply(bool_array, c(1,2), all)
  
  # reorder all the indices by their angle
  true_idx <- which(bool_mat)
  if(length(true_idx) == 0){
    # emergency case
    return(c((vec1[1]+vec2[1])/2, (vec1[2]+vec2[2])/2))
  }
  grid_true <- grid_mat[true_idx,,drop = F]
  if(nrow(grid_true) > 1){
    ang_vec <- sapply(1:nrow(grid_true), function(i){
      .angle_between_vectors(grid_true[i,], c(1,0))
    })
    grid_true <- grid_true[order(ang_vec),,drop = F]
  }
  
  # compute the common percentages among the selected values
  .max_common_cor(vec1, vec2, grid_true)
}

.search_min_along_angle <- function(mat, ang, gridsize = 100){
  rad <- ang * pi/180
  unit_vec <- c(cos(rad), sin(rad))
  xmax <- max(as.numeric(rownames(mat)))
  ymax <- max(as.numeric(colnames(mat)))
  
  unit_vec1 <- unit_vec * (xmax/unit_vec[1])
  unit_vec2 <- unit_vec * (ymax/unit_vec[2])
  
  if(.l2norm(unit_vec1) < .l2norm(unit_vec2)) search_vec <- unit_vec1 else search_vec <- unit_vec2
  perc_seq <- seq(0,1,length.out = gridsize)
  
  values <- sapply(perc_seq, function(perc){
    .search_min_in_direction_from_mat(search_vec*perc, mat)
  })

  values[.search_first_min(values)]
}

.search_min_in_direction_from_mat <- function(vec, mat){
  x_seq <- as.numeric(rownames(mat))
  y_seq <- as.numeric(colnames(mat))

  idx1 <- which.min(abs(vec[1] - x_seq)); idx2 <- which.min(abs(vec[2] - y_seq))
  abs(mat[idx1, idx2])
}

.search_first_min <- function(vec){
  lower_quant <- stats::quantile(vec, prob = 0.05, na.rm = T)
  idx <- which(vec <= lower_quant)

  # find all the first contiguous block
  if(length(idx) > 1){
    tmp <- idx[1]
    for(i in 2:length(idx)){
      if(idx[i] == idx[i-1]+1){
        tmp <- c(tmp, idx[i])
      } else {
        break()
      }
    }

    tmp[which.min(vec[tmp])]
  } else {
    idx[1]
  }
}

.max_common_cor <- function(vec1, vec2, grid_true){
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  stopifnot(len1 <= len2)
  ratio <- len1/len2
  
  common_cor <- sapply(1:nrow(grid_true), function(i){
    common_vec2 <- grid_true[i,]
    common_vec1 <- ratio*common_vec2
    
    (len1/(len1+len2))*.cor_vectors(common_vec1, vec1)^2 + (len2/(len1+len2))*.cor_vectors(common_vec2, vec2)^2
  })
  
  grid_true[which.max(common_cor),]
}

.optimize_common_vector <- function(vec1, vec2, common_vec, 
                                    min_val = 0.9, max_val = 1.1){
  len1 <- .l2norm(vec1); len2 <- .l2norm(vec2)
  stopifnot(len1 <= len2)
  ratio <- len1/len2
  
  fun <- function(x){
    common_vec2 <- x*common_vec
    common_vec1 <- ratio*common_vec2
    
    distinct_vec1 <- vec1 - common_vec1
    distinct_vec2 <- vec2 - common_vec2
    
    abs(.angle_between_vectors(distinct_vec1, distinct_vec2)-90)
  }
  
  res <- stats::optimize(fun, interval = c(min_val, max_val),
                         maximum = F)
  
  res$minimum * common_vec
}