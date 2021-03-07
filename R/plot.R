# if reserve_zero = T, then reserve zero for white
# if reserve_zero = F, then all the greens are negative and all the reds are positive
plot_heatmat <- function(dat, luminosity = T, asp = nrow(dat)/ncol(dat), 
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

plot_scores2 <- function(score_1, score_2, membership_vec){
  graphics::par(mfrow = c(1,2))
  graphics::plot(NA, xlim = range(score_1), ylim = c(0.5,3.5))
  n <- nrow(score_1)
  for(i in 1:ncol(score_1)){
    graphics::points(x = score_1[,i], y = stats::runif(n, min = i-.2, max = i+.2), col = membership_vec, pch = 16)
  }
  
  graphics::plot(NA, xlim = range(score_2), ylim = c(0.5,3.5))
  n <- nrow(score_2)
  for(i in 1:ncol(score_2)){
    graphics::points(x = score_2[,i], y = stats::runif(n, min = i-.2, max = i+.2), col = membership_vec, pch = 16)
  }
  
  invisible()
}

#######################################

.colorRamp_custom <- function(vec1, vec2, length, luminosity = T){
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
