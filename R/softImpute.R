#' softImpute diagnostic
#'
#' @param mat data matrix
#' @param K desired rank
#' @param num_val number to omit per row or per column
#' @param lambda numeric
#'
#' @return list
#' @export
### NOTE TO SELF: Honestly, if we replace all the missing values with their row-column mean and average, 
### followed by an SVD, I'm sure it'll be fine
softImpute_diagnostic <- function(mat, K, num_val = ceiling(min(c(4, dim(mat)/10))),
                                  lambda = NA){
  stopifnot(K > 1)
  
  n <- nrow(mat); p <- ncol(mat)
  idx <- .construct_missing_values(n, p, num_val = num_val)
  mat2 <- mat
  mat2[idx] <- NA
  
  if(is.na(lambda)){
    lambda0_val <- softImpute::lambda0(mat)
    lambda <- min(30, lambda0_val/100)
  }
  res <- softImpute::softImpute(mat2, rank.max = K, lambda = lambda)
  
  pred_val <- res$u %*% diag(res$d) %*% t(res$v)
  
  training_mat <- cbind(mat[-idx], pred_val[-idx])
  testing_mat <- cbind(mat[idx], pred_val[idx])
  
  colnames(training_mat) <- c("observed_val", "predicted_val")
  colnames(testing_mat) <- c("observed_val", "predicted_val")
  
  list(training_mat = training_mat, testing_mat = testing_mat)
}

#' Plotting diagnostic to determine goodness of fit
#'
#' @param diagnostic_mat one of the output matarices for \code{softImpute_diagnostic}
#' @param width parameter, controlling quantile of prediction region. The default is \code{0.8}, meaning
#' that the prediction region is by default from the 10th to 90th quantile.
#' @param scalar 2 standard deviations for a Gaussian
#' @param max_points maximum number of points to be shown in the scatterplot, purely for visualization purposes only
#' @param tol parameter between \code{0} and \code{1} for how strict (\code{1} being the strictest) to measure
#' if the principal angle falls within the prediction region
#' @param xlim plotting parameter
#' @param ylim plotting parameter
#' @param transparency plotting parameter
#' @param cex_text plotting parameter
#' @param ... additional plotting parameters
#'
#' @return a plot but nothing is returned to console
#' @export
plot_prediction_against_observed <- function(diagnostic_mat, 
                                             width = 0.8, scalar = NA, 
                                             max_points = 500000, tol = 0.95, xlim = NA,
                                             ylim = NA, transparency = 0.2, cex_text = 1, ...){
  # compute the principal angle and
  angle_val <- .compute_principal_angle(diagnostic_mat)
  
  res <- .within_prediction_region(1.5*max(diagnostic_mat), width = width,
                                   scalar = scalar, angle_val = angle_val, tol = tol,
                                   effective_max = max(diagnostic_mat[,"predicted_val"]))
  
  if(nrow(diagnostic_mat) > max_points){
    diagnostic_mat <- diagnostic_mat[sample(1:nrow(diagnostic_mat), max_points),]
  }
  
  .plot_pca_diagnostic(diagnostic_mat, seq_vec = res$seq_vec, interval_mat = res$interval_mat,
                         principal_line = res$principal_line, angle_val = angle_val,
                         xlim = xlim, ylim = ylim, transparency = transparency,
                         cex_text = cex_text, ...)
  
  invisible()
}

###########

#' Construct indices for missing values
#'
#' @param n number of rows of intended matrix
#' @param p number of columns of intended matrix
#' @param num_val number of values to make missing in each row/column
#'
#' @return vector of indices between 1 and \code{n*p}
.construct_missing_values <- function(n, p, num_val = ceiling(min(c(4, n/10, p/10)))){
  missing_mat <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:p, num_val))
  }))), do.call(rbind, (lapply(1:p, function(x){
    cbind(sample(1:n, num_val), x)
  }))))
  
  unique(apply(missing_mat, 1, function(x){
    (x[2]-1)*n + x[1]
  }))
}

#' Compute principal angle
#'
#' @param tmp_mat a matrix with \code{n} rows (for \code{n} samples) and \code{2} columns,
#' where the first column represents the observed data and the second column represents its
#' corresponding predicted values
#'
#' @return numeric
.compute_principal_angle <- function(tmp_mat){
  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  vec <- pca_res$rotation[,1]; vec <- vec/.l2norm(vec)
  if(sign(vec[1]) < 0)  vec <- -1*vec
  angle_val <- as.numeric(acos(as.numeric(c(0,1) %*% vec)))
  angle_val * 180/pi
}

.within_prediction_region <- function(max_val, width, scalar, angle_val, tol = 0.95,
                                      effective_max = max_val){
  seq_vec <- seq(0, max_val, length.out = 500)
  stopifnot(any(seq_vec <= effective_max))
  
  interval_mat <- rbind(seq_vec - scalar, seq_vec + scalar)
  rownames(interval_mat) <- c("lower", "upper")
  
  principal_line <- seq_vec * tan(angle_val*pi/180)
  idx <- which(seq_vec <= effective_max)
  
  list(seq_vec = seq_vec, interval_mat = interval_mat, principal_line = principal_line)
}

.plot_pca_diagnostic <- function(tmp_mat, seq_vec, interval_mat, principal_line, angle_val,
                                 xlim = NA, ylim = NA, transparency = 0.2, cex_text = 1, ...){
  stopifnot(ncol(interval_mat) == length(principal_line))
  
  rad <- 2/5*max(tmp_mat)
  seq_max <- 2*max(tmp_mat)
  lim_vec <- range(c(0,tmp_mat))
  if(all(is.na(xlim))) xlim <- lim_vec
  if(all(is.na(ylim))) ylim <- lim_vec
  
  graphics::plot(NA, asp = T, xlim = xlim, ylim = ylim,
                 xlab = "Predicted value", ylab = "Observed value", ...)
  graphics::polygon(c(seq_vec, rev(seq_vec)), c(interval_mat["upper",], rev(interval_mat["lower",])), col = grDevices::rgb(1,0,0,0.2),
                    border = NA, density = 30, angle = -45)
  graphics::points(tmp_mat[,"predicted_val"], tmp_mat[,"observed_val"], pch = 16, col = grDevices::rgb(0,0,0,transparency))
  
  graphics::lines(rep(0, 2), c(-2*seq_max, 2*seq_max), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), rep(0, 2), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), c(-2*seq_max, 2*seq_max), col = "red", lwd = 2)
  
  graphics::lines(seq_vec, interval_mat["lower",], col = "red", lty = 2, lwd = 2)
  graphics::lines(seq_vec, interval_mat["upper",], col = "red", lty = 2, lwd = 2)
  
  graphics::lines(seq_vec, principal_line, col = "blue", lwd = 2, lty = 2)
  
  radian_seq <- seq(0, angle_val*pi/180, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  graphics::lines(x_circ, y_circ, lwd = 2, col = "white")
  graphics::lines(x_circ, y_circ, lty = 2)
  graphics::text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(angle_val, 1), " degrees"),
                 cex = cex_text)
  
  invisible()
}