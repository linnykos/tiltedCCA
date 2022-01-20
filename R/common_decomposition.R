# a tilt_perc near-0 means that most of the distinct information
# is in modality 1, not modality 2
.common_decomposition <- function(discretization_gridsize,
                                  enforce_boundary,
                                  fix_tilt_perc,
                                  score_1,
                                  score_2,
                                  svd_1, 
                                  svd_2,
                                  target_dimred,
                                  tol = 1e-6, verbose = F){
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  circle_list <- lapply(1:rank_c, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  if(verbose) print(paste0(Sys.time(),": D-CCA: (Inner) Computing distinct percentage"))
  if(is.logical(fix_tilt_perc) && !fix_tilt_perc){
    tmp <- .search_tilt_perc(
      basis_list = basis_list,
      circle_list = circle_list,
      enforce_boundary = enforce_boundary,
      discretization_gridsize = discretization_gridsize,
      score_1 = score_1,
      score_2 = score_2,
      svd_1 = svd_1,
      svd_2 = svd_2, 
      target_dimred = target_dimred,
      verbose = verbose
    )
    tilt_perc <- tmp$percentage
    df_percentage <- tmp$df
  } else {
    if(is.logical(fix_tilt_perc) && fix_tilt_perc){
      tilt_perc <- 0.5; df_percentage <- NA
    } else if(is.numeric(fix_tilt_perc) & fix_tilt_perc >= 0 & fix_tilt_perc <= 1){
      tilt_perc <- fix_tilt_perc; df_percentage <- NA
    } else {
      stop("Invalid value of tilt_perc")
    }
  }
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing common score"))
  
  common_score <- .evaluate_radian(
    basis_list = basis_list, 
    circle_list = circle_list,
    enforce_boundary = enforce_boundary,
    percentage = tilt_perc,
    return_common_score = T,
    score_1 = score_1,
    score_2 = score_2,
    svd_1 = svd_1, 
    svd_2 = svd_2,
    target_dimred = target_dimred
  )
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, 
       tilt_perc = tilt_perc,
       df_percentage = df_percentage)
}

#####################

.search_tilt_perc <- function(basis_list,
                              circle_list,
                              discretization_gridsize,
                              enforce_boundary,
                              score_1,
                              score_2,
                              svd_1,
                              svd_2,
                              target_dimred,
                              tol = 1e-3,
                              verbose = F){
  r <- length(basis_list)
  
  # handle corner case
  if(sum(sapply(1:r, function(k){
    sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2))
  })) <= tol) return(0.5)
  
  # initialize values
  percentage_grid <- seq(0, 1, length.out = discretization_gridsize)
  value_vec <- rep(NA, discretization_gridsize)
  
  for(i in which(is.na(value_vec))){
    if(verbose) print(paste0(Sys.time(),": D-CCA : Evaluating percentage ", percentage_grid[i]))
    value_vec[i] <- .evaluate_radian(basis_list = basis_list, 
                                     circle_list = circle_list,
                                     enforce_boundary = enforce_boundary,
                                     percentage = percentage_grid[i],
                                     return_common_score = F,
                                     score_1 = score_1,
                                     score_2 = score_2,
                                     svd_1 = svd_1, 
                                     svd_2 = svd_2,
                                     target_dimred = target_dimred)
  }
  
  df <- data.frame(percentage = percentage_grid,
                   ratio_val = value_vec)
  idx_min <- .select_minimum(minimum = T,
                             x_val = percentage_grid,
                             y_val = value_vec)
  if(verbose) print(paste0(Sys.time(),": D-CCA : Selected tilt-percentage to be: ", percentage_grid[idx_min]))
  
  list(df = df, percentage = percentage_grid[idx_min])
}

############################################

.evaluate_radian <- function(basis_list, 
                             circle_list,
                             enforce_boundary,
                             percentage,
                             return_common_score,
                             score_1,
                             score_2,
                             svd_1, 
                             svd_2,
                             target_dimred){
  r <- length(basis_list)
  
  radian_vec <- sapply(1:r, function(k){
    .compute_radian(circle = circle_list[[k]],
                    enforce_boundary = enforce_boundary,
                    percentage_val = percentage, 
                    vec1 = basis_list[[k]]$rep1,
                    vec2 = basis_list[[k]]$rep2)
  })
  
  common_representation <- sapply(1:r, function(k){
    .position_from_circle(circle_list[[k]], radian_vec[k])
  })
  
  common_score <- sapply(1:r, function(k){
    basis_list[[k]]$basis_mat %*% common_representation[,k]
  })
  if(return_common_score) return(common_score)
  
  common_mat <- .convert_common_score_to_mat(common_score,
                                             score_1,
                                             score_2,
                                             svd_1, 
                                             svd_2)
  
  .grassmann_distance(orthonormal_1 = common_mat, 
                      orthonormal_2 = target_dimred)
}

.select_minimum <- function(minimum, x_val, y_val){
  stopifnot(length(x_val) == length(y_val))
  if(!minimum) y_val <- -y_val
  
  min_val <- min(y_val)
  min_idx <- which(y_val == min_val)
  if(length(min_idx) == 1){
    return(min_idx)
  }
  
  len <- length(x_val)
  extreme_pos <- sapply(min_idx, function(x){
    min(x_val[x]-min(x_val), max(x_val) - x_val[x])
  })
  min_pos <- which(extreme_pos == min(extreme_pos))
  if(length(min_pos) > 1){
    warning("Unable to determine best tilt. Returning the median.")
    mid <- stats::median(range(x_val))
    return(which.min(abs(x_val - mid)))
  } else {
    return(min_idx[min_pos])
  }
}

############################################

.compute_radian <- function(circle,
                            enforce_boundary,
                            percentage_val, 
                            vec1,
                            vec2){
  stopifnot(percentage_val >= 0, percentage_val <= 1,
            is.list(circle),
            all(sort(names(circle)) == sort(c("center", "radius"))),
            abs(.l2norm(circle$center - vec1) - circle$radius) <= 1e-6,
            abs(.l2norm(circle$center - vec2) - circle$radius) <= 1e-6)
  
  if(!enforce_boundary){
    rad1 <- .find_radian(circle, vec1)
    rad2 <- .find_radian(circle, vec2)
    stopifnot(rad1 < 0 & rad2 > 0) # must be true based on how we constructed vec1 and vec2
    rad1 <- rad1 + 2*pi # to ensure rad1 is larger than rad2
  } else {
    stopifnot(abs(vec1[2]) <= 1e-6,
              circle$radius >= circle$center[2] - 1e-6) # must be true based on how we constructed vec1
    position1 <- circle$center[1] - sqrt(circle$radius^2 - circle$center[2]^2) 
    rad1 <- .find_radian(circle, c(position1, 0))
    if(rad1 < 0) rad1 <- rad1 + 2*pi
    
    inner_rad <- atan(circle$center[1]/circle$center[2])
    if(inner_rad < 0) inner_rad <- inner_rad + 2*pi
    radmid <- 3*pi/2 - inner_rad
    if(abs(rad1 - radmid) <= 1e-6) return(rad1)
    
    stopifnot(rad1 >= radmid)
    rad2 <- radmid - (rad1 - radmid)
    
    stopifnot(rad1 >= rad2)
  }
  
  tmp <- rad2 + (rad1 - rad2)*percentage_val
  tmp
}

.convert_common_score_to_mat <- function(common_score,
                                         score_1,
                                         score_2,
                                         svd_1, 
                                         svd_2){
  rescaling_factor <- max(c(svd_1$d, svd_2$d))
  dimred_1 <- .mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*rescaling_factor)
  dimred_2 <- .mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*rescaling_factor)
  n <- nrow(dimred_1)
  
  if(ncol(common_score) < ncol(score_1)) {
    common_score1 <- cbind(common_score, matrix(0, nrow = n, ncol = ncol(score_1)-ncol(common_score)))
  } else common_score1 <- common_score
  if(ncol(common_score) < ncol(score_2)) {
    common_score2 <- cbind(common_score, matrix(0, nrow = n, ncol = ncol(score_2)-ncol(common_score)))
  } else common_score2 <- common_score
  common_1 <- common_score1 %*% crossprod(score_1, dimred_1)/n
  common_2 <- common_score2 %*% crossprod(score_2, dimred_2)/n
  
  cbind(common_1, common_2)
}
