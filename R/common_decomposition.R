.common_decomposition <- function(frnn_union,
                                  num_neigh,
                                  score_1,
                                  score_2,
                                  svd_1, 
                                  svd_2,
                                  fix_distinct_perc,
                                  radius_quantile,
                                  discretization_gridsize,
                                  iterations,
                                  tol = 1e-6, verbose = F){
  stopifnot(!any(is.na(frnn_union)) & !fix_distinct_perc | 
              all(is.na(frnn_union)) & fix_distinct_perc)
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing distinct percentage"))
  if(fix_distinct_perc){
    distinct_perc <- 0.5
  } else {
    tmp <- .search_distinct_perc(
      score_1 = score_1,
      score_2 = score_2,
      frnn_union = frnn_union,
      discretization_gridsize = discretization_gridsize,
      iterations = iterations,
      basis_list = basis_list
    )
    distinct_perc <- tmp$percentage
    df_percentage <- tmp$df
  }
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing common score"))
  # compute common score based on percentage_grid_all[idx_max]
  circle_list <- lapply(1:r, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  common_score <- .evaluate_radian(
    distinct_perc,
    frnn_union = frnn_union,
    num_neigh = num_neigh,
    score_1 = score_1,
    score_2 = score_2,
    svd_1 = svd_1, 
    svd_2 = svd_2,
    radius_quantile = radius_quantile,
    basis_list = basis_list, 
    circle_list = circle_list,
    return_common_score = T
  )
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, 
       distinct_perc = distinct_perc,
       df_percentage = df_percentage)
}

#####################

# here, a distinct_perc of 0 means that the common space is
# aligned with modality 2, meaning there is no distinct information
# in modality 2. 
# Hence, favor_modality_1 is boolean that corresponds to
# "favor_start" -- you want to find the maximizing overlap
# with as close to distinct_perc 0 as possible.
.search_distinct_perc <- function(num_neigh,
                                  score_1,
                                  score_2,
                                  svd_1,
                                  svd_2,
                                  radius_quantile,
                                  frnn_union,
                                  discretization_gridsize,
                                  iterations,
                                  basis_list,
                                  favor_modality_1,
                                  tol = 1e-3){
  stopifnot(is.logical(favor_modality_1))
  
  r <- length(basis_list)
  
  # handle corner case
  if(sum(sapply(1:r, function(k){
    sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2))
  })) <= tol) return(0.5)
  
  # initialize values
  circle_list <- lapply(1:r, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  percentage_grid_all <- numeric(0)
  value_vec_all <- numeric(0)
  
  iter <- 1
  
  while(iter <= iterations){
    if(iter == 1){
      percentage_grid <- seq(0, 1, length.out = discretization_gridsize)
      value_vec <- rep(NA, discretization_gridsize)
    } else {
      # determine grid and grab old values
      percentage_grid <- seq(percentage_grid_all[idx_1],
                             percentage_grid_all[idx_2],
                             length.out = discretization_gridsize)
      value_vec <- .grab_previous_values(percentage_grid,
                                         percentage_grid_all,
                                         value_vec_all)
    }
    
    # evaluate
    for(i in which(is.na(value_vec))){
      value_vec[i] <- .evaluate_radian(percentage_grid[i],
                                       basis_list = basis_list, 
                                       frnn_union = frnn_union,
                                       num_neigh = num_neigh,
                                       score_1 = score_1,
                                       score_2 = score_2,
                                       svd_1 = svd_1, 
                                       svd_2 = svd_2,
                                       radius_quantile = radius_quantile,
                                       circle_list = circle_list,
                                       return_common_score = F,
                                       return_frnn = F)
    }
    
    # update_values
    tmp <- .update_values(percentage_grid, percentage_grid_all,
                          value_vec, value_vec_all)
    percentage_grid_all <- tmp$percentage_grid_all
    value_vec_all <- tmp$value_vec_all
    perc_max <- .picking_maximizing_value(x_val = percentage_grid_all,
                                          y_val = value_vec_all,
                                          favor_start = favor_modality_1)
    idx_max <- which(percentage_grid_all == perc_max)
    idx_1 <- max(idx_max - 1, 1)
    idx_2 <- min(idx_max + 1, length(percentage_grid_all))
  
    iter <- iter + 1
  }
  
  df <- data.frame(percentage = percentage_grid_all,
                   overlap = value_vec_all)
  list(df = df, percentage = percentage_grid_all[idx_max])
}

############################################

.grab_previous_values <- function(percentage_grid,
                                  percentage_grid_all,
                                  value_vec_all,
                                  tol = 1e-6){
  len <- length(percentage_grid)
  value_vec <- rep(NA, len)
  
  for(i in 1:len){
    idx <- which(abs(percentage_grid[i] - percentage_grid_all) <= tol)
    if(length(idx) > 0){
      value_vec[i] <- value_vec_all[idx]
    }
  }
  
  value_vec
}

.evaluate_radian <- function(percentage_val,
                             frnn_union,
                             num_neigh,
                             score_1,
                             score_2,
                             svd_1, 
                             svd_2,
                             radius_quantile,
                             basis_list,
                             circle_list,
                             return_common_score,
                             return_frnn){
  r <- length(basis_list)
  
  radian_vec <- sapply(1:r, function(k){
    .compute_radian(percentage_val, 
                    vec1 = basis_list[[k]]$rep1,
                    vec2 = basis_list[[k]]$rep2,
                    circle = circle_list[[k]])
  })
  
  common_representation <- sapply(1:r, function(k){
    .position_from_circle(circle_list[[k]], radian_vec[k])
  })
  
  common_score <- sapply(1:r, function(k){
    basis_list[[k]]$basis_mat %*% common_representation[,k]
  })
  
  if(return_common_score){
    return(common_score)
  }
  
  common_mat <- .convert_common_score_to_mat(common_score,
                                             score_1,
                                             score_2,
                                             svd_1, 
                                             svd_2)
  
 
  
  # compute nn's
  frnn_common <- .nnlist_to_matrix(
    .construct_frnn(common_mat, 
                    radius = NA,
                    nn = num_neigh, 
                    frnn_approx = 0, 
                    resolve_isolated_nodes = T,
                    radius_quantile = 0.5,
                    verbose = F), set_to_one = T)
  
  if(return_frnn){
    return(frnn_common)
  }
  
  # compute intersection
  .computer_overlap(snn_target = frnn_union,
                    snn_query = frnn_common)
}

.update_values <- function(percentage_grid, percentage_grid_all,
                           value_vec, value_vec_all){
  stopifnot(all(!is.na(value_vec)), all(!is.na(value_vec_all)),
            length(percentage_grid) == length(value_vec),
            length(percentage_grid_all) == length(value_vec_all))
  
  vec1 <- c(percentage_grid, percentage_grid_all)
  vec2 <- c(value_vec, value_vec_all)
  
  idx <- order(vec1, decreasing = F)
  percentage_grid_all <- vec1[idx]
  value_vec_all <- vec2[idx]
  
  idx <- which(duplicated(percentage_grid_all))
  if(length(idx) > 0){
    percentage_grid_all <- percentage_grid_all[-idx]
    value_vec_all <- value_vec_all[-idx]
  }
  
  list(percentage_grid_all = percentage_grid_all,
       value_vec_all = value_vec_all)
}

.picking_maximizing_value <- function(x_val, y_val, 
                                      favor_start){
  stopifnot(is.logical(favor_start), 
            length(x_val) == length(y_val))
  
  if(!favor_start){
    x_val <- rev(x_val); y_val <- rev(y_val)
  }
  
  diff_vec <- diff(y_val)
  idx <- which(diff_vec < 0)
  if(length(idx) == 0){
    return(x_val[length(x_val)])
  }
  return(x_val[idx[1]])
}

############################################

.compute_radian <- function(percentage_val, 
                            vec1,
                            vec2,
                            circle){
  stopifnot(percentage_val >= 0, percentage_val <= 1,
            is.list(circle),
            all(sort(names(circle)) == sort(c("center", "radius"))),
            abs(.l2norm(circle$center - vec1) - circle$radius) <= 1e-6,
            abs(.l2norm(circle$center - vec2) - circle$radius) <= 1e-6)
  
  rad1 <- .find_radian(circle, vec1)
  rad2 <- .find_radian(circle, vec2)
  stopifnot(rad1 < 0 & rad2 > 0) # must be true based on how we constructed vec1 and vec2
  
  rad1 <- rad1 + 2*pi # to ensure rad1 is larger than rad2
  rad2 + (rad1 - rad2)*percentage_val
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
  
  common_1 <- common_score %*% crossprod(score_1, dimred_1)/n
  common_2 <- common_score %*% crossprod(score_2, dimred_2)/n
  
  cbind(common_1, common_2)
}

.computer_overlap <- function(snn_target,
                              snn_query){
  snn <- snn_target + snn_query
  length(which(snn@x == 2))/length(snn_target@x)
}
