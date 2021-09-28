.common_decomposition <- function(score_1, score_2, 
                                  snn_union, 
                                  fix_distinct_perc,
                                  discretization_gridsize,
                                  iterations,
                                  tol = 1e-6, verbose = F){
  stopifnot(!any(is.na(snn_union)) & !fix_distinct_perc | 
              all(is.na(snn_union)) & fix_distinct_perc)
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing distinct percentage"))
  if(fix_distinct_perc){
    distinct_perc_2 <- 0.5
  } else {
    distinct_perc_2 <- .search_distinct_perc_2()
    
    
    sapply(1:rank_c, function(k){
      if(verbose){
        if(rank_c < 10) { cat('*') } else if(k %% floor(rank_c/10) == 0) cat('*')
      } 
      else {
        .latent_distinct_perc_2(score_1[,k], score_2[,k], nn_1, nn_2, verbose = F)
      }
    })
  }
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing common score"))
  common_score <- sapply(1:rank_c, function(k){
    if(verbose){
      if(rank_c < 10) { cat('*') } else if(k %% floor(rank_c/10) == 0) cat('*')
    } 
    
    vec1 <- basis_list[[k]]$rep1; vec2 <- basis_list[[k]]$rep2
    
    if(sum(abs(vec1 - vec2)) < tol){
      common_rep <- c(vec1 + vec2)/2
    } else {
      circle <- .construct_circle(vec1, vec2)
      stopifnot(sum(abs(vec1 - rep(1,0))) <= 1e-6) ## specific to .representation_2d
      rad1 <- .find_radian(circle, vec1); rad2 <- .find_radian(circle, vec2)
      stopifnot(rad1 < 0 & rad2 > 0) # must be true based on how we constructed vec1 and vec2
      
      stopifnot(sum(abs(.position_from_circle(circle, rad1) - vec1)) <= 1e-6,
                sum(abs(.position_from_circle(circle, rad2) - vec2)) <= 1e-6)
      
      rad1 <- rad1 + 2*pi # to ensure rad1 is larger than rad2
      common_rad <- .binary_search_radian(circle, rad2, rad1, distinct_perc_2)
      common_rep <- .position_from_circle(circle, common_rad)
    }
    
    basis_list[[k]]$basis_mat %*% common_rep
  })
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, distinct_perc_2 = distinct_perc_2)
}

#####################

.search_distinct_perc_2 <- function(score_1,
                                    score_2,
                                    snn_union,
                                    discretization_gridsize,
                                    iterations,
                                    tol = 1e-3){
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
  percentage_grid_all <- seq(0, 1, length.out = discretization_gridsize)
  value_vec_all <- numeric(0)
  
  iter <- 1
  
  while(iter < iterations){
    if(iter == 1){
      percentage_grid <- percentage_grid_all
      value_vec <- rep(NA, length(discretization_gridsize))
    } else {
      # determine grid and grab old values
      percentage_grid <- seq(percentage_grid_all[idx1],
                             percentage_grid_all[idx2],
                             length.out = length(discretization_gridsize))
      value_vec <- .grab_previous_values(percentage_grid,
                                         percentage_grid_all,
                                         value_vec_all)
    }
    
    # evaluate
    for(i in which(is.na(value_vec))){
      value_vec[i] <- .evaluate_radian(percentage_grid[i],
                                       basis_list = basis_list, 
                                       snn_union = snn_union,
                                       num_neigh = num_neigh,
                                       score_1 = score_1,
                                       score_2 = score_2,
                                       svd_1 = svd_1, 
                                       svd_2 = svd_2,
                                       circle_list = circle_list,
                                       return_common_score = F)
    }
    
    # update_values
    tmp <- .update_values(percentage_grid, percentage_grid_all,
                          value_vec, value_vec_all)
    percentage_grid_all <- tmp$percentage_grid_all
    value_vec_all <- tmp$value_vec_all
    idx_max <- which.max(percentage_grid_all)
    idx_1 <- max(idx_max - 1, 1)
    idx_2 <- min(idx_max + 1, length(percentage_grid_all))
    
    iter <- iter + 1
  }
  
  # compute common score based on percentage_grid_all[idx_max]
  .evaluate_radian(percentage_grid_all[idx_max],
                   basis_list = basis_list, 
                   snn_union = snn_union,
                   num_neigh = num_neigh,
                   score_1 = score_1,
                   score_2 = score_2,
                   svd_1 = svd_1, 
                   svd_2 = svd_2,
                   circle_list = circle_list,
                   return_common_score = T)
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
                             basis_list, 
                             snn_union,
                             num_neigh,
                             score_1,
                             score_2,
                             svd_1, 
                             svd_2,
                             circle_list,
                             return_common_score){
  r <- length(basis_list)
  
  radian_vec <- sapply(1:r, function(k){
    .compute_radian(percentage_val, 
                    basis_list[[k]],
                    circle_list[[k]])
  })
  
  common_representation <- sapply(1:r, function(k){
    .position_from_circle(circle_list[[k]], radian_vec[k])
  })
  
  common_score <- t(sapply(1:r, function(k){
    basis_list[[k]]$basis_mat %*% common_representation[,k]
  }))
  
  if(return_common_score){
    return(common_score)
  }
  
  common_mat <- .convert_common_score_to_mat(common_score,
                                             score_1,
                                             score_2,
                                             svd_1, 
                                             svd_2)
  
  # compute nn's
  snn_common <- .form_snn(common_mat, num_neigh = num_neigh)
  
  # compute intersection
  .computer_overlap(snn_target = snn_union,
                    snn_query = snn_common)
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
  
  idx <- duplicated(percentage_grid_all)
  if(length(idx) > 0){
    percentage_grid_all <- percentage_grid_all[-idx]
    value_vec_all <- value_vec_all[-idx]
  }
  
  list(percentage_grid_all = percentage_grid_all,
       value_vec_all = value_vec_all)
}

############################################

.compute_radian <- function(percentage_val, 
                            basis,
                            circle){
  stopifnot(percentage_val >= 0, percentage_val <= 1,
            is.list(basis), is.list(circle),
            all(sort(names(basis)) == sort(c("rep1", "rep2", "basis_mat"))),
            all(sort(names(circle)) == sort(c("center", "radius"))))
  
  vec1 <- basis$rep1; vec2 <- basis$rep2
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
