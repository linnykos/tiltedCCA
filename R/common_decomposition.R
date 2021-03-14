.common_decomposition <- function(score_1, score_2, nn_1 = NA, nn_2 = NA, 
                                  fix_common_perc = F, tol = 1e-6){
  stopifnot((!any(is.na(nn_1)) & !any(is.na(nn_2)) & !fix_common_perc) | 
              (all(is.na(nn_1)) & all(is.na(nn_2)) & fix_common_perc))
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  if(fix_common_perc){
    common_perc <- rep(0.5, rank_c)
  } else {
    common_perc <- sapply(1:rank_c, function(k){
      if(sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2)) < tol){
        .5
      } else {
        .latent_common_perc(score_1[,k], score_2[,k], nn_1, nn_2)
      }
    })
  }
   
  common_score <- sapply(1:rank_c, function(k){
    vec1 <- basis_list[[k]]$rep1; vec2 <-  basis_list[[k]]$rep2
    if(sum(abs(vec1 - vec2)) < tol){
      common_rep <- c(vec1 + vec2)/2
    } else {
      circle <- .construct_circle(vec1, vec2)
      rad1 <- .find_radian(circle, vec1); rad2 <- .find_radian(circle, vec2)
      stopifnot(rad1 < 0 & rad2 > 0) # must be true based on how we constructed vec1 and vec2
      
      rad1 <- rad1 + 2*pi
      stopifnot(sum(abs(.position_from_circle(circle, rad1) - vec1)) <= 1e-6,
                sum(abs(.position_from_circle(circle, rad2) - vec2)) <= 1e-6)
      
      common_rad <- .binary_search_radian(circle, rad2, rad1, common_perc[k])
      common_rep <- .position_from_circle(circle, common_rad)
    }
    
    basis_list[[k]]$basis_mat %*% common_rep
  })
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, common_perc = common_perc)
}

#####################

# thresholded to be between 0.01 and 0.99
.latent_common_perc <- function(score_vec_1, score_vec_2, nn_1, nn_2, tol = 1e-6){
  stopifnot(length(score_vec_1) == length(score_vec_2))
  
  n <- length(score_vec_1)
  mode1_common_perc <- sapply(1:n, function(i){
    val_1in1 <- stats::sd(score_vec_1[i] - score_vec_1[nn_1[i,]])
    val_1in2 <- stats::sd(score_vec_1[i] - score_vec_1[nn_2[i,]])
    val_2in1 <- stats::sd(score_vec_2[i] - score_vec_2[nn_1[i,]])
    val_2in2 <- stats::sd(score_vec_2[i] - score_vec_2[nn_2[i,]])
    
    ratio1 <- ifelse(abs(val_1in1) <= tol, 0, max(val_1in2/val_1in1, 1))
    ratio2 <- ifelse(abs(val_2in2) <= tol, 0, max(val_2in1/val_2in2, 1))
    
    .sigmoid_ratio(ratio1, ratio2)
  })
  
  min(max(1-mean(mode1_common_perc), 0.01), 0.99)
}

.sigmoid_ratio <- function(a, b, const = 5){
  max_val <- max(a,b); min_val <- min(a,b)
  tmp <- (max_val - min_val)/min_val
  if(b > a) tmp <- -tmp
  
  min(max(.5*tmp+.5, 0), 1)
}

# by construction, upper_radian is always the vector on the right. common_perc
# refers to the ratio of (distinct length to left vector)/(distinct length to right vector),
# meaning the lower common_perc is, the more the resulting vector woud lean left
# (i.e., further way from the right vector)
.binary_search_radian <- function(circle, lower_radian, upper_radian, common_perc, max_iter = 10, tol = 1e-6){
  stopifnot(lower_radian < upper_radian, 0 < common_perc, common_perc < 1)
  
  lower <- lower_radian; upper <- upper_radian
  left_vec <- .position_from_circle(circle, lower_radian)
  right_vec <- .position_from_circle(circle, upper_radian)
  iter <- 1; prev_mid <- NA
  
  while(iter < max_iter){
    mid <- (lower+upper)/2
    common_vec <- .position_from_circle(circle, mid)
    left_distinct <- left_vec - common_vec
    right_distinct <- right_vec - common_vec
    ratio <- .l2norm(left_distinct)/(.l2norm(left_distinct)+.l2norm(right_distinct))
    if(abs(ratio - common_perc) <= tol) break()
    
    if(ratio > common_perc){
      # need to make left_distinct smaller, so move radian right (i.e., smaller radian)
      upper <- mid
    } else {
      lower <- mid
    }
    if(!is.na(prev_mid) && abs(mid - prev_mid) < tol) break()
    
    prev_mid <- mid
    iter <- iter+1
  }
  
  mid
}
