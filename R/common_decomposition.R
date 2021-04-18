.common_decomposition <- function(score_1, score_2, nn_1, nn_2, 
                                  fix_distinct_perc, 
                                  tol = 1e-6, verbose = F){
  stopifnot((!any(is.na(nn_1)) & !any(is.na(nn_2)) & !fix_distinct_perc) | 
              (all(is.na(nn_1)) & all(is.na(nn_2)) & fix_distinct_perc))
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  
  if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing distinct percentage"))
  if(fix_distinct_perc){
    distinct_perc_2 <- rep(0.5, rank_c)
  } else {
    distinct_perc_2 <- sapply(1:rank_c, function(k){
      if(verbose){
        if(rank_c < 10) { cat('*') } else if(k %% floor(rank_c/10) == 0) cat('*')
      } 
      if(sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2)) < tol){
        .5
      } else {
        .latent_distinct_perc_2(score_1[,k], score_2[,k], nn_1, nn_2, verbose = verbose)
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
      common_rad <- .binary_search_radian(circle, rad2, rad1, distinct_perc_2[k])
      common_rep <- .position_from_circle(circle, common_rad)
    }
    
    basis_list[[k]]$basis_mat %*% common_rep
  })
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, distinct_perc_2 = distinct_perc_2)
}

#####################

# thresholded to be between 0.01 and 0.99
.latent_distinct_perc_2 <- function(score_vec_1, score_vec_2, nn_1, nn_2,
                                    tol = 1e-6, verbose = F){
  stopifnot(length(score_vec_1) == length(score_vec_2), nrow(nn_1) == nrow(nn_2))
  
  n <- nrow(nn_1)
  distinct_perc_2 <- sapply(1:n, function(i){
    if(verbose){
      if(n < 10) { cat('*') } else if(i %% floor(n/10) == 0) cat('*')
    } 
    
    val_1in1 <- stats::sd(score_vec_1[nn_1[i,]])
    val_1in2 <- stats::sd(score_vec_1[nn_2[i,]])
    val_2in1 <- stats::sd(score_vec_2[nn_1[i,]])
    val_2in2 <- stats::sd(score_vec_2[nn_2[i,]])
    
    ratio1 <- ifelse(abs(val_1in1) <= tol, 0, max(val_1in2/val_1in1, 1))
    ratio2 <- ifelse(abs(val_2in2) <= tol, 0, max(val_2in1/val_2in2, 1))
    
    .sigmoid_ratio(ratio1, ratio2)
  })
  
  min(max(1-mean(distinct_perc_2), 0.01), 0.99)
}

.sigmoid_ratio <- function(a, b, tol = 1e-6){
  stopifnot(a >= 0, b >= 0)
  if(a <= tol & b <= tol) return(0.5)
  if(a <= tol & b >= tol) return(0)
  if(a >= tol & b <= tol) return(1)
  
  #max_val <- max(a,b); min_val <- min(a,b)
  #tmp <- (max_val - min_val)/min_val
  #if(b > a) tmp <- -tmp
  tmp <- (a-b)/min(a,b)
  
  min(max(.5*tmp+.5, 0), 1)
}

#' Find the appropriate radian for the common vector
#' 
#' Usage note: By \code{.representation_2d} works, the vector for
#' Modality 1 is always encoded by \code{right_radian}, and the
#' vector for Modality 2 is always encoded by \code{left_radian}.
#' 
#' If \code{distinct_perc_2} is close to 1, this means all the distinct information
#' (relatively speaking) lies in Modality 2. This means the returned
#' value should be very close to \code{right_radian}, meaning the distinct
#' vector for Modality 1 is near-0 (since the common vector is approximately
#' equal to the vector of Modality 1).
#'
#' @param circle object from \code{.construct_circle}
#' @param left_radian radian (between \code{-pi} and \code{pi}) of the left-most vector
#' @param right_radian radian (between \code{-pi} and \code{pi}) of the right-most vector
#' @param distinct_perc_2 the output of \code{.latent_distinct_perc_2}
#' @param max_iter positive integer
#' @param tol numeric
#'
#' @return radian (value between \code{right_radian} and \code{upper_radian})
.binary_search_radian <- function(circle, left_radian, right_radian, 
                                  distinct_perc_2, max_iter = 10, tol = 1e-6){
  stopifnot(right_radian > left_radian) # by design, we want always find a radian in the bottom-left of the circle
  stopifnot(0 < distinct_perc_2, distinct_perc_2 < 1)
  
  upper <- right_radian; lower <- left_radian
  right_vec <- .position_from_circle(circle, right_radian)
  left_vec <- .position_from_circle(circle, left_radian)
  iter <- 1; prev_mid <- NA
  
  while(iter < max_iter){
    mid <- (upper+lower)/2
    common_vec <- .position_from_circle(circle, mid)
    left_distinct <- left_vec - common_vec
    right_distinct <- right_vec - common_vec
    ratio <- .l2norm(left_distinct)/(.l2norm(left_distinct)+.l2norm(right_distinct))
    if(abs(ratio - distinct_perc_2) <= tol) break()
    
    if(ratio > distinct_perc_2){
      # the left distinct vector is too large, so to make it smaller,
      ## move the midpoint towards the larger (i.e., left) by the making
      ## the upper larger
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
