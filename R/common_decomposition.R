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
      ang <- .angle_between_vectors(vec1, vec2)
      vec_right <- .rightmost_vector(vec1, vec2)
      if(sum(abs(vec_right$vec_right - vec1)) <= sum(abs(vec_right$vec_right - vec1))){
        common_vec <- .angle_from_vector(vec_right$vec_right, (1-common_perc[k])*ang)
      } else {
        common_vec <- .angle_from_vector(vec_right$vec_right, common_perc[k]*ang)
      }
      
      a <- .l2norm(common_vec)^2; b <- -as.numeric(common_vec %*% (vec1 + vec2)); c <- as.numeric(vec1 %*% vec2)
      common_rep <- .quadratic(a,b,c) * common_vec
    }
    
    basis_list[[k]]$basis_mat %*% common_rep
  })
  
  if(length(rownames(score_1)) != 0) rownames(common_score) <- rownames(score_1)
  
  list(common_score = common_score, common_perc = common_perc)
}

#####################

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
  
  1-mean(mode1_common_perc)
}

.sigmoid_ratio <- function(a, b, const = 5){
  max_val <- max(a,b); min_val <- min(a,b)
  tmp <- (max_val - min_val)/min_val
  if(b > a) tmp <- -tmp
  
  min(max(.5*tmp+.5, 0), 1)
}

.quadratic <- function(a, b, c){
  tmp <- (-b +c(-1,1) * sqrt(b^2 - 4*a*c))/(2*a)
  tmp <- tmp[tmp > 0]
  stopifnot(length(tmp) > 0)
  min(tmp)
}

############################
