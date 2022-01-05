fine_tuning <- function(dcca_res, 
                        max_iter = 5,
                        fix_tilt_perc = NA,
                        temp_path = NULL,
                        tol = 1e-5,
                        verbose = T){
  score_1 <- dcca_res$score_1
  score_2 <- dcca_res$score_2
  
  rank_c <- min(ncol(score_1), ncol(score_2))
  stopifnot(all(is.na(fix_tilt_perc)) || length(fix_tilt_perc) == rank_c)
  
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
  
  # [[TODO: Format metacell_clustering_1 and metacell_clustering_2 from factor to list]]
  iter <- 1
  percentage_grid <- seq(0, 1, length.out = dcca_res$param_list$discretization_gridsize)
  common_score <- dcca_res$common_score
  common_score_prev <- common_score
  n <- nrow(dcca_res$common_score)
  if(dcca_res$param_list$cell_max < n){
    n_idx <- sample(1:n, size = dcca_res$param_list$cell_max)
  } else {
    n_idx <- 1:n
  }
  if(all(is.na(fix_tilt_perc))) {
    tilt_perc <- rep(dcca_res$tilt_perc[1], rank_c)
  } else {
    tilt_perc <- fix_tilt_perc
  }
  
  while(iter <= max_iter){
    for(k in 1:rank_c){
      if(all(is.na(fix_tilt_perc))){
        tmp <- percentage_grid
      } else {
        tmp <- tilt_perc[k]
      }
      
      res <- .fine_tuning_dimension(basis_list = basis_list,
                                    circle_list = circle_list,
                                    common_score = common_score,
                                    latent_dim = k,
                                    metacell_clustering_1 = dcca_res$metacell_clustering_1,
                                    metacell_clustering_2 = dcca_res$metacell_clustering_2,
                                    n_idx = n_idx,
                                    num_neigh = dcca_res$param_list$num_neigh,
                                    percentage_grid = tmp,
                                    score_1 = score_1,
                                    score_2 = score_2,
                                    svd_1 = dcca_res$svd_1, 
                                    svd_2 = dcca_res$svd_2,
                                    target_subspace = dcca_res$target_subspace,
                                    verbose = verbose)
      if(verbose) {
        print(paste0("On iteration ", iter, " for latent dimension ", k))
        print(res$df)
        print("=======")
      }
      tilt_perc[k] <- res$percentage
      common_score <- res$common_score
    }
    
    if(iter > 2 && sum(abs(common_score - common_score_prev)) <= tol) break()
    if(!is.null(temp_path) && is.character(temp_path)){
      tmp <- .compute_distinct_score(score_1, score_2, common_score)
      distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
      res <- list(common_score = common_score, 
                  distinct_score_1 = distinct_score_1,
                  distinct_score_2 = distinct_score_2,
                  score_1 = score_1, score_2 = score_2, 
                  svd_1 = dcca_res$svd_1, svd_2 = dcca_res$svd_2, 
                  tilt_perc = tilt_perc)
      save(res, file = temp_path)
    }
    
    iter <- iter + 1
    common_score_prev <- common_score
  }
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  param_list <- dcca_res$param_list
  param_list[["fine_tuning_max_iter"]] <- max_iter
  
  res <- list(common_score = common_score, 
              distinct_score_1 = distinct_score_1,
              distinct_score_2 = distinct_score_2,
              score_1 = score_1, score_2 = score_2, 
              svd_1 = dcca_res$svd_1, svd_2 = dcca_res$svd_2, 
              cca_obj = dcca_res$cca_obj, 
              df_percentage = NA,
              metacell_clustering_1 = dcca_res$metacell_clustering_1,
              metacell_clustering_2 = dcca_res$metacell_clustering_2,
              min_mat = dcca_res$min_mat,
              target_subspace = dcca_res$target_subspace,
              param_list = dcca_res$param_list,
              tilt_perc = tilt_perc
  )
  
  class(res) <- "dcca"
  res
}

##############################

.fine_tuning_dimension <- function(basis_list,
                                   circle_list,
                                   common_score,
                                   latent_dim,
                                   metacell_clustering_1,
                                   metacell_clustering_2,
                                   n_idx,
                                   num_neigh,
                                   percentage_grid,
                                   score_1,
                                   score_2,
                                   svd_1, 
                                   svd_2,
                                   target_subspace,
                                   verbose = T){
  r <- length(basis_list)
  
  value_list <- lapply(percentage_grid, function(percentage){
    radian_val <- .compute_radian(circle = circle_list[[latent_dim]],
                                  enforce_boundary = T,
                                  percentage_val = percentage, 
                                  vec1 = basis_list[[latent_dim]]$rep1,
                                  vec2 = basis_list[[latent_dim]]$rep2)
    
    common_representation_new <- .position_from_circle(circle_list[[latent_dim]], radian_val)
    common_score[,latent_dim] <- basis_list[[latent_dim]]$basis_mat %*% common_representation_new
    
    common_mat <- .convert_common_score_to_mat(common_score,
                                               score_1,
                                               score_2,
                                               svd_1, 
                                               svd_2)
    
    value <- .determine_cluster(common_mat = common_mat, 
                                target_subspace = target_subspace)
    
    list(value = value, common_vec = common_score[,latent_dim])
  })
  names(value_list) <- percentage_grid
  value_vec <- sapply(value_list, function(x){x$value})
  
  idx_min <- .select_minimum(minimum = T,
                             x_val = percentage_grid,
                             y_val = value_vec)
  df <- data.frame(percentage = percentage_grid,
                   ratio_val = value_vec)
  prev_vec <- common_score[,latent_dim]
  common_score[,latent_dim] <- value_list[[idx_min]]$common_vec
  print(sum(abs(prev_vec - common_score[,latent_dim])))
  
  list(common_score = common_score,
       df = df, 
       percentage = percentage_grid[idx_min])
}