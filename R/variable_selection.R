variable_selection <- function(tilted_res,
                               significance_vec, #larger is more significant
                               max_variables = 10,
                               cor_threshold = 0.9,
                               verbose = 0){
  stopifnot(length(names(significance_vec)) == length(significance_vec))
  
  if(verbose) print("Extracting relevant matrices")
  reference_mat <- tilted_res$common_score
  decomp_res <- tiltedCCA_decomposition(tilted_res)
  if(length(significance_vec) == nrow(tilted_res$svd_1$v) && 
     all(sort(significance_vec) == sort(rownames(tilted_res$svd_1$v)))){
    distinct_mat <- tcrossprod(.mult_mat_vec(tilted_res$svd_1$u, tilted_res$svd_1$d), tilted_res$svd_1$v)
  } else {
    distinct_mat <- tcrossprod(.mult_mat_vec(tilted_res$svd_2$u, tilted_res$svd_2$d), tilted_res$svd_2$v)
  }
  significance_vec <- significance_vec[colnames(distinct_mat)]
  
  n <- nrow(distinct_mat)
  
  candidate_list <- vector("list", length = max_variables)
  selected_variables <- numeric(0)
  
  while(length(selected_variables) < max_variables){
    if(verbose > 0) print(paste0("On iteration: ", length(selected_variables)+1))
    
    if(verbose > 0) print("Selecting variable")
    cor_vec <- sapply(1:ncol(distinct_mat), function(i){
      .variable_correlation(x_mat = reference_mat, 
                            y_vec = distinct_mat[,i])
    })
    candidate_var <- colnames(distinct_mat)[which(cor_vec <= cor_threshold)]
    candidate_list[[length(selected_variables)+1]] <- candidate_var
    
    if(verbose > 0) print(paste0(length(candidate_var), " eligble variables"))
    if(length(candidate_var) == 0) break()
    
    idx <- candidate_var[which.max(significance_vec[candidate_var])]
    if(verbose > 1) {
      tmp_k <- min(5, length(candidate_var))
      print(paste0("Top ", tmp_k, " variables:"))
      print(candidate_var[order(significance_vec[candidate_var], decreasing = T)[1:tmp_k]])
      print(paste0("Selected the variable: ", candidate_var[idx]))
    }
    selected_variables <- c(selected_variables, idx)
    reference_mat <- cbind(reference_mat, distinct_mat[,idx])
    distinct_mat <- distinct_mat[,which(!colnames(distinct_mat) %in% selected_variables),drop = F]
    if(ncol(distinct_mat) == 0) break()
  }
  
  list(selected_variables = selected_variables,
       candidate_list = candidate_list)
}

.variable_correlation <- function(x_mat, y_vec){
  df <- data.frame(x_mat)
  df <- cbind(df, y_vec)
  colnames(df)[ncol(df)] <- "y"
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
}