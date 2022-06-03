.generic_variable_selection <- function(bool_maximizing, 
                                        cor_threshold,
                                        initial_mat, # can be NULL
                                        mat,
                                        num_variables,
                                        return_candidate_list,
                                        vec,
                                        verbose = 0){
  stopifnot(is.numeric(cor_threshold), length(cor_threshold) == 1,
            cor_threshold >= 0, cor_threshold <= 1,
            is.matrix(mat), is.numeric(vec), is.logical(bool_maximizing),
            ncol(mat) == length(vec), length(vec) > 1,
            all(colnames(mat) == names(vec)),
            length(colnames(mat)) > 0, length(names(vec)) > 0)
  
  if(num_variables >= ncol(mat)) return(colnames(mat))
  if(!bool_maximizing) vec <- -vec
  
  if(verbose > 0) print("Selecting variables")
  cor_vec_intial <- numeric(0)
  selected_variables <- numeric(0)
  candidate_variables <- colnames(mat)
  candidate_list <- vector("list", length = num_variables)
  if(!all(is.null(initial_mat))){
    selected_mat <- initial_mat
  } else {
    selected_mat <- numeric(0)
  }

  while(length(selected_variables) < num_variables & length(candidate_variables) > 0){
    
    if(length(selected_mat) > 1){
      # compute correlations
      cor_vec <- sapply(candidate_variables, function(variable){
        .linear_regression(bool_include_intercept = T,
                           bool_center_x = T,
                           bool_center_y = T,
                           bool_scale_x = T,
                           bool_scale_y = T,
                           return_type = "r_squared", 
                           x_mat = selected_mat,
                           y_vec = mat[,which(colnames(mat) == variable)])
      })
      names(cor_vec) <- candidate_variables
      if(length(cor_vec_intial) == 0) cor_vec_intial <- cor_vec
      
      # screen out variables

      candidate_variables <- names(cor_vec)[which(cor_vec <= cor_threshold)]
      vec <- vec[candidate_variables]
    }
    
    if(length(candidate_variables) == 0) break()
    if(verbose > 0) print(paste0("On iteration ", length(selected_variables)+1, 
                          ": ", length(candidate_variables), " remaining candidates"))
  
    candidate_list[[length(selected_variables)+1]] <- candidate_variables
    
    # select variable and bookkeeping
    if(length(candidate_variables) == 1){
      new_variable <- candidate_variables[1]
      selected_variables <- c(selected_variables, new_variable)
      break()
    } else {
      new_variable <- names(vec)[which.max(vec)]
      
      if(verbose > 1) {
        tmp_k <- min(5, length(candidate_variables))
        print(paste0("Top ", tmp_k, " variables:"))
        print(vec[candidate_variables[order(vec[candidate_variables], decreasing = T)[1:tmp_k]]])
        print(paste0("Selected the variable: ", new_variable))
      }
      
      selected_mat <- cbind(selected_mat, mat[,new_variable,drop = F])
      selected_variables <- c(selected_variables, new_variable)
    }
    candidate_variables <- candidate_variables[!candidate_variables %in% selected_variables]
    
  }
  
  if(return_candidate_list){
    list(candidate_list = candidate_list,
         cor_vec = cor_vec_intial,
         selected_variables = selected_variables)
  } else {
    selected_variables
  }
  
}