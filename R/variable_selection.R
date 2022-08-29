.generic_variable_selection <- function(bool_maximizing, 
                                        cor_threshold,
                                        initial_mat, # can be NULL
                                        mat_1,
                                        mat_2, # can be NULL
                                        num_variables,
                                        prediction_type = c("lm", "cor")[1],
                                        return_candidate_list,
                                        vec_1,
                                        vec_2, # can be NULL
                                        verbose = 0){
  stopifnot(is.numeric(cor_threshold), length(cor_threshold) == 1,
            cor_threshold >= 0, cor_threshold <= 1,
            is.matrix(mat_1), is.numeric(vec_1), is.logical(bool_maximizing),
            ncol(mat_1) == length(vec_1), length(vec_1) > 1,
            all(colnames(mat_1) == names(vec_1)),
            length(colnames(mat_1)) > 0, length(names(vec_1)) > 0,
            prediction_type %in% c("lm", "cor"))
  if(!all(is.null(mat_2))){
    bool_second <- T
    stopifnot(!all(is.null(vec_2)), all(dim(mat_1) == dim(mat_2)),
              length(vec_1) == length(vec_2),
              all(colnames(mat_1) == colnames(mat_2)),
              all(names(vec_1) == names(vec_2)),
              length(colnames(mat_2)) > 0, length(names(vec_2)) > 0)
  } else{
    bool_second <- F
  }
  
  if(num_variables >= ncol(mat_1)) return(colnames(mat_1))
  if(!bool_maximizing) {
    vec_1 <- -vec_1
    if(bool_second) vec_2 <- -vec_2
  }
  
  if(verbose > 0) print("Gathering ingredients")
  cor_vec_intial_1 <- numeric(0)
  cor_vec_intial_2 <- numeric(0)
  selected_variables <- numeric(0)
  candidate_variables <- colnames(mat_1)
  candidate_list <- vector("list", length = num_variables)
  if(prediction_type == "lm"){
    if(!all(is.null(initial_mat))){
      selected_mat_1 <- initial_mat; selected_mat_2 <- initial_mat
    } else {
      selected_mat_1 <- numeric(0); selected_mat_2 <- numeric(0)
    }
    
    cor_mat_1 <- NULL
    cor_mat_2 <- NULL
    
  } else {
    cor_mat_1 <- stats::cor(mat_1, method = "pearson")
    cor_mat_2 <- stats::cor(mat_2, method = "pearson")
    
    mat_1 <- NULL
    mat_2 <- NULL
    selected_mat_1 <- NULL
    selected_mat_2 <- NULL
  }
 
  
  if(verbose > 0) print("Selecting variables")
  while(length(selected_variables) < num_variables & length(candidate_variables) > 0){
    
    if(!all(is.null(initial_mat)) | length(selected_variables) >= 1){
      # compute correlations
      tmp <- .var_select_compute_predictions(candidate_variables = candidate_variables,
                                             cor_mat = cor_mat_1,
                                             cor_threshold = cor_threshold,
                                             mat = mat_1,
                                             selected_mat = selected_mat_1,
                                             type = prediction_type,
                                             vec = vec_1)
      cor_vec_1 <- tmp$cor_vec; candidate_variables_1 <- tmp$candidate_variables
      if(length(cor_vec_intial_1) == 0) cor_vec_intial_1 <- cor_vec_1
      
      if(bool_second){
        tmp <- .var_select_compute_predictions(candidate_variables = candidate_variables,
                                               cor_mat = cor_mat_2,
                                               cor_threshold = cor_threshold,
                                               mat = mat_2,
                                               selected_mat = selected_mat_2,
                                               type = prediction_type,
                                               vec = vec_2)
        cor_vec_2 <- tmp$cor_vec; candidate_variables_2 <- tmp$candidate_variables
        if(length(cor_vec_intial_2) == 0) cor_vec_intial_2 <- cor_vec_2
        stopifnot(length(cor_vec_1) == length(cor_vec_2), length(selected_mat_1) == length(selected_mat_2))
        
        candidate_variables <- sort(intersect(candidate_variables_1, candidate_variables_2))
        
      } else {
        candidate_variables <- sort(candidate_variables_1)
      }
    }
    
    ####
    
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
      if(bool_second){
        rank_1 <- rank(vec_1[candidate_variables]); rank_2 <- rank(vec_2[candidate_variables])
        rank_overall <- sapply(1:length(candidate_variables), function(i){
          min(rank_1[i], rank_2[i])
        })
        names(rank_overall) <- candidate_variables
        new_variable <- names(rank_overall)[which.max(rank_overall)]
      } else {
        new_variable <- candidate_variables[which.max(vec_1[candidate_variables])]
      }
      
      if(verbose > 1) {
        tmp_k <- min(5, length(candidate_variables))
        print(paste0("Top ", tmp_k, " variables in one modality:"))
        print(vec_1[candidate_variables[order(vec_1[candidate_variables], decreasing = T)[1:tmp_k]]])
        if(bool_second){
          print(paste0("Top ", tmp_k, " variables in the other modality:"))
          print(vec_2[candidate_variables[order(vec_2[candidate_variables], decreasing = T)[1:tmp_k]]])
        }
        print(paste0("Selected the variable: ", new_variable))
      }
      
      if(prediction_type == "lm"){
        selected_mat_1 <- cbind(selected_mat_1, mat_1[,new_variable,drop = F])
        if(bool_second) selected_mat_2 <- cbind(selected_mat_2, mat_2[,new_variable,drop = F])
      }
      
      selected_variables <- c(selected_variables, new_variable)
    }
    
    candidate_variables <- candidate_variables[!candidate_variables %in% selected_variables]
  }
  
  if(return_candidate_list){
    list(candidate_list = candidate_list,
         cor_vec_1 = cor_vec_intial_1,
         cor_vec_2 = cor_vec_intial_2,
         selected_variables = selected_variables)
  } else {
    selected_variables
  }
  
}

#############

.var_select_compute_predictions <- function(candidate_variables,
                                            cor_mat,
                                            cor_threshold,
                                            mat,
                                            selected_mat,
                                            type,
                                            vec){
  if(type == "lm"){
    stopifnot(all(is.null(cor_mat)), !all(is.null(mat)), !all(is.null(selected_mat)))
    .var_select_compute_lm(candidate_variables = candidate_variables,
                           cor_threshold = cor_threshold,
                           mat = mat,
                           selected_mat = selected_mat,
                           vec = vec)
  } else if(type == "cor"){
    stopifnot(!all(is.null(cor_mat)), all(is.null(mat)), all(is.null(selected_mat)))
    .var_select_compute_correlations(candidate_variables = candidate_variables,
                                     cor_mat = cor_mat,
                                     cor_threshold = cor_threshold,
                                     vec = vec)
  } else stop("type in .var_select_compute_predictions not found")
}

.var_select_compute_lm <- function(candidate_variables,
                                   cor_threshold,
                                   mat,
                                   selected_mat,
                                   vec){
  stopifnot(all(sort(colnames(mat)) == sort(names(vec))),
            all(candidate_variables %in% names(vec)))
  
  cor_vec <- sapply(candidate_variables, function(variable){
    .linear_regression(bool_include_intercept = T,
                       bool_center_x = T,
                       bool_center_y = T,
                       bool_scale_x = T,
                       bool_scale_y = T,
                       return_type = "r_squared", 
                       x_mat = selected_mat,
                       y_vec = mat[,variable])
  })
  names(cor_vec) <- candidate_variables
  
  # screen out variables
  candidate_variables <- names(cor_vec)[which(cor_vec <= cor_threshold)]
  
  list(candidate_variables = candidate_variables,
       cor_vec = cor_vec)
}

.var_select_compute_correlations <- function(candidate_variables,
                                             cor_mat,
                                             cor_threshold,
                                             vec){
  stopifnot(length(colnames(cor_mat)) == ncol(cor_mat), 
            all(colnames(cor_mat) == rownames(cor_mat)),
            nrow(cor_mat) == length(vec),
            all(candidate_variables %in% names(vec)))
  
  candidate_idx <- which(colnames(cor_mat) %in% candidate_variables)
  candidate_nidx <- which(!colnames(cor_mat) %in% candidate_variables)
  
  cor_vec <- sapply(candidate_idx, function(i){
    max(cor_mat[i,candidate_nidx])
  })
  names(cor_vec) <- colnames(cor_mat)[candidate_idx]
  
  # screen out variables
  candidate_variables <- names(cor_vec)[which(cor_vec <= cor_threshold)]
  
  list(candidate_variables = candidate_variables,
       cor_vec = cor_vec)
}