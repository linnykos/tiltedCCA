postprocess_variable_selection <- function(input_obj,
                                           logpval_vec, #larger is more significant
                                           cor_threshold = 0.9,
                                           max_variables = 10,
                                           input_assay = 2,
                                           verbose = 1){
  if(input_assay == 1) {
    stopifnot("distinct_mat_2" %in% names(input_obj),
              any(c("common_mat_1", "common_dimred_1") %in% names(input_obj)),
              max_variables <= nrow(input_obj$svd_2$v))
  } else if (input_assay == 2){
    stopifnot("distinct_mat_1" %in% names(input_obj),
              any(c("common_mat_2", "common_dimred_2") %in% names(input_obj)),
              max_variables <= nrow(input_obj$svd_1$v))
  } else {
    stop("assay not equal to 1 or 2")
  }
  
  stopifnot(length(names(logpval_vec)) == length(logpval_vec),
            cor_threshold >= 0, cor_threshold <= 1)
  
  input_obj <- .set_defaultAssay(input_obj, assay = -input_assay+3)
  common_mat_string <- ifelse(input_assay == 1, "common_mat_1", "common_mat_2")
  common_dimred_string <- ifelse(input_assay == 1, "common_dimred_1", "common_dimred_2")
  
  if(verbose > 0) print("Extracting relevant matrices")
  if(common_mat_string %in% names(input_obj)){
    reference_dimred <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_mat")
  } else if(common_dimred_string %in% names(input_obj)){
    reference_dimred <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_dimred")
  } else {
    stop(paste0("Cannot find the appropriate common matrix for modality ", input_assay))
  }
  
  input_obj <- .set_defaultAssay(input_obj, assay = input_assay)
  distinct_mat <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_mat")
  stopifnot(all(sort(colnames(distinct_mat)) == sort(names(logpval_vec))))
  n <- nrow(distinct_mat)
  candidate_list <- vector("list", length = max_variables)
  selected_variables <- numeric(0)
  
  while(length(selected_variables) < max_variables){
    if(verbose) print(paste0("On iteration: ", length(selected_variables)+1))
    
    if(verbose) print("Selecting variable")
    cor_vec <- sapply(1:ncol(distinct_mat), function(j){
      .linear_regression(bool_include_intercept = T,
                         bool_center_x = T,
                         bool_center_y = T,
                         bool_scale_x = T,
                         bool_scale_y = T,
                         return_type = "r_squared", 
                         x_mat = reference_dimred,
                         y_vec = distinct_mat[,j])
    })
    candidate_var <- colnames(distinct_mat)[which(cor_vec <= cor_threshold)]
    candidate_list[[length(selected_variables)+1]] <- candidate_var
    
    if(verbose > 0) print(paste0(length(candidate_var), " eligble variables"))
    if(length(candidate_var) == 0) break()
    
    idx <- candidate_var[which.max(logpval_vec[candidate_var])]
    selected_variables <- c(selected_variables, idx)
    reference_dimred <- cbind(reference_dimred, distinct_mat[,idx])
    distinct_mat <- distinct_mat[,which(!colnames(distinct_mat) %in% idx),drop = F]
    if(ncol(distinct_mat) == 0) break()
  }
  
  structure(list(selected_variables = selected_variables,
                 candidate_list = candidate_list), 
            class = "varSelect")
}