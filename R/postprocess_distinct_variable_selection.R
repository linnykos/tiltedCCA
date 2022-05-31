postprocess_distinct_variable_selection <- function(
  input_obj,
  input_mat = NULL,
  logpval_vec, #larger is more significant
  cor_threshold = 0.9,
  input_assay = 2,
  max_variables = 10,
  min_subsample_cell = NULL,
  seurat_celltype_variable = "celltype",
  seurat_obj = NULL,
  verbose = 1
){
  stopifnot(input_assay %in% c(1,2))
  
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
  
  if(verbose > 0) print("Extracting relevant matrices")
  input_obj <- .set_defaultAssay(input_obj, assay = -input_assay+3)
  common_dimred_string <- ifelse(-input_assay+3 == 1, "common_dimred_1", "common_dimred_2")
  if(!common_dimred_string %in% names(input_obj)){
    input_obj <- tiltedCCA_decomposition(input_obj, 
                                         bool_modality_1_full = F,
                                         bool_modality_2_full = F,
                                         verbose = verbose)
  } 
  reference_dimred <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_dimred")
  if(verbose > 1) print(paste0("reference_dimred of dimension ", nrow(reference_dimred), " by ", ncol(reference_dimred)))
  
  input_obj <- .set_defaultAssay(input_obj, assay = input_assay)
  common_mat <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_mat")
  distinct_mat <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_mat")
  everything_mat <- common_mat + distinct_mat
  if(!all(is.null(input_mat))){
    stopifnot(all(dim(input_mat) == dim(everything_mat)), all(sort(colnames(input_mat)) == sort(colnames(input_mat))))
    if(inherits(input_mat, "dgCMatrix")) {
      everything_mat <- as.matrix(input_mat)
    } else if (inherits(input_mat, "matrix")){
      everything_mat <- input_mat
    } else {
      stop("input_mat is not a valid input")
    }
  }
  logpval_vec <- logpval_vec[colnames(everything_mat)]
  stopifnot(all(colnames(everything_mat) == names(logpval_vec)))
  cor_vec_intial <- NULL
  
  if(!is.null(min_subsample_cell)){
    if(verbose > 1) print("Reducing the number of cells")
    stopifnot(seurat_celltype_variable %in% colnames(seurat_obj@meta.data))
    membership_vec <- seurat_obj@meta.data[,seurat_celltype_variable]
    idx <- construct_celltype_subsample(membership_vec, min_subsample_cell = min_subsample_cell)
    
    ncell_before <- nrow(reference_dimred)
    reference_dimred <- reference_dimred[idx,]
    everything_mat <- everything_mat[idx,]
    ncell_after <- nrow(reference_dimred)
    if(verbose > 1) print(paste0("Reduced the number of cells from ", ncell_before, " to ", ncell_after))
  }
  
  n <- nrow(everything_mat)
  candidate_list <- vector("list", length = max_variables)
  selected_variables <- numeric(0)
  
  while(length(selected_variables) < max_variables | ncol(everything_mat) > 0){
    if(verbose > 0) print(paste0("On iteration: ", length(selected_variables)+1))
    
    if(verbose > 0) print("Selecting variable")
    cor_vec <- sapply(1:ncol(everything_mat), function(j){
      if(verbose > 1 && ncol(everything_mat) > 10 & j %% floor(ncol(everything_mat)/10) == 0) cat('*')
      if(verbose > 2) print(paste0("Working on variable ", j, " of ", ncol(everything_mat)))
      .linear_regression(bool_include_intercept = T,
                         bool_center_x = T,
                         bool_center_y = T,
                         bool_scale_x = T,
                         bool_scale_y = T,
                         return_type = "r_squared", 
                         x_mat = reference_dimred,
                         y_vec = everything_mat[,j])
    })
    names(cor_vec) <- colnames(everything_mat)
    if(all(is.null(cor_vec_intial))) cor_vec_intial <- cor_vec
    
    candidate_var <- colnames(everything_mat)[which(cor_vec <= cor_threshold)]
    candidate_list[[length(selected_variables)+1]] <- candidate_var
    
    if(verbose > 0) print(paste0(length(candidate_var), " eligble variables"))
    if(length(candidate_var) == 0) break()
    
    idx <- candidate_var[which.max(logpval_vec[candidate_var])]
    if(verbose > 1) {
      tmp_k <- min(5, length(candidate_var))
      print(paste0("Top ", tmp_k, " variables:"))
      print(logpval_vec[candidate_var[order(logpval_vec[candidate_var], decreasing = T)[1:tmp_k]]])
      print(paste0("Selected the variable: ", idx))
    }
    
    selected_variables <- c(selected_variables, idx)
    reference_dimred <- cbind(reference_dimred, everything_mat[,idx])
    colnames(reference_dimred)[ncol(reference_dimred)] <- idx
    
    stopifnot(all(selected_variables %in% colnames(reference_dimred)))
    everything_mat <- everything_mat[,which(!colnames(everything_mat) %in% selected_variables),drop = F]
    if(ncol(everything_mat) == 0) break()
  }
  
  structure(list(candidate_list = candidate_list,
                 cor_threshold = cor_threshold,
                 cor_vec_intial = cor_vec_intial,
                 logpval_vec = logpval_vec,
                 selected_variables = selected_variables), 
            class = "varSelect")
}