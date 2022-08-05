postprocess_distinct_variable_selection <- function(
    input_obj,
    input_mat = NULL,
    logpval_vec, #larger is more significant
    cor_threshold = 0.9,
    input_assay = 2,
    min_subsample_cell = NULL,
    num_variables = 10,
    seurat_celltype_variable = "celltype",
    seurat_obj = NULL,
    verbose = 1
){
  stopifnot(input_assay %in% c(1,2))
  
  if(input_assay == 1) {
    stopifnot("distinct_mat_2" %in% names(input_obj),
              any(c("common_mat_1", "common_dimred_1") %in% names(input_obj)),
              num_variables <= nrow(input_obj$svd_2$v))
  } else if (input_assay == 2){
    stopifnot("distinct_mat_1" %in% names(input_obj),
              any(c("common_mat_2", "common_dimred_2") %in% names(input_obj)),
              num_variables <= nrow(input_obj$svd_1$v))
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
  
  if(!all(is.null(input_mat))){
    if(inherits(input_mat, "dgCMatrix")) {
      everything_mat <- as.matrix(input_mat)
    } else if (inherits(input_mat, "matrix")){
      everything_mat <- input_mat
    } else {
      stop("input_mat is not a valid input")
    }
  } else {
    input_obj <- .set_defaultAssay(input_obj, assay = input_assay)
    common_mat <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_mat")
    distinct_mat <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_mat")
    everything_mat <- common_mat + distinct_mat
  }
  stopifnot(nrow(reference_dimred) == nrow(everything_mat), length(colnames(everything_mat)) > 0)
  
  logpval_vec <- logpval_vec[colnames(everything_mat)]
  stopifnot(all(colnames(everything_mat) == names(logpval_vec)))

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
  
  res <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = cor_threshold,
    initial_mat = reference_dimred,
    mat = everything_mat,
    num_variables = num_variables,
    return_candidate_list = T,
    vec = logpval_vec,
    verbose = verbose
  )
  
  structure(list(candidate_list = res$candidate_list,
                 cor_threshold = cor_threshold,
                 cor_vec_intial = res$cor_vec,
                 logpval_vec = logpval_vec,
                 selected_variables = res$selected_variables), 
            class = "varSelect")
}