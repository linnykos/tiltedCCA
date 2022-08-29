postprocess_graph_alignment <- function(
    input_obj,
    bool_use_denoised,
    bool_include_intercept = T,
    bool_use_metacells = T,
    cell_idx = NULL,
    input_assay = 1,
    return_everything_mat = F,
    seurat_obj = NULL,
    seurat_assay = NULL,
    seurat_slot = "data",
    tol = 1e-6,
    variable_names = NULL,
    verbose = 0
){
  stopifnot(input_assay %in% c(1,2),
            seurat_slot %in% c("counts", "data", "scale.data"))
  
  if(verbose > 0) print("Gathering ingredients")
  
  if(bool_use_denoised){
    input_obj <- .set_defaultAssay(input_obj, 
                                   assay = input_assay)
    common_mat_string <- paste0("common_mat_", input_assay)
    distinct_mat_string <- paste0("distinct_mat_", input_assay)
    if(!common_mat_string %in% names(input_obj) || 
       !distinct_mat_string %in% names(input_obj)){
      input_obj <- tiltedCCA_decomposition(input_obj, 
                                           bool_modality_1_full = T,
                                           bool_modality_2_full = T,
                                           verbose = verbose)
    }
    common_mat <- .get_tCCAobj(input_obj, 
                               apply_postDimred = F,
                               what = "common_mat")
    distinct_mat <- .get_tCCAobj(input_obj, 
                                 apply_postDimred = F,
                                 what = "distinct_mat")
    everything_mat <- common_mat + distinct_mat
    
  } else {
    stopifnot(inherits(seurat_obj, "Seurat"))
    if(is.null(seurat_assay)) seurat_assay <- Seurat::DefaultAssay(seurat_obj)
    
    if(seurat_slot == "counts"){
      stopifnot(length(seurat_obj[[seurat_assay]]@var.features) > 0)
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@counts[seurat_obj[[seurat_assay]]@var.features,])
      everything_mat <- as.matrix(everything_mat)
    } else if(seurat_slot == "data"){
      stopifnot(length(seurat_obj[[seurat_assay]]@var.features) > 0)
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@data[seurat_obj[[seurat_assay]]@var.features,])
      everything_mat <- as.matrix(everything_mat)
    } else if(seurat_slot == "scale.data"){
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@scale.data)
    } else {
      stop("seurat_slot invalid")
    }
  }
  
  if(!all(is.null(variable_names))){
    everything_mat <- everything_mat[,sort(intersect(colnames(everything_mat), variable_names)),drop = F]
    stopifnot(ncol(everything_mat) >= 1)
  }
  
  if(verbose > 0) print("Computing NN graph")
  common_dimred_string1 <- "common_dimred_1"
  common_dimred_string2 <- "common_dimred_2"
  if(!common_dimred_string1 %in% names(input_obj) || 
     !common_dimred_string2 %in% names(input_obj)){
    input_obj <- tiltedCCA_decomposition(input_obj, 
                                         bool_modality_1_full = F,
                                         bool_modality_2_full = F,
                                         verbose = verbose)
  }
  
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  dimred_1 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  dimred_2 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
  dimred <- cbind(dimred_1, dimred_2)
  stopifnot(nrow(dimred) == nrow(everything_mat))
  
  # account of selection of cells
  if(all(!is.null(cell_idx))){
    stopifnot(is.numeric(cell_idx), all(cell_idx >= 1),
              all(cell_idx <= nrow(dimred)),
              all(cell_idx %% 1 == 0))
    dimred <- dimred[cell_idx,]
    everything_mat <- everything_mat[cell_idx,,drop = F]
  }
  
  param <- .get_param(input_obj)
  if(bool_use_metacells & !is.null(param$mc_num_metacells)){
    metacell_clustering_list <- .get_metacell(
      input_obj,
      resolution = "cell", 
      type = "list", 
      what = "metacell_clustering"
    )
    averaging_mat <- .generate_averaging_matrix(
      metacell_clustering_list = metacell_clustering_list,
      n = n
    )
    
    if(all(!is.null(cell_idx))){
      averaging_mat <- averaging_mat[,cell_idx]
    }
    
    everything_mat <- averaging_mat %*% everything_mat
    dimred <- averaging_mat %*% dimred
  }
  
  num_neigh <- param$snn_num_neigh
  snn_mat <- .form_snn_mat(
    mat = dimred, 
    num_neigh = param$snn_num_neigh,
    bool_cosine = param$snn_bool_cosine,
    bool_intersect = param$snn_bool_intersect,
    min_deg = param$snn_min_deg,
    verbose = verbose
  )
  
  if(verbose > 0) print("Computing Laplacian")
  laplacian_basis <- .compute_laplacian_basis(
    latent_k = param$snn_latent_k,
    sparse_mat = snn_mat,
    verbose = verbose
  )
  
  if(verbose > 0) print("Regressing variables on Laplacian basis")
  p <- ncol(everything_mat)
  alignment_vec <- sapply(1:p, function(j){
    .linear_regression(bool_include_intercept = bool_include_intercept,
                       bool_center_x = T,
                       bool_center_y = T,
                       bool_scale_x = T,
                       bool_scale_y = T,
                       return_type = "r_squared", 
                       x_mat = laplacian_basis,
                       y_vec = everything_mat[,j])
  })
  names(alignment_vec) <- colnames(everything_mat)
  
  if(return_everything_mat){
    list(alignment = alignment_vec,
         everything_mat = everything_mat)
  } else {
    alignment_vec
  }
}

######################

postprocess_smooth_variable_selection <- function(
    input_obj,
    bool_use_denoised,
    bool_include_intercept = T,
    bool_use_metacells = T,
    bool_use_both_modalities = T,
    cell_idx = NULL,
    cor_threshold = 0.8,
    num_variables = 50,
    sd_quantile = 0.75,
    seurat_obj = NULL,
    seurat_assay_1 = NULL,
    seurat_assay_2 = NULL,
    seurat_slot = "data",
    variable_names = NULL,
    verbose = 0
){
  # gather ingredients
  res <- .smooth_variable_selection_helper(
    bool_include_intercept = bool_include_intercept,
    bool_use_denoised = bool_use_denoised,
    bool_use_metacells = bool_use_metacells,
    cell_idx = cell_idx,
    input_assay = 1,
    input_obj = input_obj,
    sd_quantile = sd_quantile,
    seurat_assay = seurat_assay_1,
    seurat_obj = seurat_obj,
    seurat_slot = seurat_slot,
    variable_names = variable_names,
    verbose = verbose
  )
  alignment_vec_1 <- res$alignment_vec
  everything_mat_1 <- res$everything_mat
  sd_vec_1 <- res$sd_vec
  
  if(bool_use_both_modalities){
    res <- .smooth_variable_selection_helper(
      bool_include_intercept = bool_include_intercept,
      bool_use_denoised = bool_use_denoised,
      bool_use_metacells = bool_use_metacells,
      cell_idx = cell_idx,
      input_assay = 2,
      input_obj = input_obj,
      sd_quantile = sd_quantile,
      seurat_assay = seurat_assay_2,
      seurat_obj = seurat_obj,
      seurat_slot = seurat_slot,
      variable_names = variable_names,
      verbose = verbose
    )
    alignment_vec_2 <- res$alignment_vec
    everything_mat_2 <- res$everything_mat
    sd_vec_2 <- res$sd_vec
    stopifnot(length(names(alignment_vec_2)) > 0)
    
    shared_variables <- intersect(names(alignment_vec_1), names(alignment_vec_2))
    alignment_vec_1 <- alignment_vec_1[shared_variables]
    alignment_vec_2 <- alignment_vec_2[shared_variables]
    everything_mat_1 <- everything_mat_1[,shared_variables]
    everything_mat_2 <- everything_mat_2[,shared_variables]
    sd_vec_1 <- sd_vec_1[shared_variables]
    sd_vec_2 <- sd_vec_2[shared_variables]
    stopifnot(all(names(alignment_vec_2) == names(alignment_vec_1)))
  } else {
    alignment_vec_2 <- NULL
    everything_mat_2 <- NULL
    sd_vec_2 <- NULL
  }
  
  #########
  
  p <- length(alignment_vec_1)
  if(num_variables >= p) return(list(
    alignment_1 = alignment_vec_1,
    alignment_2 = alignment_vec_2,
    cor_threshold = cor_threshold,
    sd_quantile = sd_quantile,
    sd_vec_1 = sd_vec_1,
    sd_vec_2 = sd_vec_2,
    selected_variables = names(alignment_1)
  ))
  
  if(verbose > 0) print("Selecting variables")
  selected_variables <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = cor_threshold,
    initial_mat = NULL,
    mat_1 = everything_mat_1,
    mat_2 = everything_mat_2,
    num_variables = num_variables,
    prediction_type = "cor",
    return_candidate_list = F,
    vec_1 = alignment_vec_1,
    vec_2 = alignment_vec_2,
    verbose = verbose
  )
  
  list(alignment_1 = alignment_vec_1,
       alignment_2 = alignment_vec_2,
       cor_threshold = cor_threshold,
       sd_quantile = sd_quantile,
       sd_vec_1 = sd_vec_1,
       sd_vec_2 = sd_vec_2,
       selected_variables = selected_variables)
}

###########################################

.smooth_variable_selection_helper <- function(
    bool_include_intercept,
    bool_use_denoised,
    bool_use_metacells,
    cell_idx,
    input_assay,
    input_obj,
    sd_quantile,
    seurat_assay,
    seurat_obj,
    seurat_slot,
    variable_names,
    verbose
){
  res <- postprocess_graph_alignment(
    input_obj = input_obj,
    bool_use_denoised = bool_use_denoised,
    bool_include_intercept = bool_include_intercept,
    bool_use_metacells = bool_use_metacells,
    cell_idx = cell_idx,
    input_assay = input_assay,
    return_everything_mat = T,
    seurat_obj = seurat_obj,
    seurat_assay = seurat_assay,
    seurat_slot = seurat_slot,
    variable_names = variable_names,
    verbose = verbose
  )
  alignment_vec <- res$alignment
  everything_mat <- res$everything_mat
  if(any(is.na(alignment_vec))){
    everything_mat <- everything_mat[,which(!is.na(alignment_vec)), drop = F]
    alignment_vec <- alignment_vec[which(!is.na(alignment_vec))]
  }
  
  sd_vec <- matrixStats::colSds(everything_mat)
  names(sd_vec) <- colnames(everything_mat)
  if(!is.null(sd_quantile)){
    sd_threshold <- stats::quantile(sd_vec, probs = sd_quantile)
    everything_mat <- everything_mat[,which(sd_vec >= sd_threshold), drop = F]
    alignment_vec <- alignment_vec[which(sd_vec >= sd_threshold)]
    sd_vec <- sd_vec[which(sd_vec >= sd_threshold)]
  }
  stopifnot(all(names(alignment_vec) == names(everything_mat)))
  
  list(alignment_vec = alignment_vec,
       everything_mat = everything_mat,
       sd_vec = sd_vec)
  
}