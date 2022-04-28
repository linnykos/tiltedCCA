postprocess_alignment <- function(input_obj,
                                  bool_use_denoised,
                                  bool_center = T,
                                  bool_scale = T,
                                  bool_everything_center = T,
                                  bool_everything_scale = T,
                                  bool_regression_include_intercept = T,
                                  bool_regression_center = T,
                                  bool_regression_scale = T,
                                  input_assay = 1,
                                  min_subsample_cell = NULL,
                                  seurat_celltype_variable = "celltype",
                                  seurat_obj = NULL,
                                  seurat_assay = Seurat::DefaultAssay(seurat_obj),
                                  seurat_slot = "data",
                                  verbose = 1){
  stopifnot(input_assay %in% c(1,2),
            seurat_slot %in% c("counts", "data", "scale.data"))
  
  if(verbose > 0) print("Gathering ingredients")
  input_obj <- .set_defaultAssay(input_obj, 
                                 assay = input_assay)
  common_mat <- .get_tCCAobj(input_obj, 
                             apply_postDimred = F,
                             what = "common_mat")
  if(bool_use_denoised){
    distinct_mat <- .get_tCCAobj(input_obj, 
                                 apply_postDimred = F,
                                 what = "distinct_mat")
    everything_mat <- common_mat + distinct_mat
  } else {
    stopifnot(inherits(seurat_obj, "Seurat"))
    if(seurat_slot == "counts"){
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@counts[seurat_obj[[seurat_assay]]@var.features,])
    } else if(seurat_slot == "data"){
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@data)
    } else if(seurat_slot == "scale.data"){
      everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@scale.data)
    } else {
      stop("seurat_slot invalid")
    }
  }
  
  everything_mat <- everything_mat[,colnames(common_mat)]
  stopifnot(all(dim(common_mat) == dim(everything_mat)))
  
  if(!is.null(min_subsample_cell)){
    stopifnot(seurat_celltype_variable %in% colnames(seurat_obj@meta.data))
    membership_vec <- seurat_obj@meta.data[,seurat_celltype_variable]
    idx <- construct_celltype_subsample(membership_vec, min_subsample_cell = min_subsample_cell)
    
    common_mat <- common_mat[idx,]
    everything_mat <- everything_mat[idx,]
  }
  
  if(bool_center | bool_scale){
    common_mat <- scale(common_mat,
                        center = bool_center, 
                        scale = bool_scale)
    everything_mat <- scale(everything_mat,
                            center = bool_center, 
                            scale = bool_scale)
  }
  
  if(verbose > 0) print("Starting to perform regressions")
  
  rsquare_vec <- sapply(1:ncol(common_mat), function(j){
    if(verbose > 0 && ncol(common_mat) > 10 && j %% floor(ncol(common_mat)/10) == 0) cat('*')
    .linear_regression(
      bool_include_intercept = bool_regression_include_intercept,
      bool_center_x = bool_regression_center,
      bool_center_y = bool_regression_center,
      bool_scale_x = bool_regression_scale,
      bool_scale_y = bool_regression_scale,
      return_type = "r_squared", 
      x_mat = common_mat[,j],
      y_vec = everything_mat[,j]
    )
  })
  names(rsquare_vec) <- colnames(common_mat)
  rsquare_vec
}
