postprocess_alignment <- function(multiSVD_obj,
                                  seurat_obj,
                                  bool_common_center = T,
                                  bool_common_scale = T,
                                  bool_everything_center = T,
                                  bool_everything_scale = T,
                                  bool_regression_include_intercept = T,
                                  bool_regression_center = T,
                                  bool_regression_scale = T,
                                  multiSVD_assay = 1,
                                  seurat_assay = Seurat::DefaultAssay(seurat_obj),
                                  seurat_slot = "data",
                                  verbose = 1){
  stopifnot(multiSVD_assay %in% c(1,2),
            seurat_slot %in% c("counts", "data", "scale.data"))
  
  if(verbose > 0) print("Gathering ingredients")
  multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, 
                                    assay = multiSVD_assay)
  common_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, 
                             apply_postDimred = F,
                             what = "common_mat")
  if(seurat_slot == "counts"){
    everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@counts)
  } else if(seurat_slot == "data"){
    everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@data)
  } else if(seurat_slot == "scale.data"){
    everything_mat <- Matrix::t(seurat_obj[[seurat_assay]]@scale.data)
  } else {
    stop("seurat_slot invalid")
  }
  
  everything_mat <- everything_mat[,colnames(common_mat)]
  stopifnot(all(dim(common_mat) == dim(everything_mat)))
  
  if(bool_common_center | bool_common_scale){
    common_mat <- scale(common_mat,
                        center = bool_common_center, 
                        scale = bool_common_scale)
  }
  if(bool_everything_center | bool_everything_scale){
    everything_mat <- scale(everything_mat,
                            center = bool_everything_center, 
                            scale = bool_everything_scale)
  }
  
  if(verbose > 0) print("Starting to perform regressions")
  
  rsquare_vec <- sapply(1:ncol(common_mat), function(j){
    if(verbose > 0 && ncol(common_mat) > 10 && j %% floor(ncol(common_mat)/10) == 0) cat('*')
    .univariate_regression(
      bool_include_intercept = bool_regression_include_intercept,
      bool_center_x = bool_regression_center,
      bool_center_y = bool_regression_center,
      bool_scale_x = bool_regression_scale,
      bool_scale_y = bool_regression_scale,
      return_type = "r_squared", 
      x_vec = common_mat[,j],
      y_vec = everything_mat[,j]
    )
  })
  names(rsquare_vec) <- colnames(common_mat)
  rsquare_vec
}
