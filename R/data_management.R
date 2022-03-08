#' @export
.get_defaultAssay <- function(input_obj, ...) UseMethod(".get_defaultAssay")

#' @export
.get_defaultAssay.default <- function(input_obj, ...){
  numeric(0)
}

#' @export
.get_defaultAssay.multiSVD <- function(input_obj, ...){
  input_obj$default_assay
}

####################

#' @export
.set_defaultAssay <- function(input_obj, ...) UseMethod(".set_defaultAssay")

#' @export
.set_defaultAssay.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .set_defaultAssay")
}

#' @export
.set_defaultAssay.multiSVD <- function(input_obj, 
                                       assay, ...){
  input_obj$default_assay <- assay
  input_obj
}

####################

.get_param <- function(input_obj){
  input_obj$param
}

###############

#' @export
.get_SVD <- function(input_obj, ...) UseMethod(".get_SVD")

#' @export
.get_SVD.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .get_SVD")
}

#' @export
.get_SVD.Seurat <- function(input_obj, ...){
  stop("Not yet coded")
}

#' @export
.get_SVD.matrix <- function(input_obj,
                            center,
                            dims,
                            scale, ...){
  K <- max(dims)
  tmp <- .svd_truncated(mat = input_obj, 
                        K = K, 
                        symmetric = F, 
                        rescale = F,
                        mean_vec = center, 
                        sd_vec = scale,
                        K_full_rank = (min(dim(input_obj)) == K))
  tmp <- .check_svd(tmp, dims = dims)
  .append_rowcolnames(bool_colnames = T, bool_rownames = T,
                      source_obj = input_obj,  target_obj = tmp)
}

#' @export
.get_SVD.dgCMatrix <- function(input_obj,
                               center,
                               dims,
                               scale, ...){
  .get_SVD.matrix(input_obj = input_obj,
                  center = center, dims = dims, scale = scale,
                  ...)
}

#' @export
.get_SVD.svd <- function(input_obj, ...){input_obj}

#' @export
.get_SVD.multiSVD <- function(input_obj, ...){
  if(input_obj$default_assay == 1){
    input_obj$svd_1
  } else {
    input_obj$svd_2
  }
}

###############

#' @export
.get_Dimred <- function(input_obj, ...) UseMethod(".get_Dimred")

#' @export
.get_Dimred.default <- function(input_obj, 
                                normalize_singular_value, ...){
  svd_obj <- .get_SVD(input_obj, ...)
  n <- nrow(svd_obj$u)
  if(normalize_singular_value) svd_obj$d <- svd_obj$d*sqrt(n)/svd_obj$d[1]
  dimred <- .mult_mat_vec(svd_obj$u, svd_obj$d)
  
  .append_rowcolnames(bool_colnames = F, bool_rownames = T,
                      source_obj = input_obj,  target_obj = dimred)
}

###############

#' @export
.get_postDimred <- function(input_obj, ...) UseMethod(".get_postDimred")

#' @export
.get_postDimred.default <- function(input_obj, ...){
  warning("Class of input_obj not found, using .get_Dimred")
  .get_Dimred(input_obj, ...)
}

#' @export
.get_postDimred.multiSVD <- function(input_obj, 
                                     averaging_mat, ...){
  if(input_obj$default_assay == 1){
    svd_obj <- input_obj$svd_1
  } else {
    svd_obj <- input_obj$svd_2
  }
  
  default_assay <- .get_defaultAssay(input_obj)
  param <- .get_param(input_obj)
  if(default_assay == 1){
    recenter <- param$svd_recenter_1; rescale <- param$svd_rescale_1
  } else if(default_assay == 2){
    recenter <- param$svd_recenter_2; rescale <- param$svd_rescale_2
  } else {
    stop("Invalid default assay in .get_postDimred.multiSVD")
  }
  
  .normalize_svd(input_obj = svd_obj,
                 averaging_mat = averaging_mat,
                 normalize_row = param$svd_normalize_row,
                 normalize_singular_value = param$svd_normalize_singular_value,
                 recenter = recenter, rescale = rescale)
}

#############################

#' @export
.get_metacell <- function(input_obj, ...) UseMethod(".get_metacell")

#' @export
.get_metacell.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .get_metacell")
}

#' @export
.get_metacell.meta <- function(input_obj,
                               resolution, type, what, ...){
  stopifnot(what %in% c("large_clustering_1", "large_clustering_2", "metacell_clustering"),
            type %in% c("list", "factor"),
            resolution %in% c("cell", "metacell"))
  if(what %in% c("large_clustering_1", "large_clustering_2") & 
     resolution == "cell"){
    if(what == "large_clustering_1"){
      vec <- input_obj$large_clustering_1
    } else {
      vec <- input_obj$large_clustering_2
    }
    
    if(type == "factor") return(vec) else return(.convert_factor2list(vec))
    
  } else {
    lis <- input_obj$metacell_clustering
    if(what %in% c("large_clustering_1", "large_clustering_2")){
      if(what == "large_clustering_1"){
        source_vec <- input_obj$large_clustering_1
      } else {
        source_vec <- input_obj$large_clustering_2
      }
      target_vec <- rep(NA, length(lis))
      for(i in 1:length(target_vec)){
        tab <- table(source_vec[lis[[i]]])
        target_vec[i] <- names(tab)[which.max(tab)]
      }
      target_vec <- as.factor(target_vec)
      if(type == "factor") return(target_vec) else return(.convert_factor2list(target_vec))
      
    } else {
      if(resolution == "metacell") {
        # trivial output
        if(type == "factor") return(factor(1:length(lis))) else return(lapply(1:length(lis),function(x){x}))
      } else {
        if(type == "factor") return(.convert_list2factor(lis)) else return(lis)
      }
    }
  }
  
  stop("Not able to handle .get_metacell")
}

#' @export
.get_metacell.multiSVD <- function(input_obj,
                                   resolution, type, what, ...){
  stopifnot("metacell_obj" %in% names(input_obj))
  .get_metacell(input_obj$metacell_obj)
}

#############################

#' @export
.get_snn <- function(input_obj, ...) UseMethod(".get_snn")

#' @export
.get_snn.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .get_snn")
}

#' @export
.get_snn.snn <- function(input_obj, 
                         bool_common, ...){
  if(bool_common){
    input_obj$snn_list[["common_snn"]]
  } else {
    if(input_obj$default_assay == 1){
      input_obj$snn_list[["snn_1"]]
    } else {
      input_obj$snn_list[["snn_2"]]
    }
  }
}

############################################

#' @export
.get_laplacian <- function(input_obj, ...) UseMethod(".get_laplacian")

#' @export
.get_laplacian.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .get_laplacian")
}

#' @export
.get_laplacian.snn <- function(input_obj, 
                               bool_common, ...){
  if(bool_common){
    input_obj$laplacian_list[["common_snn"]]
  } else {
    if(input_obj$default_assay == 1){
      input_obj$laplacian_list[["laplacian_1"]]
    } else {
      input_obj$laplacian_list[["laplacian_2"]]
    }
  }
}


