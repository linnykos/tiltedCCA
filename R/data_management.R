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
  .mult_mat_vec(svd_obj$u, svd_obj$d)
}

#' @export
.get_Dimred.multiSVD <- function(input_obj, ...){
  svd_obj <- .get_SVD(input_obj, ...)
  n <- nrow(svd_obj$u)
  param <- .get_param(input_obj)
  if(param$svd_normalize_singular_value) svd_obj$d <- svd_obj$d*sqrt(n)/svd_obj$d[1]
  dimred <- .mult_mat_vec(svd_obj$u, svd_obj$d)

  .append_rowcolnames(bool_colnames = F, bool_rownames = T,
                      source_obj = input_obj,  target_obj = dimred)
}

###############

.get_postDimred <- function(input_obj, 
                            averaging_mat, ...){
  stopifnot(inherits(input_obj, "multiSVD"))
  
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
                 recenter = recenter, rescale = rescale,
                 ...)
}

#############################

#' @export
.get_metacell <- function(input_obj, ...) UseMethod(".get_metacell")

#' @export
.get_metacell.default <- function(input_obj, ...){
  stop("Class of input_obj not found, using .get_metacell")
}

#' @export
.get_metacell.metacell <- function(input_obj,
                                   resolution, type, what, ...){
  stopifnot(what %in% c("large_clustering_1", "large_clustering_2", "metacell_clustering"),
            type %in% c("list", "factor"),
            resolution %in% c("cell", "metacell"))
  n <- length(input_obj$large_clustering_1)
  if(what %in% c("large_clustering_1", "large_clustering_2") & 
     resolution == "cell"){
    if(what == "large_clustering_1"){
      vec <- input_obj$large_clustering_1
    } else {
      vec <- input_obj$large_clustering_2
    }
    
    if(type == "factor") return(vec) else return(.convert_factor2list(vec))
    
  } else {
    lis <- input_obj$metacell_clustering_list
    if(all(is.null(lis))) {
      if(what %in% c("large_clustering_1", "large_clustering_2") & 
         resolution == "metacell"){
        if(what == "large_clustering_1"){
          vec <- input_obj$large_clustering_1
        } else {
          vec <- input_obj$large_clustering_2
        }
        
        if(type == "factor") return(vec) else return(.convert_factor2list(vec))
      } else {
        return(NULL)
      }
    }
    
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
        if(type == "factor") return(factor(names(lis))) else return(lapply(names(lis),function(x){x}))
      } else {
        if(type == "factor") return(.convert_list2factor(lis, n = n)) else return(lis)
      }
    }
  }
  
  stop("Not able to handle .get_metacell")
}

#' @export
.get_metacell.multiSVD <- function(input_obj, ...){
  stopifnot("metacell_obj" %in% names(input_obj))
  .get_metacell(input_obj$metacell_obj, ...)
}

#############################

.get_SNN <- function(input_obj, 
                     bool_common){
  stopifnot(inherits(input_obj, "multiSVD"))
  
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

.get_Laplacian <- function(input_obj, 
                           bool_common){
  stopifnot(inherits(input_obj, "multiSVD"))
  
  if(bool_common){
    input_obj$laplacian_list[["common_laplacian"]]
  } else {
    if(input_obj$default_assay == 1){
      input_obj$laplacian_list[["laplacian_1"]]
    } else {
      input_obj$laplacian_list[["laplacian_2"]]
    }
  }
}

############################################

.get_tCCAobj <- function(input_obj, 
                         apply_postDimred,
                         what){
  stopifnot(inherits(input_obj, "multiSVD"),
            is.logical(apply_postDimred),
            what %in% c("score", "common_basis", "common_score",
                        "distinct_score", "common_mat", "distinct_mat",
                        "common_dimred", "distinct_dimred"))
  if(what == "common_basis" & apply_postDimred){
    warning(paste0("apply_postDimred=T and what=", what," is possibly nonsensical"))
  }
  default_assay <- .get_defaultAssay(input_obj)
  
  if(what == "score"){
    stopifnot("cca_obj" %in% names(input_obj))
    if(default_assay == 1){
      tmp <- input_obj$cca_obj$score_1
    } else {
      tmp <- input_obj$cca_obj$score_2
    }
    
  } else if(what == "common_basis"){
    stopifnot("tcca_obj" %in% names(input_obj))
    tmp <- input_obj$tcca_obj$common_basis
    
  } else if(what == "common_score"){
    stopifnot("tcca_obj" %in% names(input_obj))
    tmp <- input_obj$tcca_obj$common_score
    
  } else if(what == "distinct_score"){
    stopifnot("tcca_obj" %in% names(input_obj))
    if(default_assay == 1){
      tmp <- input_obj$tcca_obj$distinct_score_1
    } else {
      tmp <- input_obj$tcca_obj$distinct_score_2
    }
    
  } else if(what == "common_mat"){
    if(default_assay == 1){
      stopifnot("common_mat_1" %in% names(input_obj))
      tmp <- input_obj$common_mat_1
    } else {
      stopifnot("common_mat_2" %in% names(input_obj))
      tmp <- input_obj$common_mat_2
    }
    
  } else if(what == "distinct_mat"){
    if(default_assay == 1){
      stopifnot("distinct_mat_1" %in% names(input_obj))
      tmp <- input_obj$distinct_mat_1
    } else {
      stopifnot("distinct_mat_2" %in% names(input_obj))
      tmp <- input_obj$distinct_mat_2
    }
  } else if(what == "common_dimred"){
    if(default_assay == 1){
      stopifnot("common_dimred_1" %in% names(input_obj))
      tmp <- input_obj$common_dimred_1
    } else {
      stopifnot("common_dimred_2" %in% names(input_obj))
      tmp <- input_obj$common_dimred_2
    }
    
  } else if(what == "distinct_dimred"){
    if(default_assay == 1){
      stopifnot("distinct_dimred_1" %in% names(input_obj))
      tmp <- input_obj$distinct_dimred_1
    } else {
      stopifnot("distinct_dimred_2" %in% names(input_obj))
      tmp <- input_obj$distinct_dimred_2
    }
  } else {
    stop(".get_tCCAobj does not have a valid argument")
  }
  
  if(apply_postDimred){
    param <- .get_param(input_obj)
    normalize_row <- param$svd_normalize_row
    normalize_singular_value <- param$svd_normalize_singular_value
    if(default_assay == 1){
      recenter <- param$svd_recenter_1
      rescale <- param$svd_rescale_1
      center <- param$svd_center_1
      dims <- param$svd_dims_1; dims <- dims - min(dims) + 1
      dims <- pmin(dims, ncol(tmp))
      scale <- param$svd_scale_1
    } else {
      recenter <- param$svd_recenter_2
      rescale <- param$svd_rescale_2
      center <- param$svd_center_2
      dims <- param$svd_dims_2; dims <- dims - min(dims) + 1
      dims <- pmin(dims, ncol(tmp))
      scale <- param$svd_scale_2
    }
    
    tmp <- .normalize_svd(input_obj = tmp,
                          averaging_mat = NULL,
                          center = center,
                          dims = dims,
                          normalize_row = normalize_row,
                          normalize_singular_value = normalize_singular_value,
                          recenter = recenter,
                          rescale = rescale,
                          scale = scale)
  } 
  if(any(dim(tmp) == 0)){
    warning("Output incomplete, possibly since the requested matrix is all-0's")
  }
  
  tmp
}

