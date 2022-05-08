#' @export
postprocess_cell_enrichment <- function(input_obj, ...) UseMethod("postprocess_cell_enrichment")

#' @export
postprocess_cell_enrichment.default <- function(input_obj, 
                                                membership_vec, 
                                                num_neigh,
                                                bool_cosine = T,
                                                bool_intersect = T,
                                                max_subsample = min(1000, length(membership_vec)),
                                                min_deg = 1,
                                                verbose = 0,
                                                ...){
  stopifnot(inherits(input_obj, "matrix"),
            nrow(input_obj) == length(membership_vec), is.factor(membership_vec))
  
  snn <- .form_snn_mat(mat = input_obj,
                       num_neigh = num_neigh,
                       bool_cosine = bool_cosine, 
                       bool_intersect = bool_intersect, 
                       min_deg = min_deg,
                       tol = 1e-4,
                       verbose = verbose)
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample)
  enrichment <- .enrichment(cell_subidx = cell_subidx,
                            g = snn, 
                            membership_vec = membership_vec, 
                            verbose = verbose)
  
  param <- list(num_neigh = num_neigh,
                bool_cosine = bool_cosine,
                bool_intersect = bool_intersect,
                max_subsample = max_subsample,
                min_deg = min_deg)
  
  list(enrichment = enrichment,
       param = param)
}

#' @export
postprocess_cell_enrichment.multiSVD <- function(input_obj, 
                                                 membership_vec, 
                                                 max_subsample = min(1000, length(membership_vec)),
                                                 verbose = 0, ...){
  stopifnot(inherits(input_obj, "multiSVD"),
            all("tcca_obj" %in% names(input_obj)),
            is.factor(membership_vec), length(membership_vec) == nrow(input_obj$svd_1$u))
  
  res <- .construct_snn_from_tcca(input_obj, verbose = verbose)
  cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample)
  
  # compute enrichment scores
  if(verbose > 0) print(paste0(Sys.time(),": Enrichment: Compute common enrichment"))
  enrichment_common <- .enrichment(cell_subidx = cell_subidx,
                                   g = res$snn_common, 
                                   membership_vec = membership_vec, 
                                   verbose = verbose)
  
  if(verbose > 0) print(paste0(Sys.time(),": Enrichment: Compute distinct 1 enrichment"))
  enrichment_distinct_1 <- .enrichment(cell_subidx = cell_subidx,
                                       g = res$snn_distinct_1, 
                                       membership_vec = membership_vec, 
                                       verbose = verbose)
  
  if(verbose > 0) print(paste0(Sys.time(),": Enrichment: Compute distinct 2 enrichment"))
  enrichment_distinct_2 <- .enrichment(cell_subidx = cell_subidx,
                                       g = res$snn_distinct_2, 
                                       membership_vec = membership_vec, 
                                       verbose = verbose)
  
  structure(list(enrichment_common = enrichment_common, 
                 enrichment_distinct_1 = enrichment_distinct_1,
                 enrichment_distinct_2 = enrichment_distinct_2), class = "enrichment")
}

############

.construct_snn_from_tcca <- function(input_obj,
                                     tol = 1e-4,
                                     verbose = 0){
  param <- .get_param(input_obj)
  num_neigh <- param$snn_num_neigh
  bool_cosine <- param$snn_bool_cosine
  bool_intersect <- param$snn_bool_intersect
  min_deg <- param$snn_min_deg
  
  # first construct the common graph
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  if("common_mat_1" %in% names(input_obj)){
    dimred_1 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_mat")
  } else if("common_dimred_1" %in% names(input_obj)){
    dimred_1 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
  } else {
    stop("Cannot find the appropriate common matrix for modality 1")
  }
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  if("common_mat_2" %in% names(input_obj)){
    dimred_2 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_mat")
  } else if("common_dimred_2" %in% names(input_obj)){
    dimred_2 <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "common_dimred")
  } else {
    stop("Cannot find the appropriate common matrix for modality 2")
  }
  dimred <- cbind(dimred_1, dimred_2)
  
  snn_common <- .form_snn_mat(mat = dimred,
                              num_neigh = num_neigh,
                              bool_cosine = bool_cosine, 
                              bool_intersect = bool_intersect, 
                              min_deg = min_deg,
                              tol = tol,
                              verbose = verbose)
  
  # next the distinct graphs
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  if("distinct_mat_1" %in% names(input_obj)){
    dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
  } else if("distinct_dimred_1" %in% names(input_obj)){
    dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_dimred")
  } else {
    stop("Cannot find the appropriate distinct matrix for modality 1")
  }
  snn_distinct_1 <- .form_snn_mat(mat = dimred,
                                  num_neigh = num_neigh,
                                  bool_cosine = bool_cosine, 
                                  bool_intersect = bool_intersect, 
                                  min_deg = min_deg,
                                  tol = tol,
                                  verbose = verbose)
  
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  if("distinct_mat_2" %in% names(input_obj)){
    dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
  } else if("distinct_dimred_2" %in% names(input_obj)){
    dimred <- .get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_dimred")
  } else {
    stop("Cannot find the appropriate distinct matrix for modality 2")
  }
  snn_distinct_2 <- .form_snn_mat(mat = dimred,
                                  num_neigh = num_neigh,
                                  bool_cosine = bool_cosine, 
                                  bool_intersect = bool_intersect, 
                                  min_deg = min_deg,
                                  tol = tol,
                                  verbose = verbose)
  
  list(snn_common = snn_common,
       snn_distinct_1 = snn_distinct_1,
       snn_distinct_2 = snn_distinct_2)
}

.construct_celltype_subsample <- function(membership_vec, max_subsample_cell){
  res <- lapply(levels(membership_vec), function(x){
    idx <- which(membership_vec == x)
    stopifnot(length(idx) > 2)
    
    if(length(idx) <= max_subsample_cell) return(idx)
    
    sample(idx, max_subsample_cell, replace = F)
  })
  
  sort(unlist(res))
}

.enrichment <- function(cell_subidx, 
                        g, 
                        membership_vec, 
                        tol = 1e-3, 
                        verbose = 0){
  stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(g))
  stopifnot(all(cell_subidx %% 1 == 0), all(cell_subidx > 0), 
            length(cell_subidx) == length(unique(cell_subidx)),
            length(cell_subidx) <= length(membership_vec))
  
  if(verbose > 0) print(paste0(Sys.time(),": Computing cell-wise enrichment"))
  
  enrichment_cell_mat <- sapply(1:length(cell_subidx), function(i){
    if(verbose > 1 && length(cell_subidx) > 10 && i %% floor(length(cell_subidx)/10) == 0) cat('*')
    
    .enrichment_cell(g, membership_vec, position = cell_subidx[i], tol = tol)
  })
  rownames(enrichment_cell_mat) <- levels(membership_vec)
  colnames(enrichment_cell_mat) <- rownames(g)[cell_subidx] 
  
  if(verbose > 0) print(paste0(Sys.time(),": Computing cell-type enrichment"))
  res <- lapply(levels(membership_vec), function(celltype){
    idx <- which(membership_vec[cell_subidx] == celltype)
    tmp_mat <- enrichment_cell_mat[,idx,drop = F]
    mean_vec <- matrixStats::rowMeans2(tmp_mat)
    sd_vec <- matrixStats::rowSds(tmp_mat)
    names(mean_vec) <- rownames(tmp_mat)
    
    mean_val <- mean_vec[which(levels(membership_vec) == celltype)]
    sd_val <- sd_vec[which(levels(membership_vec) == celltype)]
    
    list(vec = mean_vec, 
         value = mean_val,
         sd = sd_val)
  })
  names(res) <- levels(membership_vec)
  
  # reorganize values
  df <- data.frame(celltype = levels(membership_vec), 
                   value = sapply(res, function(x){x$value}),
                   sd = sapply(res, function(x){x$sd}))
  rownames(df) <- NULL
  mat <- sapply(res, function(x){x$vec})
  colnames(mat) <- paste0("from_", levels(membership_vec))
  rownames(mat) <- paste0("to_", levels(membership_vec))
  
  list(df = df, enrichment_mat = mat, enrichment_cell_mat = enrichment_cell_mat)
}

.enrichment_cell <- function(g, 
                             membership_vec, 
                             position, 
                             tol = 1e-3){
  stopifnot(position > 0, position <= nrow(g), position %% 1 == 0,
            ncol(g) == nrow(g), is.factor(membership_vec), 
            length(membership_vec) == nrow(g))
  
  n <- nrow(g)
  celltype_vec <- levels(membership_vec)
  k <- length(celltype_vec)
  target_bg <- table(membership_vec)/n
  neigh <- .nonzero_col(g, position, bool_value = F)
  len <- length(neigh)
  if(len == 0){
    return(rep(0, k))
  }
  in_len <- sapply(celltype_vec, function(celltype){
    length(which(membership_vec[neigh] == celltype))
  })
  
  enrichment_vec <- pmax((in_len/len-target_bg+tol) / (1-target_bg+tol), 0)
  names(enrichment_vec) <- celltype_vec
  
  enrichment_vec
}