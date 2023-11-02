#' Compute the negative log-10 p-values across all cell types
#' 
#' This function does the following: There are varible cell types
#' stored in \code{de_list$level_vec}. For each cell type, enumerate
#' the p-value (across different features, i.e., genes or proteins) 
#' when comparing that cell type against all other cell types.
#' Return the p-value at the quantile of \code{quantile_value}
#' among all these comparisons for this particular cell type.
#' Then, repeat for all cell types, resulting in a p-value for each
#' cell type and each feature. 
#' Lastly, for each feature, take the maximum of the -log10 p-value (across
#' all the cell types).
#'
#' @param de_list output of the \code{tiltedCCA::differential_expression} function
#' @param exclude_celltypes a vector of strings containing cell types to be 
#' excluded in this function's calculations. All the cell types should be in \code{de_list$level_vec}
#' @param maximum_output any output (i.e., negative log10 pvalue) larger than 
#' \code{maximum_output} are set to be \code{maximum_output}
#' @param quantile_value for a particular variable, when comparing one cell type
#' to all other cell types, output the p-value at the quantile of \code{quantile_value} . 
#' (prior to taking the maximum of the -log10 across all cell types)
#' @param verbose non-negative integer
#'
#' @return a vector of the -log10 p-value for each variable in \code{de_list}
#' @export
postprocess_depvalue <- function(de_list,
                                 exclude_celltypes = NULL,
                                 maximum_output = NULL,
                                 quantile_value = 0.75,
                                 verbose = 1){
  stopifnot(all(sort(names(de_list)) == sort(c("combn_mat", "de_list", "level_vec"))))
  if(!all(is.null(exclude_celltypes))){
    stopifnot(all(exclude_celltypes %in% de_list$level_vec))
  }
  
  celltype_names <- de_list$level_vec
  var_names <- sort(unique(unlist(lapply(de_list$de_list, function(x){rownames(x)}))))
  if(length(var_names) < 10) verbose <- 0
  
  logpval_vec <- sapply(1:length(var_names), function(kk){
    if(verbose > 0 & kk %% floor(length(var_names)/10) == 0) cat('*')
    var_name <- var_names[kk]
    
    # cycle through all the celltypes
    celltype_vec <- sapply(1:length(de_list$level_vec), function(i){
      idx <- which(de_list$combn_mat == i, arr.ind = T)[,2]
      vec <-  sapply(idx, function(j){
        idx <- which(rownames(de_list$de_list[[j]]) == var_name)
        if(length(idx) == 0) return(1)
        de_list$de_list[[j]][idx, "p_val"]
      })
      stats::quantile(vec, probs = quantile_value)
    })
    names(celltype_vec) <- celltype_names
    if(!all(is.null(exclude_celltypes))) celltype_vec <- celltype_vec[!names(celltype_vec) %in% exclude_celltypes]
    
    max(-log10(celltype_vec))
  })
  names(logpval_vec) <- var_names
  if(!is.null(maximum_output)) logpval_vec <- pmin(logpval_vec, maximum_output)
  
  logpval_vec
}


#' Find the marker genes the separate one cell type from all other cell types
#' 
#' Based on the output of \code{tiltedCCA::differential_expression}, for a
#' particular cell type, 1) when comparing against every other cell type,
#' find all the variables that have an adjusted p-value
#' smaller than \code{p_val_thres} and also have a log fold change larger than
#' \code{log_thres_name}, and 2) if the variable was selected for more than
#' \code{num_quantile_threshold} percentage of the comparisons with other cell 
#' types, we deem it as a marker gene for this cell type.
#'
#' @param de_list output of the \code{tiltedCCA::differential_expression} function
#' @param log_thres_name character of the column in each of \code{de_list$de_list} for the log fold change
#' @param log_thres positive value, for selecting variables with log fold change larger than \code{log_thres}
#' @param num_quantile_threshold value between 0 and 1
#' @param p_val_name character of the column in each of \code{de_list$de_list} for the p-value
#' @param p_val_adj_name character of the column in each of \code{de_list$de_list} for the multiple-testing adjusted p-value
#' @param p_val_thres small positive value between 0 and 1, for selecting variables with p-values smaller than \code{p_val_thres}
#'
#' @return a list of variable vectors, one vector for each element in \code{de_list$level_vec}
#' @export
postprocess_marker_variables <- function(de_list,
                                         log_thres_name = "avg_log2FC",
                                         log_thres = 2,
                                         num_quantile_threshold = 0.25,
                                         p_val_name = "p_val",
                                         p_val_adj_name = "p_val_adj",
                                         p_val_thres = 1e-4){
  stopifnot(all(sort(names(de_list)) == sort(c("combn_mat", "de_list", "level_vec"))))
  level_vec <- de_list$level_vec
  
  num_instances <- round(length(level_vec)*num_quantile_threshold)
  variable_list <- lapply(1:length(level_vec), function(i){
    idx <- which(de_list$combn_mat == i, arr.ind = T)[,2]
    lis <- lapply(de_list$de_list[idx], function(mat){
      idx1 <- which(abs(mat[,log_thres_name]) >= log_thres)
      idx2 <- which(abs(mat[,p_val_adj_name]) <= p_val_thres)
      rownames(mat)[intersect(idx1, idx2)]
    })
    vec <- unlist(lis)
    tab <- table(vec)
    name_vec <- names(tab)[which(tab >= num_instances)]
    
    # find the genes that have the smallest median p-value, if there's less than 5
    if(length(name_vec) <= 5){
      all_variables <- unique(unlist(lapply(de_list$de_list[idx], rownames)))
      all_variables <- all_variables[!all_variables %in% name_vec]
      
      p_val_vec <- sapply(all_variables, function(gene){
        val <- stats::median(sapply(idx, function(j){
          row_idx <- which(rownames(de_list$de_list[[j]]) == gene)
          if(length(row_idx) == 0) return(1)
          de_list$de_list[[j]][row_idx, p_val_name]
        }))
      })
      
      name_vec <- c(name_vec, all_variables[order(p_val_vec, decreasing = F)][1:(5-length(name_vec))])
    }
    
    name_vec
  })
  
  names(variable_list) <- level_vec
  variable_list
}

