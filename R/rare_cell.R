#' Detect rare cells
#'
#' @param c_g sparse matrix of class \code{dgCMatrix}, ideally from \code{combine_frnn}
#' representing the combined common embedding, where the non-zero entries represent distances
#' @param d_g_1 sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding in one modality, where the non-zero entries represent distances
#' @param d_g_2 sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding in the other modality, where the non-zero entries represent distances
#' @param idx initial indicies of cells-of-interest
#' @param common_enrich boolean, where \code{TRUE} means you want to find a set of cells
#' with a high common enrichment
#' @param distinct_enrich_1 boolean, where \code{TRUE} means you want to find a set of cells
#' with a high distinct enrichment in one modality
#' @param distinct_enrich_2 boolean, where \code{TRUE} means you want to find a set of cells
#' with a high distinct enrichment in the other modality
#' @param custom_threshold vector of length 3, where each entry is either \code{NA}
#' or a value between \code{0} and \code{1} (inclusive)
#' @param tol small positive integer
#' @param max_tries positive integer
#' @param verbose boolean
#'
#' @return list
#' @export
detect_rare_cell <- function(c_g, d_g_1, d_g_2, idx, 
                             common_enrich, distinct_enrich_1,
                             distinct_enrich_2,
                             custom_threshold = rep(NA, 3),
                             tol = 0.02, max_tries = 10, verbose = F){
  bool_vec <- c(common_enrich, distinct_enrich_1, distinct_enrich_2)
  stopifnot(any(bool_vec), any(!bool_vec))
  stopifnot(inherits(c_g, "dgCMatrix"), inherits(d_g_1, "dgCMatrix"), inherits(d_g_2, "dgCMatrix"))
  stopifnot(all(idx > 0), all(idx <= nrow(c_g)), all(idx %% 1 == 0))
  stopifnot(nrow(c_g) == ncol(c_g), all(dim(c_g) == dim(d_g_1)), all(dim(c_g) == dim(d_g_2)))
  
  baseline_scores <- compute_enrichment_scores(c_g, d_g_1, d_g_2, idx)
  
  iter <- 1
  baseline_idx <- idx
  while(TRUE){
    if(verbose) print(paste0("Iteration ", iter, ": Length of ", length(idx)))
    candidates <- .find_enrichment_candidates(c_g, d_g_1, d_g_2, idx, 
                                              common_enrich, distinct_enrich_1,
                                              distinct_enrich_2, max_tries)
    if(length(candidates) == 0) break()
    try_idx <- 1
    
    while(try_idx <= max_tries){
      bool_continue <- FALSE
      tmp_idx <- c(idx, candidates[try_idx])
      new_scores <- compute_enrichment_scores(c_g, d_g_1, d_g_2, tmp_idx)
      
      bool1 <- .compare_enrichment_scores(baseline_scores[1], new_scores[1],
                                          enrich_bool = common_enrich, 
                                          custom_threshold = custom_threshold[1], tol = tol)
      bool2 <- .compare_enrichment_scores(baseline_scores[2], new_scores[2],
                                          enrich_bool = distinct_enrich_1, 
                                          custom_threshold = custom_threshold[2], tol = tol)
      bool3 <- .compare_enrichment_scores(baseline_scores[3], new_scores[3],
                                          enrich_bool = distinct_enrich_2, 
                                          custom_threshold = custom_threshold[2], tol = tol)
      
      if(bool1 & bool2 & bool3){
        idx <- unique(tmp_idx); bool_continue <- TRUE; break()
      } else {
        try_idx <- try_idx + 1
      }
    }
    
    if(!bool_continue) break()
    if(verbose) print(round(new_scores, 2))
    iter <- iter + 1
  }
  
  list(idx = idx, scores = new_scores,
       baseline_idx = baseline_idx, baseline_scores = baseline_scores)
}

#' Compute enrichment scores
#'
#' @param c_g sparse matrix of class \code{dgCMatrix}, ideally from \code{combine_frnn}
#' representing the combined common embedding, where the non-zero entries represent distances
#' @param d_g_1 sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding in one modality, where the non-zero entries represent distances
#' @param d_g_2 sparse matrix of class \code{dgCMatrix} from \code{construct_frnn}
#' representing the distinct embedding in the other modality, where the non-zero entries represent distances
#' @param idx initial indicies of cells-of-interest
#'
#' @return vector
#' @export
compute_enrichment_scores <- function(c_g, d_g_1, d_g_2, idx){
  clisi_info1 <- sapply(1:length(idx), function(i){
    .clisi_cell(c_g, idx, position = i)
  })
  clisi_info2 <- sapply(1:length(idx), function(i){
    .clisi_cell(d_g_1, idx, position = i)
  })
  clisi_info3 <- sapply(1:length(idx), function(i){
    .clisi_cell(d_g_2, idx, position = i)
  })
  
  c(stats::median(clisi_info1["clisi_score",]), 
    stats::median(clisi_info2["clisi_score",]),
    stats::median(clisi_info3["clisi_score",]))
}

#################

.compare_enrichment_scores <- function(baseline_score, new_score,
                                       enrich_bool, custom_threshold, tol){
  stopifnot(length(custom_threshold) == 1, is.logical(enrich_bool))
  
  if(enrich_bool){
    if(is.na(custom_threshold)) {
      new_score >= baseline_score - tol
    } else {
      new_score >= custom_threshold
    }
  } else {
    if(is.na(custom_threshold)) {
      new_score <= baseline_score + tol
    } else {
      new_score <= custom_threshold
    }
  }
}

.find_enrichment_candidates <- function(c_g, d_g_1, d_g_2, idx, 
                                        common_enrich, distinct_enrich_1,
                                        distinct_enrich_2, max_tries){
  bool_vec <- c(common_enrich, distinct_enrich_1, distinct_enrich_2)
  bool_idx <- which(bool_vec)
  graph_list <- list(c_g, d_g_1, d_g_2)
  
  neigh_list <- lapply(bool_idx, function(k){
    neigh_vec <- unlist(lapply(idx, function(i){.nonzero_col(graph_list[[k]], i, bool_value = F)}))
    neigh_vec <- setdiff(neigh_vec, idx)
  })
  
  if(length(neigh_list) == 1){
    if(length(neigh_list[[1]]) == 0) return(numeric(0))
    return(as.numeric(names(sort(table(neigh_list[[1]]), decreasing = T))[1:max_tries]))
  } else{
    all_idx <- sort(.common_intersection(neigh_list))
    if(length(all_idx) == 0) return(numeric(0))
    
    df <- sapply(1:length(neigh_list), function(k){
      table(intersect(neigh_list[[k]], all_idx))
    })
    min_neigh <- apply(df, 1, min)
    all_idx[order(min_neigh, decreasing = T)][1:max_tries]
  }
}

.common_intersection <- function(lis){
  Reduce(intersect, lis)
}