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
#' @param deg_threshold for a cell to be included into the set, it must be
#' neighbors with at least \code{deg_threshold} of the cells 
#' already in the set
#' @param max_size maximum size of the desired set
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
                             deg_threshold = 0, max_size = NA,
                             tol = 0.02, max_tries = 10, verbose = F){
  bool_vec <- c(common_enrich, distinct_enrich_1, distinct_enrich_2)
  stopifnot(any(bool_vec), any(!bool_vec))
  stopifnot(inherits(c_g, "dgCMatrix"), inherits(d_g_1, "dgCMatrix"), inherits(d_g_2, "dgCMatrix"))
  stopifnot(all(idx > 0), all(idx <= nrow(c_g)), all(idx %% 1 == 0))
  stopifnot(nrow(c_g) == ncol(c_g), all(dim(c_g) == dim(d_g_1)), all(dim(c_g) == dim(d_g_2)))
  
  baseline_scores <- compute_enrichment_scores(c_g, d_g_1, d_g_2, idx)
  new_scores <- baseline_scores
  
  iter <- 1
  baseline_idx <- idx
  while(TRUE){
    if(verbose) print(paste0("Iteration ", iter, ": Length of ", length(idx)))
    candidates <- .find_enrichment_candidates(c_g, d_g_1, d_g_2, idx, 
                                              common_enrich, distinct_enrich_1,
                                              distinct_enrich_2, 
                                              deg_threshold, max_tries)
    if(length(candidates) == 0) break()
    try_idx <- 1
    
    while(try_idx <= min(max_tries, length(candidates))){
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
    
    if(!is.na(max_size) && length(idx) > max_size) break()
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
                                        distinct_enrich_2, deg_threshold,
                                        max_tries){
  bool_vec <- c(common_enrich, distinct_enrich_1, distinct_enrich_2)
  bool_idx <- which(bool_vec)
  len <- length(bool_idx)
  graph_list <- list(c_g, d_g_1, d_g_2)
  
  # extract neighbors
  neigh_list <- lapply(bool_idx, function(k){
    unlist(lapply(idx, function(i){.nonzero_col(graph_list[[k]], i, bool_value = F)}))
  })
  rank_list <- lapply(bool_idx, function(k){
    unlist(lapply(idx, function(i){rank(.nonzero_col(graph_list[[k]], i, bool_value = T))}))
  })
  
  # remove elements already in the list
  for(k in 1:len){
    tmp <- which(!neigh_list[[k]] %in% idx)
    if(length(tmp) == 0) return(numeric(0))
    neigh_list[[k]] <- neigh_list[[k]][tmp]
    rank_list[[k]] <- rank_list[[k]][tmp]
  }
  
  if(len == 1){
    if(length(neigh_list[[1]]) == 0) return(numeric(0))
    summary_df <- .organize_candidate_df(neigh_list[[1]], rank_list[[1]], sort = T)
    summary_df <- summary_df[summary_df$count >= deg_threshold*length(idx),]
    if(nrow(summary_df) == 0) return(numeric(0))
    summary_df$idx[1:min(nrow(summary_df), max_tries)]
    
  } else{
    all_idx <- sort(.common_intersection(neigh_list))
    if(length(all_idx) == 0) return(numeric(0))
 
    # keep only the element in the intersection
    for(k in 1:length(neigh_list)){
      tmp <- which(neigh_list[[k]] %in% all_idx)
      neigh_list[[k]] <- neigh_list[[k]][tmp]
      rank_list[[k]] <- rank_list[[k]][tmp]
    }
    
    summary_list <- lapply(1:len, function(k){
      .organize_candidate_df(neigh_list[[k]], rank_list[[k]], sort = F)
    })
    
    summary_df <- .merge_summary_dfs(summary_list, sort = T)
    summary_df <- summary_df[summary_df$count >= deg_threshold*length(idx),]
    if(nrow(summary_df) == 0) return(numeric(0))
    summary_df$idx[1:min(nrow(summary_df), max_tries)]
  }
}

.common_intersection <- function(lis){
  Reduce(intersect, lis)
}

# remember the rank here is negative (so larger number is closer points)
.organize_candidate_df <- function(neigh_vec, rank_vec, sort){
  stopifnot(length(neigh_vec) == length(rank_vec))
  
  uniq_val <- sort(unique(neigh_vec))
  summary_mat <- sapply(uniq_val, function(x){
    tmp <- which(neigh_vec == x)
    c(idx = x, count = length(tmp), rank = -stats::median(rank_vec[tmp]))
  })
  summary_df <- as.data.frame(t(summary_mat))
  
  if(sort){
    summary_df[order(summary_df$count, summary_df$rank, decreasing = T),]
  } else {
    summary_df
  }
}

.merge_summary_dfs <- function(summary_list, sort){
  len <- length(summary_list)
  count_vec <- apply(sapply(1:len, function(k){summary_list[[k]]$count}), 1, min)
  rank_vec <- apply(sapply(1:len, function(k){summary_list[[k]]$rank}), 1, min)
  
  summary_df <- data.frame(idx = summary_list[[1]]$idx, count = count_vec, rank = rank_vec)
  if(sort){
    summary_df[order(summary_df$count, summary_df$rank, decreasing = T),]
  } else {
    summary_df
  }
}