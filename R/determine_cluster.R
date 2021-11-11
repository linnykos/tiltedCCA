.determine_cluster <- function(mat,
                               metacell_clustering_1,
                               metacell_clustering_2,
                               n_idx,
                               num_neigh,
                               tol = 1e-3){
  stopifnot((is.list(metacell_clustering_1) & is.list(metacell_clustering_2)) ||
              (is.factor(metacell_clustering_1) & is.factor(metacell_clustering_2)),
            length(metacell_clustering_1) == nrow(mat),
            length(metacell_clustering_2) == nrow(mat))
  
  list_mode <- is.list(metacell_clustering_1)
  
  if(length(n_idx) < nrow(mat)){
    mat <- mat[n_idx,,drop = F]
    metacell_clustering_1 <- metacell_clustering_1[n_idx]
    metacell_clustering_2 <- metacell_clustering_2[n_idx]
  }
  
  n <- nrow(mat)
  snn_mat <- .form_snn_mat(bool_intersect = T,
                           mat = mat, 
                           num_neigh = num_neigh)
  nn_idx_list <- lapply(1:nrow(snn_mat), function(i){
    .nonzero_col(snn_mat, col_idx = i, bool_value = F)
  })
  
  if(list_mode){
    val_1 <- .compute_set_inclusion(list_1 = metacell_clustering_1,
                                        list_2 = nn_idx_list)
    val_2 <- .compute_set_inclusion(list_1 = metacell_clustering_2,
                                        list_2 = nn_idx_list)
  } else {
    val_1 <- .compute_factor_diversity(factor_anchor = metacell_clustering_1,
                                       factor_other = metacell_clustering_2,
                                       list_nn = nn_idx_list)
    val_2 <- .compute_factor_diversity(factor_anchor = metacell_clustering_2,
                                       factor_other = metacell_clustering_1,
                                       list_nn = nn_idx_list)
  }
  
  min(val_1, val_2)
}

.compute_set_inclusion <- function(list_1, list_2){
  stopifnot(is.list(list_1), is.list(list_2), length(list_1) == length(list_2))
  
  n <- length(list_1)
  overlap_vec <- sapply(1:n, function(i){
    # [[note to self: add a test for this corner case]]
    if(length(list_1[[i]]) == 0 && length(list_2[[i]]) == 0) return(0)
    
    length(intersect(list_1[[i]], list_2[[i]]))/length(unique(c(list_1[[i]], list_2[[i]])))
  })
  
  mean(overlap_vec)
}

.compute_factor_diversity <- function(factor_anchor, 
                                      factor_other,
                                      list_nn){
  stopifnot(is.factor(factor_anchor), is.factor(factor_other), 
            is.list(list_nn),
            length(factor_anchor) == length(factor_other),
            length(factor_anchor) == length(list_nn))
  
  level_anchor_vec <- levels(factor_anchor)
  level_other_vec <- levels(factor_other)
  tabulate_mat <- sapply(level_anchor_vec, function(level_val){
    idx <- which(factor_anchor == level_val)
    vec <- table(factor_other[idx])
    vec/sum(vec)
  })
  if(!is.matrix(tabulate_mat)) tabulate_mat <- matrix(tabulate_mat, nrow = 1, ncol = length(tabulate_mat))
  colnames(tabulate_mat) <- level_anchor_vec
  
  # for each cell in list_nn, look up its celltype in factor_anchor
  # then, see if its neighbor look similar to the intended distribution in tabulate_mat
  # (measured via KL divergence)
  # average this over all celltypes
  avg_anchor_vec <- sapply(level_anchor_vec, function(level_val){
    idx <- which(factor_anchor == level_val)
    kl_vec <- sapply(idx, function(i){
      if(length(list_nn[[i]]) == 0) return(NA)
        
      cells <- factor_other[list_nn[[i]]]
      cells <- cells[!is.na(cells)]
      vec <- table(cells)
      vec <- vec/sum(vec)
      
      -.kl_divergence(query_dist = vec,
                      reference_dist = tabulate_mat[,level_val])
    })
    
    mean(kl_vec, na.rm = T)
  })
  
  mean(avg_anchor_vec)
}

.kl_divergence <- function(query_dist, reference_dist, tol = 1e-6){
  stopifnot(abs(sum(query_dist) - 1) <= tol,
            abs(sum(reference_dist) - 1) <= tol,
            all(query_dist >= 0),
            all(reference_dist >= 0))
  
  idx <- intersect(which(reference_dist > tol), which(query_dist > tol))
  if(length(idx) == 0) return(0)
  sum(sapply(idx, function(i){
    query_dist[i]*log(query_dist[i]/reference_dist[i])
  }))
}
