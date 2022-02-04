intersect_clustering <- function(large_clustering_1, 
                                 large_clustering_2,
                                 min_size = 5,
                                 n){
  stopifnot(n >= max(unlist(large_clustering_1)),
            n >= max(unlist(large_clustering_2)))
  total_1 <- length(large_clustering_1)
  total_2 <- length(large_clustering_2)
  
  idx_df <- data.frame(modality_1 = rep(1:total_1, each = total_2),
                       modality_2 = rep(1:total_2, times = total_1))
  idx_df$total_idx <- (idx_df$modality_1 - 1)*total_2 + idx_df$modality_2
  
  idx <- matrix(NA, nrow = n, ncol = 2)
  rownames(idx) <- 1:n
  for(i in 1:total_1){
    idx[large_clustering_1[[i]],1] <- i
  }
  for(i in 1:total_2){
    idx[large_clustering_2[[i]],2] <- i
  }
  na_idx <- which(apply(idx, 1, function(vec){any(is.na(vec))}))
  if(length(na_idx) > 0) idx <- idx[-na_idx,]
  
  total_idx <- (idx[,1]-1)*total_2 + idx[,2]
  uniq_vec <- sort(unique(total_idx))
  idx_df <- idx_df[idx_df$total_idx %in% uniq_vec,]
  stopifnot(length(idx_df$total_idx) == length(uniq_vec))
  uniq_vec <- idx_df$total_idx
  
  large_clustering <- lapply(uniq_vec, function(i){
    as.numeric(which(total_idx == i))
  })
  names(large_clustering) <- paste0("large_", uniq_vec)
  large_clustering <- large_clustering[sapply(large_clustering, length) >= min_size]
  uniq_vec <- names(large_clustering)
  total_vec <- paste0("large_", idx_df$total_idx)
  
  clustering_hierarchy_1 <- sapply(1:total_1, function(i){
    idx <- which(idx_df$modality_1 == i)
    if(length(idx) > 0){
      uniq_vec[which(uniq_vec %in% total_vec[idx])]
    } else NA
  })
  names(clustering_hierarchy_1) <- names(large_clustering_1)
  clustering_hierarchy_2 <- sapply(1:total_2, function(i){
    idx <- which(idx_df$modality_2 == i)
    if(length(idx) > 0){
      uniq_vec[which(uniq_vec %in% total_vec[idx])]
    } else NA
  })
  names(clustering_hierarchy_2) <- names(large_clustering_2)
  
  list(large_clustering = large_clustering,
       clustering_hierarchy_1 = clustering_hierarchy_1,
       clustering_hierarchy_2 = clustering_hierarchy_2)
}

form_subclusters <- function(mat_1, 
                             mat_2,
                             large_clustering, 
                             target_k, 
                             verbose = T){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  df <- .compute_subcluster_splits(clustering = large_clustering, 
                                   n = nrow(mat_1), 
                                   target_k = target_k)
  
  tmp <- lapply(1:length(large_clustering), function(k){
    if(verbose) print(paste0("On cluster ", k))
    lis <- .compute_metacells(k = df$num[k], 
                              mat_1 = mat_1[large_clustering[[k]],,drop = F],
                              mat_2 = mat_2[large_clustering[[k]],,drop = F],
                              row_indices = large_clustering[[k]])
    
    num <- strsplit(names(large_clustering)[k], split = "_")[[1]][2]
    names(lis) <- paste0("metacell_", num, "_", 1:length(lis))
    lis
  })
  
  metacell_clustering <- do.call(c, tmp)
}


#######################

.compute_subcluster_splits <- function(clustering, n, target_k){
  stopifnot(is.list(clustering), 
            sum(sapply(clustering, length)) <= n,
            length(clustering) <= target_k,
            target_k > 0, target_k %% 1 == 0)
  tmp <- unlist(clustering)
  stopifnot(all(tmp %% 1 == 0), table(tmp) == 1, all(tmp > 0))
  
  df <- data.frame(total_size = sapply(clustering, length),
                   size = sapply(clustering, length), 
                   num = rep(1, length(clustering)))
  if(all(is.null(names(clustering)))){
    rownames(df) <- 1:length(clustering)
  } else {
    rownames(df) <- names(clustering)
  }
  
  while(sum(df$num) < target_k){
    idx <- which.max(df$size)
    df$num[idx] <- df$num[idx]+1
    df$size[idx] <- df$total_size[idx]/df$num[idx]
  }
  
  df
}

.compute_metacells <- function(k, 
                               mat_1, mat_2,
                               row_indices,
                               max_dim = 20){
  stopifnot(length(row_indices) == nrow(mat_1))
  if(k == 1) return(list(row_indices))
  
  svd_res_1 <- .svd_truncated(mat_1, K = min(c(k, ncol(mat_1), max_dim)), 
                              symmetric = F, 
                              rescale = F, 
                              mean_vec = T, 
                              sd_vec = F, 
                              K_full_rank = F)
  dimred_1 <- .mult_mat_vec(svd_res_1$u, svd_res_1$d/svd_res_1$d[1])
  
  svd_res_2 <- .svd_truncated(mat_2, K = min(c(k, ncol(mat_2), max_dim)), 
                              symmetric = F, 
                              rescale = F, 
                              mean_vec = T, 
                              sd_vec = F, 
                              K_full_rank = F)
  dimred_2 <- .mult_mat_vec(svd_res_2$u, svd_res_2$d/svd_res_2$d[1])
  
  dimred <- cbind(dimred_1, dimred_2)
  
  kmeans_res <- stats::kmeans(dimred, centers = k)
  metacell_clustering <- lapply(1:k, function(kk){
    row_indices[which(kmeans_res$cluster == kk)]
  })
  
  metacell_clustering
}
