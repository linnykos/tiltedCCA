form_metacells <- function(input_obj,
                           large_clustering_1, 
                           large_clustering_2,
                           num_metacells = NULL,
                           min_size = 5,
                           verbose = 0){
  stopifnot(inherits(input_obj, "multiSVD"))
  
  if(verbose >= 1) print("Extracting SVD")
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  dimred_1 <- .get_postDimred(input_obj, averaging_mat = NULL)
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  dimred_2 <- .get_postDimred(input_obj, averaging_mat = NULL)
  
  stopifnot(nrow(dimred_1) == nrow(dimred_2))
  n <- nrow(dimred_1)
  
  if(verbose >= 1) print("Computing intersection of clusterings")
  res <- .intersect_clustering(large_clustering_1 = large_clustering_1, 
                               large_clustering_2 = large_clustering_2,
                               min_size = min_size)
  
  if(verbose >= 1) print("Computing metacells")
  metacell_clustering_list <- .form_metacells(dimred_1 = dimred_1, 
                                         dimred_2 = dimred_2,
                                         large_clustering_list = res$large_clustering_list,
                                         num_metacells = num_metacells, 
                                         verbose = verbose)
  
  metacell_obj <- .create_metacell_obj(large_clustering_1 = large_clustering_1, 
                                       large_clustering_2 = large_clustering_2,
                                       metacell_clustering_list = metacell_clustering_list)
  input_obj$metacell_obj <- metacell_obj
  param <- .form_metacells_param(min_size = min_size,
                                 num_metacells = num_metacells)
  param <- .combine_two_named_lists(.get_param(input_obj), param)
  input_obj$param <- param
  
  input_obj
}


##############################

.intersect_clustering <- function(large_clustering_1, 
                                  large_clustering_2,
                                  min_size = 5){
  stopifnot(is.factor(large_clustering_1), is.factor(large_clustering_2),
            length(large_clustering_1) == length(large_clustering_2))
  
  n <- length(large_clustering_1)
  large_clustering_1_list <- .convert_factor2list(large_clustering_1)
  large_clustering_2_list <- .convert_factor2list(large_clustering_2)
  
  total_1 <- length(large_clustering_1_list)
  total_2 <- length(large_clustering_2_list)
  
  idx_df <- data.frame(modality_1 = rep(1:total_1, each = total_2),
                       modality_2 = rep(1:total_2, times = total_1))
  idx_df$total_idx <- (idx_df$modality_1 - 1)*total_2 + idx_df$modality_2
  
  idx <- matrix(NA, nrow = n, ncol = 2)
  rownames(idx) <- 1:n
  for(i in 1:total_1){
    idx[large_clustering_1_list[[i]],1] <- i
  }
  for(i in 1:total_2){
    idx[large_clustering_2_list[[i]],2] <- i
  }
  na_idx <- which(apply(idx, 1, function(vec){any(is.na(vec))}))
  if(length(na_idx) > 0) idx <- idx[-na_idx,]
  
  total_idx <- (idx[,1]-1)*total_2 + idx[,2]
  uniq_vec <- sort(unique(total_idx))
  idx_df <- idx_df[idx_df$total_idx %in% uniq_vec,]
  stopifnot(length(idx_df$total_idx) == length(uniq_vec))
  uniq_vec <- idx_df$total_idx
  
  large_clustering_list <- lapply(uniq_vec, function(i){
    as.numeric(which(total_idx == i))
  })
  names(large_clustering_list) <- paste0("large_", uniq_vec)
  large_clustering_list <- large_clustering_list[sapply(large_clustering_list, length) >= min_size]
  uniq_vec <- names(large_clustering_list)
  total_vec <- paste0("large_", idx_df$total_idx)
  
  clustering_hierarchy_1 <- lapply(1:total_1, function(i){
    idx <- which(idx_df$modality_1 == i)
    if(length(idx) > 0){
      uniq_vec[which(uniq_vec %in% total_vec[idx])]
    } else NA
  })
  names(clustering_hierarchy_1) <- names(large_clustering_1_list)
  clustering_hierarchy_2 <- lapply(1:total_2, function(i){
    idx <- which(idx_df$modality_2 == i)
    if(length(idx) > 0){
      uniq_vec[which(uniq_vec %in% total_vec[idx])]
    } else NA
  })
  names(clustering_hierarchy_2) <- names(large_clustering_2_list)
  
  list(large_clustering_list = large_clustering_list,
       clustering_hierarchy_1 = clustering_hierarchy_1,
       clustering_hierarchy_2 = clustering_hierarchy_2)
}

##############################

.form_metacells <- function(dimred_1, dimred_2,
                            large_clustering_list,
                            num_metacells, 
                            verbose = 0){
  if(is.null(num_metacells)) return(NULL)
  stopifnot(nrow(dimred_1) == nrow(dimred_1))
  n <- nrow(dimred_1)
  df <- .compute_metacell_splits(clustering = large_clustering_list, 
                                 n = n, 
                                 num_metacells = num_metacells)
  
  tmp <- lapply(1:length(large_clustering_list), function(k){
    if(verbose >= 1) print(paste0("On cluster ", k))
    lis <- .compute_metacells(dimred_1 = dimred_1[large_clustering_list[[k]],,drop = F],
                              dimred_2 = dimred_2[large_clustering_list[[k]],,drop = F],
                              k = df$num[k], 
                              row_indices = large_clustering_list[[k]])
    
    num <- strsplit(names(large_clustering_list)[k], split = "_")[[1]][2]
    names(lis) <- paste0("metacell_", num, "_", 1:length(lis))
    lis
  })
  
  do.call(c, tmp)
}

.compute_metacell_splits <- function(clustering, n, num_metacells){
  stopifnot(is.list(clustering), 
            sum(sapply(clustering, length)) <= n,
            length(clustering) <= num_metacells,
            num_metacells > 0, num_metacells %% 1 == 0)
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
  
  while(sum(df$num) < num_metacells){
    idx <- which.max(df$size)
    df$num[idx] <- df$num[idx]+1
    df$size[idx] <- df$total_size[idx]/df$num[idx]
  }
  
  df
}

.compute_metacells <- function(dimred_1, dimred_2,
                               k,
                               row_indices){
  stopifnot(length(row_indices) == nrow(dimred_1),
            nrow(dimred_1) == nrow(dimred_2))
  if(k == 1) return(list(row_indices))
  
  dimred <- cbind(dimred_1, dimred_2)
  
  kmeans_res <- stats::kmeans(dimred, centers = k)
  metacell_clustering <- lapply(1:k, function(kk){
    row_indices[which(kmeans_res$cluster == kk)]
  })
  
  metacell_clustering
}

###########################

.create_metacell_obj <- function(large_clustering_1, 
                                 large_clustering_2,
                                 metacell_clustering_list){
  structure(list(large_clustering_1 = large_clustering_1, 
                 large_clustering_2 = large_clustering_2,
                 metacell_clustering_list = metacell_clustering_list),
            class = "metacell")
}

.form_metacells_param <- function(min_size,
                                  num_metacells){
  list(mc_min_size = min_size, 
       mc_num_metacells = num_metacells)
}

###########################

.generate_averaging_matrix <- function(metacell_clustering_list, n){
  stopifnot(max(unlist(metacell_clustering_list)) <= n)
  
  k <- length(metacell_clustering_list)
  i_vec <- unlist(lapply(1:k, function(i){
    len <- length(metacell_clustering_list[[i]]); rep(i, len)
  }))
  j_vec <- unlist(metacell_clustering_list)
  x_vec <- unlist(lapply(metacell_clustering_list, function(vec){
    len <- length(vec); rep(1/len, len)
  }))
  
  Matrix::sparseMatrix(i = i_vec,
                       j = j_vec,
                       x = x_vec,
                       dims = c(k,n),
                       repr = "C")
}

