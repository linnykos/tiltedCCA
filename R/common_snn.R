.compute_common_snn_softclustering <- function(snn_1, snn_2,
                                               num_neigh,
                                               verbose = 0){
  stopifnot(all(dim(snn_1) == dim(snn_2)),
            ncol(snn_1) == nrow(snn_1))
  
  n <- nrow(snn_1)
  if(verbose >= 1) print("Getting NN information")
  nn_list <- lapply(1:n, function(i){
    if(verbose >= 2 && n > 10 && i %% floor(n/10) == 0) cat('*')
    nn_1 <- .nonzero_col(snn_1, 
                         col_idx = i,
                         bool_value = F)
    nn_2 <- .nonzero_col(snn_2, 
                         col_idx = i,
                         bool_value = F)
    
    list(common = intersect(nn_1, nn_2),
         distinct_1 = setdiff(nn_1, nn_2),
         distinct_2 = setdiff(nn_2, nn_1))
  })
  
  if(verbose >= 1) print("Determining common NN's")
  nn_common_list <- lapply(1:n, function(i){
    beta <- .compute_quadratic(common_val = length(nn_list[[i]]$common),
                               distinct_val_1 = length(nn_list[[i]]$distinct_1),
                               distinct_val_2 = length(nn_list[[i]]$distinct_2))
    if(length(nn_list[[i]]$distinct_1) <= length(nn_list[[i]]$distinct_2)){
      tmp <- nn_list[[i]]$distinct_2
      return(c(nn_list[[i]]$common, nn_list[[i]]$distinct_1, sample(tmp, size = ceiling(beta*length(tmp)))))
    } else {
      tmp <- nn_list[[i]]$distinct_1
      return(c(nn_list[[i]]$common, nn_list[[i]]$distinct_2, sample(tmp, size = ceiling(beta*length(tmp)))))
    }
  })
  
  .convert_list2sparse(n = n,
                       nn_list = nn_common_list,
                       rowname_vec = rownames(snn_1))
}

.compute_common_snn_hardclustering <- function(snn_1, snn_2,
                                               clustering_1, clustering_2,
                                               num_neigh,
                                               verbose = 0){
  stopifnot(is.factor(clustering_1), is.factor(clustering_2),
            length(clustering_1) == nrow(snn_1),
            length(clustering_1) == nrow(snn_2),
            all(dim(snn_1) == dim(snn_2)),
            ncol(snn_1) == nrow(snn_1))
  
  if(any(table(clustering_1) == 0)) clustering_1 <- droplevels(clustering_1)
  if(any(table(clustering_2) == 0)) clustering_2 <- droplevels(clustering_2)
  
  n <- nrow(snn_1)
  
  nn_list <- .l2_selection_nn(clustering_1 = clustering_1, 
                              clustering_2 = clustering_2,
                              num_neigh = num_neigh,
                              snn_1 = snn_1, 
                              snn_2 = snn_2,
                              verbose = verbose)
  
  .convert_list2sparse(n = n,
                       nn_list = nn_list,
                       rowname_vec = rownames(snn_1))
}

#######################

.compute_quadratic <- function(common_val,
                               distinct_val_1,
                               distinct_val_2){
  if(distinct_val_1 == distinct_val_2) return(1)
  if(distinct_val_1 <= distinct_val_2){
    d1 <- distinct_val_1; d2 <- distinct_val_2
  } else {
    d2 <- distinct_val_1; d1 <- distinct_val_2
  }
  c <- common_val
  
  za <- d2^2
  zb <- d1*d2 + 2*d2*c
  zc <- -(d1^2 + d1*d2 + d1*c + d2*c)
  
  return(max(c(min(c(-zb + sqrt(zb^2-4*za*zc))/(2*za), 1)), 0))
}

.l2_selection_nn <- function(clustering_1, clustering_2,
                             num_neigh,
                             snn_1, snn_2,
                             verbose = 0){
  n <- nrow(snn_1)
  
  nn_list <- lapply(1:n, function(i){
    if(verbose == 2 && n > 10 && i %% floor(n/10) == 0) cat('*')
    if(verbose == 3) print(paste0("On node ", i))
    
    nn_1 <- .nonzero_col(snn_1, 
                         col_idx = i,
                         bool_value = F)
    nn_2 <- .nonzero_col(snn_2, 
                         col_idx = i,
                         bool_value = F)
    
    prior_1 <- table(clustering_1[nn_2]); prior_1 <- prior_1/sum(prior_1)
    prior_2 <- table(clustering_2[nn_1]); prior_2 <- prior_2/sum(prior_2)
    
    idx_all <- sort(unique(c(nn_1, nn_2)))
    idx_df <- data.frame(idx = idx_all, 
                         clustering_1 = clustering_1[idx_all],
                         clustering_2 = clustering_2[idx_all])
    na_idx <- unique(c(which(is.na(idx_df$clustering_1)), which(is.na(idx_df$clustering_2))))
    if(length(na_idx) > 0){
      idx_df <- idx_df[-na_idx,]
    }
    idx_all <- idx_df$idx
    if(length(idx_all) < num_neigh) return(idx_df$idx)
    
    obs_tab <- table(clustering_1[idx_all], clustering_2[idx_all])
    obs_tab <- .remove_all_zeros_rowcol(obs_tab)
    
    if(all(dim(obs_tab) == c(1,1))){
      sample(idx_df$idx, num_neigh)
    } else {
      desired_tab <- .l2_selection_qp(num_neigh = num_neigh,
                                      obs_tab = obs_tab,
                                      prior_1 = prior_1[names(prior_1) %in% rownames(obs_tab)],
                                      prior_2 = prior_2[names(prior_2) %in% colnames(obs_tab)])
      
      .select_l2_cells(desired_tab = desired_tab,
                       idx_df = idx_df)
    }
  })
  
  names(nn_list) <- rownames(snn_1)
  nn_list
}

.remove_all_zeros_rowcol <- function(mat){
  stopifnot(is.matrix(mat), all(mat >= 0))
  
  row_idx <- which(matrixStats::rowSums2(mat) > 0)
  col_idx <- which(matrixStats::colSums2(mat) > 0)
  mat[row_idx, col_idx, drop = F]
}

.l2_selection_qp <- function(num_neigh,
                             obs_tab,
                             prior_1,
                             prior_2,
                             tol = 1e-3,
                             tol_diag = 1){
  stopifnot(nrow(obs_tab) == length(prior_1),
            ncol(obs_tab) == length(prior_2),
            abs(sum(prior_1) - 1) <= tol,
            abs(sum(prior_2) - 1) <= tol)
  p1 <- nrow(obs_tab); p2 <- ncol(obs_tab)
  if(sum(obs_tab) <= num_neigh) return(obs_tab)
  
  obj_mat <- .generate_objmat_l2(obs_tab)
  diag(obj_mat) <- diag(obj_mat) + tol_diag 
  #[[note to self: see https://stats.stackexchange.com/questions/138730/what-are-the-differences-between-various-r-quadratic-programming-solvers, switch solvers maybe?]]
  obj_vec <- .generate_objvec_l2(num_neigh = num_neigh,
                                 obs_tab = obs_tab,
                                 prior_1 = prior_1,
                                 prior_2 = prior_2)
  
  tmp <- .generate_constr_l2(obs_tab,  
                             target_count = num_neigh)
  constr_mat <- tmp$mat; constr_vec = tmp$vec
  
  res <- quadprog::solve.QP(Dmat = obj_mat, 
                            dvec = obj_vec,
                            Amat = constr_mat,
                            bvec = constr_vec,
                            meq = 1)
  
  round(res$solution*obs_tab)
}

.select_l2_cells <- function(desired_tab, idx_df){
  p1 <- nrow(desired_tab); p2 <- ncol(desired_tab)
  combn_mat <- cbind(rep(1:p1, each = p2), rep(1:p2, times = p1))
  lis <- lapply(1:nrow(combn_mat), function(k){
    i <- combn_mat[k,1]; j <- combn_mat[k,2]
    cluster_1 <- rownames(desired_tab)[i]
    cluster_2 <- colnames(desired_tab)[j]
    if(desired_tab[i,j] == 0) return(numeric(0))
    
    potential_idx <- intersect(which(idx_df$clustering_1 == cluster_1),
                               which(idx_df$clustering_2 == cluster_2))
    stopifnot(length(potential_idx) >= desired_tab[i,j])
    idx_df$idx[sample(potential_idx, size = desired_tab[i,j])]
  })
  
  sort(unlist(lis))
}

###################

.generate_indices_l2 <- function(bool_row,
                                 idx,
                                 p1, p2){
  if(bool_row) {
    idx + c(0:(p2-1))*p1
  } else {
    ((idx-1)*p1+1)+c(0:(p1-1))
  }
}

.generate_objmat_l2 <- function(obs_tab){
  p1 <- nrow(obs_tab); p2 <- ncol(obs_tab)
  
  row_objmat_list <- lapply(1:p1, function(i){
    idx_vec <- .generate_indices_l2(bool_row = T,
                                    idx = i,
                                    p1 = p1, p2 = p2)
    .generate_objmat_l2_helper(idx_vec = idx_vec,
                               num_vec = obs_tab[idx_vec],
                               total_len = p1*p2)
  })
  col_objmat_list <- lapply(1:p2, function(j){
    idx_vec <- .generate_indices_l2(bool_row = F,
                                    idx = j,
                                    p1 = p1, p2 = p2)
    .generate_objmat_l2_helper(idx_vec = idx_vec,
                               num_vec = obs_tab[idx_vec],
                               total_len = p1*p2)
  })
  Reduce("+", row_objmat_list) + Reduce("+", col_objmat_list)
}

.generate_objvec_l2 <- function(num_neigh,
                                obs_tab,
                                prior_1,
                                prior_2){
  p1 <- nrow(obs_tab); p2 <- ncol(obs_tab)
  
  row_objvec_mat <- sapply(1:p1, function(i){
    idx_vec <- .generate_indices_l2(bool_row = T,
                                    idx = i,
                                    p1 = p1, p2 = p2)
    .generate_objvec_l2_helper(idx_vec = idx_vec,
                               num_vec = obs_tab[idx_vec],
                               target_count = prior_1[i]*num_neigh,
                               total_len = p1*p2)
  })
  col_objvec_mat <- sapply(1:p2, function(j){
    idx_vec <- .generate_indices_l2(bool_row = F,
                                    idx = j,
                                    p1 = p1, p2 = p2)
    .generate_objvec_l2_helper(idx_vec = idx_vec,
                               num_vec = obs_tab[idx_vec],
                               target_count = prior_2[j]*num_neigh,
                               total_len = p1*p2)
  })
  matrixStats::rowSums2(row_objvec_mat) + matrixStats::rowSums2(col_objvec_mat)
}

################3

.generate_objmat_l2_helper <- function(idx_vec,
                                       num_vec,
                                       total_len){
  len <- length(idx_vec)
  mat <- matrix(0, nrow = total_len, ncol = total_len)
  mat[idx_vec,idx_vec] <- rep(num_vec, each = len)
  tmp_vec <- rep(0, total_len); tmp_vec[idx_vec] <- num_vec
  mat <- .mult_vec_mat(tmp_vec, mat)
  2*mat
}

.generate_objvec_l2_helper <- function(idx_vec,
                                       num_vec,
                                       target_count,
                                       total_len){
  vec <- rep(0, total_len)
  vec[idx_vec] <- 2*num_vec*target_count
  vec
}

.generate_constr_l2 <- function(obs_tab,  
                                target_count){
  p1 <- nrow(obs_tab); p2 <- ncol(obs_tab)
  total_len <- p1*p2
  mat <- matrix(0, nrow = total_len, ncol = 2*total_len+1)
  vec <- c(target_count, rep(0, total_len), rep(-1, total_len))
  
  # sums to target_count
  mat[,1] <- as.numeric(obs_tab)
  
  # at least 0
  mat[,(1:total_len)+1] <- diag(total_len)
  
  # at most 1
  mat[,(total_len+2):(2*total_len+1)] <- -diag(total_len)
  
  list(mat = mat, vec = vec)
}

###################################


.convert_list2sparse <- function(n, nn_list,
                                 rowname_vec){
  if(length(rowname_vec) > 0) stopifnot(length(rowname_vec) == n)
  
  i_vec <- rep(1:n, times = sapply(nn_list, length))
  j_vec <- unlist(nn_list)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                     j = j_vec,
                                     x = rep(1, length(i_vec)),
                                     dims = c(n,n),
                                     repr = "C")
  sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
  sparse_mat@x <- rep(1, length(sparse_mat@x))
  
  if(length(rowname_vec) > 0){
    rownames(sparse_mat) <- rowname_vec
    colnames(sparse_mat) <- rowname_vec
  }
  sparse_mat
}

