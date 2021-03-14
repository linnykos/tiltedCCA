#' Weights for cells and variables
#'
#' @param dcca_decomp Output of \code{dcca_variance_decomposition}
#' @param membership_vec integer vector
#'
#' @return 4 vectors as a list
#' @export
explained_variance <- function(dcca_decomp, membership_vec, num_pairs = 200){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  n <- nrow(dcca_decomp$score_1)
  r1 <- ncol(dcca_decomp$score_1); r2 <- ncol(dcca_decomp$score_2)
  
  idx_pair_list <- .construct_cell_pairs(membership_vec, num_pairs = num_pairs)
  modal_1 <- .explained_variance_single(dcca_decomp$score_1, dcca_decomp$common_score, dcca_decomp$distinct_score_1,
                                        dcca_decomp$svd_1, membership_vec, idx_pair_list)
  modal_2 <- .explained_variance_single(dcca_decomp$score_2, dcca_decomp$common_score, dcca_decomp$distinct_score_2,
                                        dcca_decomp$svd_2, membership_vec, idx_pair_list)
  
  list(modal_1 = modal_1, modal_2 = modal_2)
}

#####################################

.explained_variance_single <- function(score, common_score, distinct_score, svd_obj, 
                                       membership_vec, idx_pair_list = NA){
  idx <- which(!is.na(membership_vec))
  stopifnot(all(membership_vec[idx] > 0), all(membership_vec[idx] %% 1 == 0), max(membership_vec[idx]) == length(unique(membership_vec[idx])),
            all(table(membership_vec[idx]) >= 3))
  K <- max(membership_vec)
  r <- ncol(score); r1 <- ncol(common_score); n <- nrow(score)
  
  # compute the relative weights for each modality
  embedding_ev <- .construct_canonical_embedding(score, svd_obj, weight_mat = NA)
  tmp <- .project_mat2mat_colwise(common_score, score, orthogonal = F)
  embedding_co <-.construct_canonical_embedding(score, svd_obj, weight_mat = tmp)
  tmp <- .project_mat2mat_colwise(distinct_score, score, orthogonal = F)
  embedding_di <-.construct_canonical_embedding(score, svd_obj, weight_mat = tmp)
 
  weight_ev <- .l2norm(embedding_ev)^2/(n^(3/2))
  weight_co <- .l2norm(embedding_co)^2/(n^(3/2))
  weight_di <- .l2norm(embedding_di)^2/(n^(3/2))
  
  # construct the 3 nearest neighbor graphs
  g_ev <- .construct_nn_graph(embedding_ev)
  g_co <- .construct_nn_graph(embedding_co)
  g_di <- .construct_nn_graph(embedding_di)
  
  if(all(is.na(idx_pair_list))) idx_pair_list <- .construct_cell_pairs(membership_vec)
  stopifnot(length(idx_pair_list) == K)
  
  value_mat <- t(sapply(1:K, function(k){
    purity_ev <- mean(sapply(1:nrow(idx_pair_list[[k]]), function(i){
      .compute_pairwise_purity(g_ev, idx_pair_list[[k]][i,1],
                               idx_pair_list[[k]][i,2], membership_vec)
    }))
    
    purity_co <- mean(sapply(1:nrow(idx_pair_list[[k]]), function(i){
      .compute_pairwise_purity(g_co, idx_pair_list[[k]][i,1],
                               idx_pair_list[[k]][i,2], membership_vec)
    }))
    
    purity_di <- mean(sapply(1:nrow(idx_pair_list[[k]]), function(i){
      .compute_pairwise_purity(g_di, idx_pair_list[[k]][i,1],
                               idx_pair_list[[k]][i,2], membership_vec)
    }))
    
    c(purity_ev, weight_ev, purity_co, weight_co, purity_di, weight_di)
  }))
  
  rownames(value_mat) <- 1:K
  # colnames(value_mat) <- c("everything", "common", "distinct")
  
  value_mat
}

#################################

.construct_canonical_embedding <- function(score, svd_obj, weight_mat = NA){
  if(all(is.na(weight_mat))){
    score %*% .mult_mat_vec(crossprod(score, svd_obj$u), svd_obj$d)
  } else {
    weight_mat %*% .mult_mat_vec(crossprod(score, svd_obj$u), svd_obj$d)
  }
}

.construct_cell_pairs <- function(membership_vec, num_pairs = 200){
  K <- max(membership_vec, na.rm = T)
  
  lapply(1:K, function(k){
    idx <- which(membership_vec == k)
    len <- length(idx)
    if(num_pairs <= len*(len-1)/2){
      tmp <- cbind(sample(idx, size = 200, replace = T), sample(idx, size = 200, replace = T))
      tmp[-which(tmp[,1] == tmp[,2]),]
    } else {
      tmp <- t(utils::combn(len, 2))
      cbind(idx[tmp[1,]], idx[tmp[2,]])
    }
  })
}

.construct_nn_graph <- function(mat, min_neigh = 10, max_neigh = 25){
  n <- nrow(mat)
  res <- RANN::nn2(mat, k = max_neigh+1)$nn.idx[,-1]
  # -1 since RANN::nn2 counts itself as a nearest neighbor
  
  k <- min_neigh 
  prev_g <- NA; prev_num_con <- NA; edge_mat <- NA
  
  while(TRUE & k <= max_neigh){
    if(all(is.na(edge_mat))){
      edge_mat <- do.call(cbind, lapply(1:n, function(i){
        rbind(rep(i, k), res[i,1:k])
      }))
      
      g <- igraph::graph.empty(n = n, directed = F)
      g <- igraph::add_edges(g, edges = edge_mat)
      g <- igraph::simplify(g)
      
    } else {
      edge_mat <- rbind(1:n, res[,k])
      g <- igraph::add_edges(prev_g, edges = edge_mat)
      g <- igraph::simplify(g)
    }
   
    num_con <- igraph::components(g)$no
    if(!is.na(prev_num_con) && prev_num_con == num_con){
      break()
    } else {
      prev_num_con <- num_con
      prev_g <- g
    }
    
    k <- k+1
  }
 
  prev_g
}

#' Compute pairwise purity
#'
#' For two indices \code{idx1} and \code{idx2} (two nodes in the \code{igraph}
#' object \code{g}, so both \code{idx1} and \code{idx2} need to be less than
#' \code{igraph::vcount(g)}), compute the shortest path from node \code{idx1}
#' to node \code{idx2} and determine what percentage of nodes in that
#' path also share the same cluster label as node \code{idx1}
#' and node \code{idx2}.
#'
#' @param g \code{igraph} object
#' @param idx1 index between 1 and \code{igraph::vcount(g)}
#' @param idx2 index between 1 and \code{igraph::vcount(g)} not equal to \code{idx1}
#' @param membership_vec integer vector
#'
#' @return numeric
.compute_pairwise_purity <- function(g, idx1, idx2, membership_vec){
  stopifnot(membership_vec[idx1] == membership_vec[idx2])
  k <- membership_vec[idx1]
  path <- igraph::shortest_paths(g, from = idx1, to = idx2)$vpath
  
  if(length(path) == 1 && all(sapply(path, length) == 0)) return(0)
    
  len <- length(path)
  val_vec <- sapply(1:len, function(i){
    vec <- path[[i]]
    if(length(vec) == 2) return(1)
    
    length(which(membership_vec[vec[-c(1,length(vec))]] == k))/(length(vec)-2)
  })
  
  max(val_vec)
}