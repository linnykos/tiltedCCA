#' Title
#'
#' @param mat matrix
#' @param trials numeric
#'
#' @return numeric
#' @import mclust
.determine_cluster <- function(mat,
                               trials = 20){
  n <- nrow(mat)
  
  num_clusters <- sapply(1:trials, function(trial){
    vec <- .random_binary_projection(mat)
    
    cluster_res <- mclust::Mclust(vec, 
                                  G = 1:9,
                                  modelNames = "V",
                                  verbose = F)
    cluster_res$G
  })
  
  mean(num_clusters)
}

.random_binary_projection <- function(mat){
  n <- nrow(mat); p <- ncol(mat)
  sign_vec <- sample(c(-2,-1,0,1,2), size = p, replace = T,
                     prob = c(1/6, 1/6, 1/3, 1/6, 1/6))
  if(all(sign_vec == 0)){
    sign_vec <- sample(c(-2,-1,1,2), size = p, replace = T)
  }
  
  idx <- which(sign_vec == 2)
  if(length(idx) > 0){
    pos_vec2 <- 2*matrixStats::rowSums2(mat[,idx,drop = F])
  } else pos_vec2 <- rep(0, n)
  
  idx <- which(sign_vec == 1)
  if(length(idx) > 0){
    pos_vec <- matrixStats::rowSums2(mat[,idx,drop = F])
  } else pos_vec <- rep(0, n)
  
  idx <- which(sign_vec == -1)
  if(length(idx) > 0){
    neg_vec <- matrixStats::rowSums2(mat[,idx,drop = F])
  } else neg_vec <- rep(0, n)
  
  idx <- which(sign_vec == -2)
  if(length(idx) > 0){
    neg_vec2 <- 2*matrixStats::rowSums2(mat[,idx,drop = F])
  } else neg_vec2 <- rep(0, n)
  
  pos_vec2 + pos_vec - neg_vec - neg_vec2
}