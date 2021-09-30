.determine_cluster <- function(mat,
                               trials = 100){
  n <- nrow(mat)
  
  num_clusters <- sapply(1:trials, function(trial){
    vec <- .random_binary_projection(mat)

    cluster_res <- Ckmeans.1d.dp::Ckmeans.1d.dp(vec)
    max(cluster_res$cluster)
  })
  
  mean(num_clusters)
}

.random_binary_projection <- function(mat){
  n <- nrow(mat); p <- ncol(mat)
  sign_vec <- sample(c(-2,-1,0,1,2), size = p, replace = T)
  if(all(sign_vec == 0)){
    sign_vec <- sample(c(-1,1), size = p, replace = T)
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