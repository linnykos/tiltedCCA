#' Weights for cells and variables
#'
#' @param dcca_decomp Output of \code{dcca_variance_decomposition}
#' @param membership_vec membership vector
#'
#' @return 4 vectors as a list
#' @export
explained_variance <- function(dcca_decomp, membership_vec){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  n <- nrow(dcca_decomp$score_1)
  r1 <- ncol(dcca_decomp$score_1); r2 <- ncol(dcca_decomp$score_2)
  
  modality_1 <- .explained_variance_single(dcca_decomp$score_1, dcca_decomp$common_score, dcca_decomp$distinct_score_1,
                                           dcca_decomp$svd_1, membership_vec)
  modality_2 <- .explained_variance_single(dcca_decomp$score_2, dcca_decomp$common_score, dcca_decomp$distinct_score_2,
                                           dcca_decomp$svd_2, membership_vec)
  
  modality_1 <- t(apply(modality_1, 1, function(x){x/sum(x)}))
  modality_2 <- t(apply(modality_2, 1, function(x){x/sum(x)}))
  
  list(modality_1 = modality_1, modality_2 = modality_2)
}

#####################################

.explained_variance_single <- function(score, common_score, distinct_score, svd_obj, 
                                       membership_vec){
  stopifnot(all(membership_vec > 0), all(membership_vec %% 1 == 0), max(membership_vec) == length(unique(membership_vec)),
            all(table(membership_vec) >= 3))
  K <- max(membership_vec); r <- ncol(score); r1 <- ncol(common_score)
  
  weights <- sapply(1:ncol(score), function(j){
    abs(score[,j] %*% svd_obj$u) %*% svd_obj$d/sum(svd_obj$d)
  })
  
  score_var <- apply(score, 2, stats::sd)
  distinct_var <- apply(distinct_score, 2, stats::sd)
  common_var <- apply(common_score, 2, stats::sd)
  if(length(common_var) <= length(score_var)) common_var <- c(common_var, rep(0,length(score_var) - length(common_var)))
  
  cell_mat <- t(sapply(1:K, function(k){
    idx <- which(membership_vec == k)
    
    common_cell <- sum(sapply(1:r, function(j){
      if(j <= r){(1 - min(stats::sd(common_score[idx,j])/common_var[j], 1)) * common_var[j]/score_var[j] * weights[j]} else 0
    }))
    
    distinct_cell <- sum(sapply(1:r, function(j){
      (1 - min(stats::sd(distinct_score[idx,j])/distinct_var[j], 1)) * distinct_var[j]/score_var[j] * weights[j]
    }))
    
    c(common_cell, distinct_cell)
  }))
  colnames(cell_mat) <- c("common", "distinct")
  rownames(cell_mat) <- 1:K
  
  cell_mat
}