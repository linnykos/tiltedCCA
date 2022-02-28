context("Test tiltedCCA_decomposition")

## tiltedCCA_decomposition is correct

test_that("(Basic) tiltedCCA_decomposition works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
  
  n <- nrow(mat_1)
  expect_true(class(res) == "tiltedCCA_decomp")
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("common_score", 
                                             "cca_obj", 
                                             "distinct_score_1", 
                                             "distinct_score_2",
                                             "common_mat_1", 
                                             "common_mat_2",
                                             "distinct_mat_1", 
                                             "distinct_mat_2", 
                                             "tilt_perc",
                                             "df_percentage",
                                             "svd_1", "svd_2", 
                                             "score_1", "score_2",
                                             "common_basis",
                                             "target_dimred", "param_list"))))
  expect_true(all(dim(res$common_score) == c(n, K)))
})

test_that("(Coding) tiltedCCA_decomposition preserves rownames and colnames", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
  
  expect_true(length(rownames(res$common_score)) > 1)
  expect_true(length(rownames(res$distinct_score_1)) > 1)
  expect_true(length(rownames(res$distinct_score_2)) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_score_2)))
  
  expect_true(length(rownames(res$common_mat_1)) > 1)
  expect_true(length(rownames(res$common_mat_2)) > 1)
  expect_true(length(rownames(res$distinct_mat_1)) > 1)
  expect_true(length(rownames(res$distinct_mat_2)) > 1)
  expect_true(length(colnames(res$common_mat_1)) > 1)
  expect_true(length(colnames(res$common_mat_2)) > 1)
  expect_true(length(colnames(res$distinct_mat_1)) > 1)
  expect_true(length(colnames(res$distinct_mat_2)) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$common_mat_2)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_1)))
  expect_true(all(rownames(mat_1) == rownames(res$distinct_mat_2)))
  expect_true(all(colnames(mat_1) == colnames(res$common_mat_1)))
  expect_true(all(colnames(mat_1) == colnames(res$distinct_mat_1)))
  expect_true(all(colnames(mat_2) == colnames(res$common_mat_2)))
  expect_true(all(colnames(mat_2) == colnames(res$distinct_mat_2)))
})

test_that("(Math) tiltedCCA_decomposition yields uncorrelated distinct matrices", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    target_dimred <- test_data$target_dimred
    K <- 2
    
    tilted_res <- tiltedCCA(mat_1, mat_2, 
                            dims_1 = 1:K, dims_2 = 1:K, 
                            target_dimred = target_dimred,
                            snn_k = 2,
                            snn_min_deg = 1,
                            snn_num_neigh = 10,
                            verbose = F)
    res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
    
    
    tmp <- crossprod(res$distinct_mat_1, res$distinct_mat_2)
    
    sum(abs(tmp)) <= 1e-4
  })

  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition yields a low-rank matrix", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    target_dimred <- test_data$target_dimred
    K <- 2
    
    tilted_res <- tiltedCCA(mat_1, mat_2, 
                            dims_1 = 1:K, dims_2 = 1:K, 
                            target_dimred = target_dimred,
                            snn_k = 2,
                            snn_min_deg = 1,
                            snn_num_neigh = 10,
                            verbose = F)
    res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
    
    bool1 <- Matrix::rankMatrix(res$common_mat_1) == K
    bool2 <- Matrix::rankMatrix(res$common_mat_2) == K
    bool3 <- Matrix::rankMatrix(res$distinct_mat_1) == K
    bool4 <- Matrix::rankMatrix(res$distinct_mat_2) == K
    bool5 <- Matrix::rankMatrix(res$common_mat_1 + res$distinct_mat_1) == K
    bool6 <- Matrix::rankMatrix(res$common_mat_2 + res$distinct_mat_2) == K
    
    bool1 & bool2 & bool3 & bool4 & bool5 & bool6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition yields common matrices with the same column space", {
  # load("tests/assets/test_data1.RData")
  
  bool_vec <- sapply(1:4, function(i){
    load(paste0("../assets/test_data", i, ".RData"))
    mat_1 <- test_data$mat_1
    mat_2 <- test_data$mat_2
    target_dimred <- test_data$target_dimred
    K <- 2
    
    tilted_res <- tiltedCCA(mat_1, mat_2, 
                            dims_1 = 1:K, dims_2 = 1:K, 
                            target_dimred = target_dimred,
                            snn_k = 2,
                            snn_min_deg = 1,
                            snn_num_neigh = 10,
                            verbose = F)
    res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
    
    svd_1 <- svd(res$common_mat_1)$u[,1:K]
    svd_2 <- svd(res$common_mat_2)$u[,1:K]
    
    sum(abs(tcrossprod(svd_1) - tcrossprod(svd_2))) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

test_that("(Math) tiltedCCA_decomposition can obtain the same result when fed into itself (i.e., stability/identifiability)", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  K <- 2
  
  tilted_res <- tiltedCCA(mat_1, mat_2, 
                          dims_1 = 1:K, dims_2 = 1:K, 
                          target_dimred = target_dimred,
                          snn_k = 2,
                          snn_min_deg = 1,
                          snn_num_neigh = 10,
                          verbose = F)
  res <- tiltedCCA_decomposition(tilted_res, rank_c = K, verbose = F)
  
  tilted_res2 <- tiltedCCA(res$common_mat_1 + res$distinct_mat_1, 
                           res$common_mat_2 + res$distinct_mat_2,
                         dims_1 = 1:K, dims_2 = 1:K, 
                         target_dimred = target_dimred,
                         snn_k = 2,
                         snn_min_deg = 1,
                         snn_num_neigh = 10,
                         verbose = F)
  res2 <- tiltedCCA_decomposition(tilted_res2, rank_c = K, verbose = F)
  
  expect_true(sum(abs(res$common_mat_1 - res2$common_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$common_mat_2 - res2$common_mat_2)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_1 - res2$distinct_mat_1)) <= 1e-6)
  expect_true(sum(abs(res$distinct_mat_2 - res2$distinct_mat_2)) <= 1e-6)
})
