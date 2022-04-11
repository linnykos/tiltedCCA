context("Test utility")

## .mult_vec_mat is correct

test_that(".mult_vec_mat works", {
  trials <- 100
  n <- 10
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- round(100*runif(n))
    mat <- matrix(100*runif(n^2), n, n)
    res1 <- .mult_vec_mat(vec, mat)
    res2 <- diag(vec) %*% mat
    
    sum(abs(res1 - res2)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

########################

## .mult_mat_vec is correct

test_that(".mult_mat_vec works", {
  trials <- 100
  n <- 10
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- round(100*runif(n))
    mat <- matrix(100*runif(n^2), n, n)
    res1 <- .mult_mat_vec(mat, vec)
    res2 <- mat %*% diag(vec)
    
    sum(abs(res1 - res2)) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

test_that(".mult_mat_vec can be used to compute inv_onehalf", {
  trials <- 100
  p <- 6
  cov_mat <- matrix(0, p, p)
  cov_mat[1:(p/2), 1:(p/2)] <- 2
  cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
  diag(cov_mat) <- c(rep(5, p/2), rep(1, p/2))
  n <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- MASS::mvrnorm(n = 100, mu = rep(0, p), Sigma = cov_mat)
    mat <- stats::cov(dat)
    eigen_res <- eigen(mat)
    res <- .mult_mat_vec(eigen_res$vectors, 1/sqrt(eigen_res$values))
    
    sum(abs(crossprod(res, mat) %*% res - diag(p))) <= 1e-4
  })
  
  expect_true(all(bool_vec))
})

#############################

# .matching_idx is correct

test_that(".matching_idx is correct", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    
    vec1 <- c(11:30)[sample(1:20)]
    vec2 <- c(11:30)[sample(1:20)]
    
    res <- .matching_idx(vec1, vec2)
    
    abs(vec1[res] - vec2) == 0
  })
  
  expect_true(all(bool_vec))
})

###############################

## .append_rowcolnames is correct

test_that(".append_rowcolnames works for two matrices", {
  mat1 <- matrix(1:30, nrow = 6, ncol = 5)
  mat2 <- matrix(1:30, nrow = 6, ncol = 5)
  rownames(mat2) <- paste0("c", 1:6)
  colnames(mat2) <- paste0("g", 1:5)
  
  res <- .append_rowcolnames(bool_colnames = T,
                             bool_rownames = T,
                             source_obj = mat2,
                             target_obj = mat1)
  
  expect_true(all(rownames(res) == rownames(mat2)) & length(rownames(res)) == 6)
  expect_true(all(colnames(res) == colnames(mat2)) & length(colnames(res)) == 5)
})


test_that(".append_rowcolnames works for svd from matrix", {
  mat1 <- matrix(1:30, nrow = 6, ncol = 5)
  svd1 <- irlba::irlba(mat1, nv = 2)
  class(svd1) <- "svd"
  mat2 <- matrix(1:30, nrow = 6, ncol = 5)
  rownames(mat2) <- paste0("c", 1:6)
  colnames(mat2) <- paste0("g", 1:5)
  
  res <- .append_rowcolnames(bool_colnames = T,
                             bool_rownames = T,
                             source_obj = mat2,
                             target_obj = svd1)
  
  expect_true(all(rownames(res$u) == rownames(mat2)) & length(rownames(res$u)) == 6)
  expect_true(all(rownames(res$v) == colnames(mat2)) & length(rownames(res$v)) == 5)
})


test_that(".append_rowcolnames works for matrix from svd", {
  mat1 <- matrix(1:30, nrow = 6, ncol = 5)
  mat2 <- matrix(1:30, nrow = 6, ncol = 5)
  svd2 <- irlba::irlba(mat1, nv = 2)
  class(svd2) <- "svd"
  rownames(svd2$u) <- paste0("c", 1:6)
  rownames(svd2$v) <- paste0("g", 1:5)
  
  res <- .append_rowcolnames(bool_colnames = T,
                             bool_rownames = T,
                             source_obj = svd2,
                             target_obj = mat1)
  
  expect_true(all(rownames(res) == rownames(svd2$u)) & length(rownames(res)) == 6)
  expect_true(all(colnames(res) == rownames(svd2$v)) & length(colnames(res)) == 5)
})


test_that(".append_rowcolnames works for two svds", {
  mat1 <- matrix(1:30, nrow = 6, ncol = 5)
  svd1 <- irlba::irlba(mat1, nv = 2)
  class(svd1) <- "svd"
  mat2 <- matrix(1:30, nrow = 6, ncol = 5)
  svd2 <- irlba::irlba(mat1, nv = 2)
  class(svd2) <- "svd"
  rownames(svd2$u) <- paste0("c", 1:6)
  rownames(svd2$v) <- paste0("g", 1:5)
  
  res <- .append_rowcolnames(bool_colnames = T,
                             bool_rownames = T,
                             source_obj = svd2,
                             target_obj = svd1)
  
  expect_true(all(rownames(res$u) == rownames(svd2$u)) & length(rownames(res$u)) == 6)
  expect_true(all(rownames(res$v) == rownames(svd2$v)) & length(rownames(res$v)) == 5)
})

##############################

## .combine_two_named_lists is correct

test_that(".combine_two_named_lists works", {
  list1 <- list(a = 1, b = 1:5, c = 3)
  list2 <- list(d = 2, e = 1:10)
  res <- .combine_two_named_lists(list1, list2)
  
  expect_true(all(names(res) == c("a", "b", "c", "d", "e")))
})

test_that(".combine_two_named_lists respects the first named conflict", {
  list1 <- list(a = 1, b = 1:5, c = 3)
  list2 <- list(d = 2, a = 1:10)
  res <- .combine_two_named_lists(list1, list2)
  
  expect_true(all(names(res) == c("a", "b", "c", "d")))
  expect_true(res$a == 1)
})

test_that(".combine_two_named_lists works with there are NULLs", {
  list1 <- list(a=5, b=2, c=NULL)
  list2 <- list(d=3, e=NULL, f=10)
  res <- .combine_two_named_lists(list1, list2)
  
  expect_true(all(names(res) == c("a", "b", "c", "d", "e", "f")))
  expect_true(is.null(list2[["c"]]))
  expect_true(is.null(list2[["e"]]))
})

##########################

## .convert_factor2list is correct

test_that(".convert_factor2list works", {
  vec <- factor(rep(c("a","b","c","d","e"), each = 5))
  res <- .convert_factor2list(vec)
  
  expect_true(all(names(res) == c("a","b","c","d","e")))
  expect_true(is.list(res))
  for(i in 1:length(res)){
    all(sort(which(vec == names(res)[i])) == sort(res[[i]]))
  }
})

test_that(".convert_factor2list works with NAs", {
  vec <- factor(rep(c("a","b","c","d","e"), each = 5))
  vec[sample(length(vec), 10)] <- NA
  vec <- droplevels(vec)
  res <- .convert_factor2list(vec)
  
  expect_true(all(sort(names(res)) == sort(levels(vec))))
  expect_true(is.list(res))
  for(i in 1:length(res)){
    all(sort(which(vec == names(res)[i])) == sort(res[[i]]))
  }
})

#######################

## .convert_list2factor is correct

test_that(".convert_list2factor works", {
  vec <- factor(rep(c("a","b","c","d","e"), each = 5))
  vec[sample(length(vec), 10)] <- NA
  vec <- droplevels(vec)
  lis <- .convert_factor2list(vec)
  res <- .convert_list2factor(lis, n = length(vec))
  
  expect_true(all(vec[!is.na(vec)] == res[!is.na(res)]))
  expect_true(all(is.na(vec) == is.na(res)))
})
