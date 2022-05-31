context("Test variable selection")

## .generic_variable_selection is correct

test_that(".generic_variable_selection works", {
  set.seed(10)
  n <- 100; p <- 50
  mat <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
  vec <- runif(p)
  colnames(mat) <- paste0("g_", 1:p)
  rownames(mat) <- paste0("c_", 1:n)
  names(vec) <- colnames(mat)
  
  res <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.8,
    initial_mat = NULL,
    mat = mat,
    num_variables = 10,
    return_candidate_list = T,
    vec = vec
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("candidate_list", "variables"))))
  expect_true(names(vec)[which.max(vec)] == res$variables[1])
  expect_true(length(res$variables) <= 10)
  expect_true(length(res$candidate_list) <= length(res$variables))
  for(i in 2:length(res$candidate_list)){
    expect_true(all(res$candidate_list[[i]] %in% res$candidate_list[[i-1]]))
  }
})

test_that(".generic_variable_selection works with an initial mat", {
  set.seed(10)
  n <- 100; p <- 50
  mat <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
  vec <- c(1,0.99,0.999, runif(p-3))
  colnames(mat) <- paste0("g_", 1:p)
  rownames(mat) <- paste0("c_", 1:n)
  names(vec) <- colnames(mat)
  initial_mat <- mat[,1:3] + MASS::mvrnorm(n = n, mu = rep(0, 3), Sigma = 0.01*diag(3))
  
  res1 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.8,
    initial_mat = NULL,
    mat = mat,
    num_variables = 10,
    return_candidate_list = T,
    vec = vec
  )
  
  res2 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.8,
    initial_mat = initial_mat,
    mat = mat,
    num_variables = 10,
    return_candidate_list = T,
    vec = vec
  )
  
  expect_true(is.list(res2))
  expect_true(all(sort(names(res2)) == sort(c("candidate_list", "variables"))))
  expect_true(names(vec)[which.max(vec)] == res1$variables[1])
  expect_true(!any(names(vec)[1:3] %in% res2$variables))
})

test_that(".generic_variable_selection selects uncorrelated variables", {
  set.seed(10)
  n <- 100; p <- 30
  block_mat <- matrix(0, p, p)
  for(i in 1:3){
    block_mat[((i-1)*10+1):(i*10), ((i-1)*10+1):(i*10)] <- 0.9
  }
  diag(block_mat) <- 1
  
  mat <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = block_mat)
  vec <- seq(0, 1, length.out = p)
  colnames(mat) <- paste0("g_", 1:p)
  rownames(mat) <- paste0("c_", 1:n)
  names(vec) <- colnames(mat)
  
  res1 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.6,
    initial_mat = NULL,
    mat = mat,
    num_variables = 3,
    return_candidate_list = T,
    vec = vec
  )
  
  res2 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.6,
    initial_mat = NULL,
    mat = mat,
    num_variables = 4,
    return_candidate_list = T,
    vec = vec
  )
  
  expect_true(all(sort(res1$variables) == sort(res2$variables)))
})