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
    mat_1 = mat,
    mat_2 = NULL,
    num_variables = 10,
    return_candidate_list = T,
    vec_1 = vec,
    vec_2 = NULL
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("candidate_list", "selected_variables", "cor_vec"))))
  expect_true(names(vec)[which.max(vec)] == res$selected_variables[1])
  expect_true(length(res$variables) <= 10)
  expect_true(length(res$candidate_list) <= length(res$selected_variables))
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
    mat_1 = mat,
    mat_2 = NULL,
    num_variables = 10,
    return_candidate_list = T,
    vec_1 = vec,
    vec_2 = NULL
  )
  
  res2 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.8,
    initial_mat = initial_mat,
    mat_1 = mat,
    mat_2 = NULL,
    num_variables = 10,
    return_candidate_list = T,
    vec_1 = vec,
    vec_2 = NULL
  )
  
  expect_true(is.list(res2))
  expect_true(all(sort(names(res2)) == sort(c("candidate_list", "selected_variables", "cor_vec"))))
  expect_true(length(res2$cor_vec_1) == p)
  expect_true(min(res2$cor_vec_1[1:3]) > max(res2$cor_vec_1[-c(1:3)]))
  expect_true(names(vec)[which.max(vec)] == res1$selected_variables[1])
  expect_true(all(names(vec)[1:3] %in% res2$selected_variables))
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
    mat_1 = mat,
    mat_2 = NULL,
    num_variables = 3,
    return_candidate_list = T,
    vec_1 = vec,
    vec_2 = NULL,
    verbose = 0
  )
  
  res2 <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.6,
    initial_mat = NULL,
    mat_1 = mat,
    mat_2 = NULL,
    num_variables = 4,
    return_candidate_list = T,
    vec_1 = vec,
    vec_2 = NULL
  )
  
  selected_idx <- sapply(res1$selected_variables, function(x){
    as.numeric(strsplit(x, split = "_")[[1]][2])
  })
  expect_true(all(table(floor(selected_idx/10)) == 1))
  expect_true(all(res1$selected_variables %in% res2$selected_variables))
})


test_that(".generic_variable_selection works with two modalities", {
  set.seed(10)
  n <- 100; p <- 30
  block_mat <- matrix(0, p, p)
  for(i in 1:3){
    block_mat[((i-1)*10+1):(i*10), ((i-1)*10+1):(i*10)] <- 0.9
  }
  diag(block_mat) <- 1
  mat_2 <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = block_mat)
  mat_1 <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
  
  vec_1 <- seq(0, 1, length.out = p)
  vec_2 <- seq(0, 1, length.out = p) + stats::runif(p,max=0.1)
  
  colnames(mat_1) <- paste0("g_", 1:p)
  rownames(mat_1) <- paste0("c_", 1:n)
  names(vec_1) <- colnames(mat_1)
  colnames(mat_2) <- colnames(mat_1)
  rownames(mat_2) <- rownames(mat_1)
  names(vec_2) <- names(vec_1)
  
  res <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.6,
    initial_mat = NULL,
    mat_1 = mat_1,
    mat_2 = mat_2,
    num_variables = 3,
    return_candidate_list = T,
    vec_1 = vec_1,
    vec_2 = vec_2,
    verbose = 0
  )
  
  selected_idx <- sapply(res$selected_variables, function(x){
    as.numeric(strsplit(x, split = "_")[[1]][2])
  })
  expect_true(all(table(floor(selected_idx/10)) == 1))
})


test_that(".generic_variable_selection works with two modalities and type = 'cor'", {
  set.seed(10)
  n <- 100; p <- 30
  block_mat <- matrix(0, p, p)
  for(i in 1:3){
    block_mat[((i-1)*10+1):(i*10), ((i-1)*10+1):(i*10)] <- 0.9
  }
  diag(block_mat) <- 1
  mat_2 <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = block_mat)
  mat_1 <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
  
  vec_1 <- seq(0, 1, length.out = p)
  vec_2 <- seq(0, 1, length.out = p) + stats::runif(p,max=0.1)
  
  colnames(mat_1) <- paste0("g_", 1:p)
  rownames(mat_1) <- paste0("c_", 1:n)
  names(vec_1) <- colnames(mat_1)
  colnames(mat_2) <- colnames(mat_1)
  rownames(mat_2) <- rownames(mat_1)
  names(vec_2) <- names(vec_1)
  
  res <- .generic_variable_selection(
    bool_maximizing = T, 
    cor_threshold = 0.6,
    initial_mat = NULL,
    mat_1 = mat_1,
    mat_2 = mat_2,
    num_variables = 3,
    prediction_type = "cor",
    return_candidate_list = T,
    vec_1 = vec_1,
    vec_2 = vec_2,
    verbose = 0
  )
  
  selected_idx <- sapply(res$selected_variables, function(x){
    as.numeric(strsplit(x, split = "_")[[1]][2])
  })
  expect_true(all(table(floor(selected_idx/10)) == 1))
})

