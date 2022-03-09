context("Test tiltedCCA")

## tiltedCCA is correct

test_that("(Basic) tiltedCCA works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  
  res <- tiltedCCA(input_obj = multiSVD_obj)
  
  expect_true(is.list(res))
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("tcca_obj", "cca_obj") %in% names(res)))
  expect_true(inherits(res$cca_obj, "cca"))
  expect_true(inherits(res$tcca_obj, "tcca"))
  expect_true(all(sort(names(res$cca_obj)) == sort(c("score_1", "score_2", "cca_obj"))))
  expect_true(all(sort(names(res$tcca_obj)) == sort(c("common_basis",
                                                      "common_score",
                                                      "distinct_score_1",
                                                      "distinct_score_2",
                                                      "df_percentage",
                                                      "tilt_perc"))))
  expect_true(length(grep("^tcca*", names(.get_param(res)))) == 3)
  expect_true(all(dim(res$tcca_obj$common_score) == c(n,2)))
})


test_that("(Basic) tiltedCCA works with variable dimensions", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  target_dimred <- test_data$target_dimred
  
  set.seed(10)
  mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1))), nrow = nrow(mat_1), ncol = ncol(mat_1))
  mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2))), nrow = nrow(mat_2), ncol = ncol(mat_2))
  
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:4, dims_2 = 2:3,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  
  res <- tiltedCCA(input_obj = multiSVD_obj)
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("tcca_obj", "cca_obj") %in% names(res)))
  expect_true(all(dim(res$tcca_obj$common_score) == c(n,2)))
  expect_true(all(dim(res$tcca_obj$distinct_score_1) == c(n,4)))
  expect_true(all(dim(res$tcca_obj$distinct_score_2) == c(n,2)))
})

test_that("(Basic) tiltedCCA works with a sparse matrix", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  
  n <- nrow(mat_1)
  mat_1[sample(1:prod(dim(mat_1)), round(prod(dim(mat_1))/2))] <- 0
  mat_1 <- Matrix::Matrix(mat_1, sparse = T)
  mat_2[sample(1:prod(dim(mat_2)), round(prod(dim(mat_1))/2))] <- 0
  mat_2 <- Matrix::Matrix(mat_2, sparse = T)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  
  res <- tiltedCCA(input_obj = multiSVD_obj)
  
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("tcca_obj", "cca_obj") %in% names(res)))
  expect_true(all(dim(res$tcca_obj$common_score) == c(n,2)))
  expect_true(all(dim(res$tcca_obj$distinct_score_1) == c(n,2)))
  expect_true(all(dim(res$tcca_obj$distinct_score_2) == c(n,2)))
})

test_that("(Coding) tiltedCCA preserves rownames and colnames", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  
  res <- tiltedCCA(input_obj = multiSVD_obj)
  
  expect_true(length(res$tcca_obj$common_score) > 1)
  expect_true(length(res$tcca_obj$distinct_score_1) > 1)
  expect_true(length(res$tcca_obj$distinct_score_2) > 1)
  expect_true(length(res$cca_obj$score_1) > 1)
  expect_true(length(res$cca_obj$score_2) > 1)
  expect_true(all(rownames(mat_1) == rownames(res$tcca_obj$common_score)))
  expect_true(all(rownames(mat_1) == rownames(res$tcca_obj$distinct_score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$tcca_obj$distinct_score_2)))
  expect_true(all(rownames(mat_1) == rownames(res$cca_obj$score_1)))
  expect_true(all(rownames(mat_1) == rownames(res$cca_obj$score_2)))
})

test_that("(Math) tiltedCCA is symmetric if the arguments are flipped", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj1 <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = F, center_2 = F,
                                   normalize_row = T,
                                   normalize_singular_value = F,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = F, scale_2 = F)
  multiSVD_obj1 <- form_metacells(input_obj = multiSVD_obj1,
                                  large_clustering_1 = large_clustering_1, 
                                  large_clustering_2 = large_clustering_2,
                                  num_metacells = NULL)
  multiSVD_obj1 <- compute_snns(input_obj = multiSVD_obj1,
                                latent_k = 2,
                                num_neigh = 10,
                                bool_cosine = T,
                                bool_intersect = T,
                                min_deg = 1)
  res1 <- tiltedCCA(input_obj = multiSVD_obj1)
  
  ####
  
  multiSVD_obj2 <- create_multiSVD(mat_1 = mat_2, mat_2 = mat_1,
                                   dims_1 = 1:2, dims_2 = 1:2,
                                   center_1 = F, center_2 = F,
                                   normalize_row = T,
                                   normalize_singular_value = F,
                                   recenter_1 = F, recenter_2 = F,
                                   rescale_1 = F, rescale_2 = F,
                                   scale_1 = F, scale_2 = F)
  multiSVD_obj2 <- form_metacells(input_obj = multiSVD_obj2,
                                  large_clustering_1 = large_clustering_2, 
                                  large_clustering_2 = large_clustering_1,
                                  num_metacells = NULL)
  multiSVD_obj2 <- compute_snns(input_obj = multiSVD_obj2,
                                latent_k = 2,
                                num_neigh = 10,
                                bool_cosine = T,
                                bool_intersect = T,
                                min_deg = 1)
  res2 <- tiltedCCA(input_obj = multiSVD_obj2)
  
  expect_true(abs(res1$tcca_obj$tilt_perc - (1-res2$tcca_obj$tilt_perc)) <= 1e-6)
  
  tmp1 <- res1$tcca_obj$common_score
  tmp2 <- res2$tcca_obj$common_score
  for(j in 1:ncol(tmp1)){
    if(sign(sum(tmp1[1:10,j])) != sign(sum(tmp2[1:10,j]))) {
      tmp1[,j] <- -tmp1[,j]
      res1$tcca_obj$distinct_score_1[,j] <- -res1$tcca_obj$distinct_score_1[,j]
      res1$tcca_obj$distinct_score_2[,j] <- -res1$tcca_obj$distinct_score_2[,j]
    }
  }
  
  expect_true(sum(abs(tmp1 - tmp2)) <= 1e-6)
  expect_true(sum(abs(res1$distinct_score_1 - res2$distinct_score_2)) <= 1e-6)
  expect_true(sum(abs(res1$distinct_score_2 - res2$distinct_score_1)) <= 1e-6)
})

test_that("(Basic) tiltedCCA works with num_metacells", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = F,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = 100)
  multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                               latent_k = 2,
                               num_neigh = 10,
                               bool_cosine = T,
                               bool_intersect = T,
                               min_deg = 1)
  res <- tiltedCCA(input_obj = multiSVD_obj)
  
  expect_true(inherits(res, "multiSVD"))
  expect_true(all(names(multiSVD_obj) %in% names(res)))
  expect_true(all(c("tcca_obj", "cca_obj") %in% names(res)))
  expect_true(all(dim(res$tcca_obj$common_score) == c(n,2)))
  expect_true(all(dim(res$tcca_obj$distinct_score_1) == c(n,2)))
  expect_true(all(dim(res$tcca_obj$distinct_score_2) == c(n,2)))
  expect_true(all(rownames(res$tcca_obj$common_score) == rownames(mat_1)) & length(rownames(res$tcca_obj$common_score)) > 0)
  expect_true(length(rownames(res$tcca_obj$common_basis)) == res$param$mc_num_metacells)
})
