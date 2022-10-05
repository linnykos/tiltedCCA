context("Test graph alignment")

## postprocess_graph_alignment is correct

test_that("postprocess_graph_alignment works", {
  # load("tests/assets/test_data2.RData")
  load("../assets/test_data2.RData")
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
                               min_deg = 10)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  res <- postprocess_graph_alignment(
    input_obj = multiSVD_obj,
    bool_use_denoised = T,
    bool_use_metacells = F,
    input_assay = 1
  )
  expect_true(is.numeric(res))
  expect_true(all(sort(names(res)) == sort(colnames(mat_1))))
  
  res <- postprocess_graph_alignment(
    input_obj = multiSVD_obj,
    bool_use_denoised = T,
    bool_use_metacells = F,
    input_assay = 2
  )
  expect_true(is.numeric(res))
  expect_true(all(sort(names(res)) == sort(colnames(mat_2))))
})

test_that("postprocess_graph_alignment works with singular variables", {
  # load("tests/assets/test_data2.RData")
  load("../assets/test_data2.RData")
  mat_1 <- test_data$mat_1
  mat_1 <- cbind(mat_1, matrix(0, nrow = nrow(mat_1), ncol = 2))
  colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
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
                               min_deg = 10)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  res <- postprocess_graph_alignment(
    input_obj = multiSVD_obj,
    bool_use_denoised = T,
    bool_use_metacells = F,
    input_assay = 1
  )
  expect_true(length(res) == ncol(mat_1))
  expect_true(all(is.na(res[(length(res)-1):length(res)])))
  expect_true(all(names(res) == colnames(mat_1)))
})

###############################

## postprocess_smooth_variable_selection is correct

test_that("postprocess_smooth_variable_selection works", {
  # load("tests/assets/test_data1.RData")
  load("../assets/test_data1.RData")
  mat_1 <- test_data$mat_1
  mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1))), nrow = nrow(mat_1), ncol = ncol(mat_1))
  mat_2 <- test_data$mat_2
  mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2))), nrow = nrow(mat_2), ncol = ncol(mat_2))
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["RNA"]]@var.features <- colnames(mat_1)
  
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
                               min_deg = 10)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  res <- postprocess_smooth_variable_selection(
    input_obj = multiSVD_obj,
    bool_use_denoised = F,
    bool_use_metacells = F,
    num_variables = 5,
    sd_quantile = 0,
    seurat_obj = seurat_obj,
    seurat_assay_1 = "RNA",
    seurat_slot = "counts"
  )
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("alignment_1", "alignment_2", "cor_threshold", "selected_variables", "sd_quantile", "sd_vec_1", "sd_vec_2"))))
  expect_true(all(names(res$alignment) == colnames(mat_1)))
  expect_true(length(res$selected_variables) <= 5)
  expect_true(all(res$selected_variables %in% colnames(mat_1)))
})

test_that("postprocess_smooth_variable_selection works with singular variables", {
  # load("tests/assets/test_data2.RData")
  load("../assets/test_data2.RData")
  mat_1 <- test_data$mat_1
  mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1))), nrow = nrow(mat_1), ncol = ncol(mat_1))
  mat_1 <- cbind(mat_1, matrix(0, nrow = nrow(mat_1), ncol = 2))
  colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
  mat_2 <- test_data$mat_2
  mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2))), nrow = nrow(mat_2), ncol = ncol(mat_2))
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["RNA"]]@var.features <- colnames(mat_1)
  
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
                               min_deg = 10)
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                            verbose = F)
  multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)
  
  res <- postprocess_smooth_variable_selection(
    input_obj = multiSVD_obj,
    bool_use_denoised = F,
    bool_use_metacells = F,
    num_variables = 5,
    sd_quantile = 0,
    seurat_obj = seurat_obj,
    seurat_assay_1 = "RNA",
    seurat_slot = "counts"
  )
  
  expect_true(all(is.na(res$alignment[(length(res$alignment)-1):length(res$alignment)])))
  expect_true(all(names(res$alignment) == colnames(mat_1)))
  na_vars <- names(res$alignment)[is.na(res$alignment)]
  if(length(na_vars) > 0) expect_true(!all(na_vars %in% res$selected_variables))
})

