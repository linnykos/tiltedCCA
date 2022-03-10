context("Test Seurat helpers")

## create_SeuratDim is correct

test_that("create_SeuratDim works", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = T,
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
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
  multiSVD_obj <- tiltedCCA_decomposition(input_obj = multiSVD_obj)
  
  res <- create_SeuratDim(input_obj = multiSVD_obj,
                          what = "common",
                          aligned_umap_assay = NULL,
                          seurat_obj = NULL)
  expect_true(inherits(res, "DimReduc"))
  expect_true(all(dim(res@cell.embeddings) == c(300,2)))
  expect_true(length(rownames(res@cell.embeddings)) > 0 & all(rownames(res@cell.embeddings) == rownames(mat_1)))
  
  res <- create_SeuratDim(input_obj = multiSVD_obj,
                          what = "distinct_1",
                          aligned_umap_assay = NULL,
                          seurat_obj = NULL)
  expect_true(inherits(res, "DimReduc"))
  expect_true(all(dim(res@cell.embeddings) == c(300,2)))
  expect_true(length(rownames(res@cell.embeddings)) > 0 & all(rownames(res@cell.embeddings) == rownames(mat_1)))
  
  res <- create_SeuratDim(input_obj = multiSVD_obj,
                          what = "distinct_2",
                          aligned_umap_assay = NULL,
                          seurat_obj = NULL)
  expect_true(inherits(res, "DimReduc"))
  expect_true(all(dim(res@cell.embeddings) == c(300,2)))
  expect_true(length(rownames(res@cell.embeddings)) > 0 & all(rownames(res@cell.embeddings) == rownames(mat_1)))
})

test_that("create_SeuratDim works with Seurat objects", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = T,
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
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
  multiSVD_obj <- tiltedCCA_decomposition(input_obj = multiSVD_obj)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["umap"]] <- Seurat::RunUMAP(.get_Dimred(multiSVD_obj),
                                          assay = "RNA",
                                          verbose = F)
  
  res <- create_SeuratDim(input_obj = multiSVD_obj,
                          what = "common",
                          aligned_umap_assay = "umap",
                          seurat_obj = seurat_obj)
  expect_true(inherits(res, "DimReduc"))
  expect_true(all(dim(res@cell.embeddings) == c(300,2)))
  expect_true(length(rownames(res@cell.embeddings)) > 0 & all(rownames(res@cell.embeddings) == rownames(mat_1)))
})

######################

## .translate_celltype is correct

test_that(".translate_celltype works", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = T,
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
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj$celltype <- large_clustering_1
  
  metacell_list <- .get_metacell(input_obj = multiSVD_obj, 
                                 resolution = "cell",
                                 type = "list",
                                 what = "metacell_clustering")
  expect_true(all(names(metacell_list) == rownames(.get_SNN(multiSVD_obj,
                                                            bool_common = F))))
  
  res <- .translate_celltype(input_obj = multiSVD_obj,
                             celltype_vec = seurat_obj$celltype,
                             metacell_list = metacell_list)
  expect_true(length(res) == length(metacell_list))
  expect_true(is.factor(res))
  expect_true(all(sort(levels(res)) == sort(levels(large_clustering_1))))
})

test_that(".translate_celltype works with NULL metacell_list", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
  mat_1 <- test_data$mat_1
  mat_2 <- test_data$mat_2
  n <- nrow(mat_1)
  large_clustering_1 <- test_data$clustering_1
  large_clustering_2 <- test_data$clustering_2
  multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                  dims_1 = 1:2, dims_2 = 1:2,
                                  center_1 = F, center_2 = F,
                                  normalize_row = T,
                                  normalize_singular_value = T,
                                  recenter_1 = F, recenter_2 = F,
                                  rescale_1 = F, rescale_2 = F,
                                  scale_1 = F, scale_2 = F)
  multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                                 large_clustering_1 = large_clustering_1, 
                                 large_clustering_2 = large_clustering_2,
                                 num_metacells = NULL)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj$celltype <- large_clustering_1
  
  metacell_list <- .get_metacell(input_obj = multiSVD_obj, 
                                 resolution = "cell",
                                 type = "list",
                                 what = "metacell_clustering")
  
  res <- .translate_celltype(input_obj = multiSVD_obj,
                             celltype_vec = seurat_obj$celltype,
                             metacell_list = metacell_list)
  
  expect_true(all(res == seurat_obj$celltype))
  expect_true(is.factor(res))
})

##########################

## create_reducedSeuratObj is correct


test_that("create_reducedSeuratObj works", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
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
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
  multiSVD_obj <- tiltedCCA_decomposition(input_obj = multiSVD_obj)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["umap"]] <- Seurat::RunUMAP(.get_Dimred(multiSVD_obj),
                                          assay = "RNA",
                                          verbose = F)
  seurat_obj$celltype <- large_clustering_1
  
  res <- create_reducedSeuratObj(input_obj = multiSVD_obj,
                                 what = "laplacian_1",
                                 aligned_umap_assay = "umap",
                                 seurat_celltype = "celltype",
                                 seurat_obj = seurat_obj)
  
  expect_true(inherits(res, "Seurat"))
  expect_true(inherits(res[["umap"]], "DimReduc"))
  expect_true(all(dim(res[["umap"]]@cell.embeddings) == c(100,2)))
})


test_that("create_reducedSeuratObj works with no metacells", {
  # load("tests/assets/test_data3.RData")
  load("../assets/test_data3.RData")
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
  multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
  multiSVD_obj <- tiltedCCA_decomposition(input_obj = multiSVD_obj)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["umap"]] <- Seurat::RunUMAP(.get_Dimred(multiSVD_obj),
                                          assay = "RNA",
                                          verbose = F)
  seurat_obj$celltype <- large_clustering_1
  
  res <- create_reducedSeuratObj(input_obj = multiSVD_obj,
                                 what = "laplacian_1",
                                 aligned_umap_assay = "umap",
                                 seurat_celltype = "celltype",
                                 seurat_obj = seurat_obj)
  
  expect_true(inherits(res, "Seurat"))
  expect_true(inherits(res[["umap"]], "DimReduc"))
  expect_true(all(dim(res[["umap"]]@cell.embeddings) == c(300,2)))
})
