plot_embeddings2 <- function(c_g, d_g, e_g, only_embedding = F,
                             verbose = T, ...){
  list_g <- list(c_g = c_g, d_g = d_g, e_g = e_g)
  list_output <- vector("list", 3)

  for(i in 1:3){
    # symmetrize
    list_g[[i]] <- .symmetrize_sparse(list_g[[i]], set_ones = F)
  
    # convert into kernels
    list_g[[i]] <- .distance_to_kernel(list_g[[i]])
    list_g[[i]] <- SeuratObject::as.Graph(list_g[[i]])
    
    # use Seurat
    tmp <- Seurat::RunUMAP(list_g[[i]], verbose = verbose, ...)
    list_output[[i]] <- tmp@cell.embeddings
    
    # # pass into python
    # # see https://rstudio.github.io/reticulate/articles/calling_python.html
    # csc_matrix <- reticulate::r_to_py(list_g[[i]])
    # if (!reticulate::py_module_available(module = 'umap')) {
    #   stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
    # }
    # umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
    # umap <- umap_import$UMAP(
    #   graph = csc_matrix,
    #   n_components = as.integer(x = 2),
    #   initial_alpha = 1,
    #   a = NULL,
    #   b = NULL,
    #   gamma = 1,
    #   negative_sample_rate = 5,
    #   n_epochs = NULL,
    #   init = "spectral",
    #   random_state = NULL,
    #   metric = "cosine",
    #   densmap = FALSE,
    #   densmap_kwds = NULL,
    #   output_dens = FALSE,
    #   output_metric = "euclidean",
    #   output_metric_kwds = NULL,
    #   euclidean_output = TRUE,
    #   parallel = FALSE,
    #   verbose = TRUE)
    # umap$simplicial_set_embedding(as.matrix(x = object))
  }
  
  if(only_embedding) return(list_output)
}

######################

.distance_to_kernel <- function(mat_g){
  stopifnot(inherits(mat_g, "dgCMatrix"))
  
  n <- nrow(mat_g)
  x_val <- unlist(lapply(1:n, function(col_idx){
    val1 <- mat_g@p[col_idx]+1
    val2 <- mat_g@p[col_idx+1]+1
    
    if(val1 == val2) return(numeric(0))
    vec <- mat_g@x[(val1+1):val2]
    min_val <- min(vec)
    max_val <- max(vec)
    exp(-(vec - min_val)/max_val)
  }))
  
  mat_g@x <- x_val
  mat_g
}