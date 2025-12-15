# Rotate embedding in a Seurat object

This method rotates the embedding in `target_embedding` to look
(visually) most similar to `source_embedding`

## Usage

``` r
rotate_seurat_embeddings(seurat_obj, source_embedding, target_embedding)
```

## Arguments

- seurat_obj:

  object of class `Seurat`

- source_embedding:

  character vector, so that `seurat_obj[[source_embedding]]` is a object
  of class `DimReduc`

- target_embedding:

  character vector, so that `seurat_obj[[target_embedding]]` is a object
  of class `DimReduc`

## Value

an updated `seurat_obj`
