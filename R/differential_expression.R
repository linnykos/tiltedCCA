differential_expression <- function(seurat_obj,
                                    assay,
                                    idents,
                                    slot = "data",
                                    test_use = "wilcox",
                                    verbose = T){
  stopifnot(assay %in% names(seurat_obj),
            idents %in% colnames(seurat_obj@meta.data),
            test_use %in% c("MAST", "wilcox"),
            slot %in% c("counts", "data", "scale.data"),
            length(seurat_obj[[assay]]@var.features) > 0)
  
  Seurat::DefaultAssay(seurat_obj) <- assay
  Seurat::Idents(seurat_obj) <- idents
  
  level_vec <- sort(levels(Seurat::Idents(seurat_obj)))
  num_celltype <- length(level_vec)
  combn_mat <- utils::combn(num_celltype, 2)
  variable_vec <- seurat_obj[[assay]]@var.features
  
  de_list <- lapply(1:ncol(combn_mat), function(i){
    if(verbose) print(paste0("Applying combination ", i, " of ", ncol(combn_mat)))
    ident_1 <- level_vec[combn_mat[1,i]]
    ident_2 <- level_vec[combn_mat[2,i]]
    
    Seurat::FindMarkers(seurat_obj,
                        features = variable_vec,
                        ident.1 = ident_1,
                        ident.2 = ident_2,
                        test.use = test_use,
                        slot = slot,
                        min.pct = 0,
                        logfc.threshold = 0,
                        only.pos = F,
                        verbose = verbose)
  })
  
  list(combn_mat = combn_mat,
       de_list = de_list,
       level_vec = level_vec)
}