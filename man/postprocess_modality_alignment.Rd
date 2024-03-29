% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocess_modality_alignment.R
\name{postprocess_modality_alignment}
\alias{postprocess_modality_alignment}
\title{Compute how aligned the common component for a feature is with that feature in the original data modality}
\usage{
postprocess_modality_alignment(
  input_obj,
  bool_use_denoised,
  input_assay,
  bool_center = T,
  bool_scale = T,
  bool_regression_include_intercept = T,
  min_subsample_cell = NULL,
  seurat_celltype_variable = "celltype",
  seurat_obj = NULL,
  seurat_assay = NULL,
  seurat_slot = "data",
  verbose = 1
)
}
\arguments{
\item{input_obj}{a \code{multiSVD_obj} that was the output of \code{tiltedCCA::tiltedCCA_decomposition}}

\item{bool_use_denoised}{boolean. If \code{TRUE}, then the common component is compared against
the common plus distinct component. If \code{FALSE}, then the common component is compared against
the original data matrix in slot \code{seurat_slot} in \code{seurat_obj[[seurat_assay]]}}

\item{input_assay}{integer of \code{1} or \code{2}, denoting which modality is being analyzed}

\item{bool_center}{boolean if all the features in the common component are centered prior to the comparison}

\item{bool_scale}{boolean if all the features in the common component are rescaled prior to the comparison}

\item{bool_regression_include_intercept}{boolean if the regression analysis}

\item{min_subsample_cell}{if not \code{NULL}, subsample \code{min_subsample_cell} number of cells of each cell type
(denoted by in \code{seurat_obj$seurat_celltype_variable})}

\item{seurat_celltype_variable}{a character where \code{seurat_obj$seurat_celltype_variable} 
denotes the cell type for each cell in \code{seurat_obj}}

\item{seurat_obj}{the \code{Seurat} object that was used to compute \code{input_obj}, the \code{multiSVD_obj}}

\item{seurat_assay}{the assay to extract the data matrix, which is relevant \code{bool_use_denoised=FALSE}}

\item{seurat_slot}{the slot to extract the data matrix, which is relevant \code{bool_use_denoised=FALSE}}

\item{verbose}{non-negative integer}
}
\value{
a vector of R-squared values for each variable
}
\description{
The common matrix for modality \code{input_assay} is extracted, and then
compared against another matrix (depending on \code{bool_use_denoised} is set).
(For example, if \code{bool_use_denoised=TRUE}, then this is the common plus distinct
matrix.) Then, a regression is performed, one per feature (i.e., gene or protein)
that regresses latter matrix onto the common matrix, and the R-squared
(one per feature) is returned.
}
