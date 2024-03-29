% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tiltedCCA_common_score.R
\name{.tiltedCCA_common_score}
\alias{.tiltedCCA_common_score}
\title{Main workhorse of dcca_factor}
\usage{
.tiltedCCA_common_score(
  averaging_mat,
  cca_res,
  discretization_gridsize,
  enforce_boundary,
  fix_tilt_perc,
  snn_bool_cosine,
  snn_bool_intersect,
  snn_k,
  snn_min_deg,
  snn_num_neigh,
  svd_1,
  svd_2,
  target_dimred,
  verbose = 0
)
}
\arguments{
\item{averaging_mat}{sparse matrix}

\item{cca_res}{returned object from \code{.cca}}

\item{discretization_gridsize}{positive integer for how many values between 0 and 1 (inclusive) to search the 
appropriate amount of tilt over}

\item{enforce_boundary}{boolean, on whether or not the tilt is required to stay between
the two canonical score vectors}

\item{fix_tilt_perc}{boolean or a numeric. If \code{FALSE}, then the tilt is adaptively
determined, and if \code{TRUE}, then the tilt is set to be equal to 
\code{0.5}. If numeric, the value should be between \code{0} and \code{1},
which the tilt will be set to.}

\item{snn_bool_cosine}{boolean}

\item{snn_bool_intersect}{boolean}

\item{snn_k}{integer}

\item{snn_min_deg}{integer}

\item{snn_num_neigh}{integer}

\item{svd_1}{SVD of the denoised variant of \code{mat_1} from \code{dcca_factor}}

\item{svd_2}{SVD of the denoised variant of \code{mat_2} from \code{dcca_factor}}

\item{target_dimred}{matrix}

\item{verbose}{non-negative integer}
}
\value{
list
}
\description{
Given the two matrices (given by \code{svd_1} and \code{svd_2}) and the
CCA solution in \code{cca_res}, compute the common scores.
This calls the functions
\code{.common_decomposition} and \code{.compute_distinct_score}.
}
