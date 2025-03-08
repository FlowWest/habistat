% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vector.R
\name{vector_import_srh2d}
\alias{vector_import_srh2d}
\title{Import tabular SRH-2D results}
\usage{
vector_import_srh2d(
  csv_tbl,
  units = "ft",
  units_xy = units,
  units_z = units,
  has_vid_col = FALSE
)
}
\arguments{
\item{csv_tbl}{A data frame or tibble containing the file path to each CSV containing SRH-2D model results in the standard SRH-2D export format (columns for vertex id \link{optional}, X coordinate, Y coordinate, Z coordinate, WSE, depth, velocity X component, velocity Y component, velocity magnitude, froude number, and shear stress). There should be a numeric column called \code{flow_cfs} and a character column called \code{filename} containing valid file paths to the CSV corresponding to each flow.}

\item{units}{A string indicating units of \code{ft} or \code{m}}

\item{units_xy}{A string indicating horizontal coordinate units of \code{ft} or \code{m}. Defaults to \code{units} if not separately specified.}

\item{units_z}{A string indicating vertical coordinate units of \code{ft} or \code{m}. Defaults to \code{units} if not separately specified.}

\item{has_vid_col}{Boolean indicating whether or not there is already a vertex id (\code{vid}) column in the CSV. If \code{FALSE}, one will be generated based on the row number.}
}
\value{
A \code{tbl_df} data frame containing all SRH-2D model results at all flows with standard column names, one row per vertex per flow.
}
\description{
Import tabular SRH-2D results
}
