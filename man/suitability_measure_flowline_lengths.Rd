% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/suitability.R
\name{suitability_measure_flowline_lengths}
\alias{suitability_measure_flowline_lengths}
\title{Flowline length helper function: Extract NHD flowlines that correspond to groups and calculate length}
\usage{
suitability_measure_flowline_lengths(group_polygons, .group_var = comid)
}
\arguments{
\item{group_polygons}{An \code{sf} polygon feature collection with one polygon for each analysis reach. These may have been delineated previously by splitting the hydraulic model domain into subsections.}

\item{.group_var}{The unquoted name of the attribute used to identify groups in the \code{group_polygons} layer. Defaults to \code{comid}}
}
\value{
A \code{tbl_df} data frame, one row per group (\code{comid}), containing the calculated length of the associated NHDPlusV2 flowline (joined by \code{comid}) clipped to the group polygon extent.
}
\description{
Flowline length helper function: Extract NHD flowlines that correspond to groups and calculate length
}
