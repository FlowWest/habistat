% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duration.R
\name{flow_threshold_remove_baseflow_apply}
\alias{flow_threshold_remove_baseflow_apply}
\title{Remove Baseflow via Flow Threshold (Multiple)}
\usage{
flow_threshold_remove_baseflow_apply(
  data,
  .flow_threshold_var,
  .flow_var = flow_cfs,
  .habitat_var = wua_per_lf
)
}
\arguments{
\item{data}{A \code{data.frame} or \code{tbl_df} containing flow-to-suitable-area curves for multiple groups. All variables except for the \code{.flow_var} and \code{.habitat_var} will be treated as groups.}

\item{.flow_threshold_var}{Unquoted column name for the flow threshold defining baseflow for each group.}

\item{.flow_var}{Unquoted column name for the flow in cfs, defaults to \code{flow_cfs}}

\item{.habitat_var}{Unquoted column name for the suitable area, defaults to \code{wua_per_lf}}
}
\value{
A modified version of the original \code{data.frame} or \code{tbl_df} with additional flow and habitat columns identified by "_nbfc" ("no baseflow channel)
}
\description{
Remove Baseflow via Flow Threshold (Multiple)
}
