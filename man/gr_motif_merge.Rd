% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gr_motif_merge.R
\name{gr_motif_merge}
\alias{gr_motif_merge}
\title{Count and merge DNA shape motif location given by Expectation–Maximization algorithm}
\usage{
gr_motif_merge(
  motif_location_array,
  gr_file,
  motifLength = NULL,
  filename = "location_merge"
)
}
\arguments{
\item{motif_location_array}{the motif location array returned by function 'SMEM'.}

\item{gr_file}{the GRange file returned by function 'bed_input'.}

\item{motifLength}{integer, length of shape motif to be discovered, default is 12.}

\item{filename}{character, name of output file, default is "location_merge".}
}
\value{
A list of motif location data frame and motif shape data
}
\description{
This function reads the location array of shape motif and output final motif location results. The motif discovery with GRanges file input can get the genome location of the motif.
}
