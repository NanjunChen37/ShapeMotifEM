% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMEM.R
\name{SMEM}
\alias{SMEM}
\title{DNA shape motif discovery with Expectation–Maximization algorithm}
\usage{
SMEM(
  dna_shape_data,
  filename = "SMEM",
  peakCount = 20,
  motifLength = 12,
  motifCount = 2,
  replicates = 50,
  tolerance = 1e-06
)
}
\arguments{
\item{dna_shape_data}{Input DNA shape data convert by function fasta_input or bed_input.}

\item{filename}{character, name of output .csv file.}

\item{peakCount}{integer, the number of peaks to retrieve per batch, default is 20.}

\item{motifLength}{integer, length of shape motif to be discovered, default is 12.}

\item{motifCount}{integer, number of motifs to be discovered for each peak, default is 2.}

\item{replicates}{integer, number of repetitions for each batch, default is 50.}

\item{tolerance}{double, tolerance threshold for EM iteration, default is 0.000001.}
}
\value{
A list of motif location arrays
}
\description{
This function reads the DNA shape data file and conduct shape motif discovery using the Expectation–Maximization algorithm.
}
