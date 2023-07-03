#' Read and process DNA sequence from a FASTA file
#'
#' This function reads a DNA sequence from a FASTA file and processes it to the DNA shape data. The DNA shape data can be used as the input of function SMEM() or SMEM_Gibbs() for DNA shape motif discovery
#'
#' @param fasta_file path to the FASTA file containing the DNA sequence.
#' @param shapeind integer, the index of shape feature to process, default is "all". 1 to 5 indicate feature "MGW", "HelT", "ProT", "Roll", and "EP" respectively.
#' @return An array of processed shape data.
#' @import Biostrings
#' @import DNAshapeR
#' @export fasta_input
fasta_input <- function(fasta_file, shapeind = "all") {

  # Check that the file exists
  if (!file.exists(fasta_file)) {
    stop(paste0("The file '", fasta_file, "' does not exist."))
  }

  # Check if shapeind is valid
  if(!is.numeric(shapeind) && shapeind != "all") {
    stop("shapeind must be an integer between 1 and 5 or 'all'")
  }

  # Load required packages
  library(Biostrings)
  library(DNAshapeR)

  # Read FASTA
  sequences = readDNAStringSet(paste0(fasta_file))
  sequences_set = unique(unlist(as.character(sequences)))

  if (length(unique(nchar(sequences_set))) == 1){
    prePeakLength = unique(nchar(sequences_set))
    peakLength = prePeakLength - 4
    print(paste0("The length of the input sequence is ", prePeakLength, ". The peakLength used for motif discovery is ", peakLength))
  } else {
    stop("The length of FASTA sequence is not uniform, pleas check.")
  }

  # Check if peakLength is valid
  if (!is.numeric(peakLength) || peakLength <= 0) {
    stop("peakLength must be a numeric value greater than 0.")
  }

  # Process DNA sequence using DNAshapeR
  # shape_data = getShape(fasta_file)
  #
  # fn_AHFile = fasta_name
  shape_data = getShape(fasta_file)
  shape_name = c('Rise', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'Shear', 'Stagger', 'Stretch')
  for (sn in shape_name){
    temp_shape = getShape(fasta_file, shapeType = sn)
    shape_data = append(shape_data, temp_shape)
  }

  # Check that the peak length is avaliable
  if (prePeakLength > ncol(shape_data[[1]])) {
    stop(paste0("peakLength too large, select a suitable peakLength"))
  }

  # Transpose and process shape data based on shapeind
  if (shapeind == "all") {
    shape_data_transpose = array(0, dim = c(length(shape_data), prePeakLength - 4, nrow(shape_data[[1]])))
    for (a in 1:length(shape_data)) {
      shape_data_transpose[a , , ] = matrix(t(shape_data[[a]][, 3:(prePeakLength - 2)]), prePeakLength - 4, nrow(shape_data[[1]]))
    }
  } else if (is.numeric(shapeind) && shapeind >= 1 && shapeind <= 14) {
    shape_data_transpose = array(0, dim = c(length(shape_data) / 14, prePeakLength - 4, nrow(shape_data[[1]])))
    shape_data_transpose[1 , ,] = matrix(t(shape_data[[shapeind]][, 3:(prePeakLength - 2)]), prePeakLength - 4, nrow(shape_data[[1]]))
  } else {
    stop("shapeind must be an integer between 1 and 14 or 'all'")
  }

  # Return processed shape data
  return(shape_data_transpose)
}

# fasta_file = "fasta.fa"
# test = fasta_input(fasta_file, shapeind = 'all')
# test
