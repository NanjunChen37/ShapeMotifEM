#' Count and merge DNA shape motif location given by Expectation–Maximization algorithm
#'
#' This function reads the location array of shape motif and output final motif location results. The motif discovery with BED file input can get the genome location of the motif. The FASTA input will get results only based on the FASTA file.
#'
#' @param motif_location_array the motif location array returned by function 'SMEM'.
#' @param fasta_file the fasta file used for motif discovery.
#' @param motifLength integer, length of shape motif to be discovered, default is 12.
#' @param filename character, name of output file, default is "location_merge".
#' @return A list of motif location data frame and motif shape data
#' @import Biostrings
#' @import GenomicRanges
#' @import DNAshapeR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom stats cor dist dnorm rnorm sd
#' @importFrom utils read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @export fasta_motif_merge
fasta_motif_merge <- function(motif_location_array, fasta_file, motifLength = NULL, filename = "location_merge"){

  # Load required packages
  library(Biostrings)
  library(DNAshapeR)
  library(GenomicRanges)

  # Read FASTA
  sequences = readDNAStringSet(paste0(fasta_file))
  sequences_set = unique(unlist(as.character(sequences)))

  if (length(unique(nchar(sequences_set))) == 1){
    prePeakLength = unique(nchar(sequences_set))
    peakLength = prePeakLength - 4
  } else {
    stop("The length of FASTA sequence is not uniform, pleas check.")
  }

  # Function of getting reduced motif location, merge-to location, and frequency on each peak
  get_motif_location <- function(location_vector, motifLength, peakLength) {

    # Variable initialization
    tempVecLength = 1
    location_frequncy = 0
    reduce_location_vector = frequency = merge_location = c()
    temp_vector = location_vector
    min_value = min(location_vector)

    # Process on each column vector
    for (vectorInd in 1:length(location_vector)) {

      if ((temp_vector[vectorInd] - min_value) == 0 | (temp_vector[vectorInd] - min_value) < motifLength) {
        # Merge location with in motifLength
        reduce_location_vector = append(reduce_location_vector, min_value)
        temp_merge_location = min_value + (motifLength - 1) + (temp_vector[vectorInd] - min_value)

        # Record frequency
        location_frequncy = location_frequncy + 1
        frequency[tempVecLength] = location_frequncy / length(location_vector)

      } else {
        # Append non-overlapping motifs
        reduce_location_vector = append(reduce_location_vector, temp_vector[vectorInd])

        min_value = temp_vector[vectorInd]
        temp_merge_location = min_value + (motifLength - 1)

        # Record frequency
        location_frequncy = 1
        tempVecLength = tempVecLength + 1
        frequency[tempVecLength] = location_frequncy / length(location_vector)
      }

      if (vectorInd > 1 &&
          temp_merge_location - max(merge_location) < motifLength &&
          temp_vector[vectorInd] - reduce_location_vector[vectorInd - 1] < motifLength) {
        merge_location[length(merge_location)] = temp_merge_location
      }

      merge_location = append(merge_location, temp_merge_location)
      merge_location = unique(merge_location)
    }

    reduce_location = unique(reduce_location_vector)

    # Filter some of the motifs for shape extraction
    filter_index = reduce_location >= 3 & reduce_location <= (peakLength - motifLength + 1)

    reduce_location_filtered = reduce_location[filter_index]
    merge_location_filtered = merge_location[filter_index]
    merge_location_filtered[merge_location_filtered > peakLength] = peakLength
    frequency_filtered = frequency[filter_index]

    if (length(reduce_location_filtered) == 0) {
      reduce_location_filtered = 1
      merge_location_filtered = 2
      frequency_filtered = 0
    }


    return(list("reduce_location" = reduce_location_filtered, "merge_location" = merge_location_filtered, "frequency" = frequency_filtered))

  }

  # Function of getting shape data of each peak
  get_peak_shape <- function(fasta_file, prePeakLength) {

    shape_data = getShape(fasta_file)
    shape_name = c('Rise', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'Shear', 'Stagger', 'Stretch')
    for (sn in shape_name){
      temp_shape = getShape(fasta_file, shapeType = sn)
      shape_data = append(shape_data, temp_shape)
    }

    shape_data_transpose = array(0, dim = c(length(shape_data), prePeakLength - 4, nrow(shape_data[[1]])))
    for (a in 1:length(shape_data)) {
      shape_data_transpose[a , , ] = matrix(t(shape_data[[a]][, 3:(prePeakLength - 2)]), prePeakLength - 4, nrow(shape_data[[1]]))
    }

    return(shape_data_transpose)
  }

  # Function of getting shape data of each motif
  get_motif_shape_seq <- function(shape_data, sequences_set, motifLength, reduce_location_vector, merge_location_vector) {

    merge_shape_set = merge_sequence_set = peak_number_vector = c()

    for (vectorLen in 1:length(reduce_location_vector)) {

      # Shape data extraction
      extendStart = reduce_location_vector[vectorLen]
      extendEnd = merge_location_vector[vectorLen]

      # Extract merged shape data
      merge_shape = shape_data[, , peakNum][, extendStart:extendEnd]
      merge_shape_set[[vectorLen]] = merge_shape

      merge_sequence = subseq(sequences_set[peakNum], start = (extendStart + 2), end = (extendEnd + 2))
      merge_sequence_set[[vectorLen]] = merge_sequence

      # Record repeat peak number
      peak_number_vector = append(peak_number_vector, peakNum)
      peak_number_set[[peakNum]] = peak_number_vector

    }

    return(list("merge_shape_set" = merge_shape_set, "merge_sequence_set" = merge_sequence_set, "peak_number_set" = peak_number_set))
  }

  ###############

  shape_data = get_peak_shape(fasta_file, prePeakLength)

  # Sort motif location matrix
  sort_location = apply(motif_location_array, 2, sort)

  # Read fasta file
  sequences = readDNAStringSet(paste0(fasta_file))
  sequences_set = unique(unlist(as.character(sequences)))

  # Initialization of global variables
  final_location_set = merge_location_set = frequency_set = c()
  peak_number_set = motif_shape_data = motif_sequence_data = c()

  pb = txtProgressBar(style = 3)
  for (peakNum in 1:ncol(motif_location_array)) {
    # Retrive each peak vector
    location_vector = sort_location[, peakNum]

    # Call location merge
    location_merge = get_motif_location(location_vector, motifLength, peakLength)
    reduce_location_vector = location_merge$reduce_location
    merge_location_vector = location_merge$merge_location
    motif_frequncy = location_merge$frequency

    # Call get_motif_shape
    get_shape_seq = get_motif_shape_seq(shape_data, sequences_set, motifLength, reduce_location_vector, merge_location_vector)
    merge_shape_set = get_shape_seq$merge_shape_set
    merge_sequence_set = get_shape_seq$merge_sequence_set
    peak_number_set = get_shape_seq$peak_number_set



    # Get final reduced location, merge-to location, and frequency
    final_location_set[[peakNum]] = reduce_location_vector
    merge_location_set[[peakNum]] = merge_location_vector
    frequency_set[[peakNum]] = motif_frequncy

    # Correspond merge shape data
    motif_shape_data[[peakNum]] = merge_shape_set
    motif_sequence_data[[peakNum]] = merge_sequence_set

    setTxtProgressBar(pb, peakNum / ncol(motif_location_array))
  }

  motif_data_frame = data.frame(
    peakNo. = unlist(peak_number_set),
    seriesNo. = unlist(lapply(peak_number_set, function(x) seq_along(x))),
    locationStrat = unlist(final_location_set),
    locationEnd = unlist(merge_location_set),
    motifLength = motifLength,
    mergerMotifLength = unlist(merge_location_set) - unlist(final_location_set) + 1,
    sequence = unlist(motif_sequence_data),
    motifFrequency = round(unlist(frequency_set), 4)
  )

  write.csv(motif_data_frame, file = paste0(filename, "_motif_data_frame.csv"), row.names = FALSE)

  return(list("motif_data_frame" = motif_data_frame, "motif_shape_data" = motif_shape_data))

}


# motif_location_array = source("test_location_new")$value
# fasta_file = "HepG2_AH22967.fa"
#
# fa_mege_test = fasta_motif_merge(em_test, fasta_file, motifLength = 10, filename = "fasta_location_merge")
