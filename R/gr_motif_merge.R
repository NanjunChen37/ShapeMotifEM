#' Count and merge DNA shape motif location given by Expectation–Maximization algorithm
#'
#' This function reads the location array of shape motif and output final motif location results. The motif discovery with GRanges file input can get the genome location of the motif.
#'
#' @param motif_location_array the motif location array returned by function 'SMEM'.
#' @param gr_file the GRange file returned by function 'bed_input'.
#' @param motifLength integer, length of shape motif to be discovered, default is 12.
#' @param filename character, name of output file, default is "location_merge".
#' @return A list of motif location data frame and motif shape data
#' @import Biostrings
#' @import GenomicRanges
#' @import DNAshapeR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom stats cor dist dnorm rnorm sd
#' @importFrom utils read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @export gr_motif_merge
gr_motif_merge <- function(motif_location_array, gr_file, motifLength = NULL, filename = "location_merge"){

  # Initialization of global parameters
  peakLength = unique(width(gr_file))
  lenGrSet = ncol(motif_location_array)
  grResize = gr_file

  # Load required packages
  library(Biostrings)
  library(DNAshapeR)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # Function of getting extended gr_file and shape data
  get_extend_gr_shape <- function(gr_file, motifLength, peakLength, filename){

    # Get GRange file with extended length
    grExtendUp = flank(grResize, motifLength)
    grExtend = resize(grExtendUp, peakLength + motifLength * 2)
    grWidth = unique(width(grExtend))

    # Get FASTA and shape data of extended GRange file
    getFasta(grExtend, Hsapiens, width = grWidth + 4, filename = paste0(filename, "_extend_fasta.fa"))
    fn_Final = paste0(filename, "_extend_fasta.fa")
    extend_shape = getShape(fn_Final)
    shape_name = c('Rise', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'Shear', 'Stagger', 'Stretch')
    for (sn in shape_name){
      temp_shape = getShape(fn_Final, shapeType = sn)
      extend_shape = append(extend_shape, temp_shape)
    }

    # Transpose shape data
    extend_shape_transpose = array(0, dim = c(length(extend_shape), grWidth, nrow(extend_shape[[1]])))
    for (a in 1:length(extend_shape)) {
      extend_shape_transpose[a , ,] = matrix(t(extend_shape[[a]][, 3:(grWidth + 2)]), grWidth, nrow(extend_shape[[1]]))
    }

    return(list("extend_gr_file" = grExtend, "extend_shape" = extend_shape_transpose))
  }

  # Function of getting reduced motif location, merge-to location, and frequency on each peak
  get_motif_location <- function(location_vector, motifLength){

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

    return(list("reduce_location" = reduce_location, "merge_location" = merge_location, "frequency" = frequency))

  }

  # Function of getting genomic motif location
  get_genomic_location <- function(grSplit, reduce_location_vector, merge_location_vector) {

    # Variable initialization
    genomic_location = peak_number_vector = c()
    grEach = rep(grSplit[[peakNum]], times = length(reduce_location_vector))
    # mergeShapeSet

    for (vectorLen in 1:length(reduce_location_vector)) {
      # genomic level location output
      temp_genomic_location = start(grResize)[[peakNum]] + reduce_location_vector[vectorLen] - 1

      # Get genomic location from the GRange file
      genomic_location = append(genomic_location, temp_genomic_location)
      genomic_location_set[[peakNum]] = genomic_location

      # Record repeat peak number
      peak_number_vector = append(peak_number_vector, peakNum)
      peak_number_set[[peakNum]] = peak_number_vector

      # Get GRange file for each motif
      if (reduce_location_vector[vectorLen] == 0) {reduce_location_vector[vectorLen] = 1}
      grEach[vectorLen] = narrow(grEach[vectorLen], reduce_location_vector[vectorLen])
      grEach[vectorLen] = resize(grEach[vectorLen], merge_location_vector[vectorLen] - reduce_location_vector[vectorLen] + 1)
      # grEach[vectorLen] = resize(grEach[vectorLen], merge_location_set[[peakNum]][[vectorLen]] - final_location_set[[peakNum]][[vectorLen]] + 1)

      # Merge motifs GRange on the same peak
      if (vectorLen == 1) {
        grMotif[[peakNum]] = grEach[vectorLen]
      } else {
        grMotif[[peakNum]] = c(grMotif[[peakNum]], grEach[vectorLen])
      }
    }

    return(list("grMotif" = grMotif[[peakNum]], "peak_number_set" = peak_number_set))

  }

  # Function of getting shape data of each motif
  get_motif_shape <- function(extend_shape, motifLength, reduce_location_vector, merge_location_vector) {

    merge_shape_set = c()

    for (vectorLen in 1:length(reduce_location_vector)) {

      # shape data extraction
      extendStart = motifLength + reduce_location_vector[vectorLen] + 1
      extendEnd = motifLength + merge_location_vector[vectorLen] + 1

      #
      merge_shape = extend_shape[, , peakNum][, extendStart:extendEnd]
      merge_shape_set[[vectorLen]] = merge_shape
    }

    return(merge_shape_set)
  }

  # Call get extended grfile and shape data
  extend = get_extend_gr_shape(gr_file, motifLength, peakLength, filename)
  extend_gr_file = extend$extend_gr_file
  extend_shape = extend$extend_shape

  # Sort motif location matrix
  sort_location = apply(motif_location_array, 2, sort)

  # Initialization of global variables
  final_location_set = merge_location_set = frequency_set = c()
  genomic_location_set = peak_number_set = grMotif = motif_genomic_location_set = motif_shape_data = c()

  grSplit = split(grResize, rep(1:length(grResize), each = 1))

  pb = txtProgressBar(style = 3)
  for (peakNum in 1:ncol(motif_location_array)) {
    # Retrive each peak vector
    location_vector = sort_location[, peakNum]

    # Call location merge
    location_merge = get_motif_location(location_vector, motifLength)
    reduce_location_vector = location_merge$reduce_location
    merge_location_vector = location_merge$merge_location
    motif_frequncy = location_merge$frequency

    # Call get_genomic_location
    genomic_location = get_genomic_location(grSplit, reduce_location_vector, merge_location_vector)
    motif_genomic_location = genomic_location$grMotif
    peak_number_set = genomic_location$peak_number_set

    # Call get_motif_shape
    merge_shape_set = get_motif_shape(extend_shape, motifLength, reduce_location_vector, merge_location_vector)


    # Get final reduced location, merge-to location, and frequency
    final_location_set[[peakNum]] = reduce_location_vector
    merge_location_set[[peakNum]] = merge_location_vector
    frequency_set[[peakNum]] = motif_frequncy

    # Output the motif location GRange file
    if (peakNum == 1) {
      motif_genomic_location_set = motif_genomic_location
    } else {
      motif_genomic_location_set = c(motif_genomic_location_set, motif_genomic_location)
    }

    # Correspond merge shape data
    motif_shape_data[[peakNum]] = merge_shape_set

    setTxtProgressBar(pb, peakNum / ncol(motif_location_array))
  }

  # Output result data frame
  write.table(motif_genomic_location_set, file = paste0(filename, "_BED_Table"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE )

  bed_table = read.table(paste0(filename, "_BED_Table"), header = FALSE)
  colnames(bed_table) = c("seqnames", "start", "end", "width", "strand")

  motif_data_frame = data.frame(
    peakNo. = unlist(peak_number_set),
    seriesNo. = unlist(lapply(peak_number_set, function(x) seq_along(x))),
    chrNo. = bed_table$seqnames,
    location = unlist(final_location_set),
    genoStart = bed_table$start + 1,
    genoEnd = bed_table$end + 1,
    motiflength = motifLength,
    mergerMotifLength = bed_table$end - bed_table$start + 1,
    motifFrequency = round(unlist(frequency_set), 4)
  )

  write.csv(motif_data_frame, file = paste0(filename, "_motif_data_frame.csv"), row.names = FALSE)

  return(list("motif_data_frame" = motif_data_frame, "motif_shape_data" = motif_shape_data))

}


# motif_location_array = source("test_location_new")$value
# gr_file = source("gr_test")$value
#
# gr_merge_test = gr_motif_merge(motif_location_array, gr_file, motifLength = 10, filename = "sortlocation")
#
# merge_test$motif_data_frame[1:3,]
# merge_test$motif_shape_data[1:3]

# gr_merge_test = gr_motif_merge(em_test, gr_test$gr_file, motifLength = 10, filename = "sortlocation")
# motif_location_array = source("test_location_new")$value
# gr_file = source("gr_test")$value
# peakCount = 20
# motifLength = 10
# peakLength = unique(width(gr_file))
# lenGrSet = ncol(motif_location_array)
# filename = "sortlocation"
