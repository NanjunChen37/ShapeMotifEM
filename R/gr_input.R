#' Read and process DNA sequence from a GRanges file
#'
#' This function reads a DNA sequence from a GRanges file and processes it to the DNA shape data. The DNA shape data can be used as the input of function SMEM() or SMEM_Gibbs() for DNA shape motif discovery.
#'
#' @param gr_file The GRanges file containing the peak information.
#' @param peakLength integer, length of center peak to be processed, default is 100.
#' @param total_peak_number integer, total peak number, default is 100.
#' @param fasta_name file name of the output FASTA file, default is "fasta.fa".
#' @param shapeind integer, the index of shape feature to process, default is "all". 1 to 14 indicate feature '1.MGW', '2.HelT', '3.ProT', '4.Roll', '5.EP', '6.Rise', '7.Shift', '8.Slide', '9.Tilt', '10.Buckle', '11.Opening', '12.Shear', '13.Stagger', '14.Stretch' respectively.
#' @param sort_by parameter used to sort the gr file, default is "signalValue".
#' @param sort_descend sort the peaks in a descending order according to 'sort_by', default is "TRUE".
#' @return An array of processed shape data
#' @import GenomicRanges
#' @import DNAshapeR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom stats cor dist dnorm rnorm sd
#' @importFrom utils read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @export gr_input
gr_input <- function(gr_file, peakLength = 100, total_peak_number = 100, fasta_name = "fasta.fa",
                     shapeind = "all", sort_by = "signalValue", sort_descend = TRUE) {

  # Load required packages
  library(GenomicRanges)
  library(DNAshapeR)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # Parameter initialization
  lenGrTop = lenGrSet = total_peak_number
  lenGrReduce = 0
  prePeakLength = peakLength + 4
  halfPeakLength = (prePeakLength - 4) / 2
  shape_data = c()

  # Determine the column name used to sort the gr file
  if (sort_by %in% c("score", "signalValue", "pValue", "qValue", "peak")) {
    if (sort_by == "signalValue") {
      grOrder = gr_file[order(gr_file$signalValue, decreasing = sort_descend)]
    } else if (sort_by == "score") {
      grOrder = gr_file[order(gr_file$score, decreasing = sort_descend)]
    } else if (sort_by == "pValue") {
      grOrder = gr_file[order(gr_file$pValue, decreasing = sort_descend)]
    } else if (sort_by == "qValue") {
      grOrder = gr_file[order(gr_file$qValue, decreasing = sort_descend)]
    } else {
      grOrder = gr_file[order(gr_file$peak, decreasing = sort_descend)]
    }
  } else if (sort_by == "none") {
    print("User defined gr file")
    grOrder = gr_file
  } else {
    stop( "Invalid value for 'sort_by'. Choose 'sort_by = 'score, signalValue, pValue, qValue, peak' or 'none''")
  }

  # Reduce overlap and extract peaks
  while (lenGrReduce < lenGrSet) {
    grTop = grOrder[1:lenGrTop]
    grReduce = reduce(grTop)
    lenGrTop = length(grTop)
    lenGrReduce = length(grReduce)
    overlap = lenGrSet - lenGrReduce
    lenGrTop = lenGrTop + overlap
    if (lenGrTop >= length(grOrder)) {
      break
    }
  }

  # Report the shorteset peak that short than peakLength
  if (min(width(grReduce)) < peakLength) {
    print(paste0(
      "The shortest peak is peak number ",
      which.min(width(grReduce)),
      ". Its length is less than peakLength, peakLength has been set to ",
      min(width(grReduce))
    ))
  }

  # Get the center base pair of the gr file
  # Split the gr file for processing
  grSplit <- split(grReduce, rep(1:length(grReduce), each = 1))

  print("Extracting peaks, please wait...")
  Sys.sleep(1)
  # pb <- txtProgressBar(style = 3)

  for (a in 1:length(grReduce)) {
    halfWidth = width(grSplit[[a]]) %/% 2
    narrowRange = halfWidth - halfPeakLength
    grSplit[[a]] = narrow(grSplit[[a]], narrowRange)
    grSplit[[a]] = resize(grSplit[[a]], (prePeakLength - 4))
    if (a == 2) {
      grResize <- c(grSplit[[1]], grSplit[[2]])
    }
    if (a > 2) {
      grResize <- c(grResize, grSplit[[a]])
    }
    # setTxtProgressBar(pb, a / length(grReduce))
  }

  # Get FASTA sequence and keep
  tryCatch({
    getFasta(grResize, Hsapiens, width = prePeakLength, filename = fasta_name)
  }, error = function(e) {
    # If an error occurs, print a prompt message and continue with the next command
    message(paste0("Error occurred: fail to get FASTA sequence \n",
                   e$message, "\n",
                   "Please check: the packages needed is install correctly \n",
                   "If error still occurred: \n",
                   "1. Use other methods to obtain FASTA \n",
                   "2. Then use 'FASTAInput()' to obtain shape information"))
  })

  # Retrive fasta file and predict the shape
  tryCatch({
    fn_AHFile = fasta_name
    shape_data = getShape(fn_AHFile)
    shape_name = c('Rise', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'Shear', 'Stagger', 'Stretch')
    for (sn in shape_name){
      temp_shape = getShape(fn_AHFile, shapeType = sn)
      shape_data = append(shape_data, temp_shape)
    }

  }, error = function(e) {
    message(paste0("Error occurred: FASTA not get \n", e$message))
  })

  # Transpose and process shape data based on shapeind
  tryCatch({
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
  }, error = function(e) {
    message(paste0("Error occurred: DNA shape data not get \n", e$message))
  }, finally = {
    if (is.null(shape_data)) {
      print("The modified gr file is returned as 'gr_file'")
      return("gr_file" = grResize)
    } else {
      print("the modified gr file and DNA shape data is returned as '$gr_file' and '$dna_shape_data'")
      return(list("gr_file" = grResize, "dna_shape_data" = shape_data_transpose))
    }
  })

}

# example_data = system.file("data", "A549_MAX_gr.RData", package = "ShapeMotifEM")
# example_gr_file = load(example_data)
# gr_file = A549_MAX_gr
#
# gr_test = gr_input(gr_file, peakLength = 100, total_peak_number = 100, fasta_name = "fasta.fa",
#                     shapeind = "all", sort_by = "signalValue", sort_descend = TRUE)
