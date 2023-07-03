#' Visualize the DNA shape motif output by the merge functions
#'
#' This function reads the motif shape data and the motif data frame
#'
#' This function output possible shape motifs shape feature and sequence logo
#'
#' @param motif_data_frame the motif motif_data_frame returned by function 'fasta/bed_motif_merge'.
#' @param motif_shape_data the motif motif_shape_data returned by function 'fasta/bed_motif_merge'.
#' @param motifFrequency integer, the threshold for selected motifs with top exsitence frequency.
#' @param shapeIndex integer, the index of shape feature to process, 1 to 5 indicate feature "MGW", "HelT", "ProT", "Roll", and "EP" respectively.
#' @param motifLength integer, length of shape motif to be aligned and visualized.
#' @param align_method character, the metric used for alignment, default is "pearson", can also choose "distance".
#' @param ref_number integer, the reference motif for alignment, default is 1.
#' @param input_form integer, '0' for 'fasta_motif_merge' or '1' for 'bed_motif_merge', user input while the function is running.
#' @param continue_alignment integer, '0' for 'No' or '1' for 'Yes', user input while the function is running.
#' @param top_motif integer, number of motifs to align, user input while the function is running.
#' @param threshold double, the threshold for align method 'pearson', user input while the function is running.
#' @return A list of motif location data frame and motif shape data
#' @import Biostrings
#' @import reshape2
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom stats cor dist dnorm rnorm sd
#' @importFrom utils read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @export motif_visualizer
motif_visualizer <- function(motif_data_frame, motif_shape_data, motifFrequency = 0.7, shapeIndex = NULL, motifLength = NULL,
                             align_method = "pearson", ref_number = 1,
                             input_form = NULL, continue_alignment = NULL, top_motif = NULL, threshold = NULL) {

  # Load required packages
  library(reshape2) # for reshape2::melt()
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # Function of motif aligner pearson
  motif_aligner_pearson <- function(aindex, extractMotifLength, threshold, input_form, sequence_data) {

    # Get compare motifs
    # a: the index of compared motifs
    # m: the total length for the original vector
    # n: total length of reduced motifs
    partialVecListAll = c()
    n = extractMotifLength
    m = length(aindex[[1]])
    for (a in 2:length(aindex)) {
      partialVecList = partialVecSet = c()

      for (i in 1:(m - n + 1)) {
        partialVec = aindex[[a]][i:(i + n - 1)]
        partialVecSet[[i]] = partialVec
      }

      partialVecListAll[[a - 1]] = partialVecSet
    }

    # Get center reference motifs and align
    compareIndexrec = nrec = targetNumrec = compareNumSetAll = vectorNumSetAll = sequenceIndexSetAll = c()
    targetVecList = compareMatrixList = targetVecSet = compareMatrixSet = c()
    outputMatrixAll = outputMatrix = distRecordSetAll = seqListAll = c()

    if (input_form == 0) {
      print("Shape motif alignment and get sequence ...")
    } else {
      print("Shape motif alignment ...")
    }

    pb <- txtProgressBar(style = 3)

    # Measure distance for each target vector
    for (i in 1:(m - n + 1)) {
      targetVec = aindex[[1]][i:(i + n - 1)]

      # Build compare matrix
      outputMatrix = targetVec

      # Index for each compare vector set
      compareNumSet = c(i)
      vectorNumSet = c(1)
      distRecordSet = c()

      if (input_form == 0) {
        sequenceSet = c(subseq(sequence_data[[1]], start = compareNumSet, end = (compareNumSet + motifLength - 1)))
      }

      for (b in 1:length(partialVecListAll)) {
        minDist = Inf
        compareMatrix = targetVec

        # Index for all compare vectors in one set
        for (j in 1:length(partialVecListAll[[b]])) {
          if (!any(is.na(partialVecListAll[[b]][[j]]))) {
            compareMatrix = unname(rbind(compareMatrix, partialVecListAll[[b]][[j]]))
          }
        }

        max_Pearson = 0
        output_vector = c()
        vectorNum = compareNum = sequenceNum = NULL
        for (k in 1:(nrow(compareMatrix) - 1)) {
          vector_Pearson = cor(compareMatrix[1, ], compareMatrix[k + 1, ])
          #   print(paste0("k =", k, ",pearson = ", vector_Pearson))

          if (vector_Pearson > max_Pearson & vector_Pearson > threshold) {
            max_Pearson = vector_Pearson
            output_vector = compareMatrix[k + 1, ]
            vectorNum = b + 1
            compareNum = k
            if (input_form == 0) {
              sequenceNum = subseq(sequence_data[[b + 1]], start = k, end = (k + motifLength - 1))
            }
          }
        }
          # paste("vectorNum = ", b+1, "compareNum =", compareNum, "max_Pearson =", max_Pearson))
        outputMatrix = rbind(outputMatrix, output_vector)
        vectorNumSet = append(vectorNumSet, vectorNum)
        compareNumSet = append(compareNumSet, compareNum) # use compare number to record each fragment
        sequenceIndexSet = compareNumSet
        distRecordSet = append(distRecordSet, max_Pearson)
        if (input_form == 0) {
          sequenceSet = append(sequenceSet, sequenceNum)
        }
      }

      outputMatrixAll[[i]] = outputMatrix
      vectorNumSetAll[[i]] = vectorNumSet
      compareNumSetAll[[i]] = compareNumSet
      sequenceIndexSetAll[[i]] = sequenceIndexSet
      distRecordSetAll[[i]] = distRecordSet
      if (input_form == 0) {
        seqListAll[[i]] = sequenceSet
      }

      setTxtProgressBar(pb, i / (m - n + 1))
    }

    # Pearson threshold to filter motifs, remove single results
    single_indices = which(lengths(vectorNumSetAll) == 1)
    if (length(single_indices) > 0) {
      outputMatrixAll = outputMatrixAll[-single_indices]
      vectorNumSetAll = vectorNumSetAll[-single_indices]
      compareNumSetAll = compareNumSetAll[-single_indices]
      sequenceIndexSetAll = sequenceIndexSetAll[-single_indices]
      distRecordSetAll = distRecordSetAll[-single_indices]
      seqListAll = seqListAll[-single_indices]
    }

    means = numeric(length(distRecordSetAll))
    for (i in seq_along(distRecordSetAll)) {
      non_zero_values = distRecordSetAll[[i]][distRecordSetAll[[i]] != 0]
      means[i] = mean(non_zero_values)
    }

    distRecordSetAllmean = means

    # Get sequence when use GRange input
    if (input_form == 1){

      getSeqData <- function(chrName, startPoint, extractMotifLength){
        seq = getSeq(Hsapiens, chrName, start = startPoint, end = startPoint + extractMotifLength - 1)
        seqString = toString(seq)
        return(seqString)
      }

      print("Get alignment sequence ...")
      pb <- txtProgressBar(style = 3)

      seqList = c()
      for (i in 1:length(vectorNumSetAll)){
        vector_index = vectorNumSetAll[[i]]
        sequence_index = sequenceIndexSetAll[[i]]
        seqList = c()
        for (j in 1:length(vector_index)){
          chrName = df_freq_order$chrNo.[vector_index[j]]
          startPoint = df_freq_order$genoStart[vector_index[j]] + sequence_index[j]
          seq = getSeqData(chrName, startPoint, extractMotifLength)
          seqList = c(seqList, seq)
        }

        seqListAll[[i]] = seqList

        setTxtProgressBar(pb, i / length(vectorNumSetAll))
      }
    }

    return(
      list(
        "alignedShapeData" = outputMatrixAll,
        "vectorNumSetAll" = vectorNumSetAll,
        "locationindexSetAll" = compareNumSetAll,
        "sequenceIndexSetAll" = sequenceIndexSetAll,
        "distRecordSetAll" = distRecordSetAllmean,
        "seqListAll" = seqListAll
      )
    )

  }

  # Function of motif aligner pearson with all features
  motif_aligner_all_pearson <- function(aindex, extractMotifLength, threshold, input_form, sequence_data) {

    # Get compare motifs
    # a: the index of compared motifs
    # m: the total length for the original vector
    # n: total length of reduced motifs
    partialVecListAll = c()
    n = extractMotifLength
    m = ncol(aindex[[1]])
    for (a in 2:length(aindex)) {
      partialVecList = partialVecSet = c()
      for (i in 1:(m - n + 1)) {
        if (i + n - 1 > ncol(aindex[[a]])) {
          next
        }
        partialVec = aindex[[a]][,i:(i + n - 1)]
        partialVecSet[[i]] = partialVec
      }
      partialVecListAll[[a - 1]] = partialVecSet
    }

    # Get center reference motifs and align
    compareIndexrec = nrec = targetNumrec = compareNumSetAll = vectorNumSetAll = sequenceIndexSetAll = c()
    targetVecList = compareMatrixList = targetVecSet = compareMatrixSet = c()
    outputMatrixAll = outputMatrix = distRecordSetAll = seqListAll = c()

    if (input_form == 0) {
      print("Shape motif alignment and get sequence ...")
    } else {
      print("Shape motif alignment ...")
    }

    pb <- txtProgressBar(style = 3)

    # Measure distance for each target vector
    for (i in 1:(m - n + 1)) {
      outputMatrix = c()

      targetVec = aindex[[1]][,i:(i + n - 1)] #

      # Build compare matrix
      outputMatrix[[1]] = targetVec #

      # Index for each compare vector set
      compareNumSet = c(i)
      vectorNumSet = c(1)
      distRecordSet = compareMatrix = c() #

      if (input_form == 0) {
        sequenceSet = c(subseq(sequence_data[[1]], start = compareNumSet, end = (compareNumSet + motifLength - 1)))
      }

      for (b in 1:length(partialVecListAll)) {
        minDist = Inf
        compareMatrix[[1]] = targetVec #

        # Index for all compare vectors in one set
        for (j in 1:length(partialVecListAll[[b]])) {
          compareMatrix[[j + 1]] = partialVecListAll[[b]][[j]]
          # compareMatrix = unname(rbind(compareMatrix, partialVecListAll[[b]][[j]]))
        }

        max_Pearson = 0
        Pearsons = -Inf
        output_vector = c()
        vectorNum = compareNum = sequenceNum = NULL
        for (k in 1:(length(compareMatrix) - 1)) { #
          # vector_Pearson = cor(compareMatrix[1, ], compareMatrix[k + 1, ])
          #   print(paste0("k =", k, ",pearson = ", vector_Pearson))
          tempPearson = 0
          for(p in 1:nrow(compareMatrix[[1]])){
            tempPearson = tempPearson + cor(compareMatrix[[k + 1]][p,], compareMatrix[[1]][p,], method='pearson')
          }
          Pearsons = tempPearson / nrow(compareMatrix[[1]])
          # print(Pearsons)

          if (Pearsons > max_Pearson & Pearsons > threshold) {
            max_Pearson = Pearsons
            output_vector = compareMatrix[[k + 1]]
            vectorNum = b + 1
            compareNum = k
            if (input_form == 0) {
              sequenceNum = subseq(sequence_data[[b + 1]], start = k, end = (k + motifLength - 1))
            }
          }
        }
        # paste("vectorNum = ", b+1, "compareNum =", compareNum, "max_Pearson =", max_Pearson))
        outputMatrix[[length(outputMatrix) + 1]] = output_vector
        vectorNumSet = append(vectorNumSet, vectorNum)
        compareNumSet = append(compareNumSet, compareNum) # use compare number to record each fragment
        sequenceIndexSet = compareNumSet
        distRecordSet = append(distRecordSet, max_Pearson)
        if (input_form == 0) {
          sequenceSet = append(sequenceSet, sequenceNum)
        }
      }

      outputMatrixAll[[i]] = outputMatrix
      vectorNumSetAll[[i]] = vectorNumSet
      compareNumSetAll[[i]] = compareNumSet
      sequenceIndexSetAll[[i]] = sequenceIndexSet
      distRecordSetAll[[i]] = distRecordSet
      if (input_form == 0) {
        seqListAll[[i]] = sequenceSet
      }

      setTxtProgressBar(pb, i / (m - n + 1))
    }

    # Pearson threshold to filter motifs, remove single results
    single_indices = which(lengths(vectorNumSetAll) == 1)
    if (length(single_indices) > 0) {
      outputMatrixAll = outputMatrixAll[-single_indices]
      vectorNumSetAll = vectorNumSetAll[-single_indices]
      compareNumSetAll = compareNumSetAll[-single_indices]
      sequenceIndexSetAll = sequenceIndexSetAll[-single_indices]
      distRecordSetAll = distRecordSetAll[-single_indices]
      seqListAll = seqListAll[-single_indices]
    }

    means = numeric(length(distRecordSetAll))
    for (i in seq_along(distRecordSetAll)) {
      non_zero_values = distRecordSetAll[[i]][distRecordSetAll[[i]] != 0]
      means[i] = mean(non_zero_values)
    }

    distRecordSetAllmean = means

    # Get sequence when use GRange input
    if (input_form == 1){

      getSeqData <- function(chrName, startPoint, extractMotifLength){
        seq = getSeq(Hsapiens, chrName, start = startPoint, end = startPoint + extractMotifLength - 1)
        seqString = toString(seq)
        return(seqString)
      }

      print("Get alignment sequence ...")
      pb <- txtProgressBar(style = 3)

      seqList = c()
      for (i in 1:length(vectorNumSetAll)){
        vector_index = vectorNumSetAll[[i]]
        sequence_index = sequenceIndexSetAll[[i]]
        seqList = c()
        for (j in 1:length(vector_index)){
          chrName = df_freq_order$chrNo.[vector_index[j]]
          startPoint = df_freq_order$genoStart[vector_index[j]] + sequence_index[j] #- 3 ##########
          seq = getSeqData(chrName, startPoint, extractMotifLength)
          seqList = c(seqList, seq)
        }

        seqListAll[[i]] = seqList

        setTxtProgressBar(pb, i / length(vectorNumSetAll))
      }
    }

    return(
      list(
        "alignedShapeData" = outputMatrixAll,
        "vectorNumSetAll" = vectorNumSetAll,
        "locationindexSetAll" = compareNumSetAll,
        "sequenceIndexSetAll" = sequenceIndexSetAll,
        "distRecordSetAll" = distRecordSetAllmean,
        "seqListAll" = seqListAll
      )
    )

  }

  # Function of getting pearson aligned data frame
  pearson_aligned_df <- function(alignedShapeAll) {

    # Draw motif line data frame
    shapeMatrixSet = areadfmnSet = shapedfmSet = c()
    for (a in 1:length(alignedShapeAll)) {
      alignedShape = alignedShapeAll[[a]]
      colindex = c()

      for (i in 1:nrow(alignedShape)) {
        colindex[i] = paste("Motif", i, sep = "")
      }

      label <- c(1:ncol(alignedShape))
      dimnames = c("label", colindex)
      shapeMatrix = matrix(nrow = nrow(alignedShape) + 1, ncol = ncol(alignedShape))
      shapeMatrix[1, ] = c(1:ncol(alignedShape))

      for (i in 1:nrow(alignedShape)) {
        shapeMatrix[(i + 1), ] = alignedShape[i, ]
      }

      shapedf = as.data.frame(t(shapeMatrix))
      names(shapedf) = dimnames
      shapedfm = reshape2::melt(shapedf, id = "label")
      colnames(shapedfm) = c("label", "Legend", "value")

      shapedfmSet[[a]] = shapedfm
    }

    # Draw motif area data frame
    for (a in 1:length(alignedShapeAll)) {
      alignedShape = alignedShapeAll[[a]]
      minv2 = maxv2 = meanv2 = c()
      meanv2 = apply(alignedShape, 2, mean)
      sd2 = round(apply(alignedShape, 2, sd), 3)
      maxv2 = meanv2 + sd2
      minv2 = meanv2 - sd2
      label = c(1:ncol(alignedShape))

      areadfmn = data.frame(label = label, maxv = maxv2, minv = minv2, meanv = meanv2)

      areadfmnSet[[a]] = areadfmn
    }

    return(list("motif_shape_df_melt" = shapedfmSet, "motif_area_df_melt" = areadfmnSet))
  }

  # Function of motif aligner Edist
  motif_aligner_dist <- function(aindex, extractMotifLength, input_form, sequence_data) {

    # Get compare motifs
    # a: the index of compared motifs
    # m: the total length for the original vector
    # n: total length of reduced motifs
    partialVecListAll = c()
    n = extractMotifLength
    m = length(aindex[[1]])
    for (a in 2:length(aindex)) {
      partialVecSet = c()

      for (i in 1:(m - n + 1)) {
        partialVec = aindex[[a]][i:(i + n - 1)]
        partialVecSet[[i]] = partialVec
      }

      partialVecListAll[[a - 1]] = partialVecSet
    }

    # Get center reference motifs and align
    compareIndexrec = nrec = targetNumrec = compareNumSetAll = sequenceIndexSetAll = seqListAll = c()
    targetVecList = compareMatrixList = targetVecSet = compareMatrixSet = c()
    outputMatrixAll = outputMatrix = c()
    distRecordSet = distRecordSetAll = c()

    if (input_form == 0) {
      print("Shape motif alignment and get sequence ...")
    } else {
      print("Shape motif alignment ...")
    }

    pb <- txtProgressBar(style = 3)

    # Index for each compare vector set
    for (i in 1:(m - n + 1)) {
      targetVec = aindex[[1]][i:(i + n - 1)]

      # Built compare matrix
      outputMatrix = targetVec

      # Index for each compare vector set
      compareNumSet = c(i)

      if (input_form == 0) {
        sequenceSet = c(subseq(sequence_data[[1]], start = compareNumSet, end = (compareNumSet + motifLength - 1)))
      }

      for (b in 1:length(partialVecListAll)) {
        minDist = Inf
        compareMatrix = targetVec

        # Index for all compare vectors in one set
        for (j in 1:length(partialVecListAll[[b]])) {
          if (!any(is.na(partialVecListAll[[b]][[j]]))){
            compareMatrix = unname(rbind(compareMatrix, partialVecListAll[[b]][[j]]))
          }
        }

        # Distance calculation
        Dist = dist(compareMatrix, p = 2)
        # print(Dist)

        for (l in 1:(nrow(compareMatrix) - 1)) {
          tempMinDist = min(Dist[1:l])
          if (tempMinDist < minDist) {
            minDist = tempMinDist
            compareNum = which.min(Dist[1:l])
            if (input_form == 0) {
              sequenceNum = subseq(sequence_data[[b + 1]], start = compareNum, end = (compareNum + motifLength - 1))
            }
          }
        }

        outputMatrix = rbind(outputMatrix, compareMatrix[(compareNum + 1),])
        compareNumSet = append(compareNumSet, compareNum) # use compare number to record each fragments
        sequenceIndexSet = compareNumSet
        distRecordSet = minDist
        if (input_form == 0) {
          sequenceSet = append(sequenceSet, sequenceNum)
        }
      }

      # print(compareNumSet)
      outputMatrixAll[[i]] = outputMatrix
      compareNumSetAll[[i]] = compareNumSet
      sequenceIndexSetAll[[i]] = sequenceIndexSet
      distRecordSetAll[[i]] = distRecordSet
      if (input_form == 0) {
        seqListAll[[i]] = sequenceSet
      }

      setTxtProgressBar(pb, i / (m - n + 1))
    }

    if (input_form == 1){
    # get sequence
      alignedShapeAll = outputMatrixAll
      sequenceIndexAll = sequenceIndexSetAll
      startindexAll = seqListAll = c()

      print("Get alignment sequence ...")
      pb <- txtProgressBar(style = 3)

      for (a in 1:length(alignedShapeAll)) {
        sequenceIndex = sequenceIndexAll[[a]]
        for (dforder in 1:top_motif) {
          chrindex[dforder] = df_freq_order$chrNo.[dforder]
          startindex[dforder] = df_freq_order$genoStart[dforder] + sequenceIndex[dforder]
        }

        locdf = data.frame(no = c(1:nrow(alignedShapeAll[[a]])), chr = chrindex, start = startindex )

        getSeqData <- function(chrName, startPoint, extractMotifLength){
          seq = getSeq(Hsapiens, chrName, start = startPoint, end = startPoint + extractMotifLength - 1)
          # print(seq)
          seqString = toString(seq)
          return(seqString)
        }

        seqList = c()
        for (seqNum in 1:nrow(locdf)) {
          chrName = locdf$chr[seqNum]
          startPoint = locdf$start[seqNum]
          sequenceOutput = getSeqData(chrName, startPoint, extractMotifLength)
          seqList[seqNum] = sequenceOutput
        }
        # print(seqList)

        seqListAll[[a]] = seqList
        startindexAll[[a]] = startindex

        setTxtProgressBar(pb, a / length(alignedShapeAll))
      }
    }

    return(
      list(
        "alignedShapeData" = outputMatrixAll,
        "locationindexSetAll" = compareNumSetAll,
        "sequenceIndexSetAll" = sequenceIndexSetAll,
        "distRecordSetAll" = distRecordSetAll,
        "seqListAll" = seqListAll
      )
    )
  }

  # Function of getting distance aligned data frame
  dist_aligned_df <- function(alignedShapeAll) {

    # Draw motif line data frame
    shapeMatrixSet = areadfmnSet = shapedfmSet = c()
    for (a in 1:length(alignedShapeAll)) {
      alignedShape = alignedShapeAll[[a]]
      colindex = c()
      for (i in 1:nrow(alignedShape)) {
        colindex[i] = paste("Motif", i, sep = "")
      }
      label <- c(1:ncol(alignedShape))
      dimnames = c("label", colindex)
      shapeMatrix <-
        matrix(nrow = nrow(alignedShape) + 1,
               ncol = ncol(alignedShape))
      shapeMatrix[1, ] = c(1:ncol(alignedShape))
      for (i in 1:nrow(alignedShape)) {
        shapeMatrix[(i + 1), ] = alignedShape[i, ]
      }
      shapedf = as.data.frame(t(shapeMatrix))
      names(shapedf) <- dimnames
      shapedfm <- reshape2::melt(shapedf, id = "label")
      colnames(shapedfm) = c("label", "Legend", "value")

      shapedfmSet[[a]] = shapedfm
    }

    # Draw motif area data frame
    for (a in 1:length(alignedShapeAll)) {
      alignedShape = alignedShapeAll[[a]]
      minv2 = maxv2 = meanv2 = c()
      meanv2 = apply(alignedShape, 2, mean)
      sd2 = round(apply(alignedShape, 2, sd), 3)
      maxv2 = meanv2 + sd2
      minv2 = meanv2 - sd2
      label <- c(1:ncol(alignedShape))

      areadfmn <- data.frame(label = label, maxv = maxv2, minv = minv2, meanv = meanv2)

      areadfmnSet[[a]] = areadfmn
    }

    return(list("motif_shape_df_melt" = shapedfmSet, "motif_area_df_melt" = areadfmnSet))
  }

  # Function of getting align score
  align_score_sort <- function(distRecordSetAll) {

    distanceVec = c()
    for (dis in 1:length(distRecordSetAll)) {
      distanceVec = append(distanceVec, distRecordSetAll[[dis]])
    }

    align_index = order(distanceVec)
    align_score = sort(distanceVec)
    cat("\n", "Alignment index for each element: \n", align_index, "\n")
    cat("Pearson/Distance value: \n", align_score, "\n")

    align_score = rbind(align_index, align_score)

    return(align_score)

  }

  if (is.null(input_form)){
    input_form = as.numeric(readline("Select your input form: Enter '0' for 'fasta_motif_merge' or '1' for 'bed_motif_merge': "))
  }

  if (!(input_form %in% c(0, 1))) {
    stop(paste0("Input of this function must be the result from 'fasta_motif_merge' or 'bed_motif_merge'"))
  }

  # Processing the motif data frame
  df_order = motif_data_frame[order(-motif_data_frame$motifFrequency),]
  df_freq_order = df_order[df_order$motifFrequency >= motifFrequency,]

  # Change reference
  a = df_order[1,]
  df_freq_order[1,] = df_freq_order[ref_number,]
  df_freq_order[ref_number,] = a

  # Confirm the parameter setting and continue
  print(paste0("There are '", nrow(df_freq_order), "' results with a frequency greater than '",
              motifFrequency, "' , the reference sequence is ", ref_number))

  if (is.null(continue_alignment)){
    continue_alignment = readline("Do you want to continue with this setting? Enter '1' for 'Yes' or '0' for 'No': ")
  }

  if (continue_alignment == "1") {

    if (is.null(top_motif)){
      top_motif = as.numeric(readline(paste0("Please enter the top_motif number between 0 and ", nrow(df_order), " for alignment : ")))
    }

    # Check if the top_motif is a valid value
    if (top_motif < 0 || top_motif > nrow(df_order)) {
      stop(paste0("Invalid value for top_motif parameter. Please enter a value between 0 and ", nrow(df_order), "."))
    }

    # Prepare motif shape data for alignment
    shape_data = sequence_data = chrindex = startindex = c()

    if (shapeIndex == "all"){

      aindex = cindex = chrindex = startindex = c()
      for (dforder in 1:topsample) {
        aindex[[dforder]] = motifShapeData[[df05order$peakNo.[dforder]]][[df05order$seriesNo.[dforder]]]

        if (input_form == 0) {
          sequence_data[[dfOrder]] = df_freq_order[dfOrder,]$sequence
        }
      }

      print("Use the align method pearson on all features for alignment")

      if (is.null(threshold)){
        threshold = readline("Select a filter threshold between 0 and 1: ")
      }

      # Call motif aligner pearson
      align_results = motif_aligner_pearson(shape_data, motifLength, threshold, input_form, sequence_data)
      alignedShapeAll = align_results$alignedShapeData
      vectorNumSetAll = align_results$vectorNumSetAll
      locationindexAll = align_results$locationindexSetAll
      sequenceIndexAll = align_results$sequenceIndexSetAll
      distRecordSetAll = align_results$distRecordSetAll
      seqListAll = align_results$seqListAll

      align_df = pearson_aligned_df(alignedShapeAll)
      align_score = align_score_sort(distRecordSetAll)

    } else {

      aindex = cindex = chrindex = startindex = c()
      for (dfOrder in 1:top_motif) {
        shape_data[[dfOrder]] = motif_shape_data[[df_freq_order$peakNo.[dfOrder]]][[df_freq_order$seriesNo.[dfOrder]]]
        shape_data[[dfOrder]] = shape_data[[dfOrder]][shapeIndex,]

        if (input_form == 0) {
          sequence_data[[dfOrder]] = df_freq_order[dfOrder,]$sequence
        }
      }

      # Call motif aligner
      if (align_method == "pearson"){

        print("Use the align method pearson for alignment")

        if (is.null(threshold)){
          threshold = readline("Select a filter threshold between 0 and 1: ")
        }

        # Call motif aligner pearson
        align_results = motif_aligner_pearson(shape_data, motifLength, threshold, input_form, sequence_data)
        alignedShapeAll = align_results$alignedShapeData
        vectorNumSetAll = align_results$vectorNumSetAll
        locationindexAll = align_results$locationindexSetAll
        sequenceIndexAll = align_results$sequenceIndexSetAll
        distRecordSetAll = align_results$distRecordSetAll
        seqListAll = align_results$seqListAll

        align_df = pearson_aligned_df(alignedShapeAll)
        align_score = align_score_sort(distRecordSetAll)

      } else if (align_method == "distance") {

        print("Use the align method distance for alignment")
        Sys.sleep(1.5)

        # Call motif aligner dist
        align_results = motif_aligner_dist(shape_data, motifLength, input_form, sequence_data)
        alignedShapeAll = align_results$alignedShapeData
        locationindexAll = align_results$locationindexSetAll
        sequenceIndexAll = align_results$sequenceIndexSetAll
        distRecordSetAll = align_results$distRecordSetAll
        seqListAll = align_results$seqListAll

        # Get aligned data frame
        align_df = dist_aligned_df(alignedShapeAll)
        align_score = align_score_sort(distRecordSetAll)

      } else {
        stop("Invalid value for 'align_method'. Please select 'pearson' or 'distance'.")
      }
    }

  } else {
    stop("Reset the 'motifFrequency' from 0 to 1.")
  }

  return(list("align_score" = align_score, "align_motif_dataframe" = align_df, "align_motif_sequence" = seqListAll))

}


# vis_test = source("merge_test")$value
# motif_data_frame = vis_test$motif_data_frame
# motif_shape_data = vis_test$motif_shape_data
#
# fasta_vis_test = source("fa_merge_test")$value
# motif_data_frame = fasta_vis_test$motif_data_frame
# motif_shape_data = fasta_vis_test$motif_shape_data
#
# test_v = motif_visualizer(motif_data_frame, motif_shape_data, motifFrequency = 0.7, shapeIndex = 1,
#                           motifLength = 10, align_method = "distance", ref_number = 1)
#
# motifFrequency = 0.7
# ref_number = 1
# shapeIndex = 1
# motifLength = extractMotifLength = 10
# top_motif = 20
# input_form = 0
#
# aindex = shape_data























