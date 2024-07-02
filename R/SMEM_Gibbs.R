#' DNA shape motif discovery with Expectation–Maximization algorithm
#'
#' This function reads the DNA shape data file and conduct shape motif discovery using the Expectation–Maximization algorithm and Gibbs sampling.
#'
#' @param dna_shape_data Input DNA shape data convert by function fasta_input or bed_input.
#' @param filename character, name of output .csv file.
#' @param peakCount integer, the number of peaks to retrieve per batch, default is 20.
#' @param motifLength integer, length of shape motif to be discovered, default is 12.
#' @param motifCount integer, number of motifs to be discovered for each peak, default is 2.
#' @param replicates integer, number of repetitions for each batch, default is 50.
#' @param tolerance double, tolerance threshold for EM iteration, default is 0.000001.
#' @return A list of motif location arrays
#' @importFrom stats cor dist dnorm rnorm sd
#' @importFrom utils read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @export SMEM_Gibbs
SMEM_Gibbs <- function(dna_shape_data, filename = "file", peakCount = 20, motifLength = 12, motifCount = 2, replicates = 50, tolerance = 0.000001) {

  # Function of the evaluation metric
  evaluateMatrices <- function(timeSeriesMotifs, finalTimeSeriesMotifs) {

    EDists = rep(Inf, length(timeSeriesMotifs))
    Pearsons = rep(-Inf, length(timeSeriesMotifs))
    MAPEs = rep(Inf, length(timeSeriesMotifs))
    for (k in 1:length(timeSeriesMotifs)) {
      for (k2 in 1:length(finalTimeSeriesMotifs)) {
        EDists[k] = min(EDists[k], sqrt(sum((timeSeriesMotifs[[k]] - finalTimeSeriesMotifs[[k2]]) ^ 2
        )) / prod(dim(timeSeriesMotifs[[k]])))

        tempPearson = 0
        for (i in 1:nrow(timeSeriesMotifs[[k]])) {
          tempPearson = tempPearson + cor(timeSeriesMotifs[[k]][i,], finalTimeSeriesMotifs[[k2]][i,],
                                          method = 'pearson')
        }
        Pearsons[k] = max(Pearsons[k], tempPearson / nrow(timeSeriesMotifs[[k]]))

        MAPEs[k] = min(MAPEs[k], sum(
          abs(timeSeriesMotifs[[k]] - finalTimeSeriesMotifs[[k2]]) / abs(timeSeriesMotifs[[k]])
        ) / prod(dim(timeSeriesMotifs[[k]])))
      }
    }
    print(sum(EDists) / length(EDists))
    print(sum(Pearsons) / length(Pearsons))
    print(sum(MAPEs) / length(MAPEs))

    return(list("EDists" = EDists, "Pearsons" = Pearsons, "MAPEs" = MAPEs))
  }

  # Function of the EM_Gbbis algorithm for single shape feature
  EMgibbs_single <- function(inputData, peakCount, K, inputMotifLength, tol) {

    numOfFeatures = 1
    shapeFeatureCount = numOfFeatures
    pseudoPosteriorCount = 0

    # Initialization
    initialTimeSeriesMotifs = list()

    for (k in 1:K) {
      randomSeqIndex = sample(1:length(inputData), 1)
      randomSeqPosition = sample(1:(length(inputData[[randomSeqIndex]]) - inputMotifLength + 1), 1)
      initialTimeSeriesMotifs[[k]] = inputData[[randomSeqIndex]][randomSeqPosition:(randomSeqPosition + inputMotifLength - 1)]
    }

    # Loop Modified
    minTol = Inf
    maxIterations = 10 * length(inputData[[1]])
    currentTol = 1
    iter = 1
    currentTimeSeriesMotifs = initialTimeSeriesMotifs

    pi = rep(1 / K, K)
    prior = matrix(0, nrow = K, ncol = peakCount)
    posterior = list()
    for (k in 1:K) {
      posterior[[k]] = list()
      for (i in 1:peakCount) {
        peakLen = length(inputData[[i]])
        posterior[[k]][[i]] = rep(1 / peakLen, peakLen)
      }
    }
    tempData = c()
    for (i in 1:peakCount) {
      tempData = cbind(tempData, inputData[[i]])
    }
    BGmeans = rowMeans(tempData)
    BGsd = apply(tempData, 1, sd)

    while (currentTol > tol && iter < maxIterations) {
      print(sprintf("Iteration = %d  | currentTol = %0.10f", iter, currentTol))

      #### E-step ######
      motifLocations = matrix(0, nrow = K, ncol = peakCount)
      for (k in 1:K) {
        for (i in 1:peakCount) {
          prior[k, i] = 0
          probb = rep(0, (length(inputData[[i]]) - inputMotifLength + 1))
          for (j in 1:(length(inputData[[i]]) - inputMotifLength + 1)) {
            tempMotif = inputData[[i]][j:(j + inputMotifLength - 1)]
            probb[j] = prod(dnorm(tempMotif, mean = currentTimeSeriesMotifs[[k]], sd = BGsd)) / prod(dnorm(tempMotif, BGmeans, BGsd))
          }
          prior[k, i] = sum(probb)
          prior[k, i] = prior[k, i] * pi[k] # *dnorm(inputData[[i]],BGmeans,BGsd)
          if (sum(probb) == 0 || is.na(sum(probb))) {
            probb = rep(1 / length(probb), length(probb))  # Numerical Underflow Handling
          }
          probb = probb / sum(probb)
          motifLocations[k, i] = sample(1:length(probb), 1, prob = probb)
        }
      }
      colSumMask = matrix(rep(colSums(prior), K), nrow = K, byrow = TRUE)
      prior = prior / colSumMask
      prior[is.na(colSumMask) | colSumMask == 0] = 1 / K  # Numerical Underflow Handling

      #### M-step ######
      pi = rowSums(prior) / peakCount
      tempTimeSeriesMotifs = list()
      for (k in 1:K) {
        tempTimeSeriesMotifs[[k]] = matrix(0, nrow = shapeFeatureCount, ncol = inputMotifLength)
        for (i in 1:peakCount) {

          currentMotif = inputData[[i]][motifLocations[k, i]:(motifLocations[k, i] + inputMotifLength - 1)]
          tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] + currentMotif * prior[k, i]
        }
        tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] / sum(prior[k, ])
      }

      currentTol = 0
      for (k in 1:K) {
        currentTol = currentTol + sum(abs(tempTimeSeriesMotifs[[k]] - currentTimeSeriesMotifs[[k]])) / sum(abs(currentTimeSeriesMotifs[[k]]))
      }
      currentTol = currentTol / K
      currentTimeSeriesMotifs = tempTimeSeriesMotifs
      iter = iter + 1

      if (is.na(currentTol)) {
        return()
      }

      if (currentTol < minTol) {
        minIter = iter
        minTol = currentTol
        minTempTimeSeriesMotifs = tempTimeSeriesMotifs
        minPrior = prior
        minMotifLocation = motifLocations
      }
    }

    if (iter == maxIterations & currentTol != minTol) {
      iter = minIter
      currentTol = minTol
      tempTimeSeriesMotifs = minTempTimeSeriesMotifs
      prior = minPrior
      motifLocations = minMotifLocation
    }

    print(sprintf("Miniteration = %d  | minTol = %0.10f", iter, currentTol))

    # The Phase Shift Problem
    finalLocationsInitial = matrix(0, nrow = K, ncol = peakCount)
    for (k in 1:K) {
      for (i in 1:peakCount) {
        finalLocationsInitial[k, i] = motifLocations[k, i]
      }
    }
    finalLocations = matrix(0, nrow = K, ncol = peakCount)
    finalTimeSeriesMotifs = list()
    for (k in 1:K) {
      minPD = Inf
      maxLProb = -Inf
      for (shift in (-inputMotifLength + 1):(inputMotifLength - 1)) {
        shifttedLocations = finalLocationsInitial[k, ] + shift
        newTimeSeriesMotif = matrix(0, nrow = shapeFeatureCount, ncol = inputMotifLength)
        actualCount = 0
        for (i in 1:peakCount) {
          if (shifttedLocations[i] >= 1 & (shifttedLocations[i] + inputMotifLength - 1) <= length(inputData[[i]])) {
            newTimeSeriesMotif = newTimeSeriesMotif + prior[k, i] * inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)]
            actualCount = actualCount + prior[k, i]
          }
        }
        newTimeSeriesMotif = newTimeSeriesMotif / actualCount

        l_totalProb = 0
        for (i in 1:peakCount) {
          allBG = sum(log10(dnorm(inputData[[i]], BGmeans, BGsd)))
          if (shifttedLocations[i] >= 1 & (shifttedLocations[i] + inputMotifLength - 1) <= length(inputData[[i]])) {
            l_totalProb = l_totalProb +  prior[k, i] * 10 ^ (sum(log10(
              dnorm(inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)], mean = newTimeSeriesMotif, sd = BGsd)))
              - sum(log10(dnorm(inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)], mean = BGmeans, sd = BGsd))))
          } else {
            l_totalProb = l_totalProb + prior[k, i]
          }
        }
        if (l_totalProb > maxLProb) {
          finalLocations[k, ] = shifttedLocations
          finalTimeSeriesMotifs[[k]] = newTimeSeriesMotif
          maxLProb = l_totalProb
        }
      }
    }

    return(list("finalTimeSeriesMotifs" = finalTimeSeriesMotifs, "finalLocations" = finalLocations))
  }

  # Function of the EM_Gbbis algorithm for multiple shape features
  EMgibbs_mod <- function(inputData, peakCount, K, inputMotifLength, tol) {

      numOfFeatures = nrow(inputData[[1]])
      shapeFeatureCount = numOfFeatures
      pseudoPosteriorCount = 0

      # Initialization
      initialTimeSeriesMotifs = list()

      for (k in 1:K) {
        randomSeqIndex = sample(1:length(inputData), 1)
        randomSeqPosition = sample(1:(ncol(inputData[[randomSeqIndex]]) - inputMotifLength + 1), 1)
        initialTimeSeriesMotifs[[k]] = inputData[[randomSeqIndex]][, randomSeqPosition:(randomSeqPosition + inputMotifLength - 1)]
      }

      # Loop Modified
      minTol = Inf
      maxIterations = 10 * ncol(inputData[[1]])
      currentTol = 1
      iter = 1
      currentTimeSeriesMotifs = initialTimeSeriesMotifs

      pi = rep(1 / K, K)
      prior = matrix(0, nrow = K, ncol = peakCount)
      posterior = list()
      for (k in 1:K) {
        posterior[[k]] = list()
        for (i in 1:peakCount) {
          peakLen = ncol(inputData[[i]])
          posterior[[k]][[i]] = rep(1 / peakLen, peakLen)
        }
      }
      tempData = c()
      for (i in 1:peakCount) {
        tempData = cbind(tempData, inputData[[i]])
      }
      BGmeans = rowMeans(tempData)
      BGsd = apply(tempData, 1, sd)

      while (currentTol > tol && iter < maxIterations) {
        print(sprintf("Iteration = %d  | currentTol = %0.10f", iter, currentTol))

        #### E-step ######
        motifLocations = matrix(0, nrow = K, ncol = peakCount)
        for (k in 1:K) {
          for (i in 1:peakCount) {
            prior[k, i] = 0
            probb = rep(0, (ncol(inputData[[i]]) - inputMotifLength + 1))
            for (j in 1:(ncol(inputData[[i]]) - inputMotifLength + 1)) {
              tempMotif = inputData[[i]][, j:(j + inputMotifLength - 1)]
              probb[j] = prod(dnorm(tempMotif, mean = currentTimeSeriesMotifs[[k]], sd = BGsd)) / prod(dnorm(tempMotif, BGmeans, BGsd))
            }
            prior[k, i] = sum(probb)
            prior[k, i] = prior[k, i] * pi[k]
            if (sum(probb) == 0 || is.na(sum(probb))) {
              probb = rep(1 / length(probb), length(probb))
            }
            probb = probb / sum(probb)
            motifLocations[k, i] = sample(1:length(probb), 1, prob = probb)
          }
        }
        colSumMask = matrix(rep(colSums(prior), K), nrow = K, byrow = TRUE)
        prior = prior / colSumMask
        prior[is.na(colSumMask) |
                colSumMask == 0] = 1 / K

        #### M-step ######
        pi = rowSums(prior) / peakCount
        tempTimeSeriesMotifs = list()
        for (k in 1:K) {
          tempTimeSeriesMotifs[[k]] = matrix(0, nrow = shapeFeatureCount, ncol = inputMotifLength)
          for (i in 1:peakCount) {
            currentMotif = inputData[[i]][, motifLocations[k, i]:(motifLocations[k, i] + inputMotifLength - 1)]
            tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] + currentMotif * prior[k, i]
          }
          tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] / sum(prior[k,])
        }

        currentTol = 0
        for (k in 1:K) {
          currentTol = currentTol + sum(abs(tempTimeSeriesMotifs[[k]] - currentTimeSeriesMotifs[[k]])) / sum(abs(currentTimeSeriesMotifs[[k]]))
        }
        currentTol = currentTol / K

        if (currentTol < minTol) {
          minIter = iter
          minTol = currentTol
          minTempTimeSeriesMotifs = tempTimeSeriesMotifs
          minPrior = prior
          minMotifLocation = motifLocations
        }

        currentTimeSeriesMotifs = tempTimeSeriesMotifs
        iter = iter + 1
        if (is.na(currentTol)) {
          return()
        }
      }

      if (iter == maxIterations & currentTol > minTol) {
        iter = minIter
        currentTol = minTol
        tempTimeSeriesMotifs = minTempTimeSeriesMotifs
        prior = minPrior
        motifLocations = minMotifLocation
      }

      print(sprintf("Miniteration = %d  | minTol = %0.10f", iter, currentTol))

      # The Phase Shift Problem
      finalLocationsInitial = matrix(0, nrow = K, ncol = peakCount)
      for (k in 1:K) {
        for (i in 1:peakCount) {
          finalLocationsInitial[k, i] = motifLocations[k, i]
        }
      }
      finalLocations = matrix(0, nrow = K, ncol = peakCount)
      finalTimeSeriesMotifs = list()
      for (k in 1:K) {
        minPD = Inf
        maxLProb = -Inf
        for (shift in (-inputMotifLength + 1):(inputMotifLength - 1)) {
          shifttedLocations = finalLocationsInitial[k,] + shift
          newTimeSeriesMotif = matrix(0, nrow = shapeFeatureCount, ncol = inputMotifLength)
          actualCount = 0
          for (i in 1:peakCount) {
            if (shifttedLocations[i] >= 1 & (shifttedLocations[i] + inputMotifLength - 1) <= ncol(inputData[[i]])) {
              newTimeSeriesMotif = newTimeSeriesMotif + prior[k, i] * inputData[[i]][, shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)]
              actualCount = actualCount + prior[k, i]
            }
          }
          newTimeSeriesMotif = newTimeSeriesMotif / actualCount

          l_totalProb = 0
          for (i in 1:peakCount) {
            allBG = sum(log10(dnorm(inputData[[i]], BGmeans, BGsd)))
            if (shifttedLocations[i] >= 1 &
                (shifttedLocations[i] + inputMotifLength - 1) <= ncol(inputData[[i]])) {
              l_totalProb = l_totalProb +  prior[k, i] * 10 ^ (sum(log10(
                dnorm(inputData[[i]][, shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)], mean = newTimeSeriesMotif, sd = BGsd)))
                - sum(log10(dnorm(inputData[[i]][, shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength - 1)], mean = BGmeans, sd = BGsd))))
            } else{
              l_totalProb = l_totalProb + prior[k, i]
            }
          }
          if (l_totalProb > maxLProb) {
            finalLocations[k,] = shifttedLocations
            finalTimeSeriesMotifs[[k]] = newTimeSeriesMotif
            maxLProb = l_totalProb
          }
        }
      }

      return(list("finalTimeSeriesMotifs" = finalTimeSeriesMotifs, "finalLocations" = finalLocations))
    }

  # Data Integration and Initialization
  cycleCount = reduce_column = 0
  experimentData = list()
  finalLocationArray = c()
  timeSeriesMotifCount = motifCount
  NumOfSDs   = dim(dna_shape_data)[1]
  peakLength = dim(dna_shape_data)[2]
  lenGrSet   = dim(dna_shape_data)[3]

  # Check that the peak count is available
  if (peakCount > lenGrSet) {
    stop(paste("peakCount", peakCount, "larger than total data length", lenGrSet))
  }

  # Check that the motif length is available
  if (motifLength > dim(dna_shape_data)[2]) {
    stop(paste("motifLength", motifLength, "should not larger than peak length", dim(dna_shape_data)[2]))
  }

  # Check that the peak count is suitable
  if (lenGrSet %% peakCount != 0) {
    message(paste0("peakCount ", peakCount, " is not divisible by total peak number ", lenGrSet,". Consider choosing a suitable peakCount."))
    Sys.sleep(3)
  }

  # Check if motif discovery is conduct on a single shape
  if (NumOfSDs > 1){
    print(paste("Start shape motif discovery on", NumOfSDs, "shape features"))
    Sys.sleep(2)

    # Discover shape motif according to peakCount
    while ((peakCount * cycleCount) < lenGrSet) {

      # Start point counter
      startPoint = cycleCount * peakCount + 1
      endPoint = cycleCount * peakCount + peakCount
      if (endPoint > lenGrSet){endPoint = lenGrSet}

      # Peak number initialization
      peakNumber = 1

      # Motif discovery according to peak count
      for (x in startPoint:endPoint) {
        experimentData[[peakNumber]] <- dna_shape_data[, , x]
        peakNumber = peakNumber + 1
      }

      # Output raw location file
      outputFinalLocation = paste("Location_", filename, "_Top_", lenGrSet, "_", startPoint, "_to_", endPoint, ".csv", sep = "")
      fileConnFinalLocation = file(outputFinalLocation, open = "wt")

      # Output evaluation metrics file
      outputFileName = paste("Experiment_", filename, "_Top_", lenGrSet, "_", startPoint, "_to_", endPoint, ".csv", sep = "")
      fileConn = file(outputFileName, open = "wt")
      writeLines(
        sprintf(
          "PeakCount, PeakLength, MotifCount, MotifLength, NumOfSDs, Method, Time, Edist, Pearson, MAPE"
        ),
        fileConn
      )

      tempLocationArray = c()
      # EM for motif discovery according to the set number of repetitions
      for (r in 1:replicates) {

        # Get the background information of current data batch
        tempPredData <- c()
        processBar <- txtProgressBar(style = 3)

        for (y in 1:peakCount) {
          tempPredData <- cbind(tempPredData, experimentData[[y]])
          setTxtProgressBar(processBar, y / peakCount)
        }

        close(processBar)
        NoiseMean <- rowMeans(tempPredData)
        NoiseSD <- apply(tempPredData, 1, sd)

        shapeFeatureCount = length(NoiseMean)
        timeSeriesMotifs = list()

        for (c in 1:timeSeriesMotifCount) {
          tempTimeSeriesMotifMatrix = c()

          for (f in 1:shapeFeatureCount) {
            tempTimeSeriesMotifMatrix = rbind(tempTimeSeriesMotifMatrix, rnorm(motifLength, mean = NoiseMean[f], sd = NoiseSD[f]))
          }

          timeSeriesMotifs[[c]] =  tempTimeSeriesMotifMatrix
        }

        # Set the tolerance threshold for iteration
        tol = tolerance

        # Apply EM algorithm for motif discovery
        possibleError <- tryCatch({
          start_time <- Sys.time()
          results = EMgibbs_mod(experimentData, peakCount, timeSeriesMotifCount, motifLength, tol)
          timeTaken = as.numeric(Sys.time()) - as.numeric(start_time)
          evaluated = evaluateMatrices(timeSeriesMotifs, results$finalTimeSeriesMotifs)

          # Output evaluation metrics
          writeLines(
            sprintf(
              "%d, %d, %d, %d, %d, EM, %f,  %f,  %f,  %f",
              peakCount, peakLength, timeSeriesMotifCount, motifLength, NumOfSDs, timeTaken,
              mean(evaluated$EDists),
              mean(evaluated$Pearsons),
              mean(evaluated$MAPEs)
            ), fileConn)
        },
        error = function(e)
          e,
        finally = {
        })

        # Record the motif location for each peak
        tempLocationArray = rbind(tempLocationArray, results$finalLocations)
      }

      # Output location
      write.csv(tempLocationArray, fileConnFinalLocation)

      # Output the motif location
      finalLocationArray[[cycleCount + 1]] = tempLocationArray

      cycleCount = cycleCount + 1

      # close file connection
      close(fileConnFinalLocation)
      close(fileConn)
    }

  }else if (NumOfSDs == 1) {
    print(paste("Start shape motif discovery on single shape features"))
    Sys.sleep(2)

    # Discover shape motif according to peakCount
    while ((peakCount * cycleCount) < lenGrSet) {

      # Start point counter
      startPoint = cycleCount * peakCount + 1
      endPoint = cycleCount * peakCount + peakCount
      if (endPoint > lenGrSet){endPoint = lenGrSet}

      # Peak number initialization
      peakNumber = 1

      # Motif discovery according to peak count
      for (x in startPoint:endPoint) {
        experimentData[[peakNumber]] <- dna_shape_data[, , x]
        peakNumber = peakNumber + 1
      }

      # Output raw location file
      outputFinalLocation = paste("Location_", filename, "_Top_", lenGrSet, "_", startPoint, "_to_", endPoint, "EM_Gibbs", ".csv", sep = "")
      fileConnFinalLocation = file(outputFinalLocation, open = "wt")

      # Output evaluation metrics file
      outputFileName = paste("Experiment_", filename, "_Top_", lenGrSet, "_", startPoint, "_to_", endPoint, "EM_Gibbs", ".csv", sep = "")
      fileConn = file(outputFileName, open = "wt")
      writeLines(
        sprintf(
          "PeakCount, PeakLength, MotifCount, MotifLength, NumOfSDs, Method, Time, Edist, Pearson, MAPE"
        ),
        fileConn
      )

      tempLocationArray = c()
      # EM for motif discovery according to the set number of repetitions
      for (r in 1:replicates) {

        # Get the background information of current data batch
        tempPredData <- c()
        processBar <- txtProgressBar(style = 3)

        for (y in 1:peakCount) {
          tempPredData <- cbind(tempPredData, experimentData[[y]])
          setTxtProgressBar(processBar, y / peakCount)
        }

        close(processBar)
        NoiseMean <- mean(tempPredData)
        NoiseSD <- sd(tempPredData)

        shapeFeatureCount = length(NoiseMean)
        timeSeriesMotifs = list()

        for (c in 1:timeSeriesMotifCount) {
          tempTimeSeriesMotifMatrix = c()

          for (f in 1:shapeFeatureCount) {
            tempTimeSeriesMotifMatrix = rbind(tempTimeSeriesMotifMatrix, rnorm(motifLength, mean = NoiseMean[f], sd = NoiseSD[f]))
          }

          timeSeriesMotifs[[c]] =  tempTimeSeriesMotifMatrix
        }

        # Set the tolerance threshold for iteration
        tol = tolerance

        # Apply EM algorithm for motif discovery
        possibleError <- tryCatch({
          start_time <- Sys.time()
          results = EMgibbs_single(experimentData, peakCount, timeSeriesMotifCount, motifLength, tol)
          timeTaken = as.numeric(Sys.time()) - as.numeric(start_time)
          evaluated = evaluateMatrices(timeSeriesMotifs, results$finalTimeSeriesMotifs)

          # Output location
          write.csv(results$finalLocations, fileConnFinalLocation)

          # Output evaluation metrics
          writeLines(
            sprintf(
              "%d, %d, %d, %d, %d, EM, %f,  %f,  %f,  %f",
              peakCount, peakLength, timeSeriesMotifCount, motifLength, NumOfSDs, timeTaken,
              mean(evaluated$EDists),
              mean(evaluated$Pearsons),
              mean(evaluated$MAPEs)
            ),
            fileConn
          )
        },
        error = function(e)
          e,
        finally = {
        })

        # Record the motif location for each peak
        tempLocationArray = rbind(tempLocationArray, results$finalLocations)
      }

      # Output the motif location
      finalLocationArray[[cycleCount + 1]] = tempLocationArray

      cycleCount = cycleCount + 1

      # close file connection
      close(fileConnFinalLocation)
      close(fileConn)
    }

  } else {
    stop("Number of shape features must be a positive integer")
  }

  # Merge all the location results and output csv
  finalLocationArray = do.call(cbind, finalLocationArray)[, 1:lenGrSet]
  LocationArrayFilename = paste("Location_", filename, "_Top_", lenGrSet, "_1_to_", endPoint, "EM_Gibbs", ".csv", sep = "")
  write.csv(finalLocationArray, file = LocationArrayFilename,
            row.names = c(paste0("rep", 1:(motifCount * replicates))))

  return("finalLocationArray" = finalLocationArray)

}


# test = source("single_test_shape")
# test = test$value
#
# em_test = SMEM(test, filename = "file", peakCount = 22, motifLength = 10, motifCount = 2, replicates = 5, tolerance = 0.000001)

# filename = "file"
# peakCount = 20
# motifLength = 12
# motifCount = 2
# replicates = 5
# tolerance = 0.000001
# dna_shape_data = test
