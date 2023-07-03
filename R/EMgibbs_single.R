EMgibbs_single <- function(inputData, peakCount, K, inputMotifLength, tol) {
  # inputData = simulatedData
  # K = timeSeriesMotifCount
  # inputMotifLength = motifLength
  # tol = 0.000001

  numOfFeatures = 1
  shapeFeatureCount = numOfFeatures
  pseudoPosteriorCount = 0

  # Initialization
  initialTimeSeriesMotifs = list()

  for (k in 1:K) {
    randomSeqIndex = sample(1:length(inputData), 1)
    randomSeqPosition = sample(1:(length(inputData[[randomSeqIndex]]) - inputMotifLength +
                                    1), 1)
    initialTimeSeriesMotifs[[k]] = inputData[[randomSeqIndex]][randomSeqPosition:(randomSeqPosition +
                                                                                      inputMotifLength - 1)]
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
          probb[j] = prod(dnorm(tempMotif, mean = currentTimeSeriesMotifs[[k]], sd =
                                  BGsd)) / prod(dnorm(tempMotif, BGmeans, BGsd))
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
    prior[is.na(colSumMask) |
            colSumMask == 0] = 1 / K  # Numerical Underflow Handling


    #### M-step ######
    pi = rowSums(prior) / peakCount
    tempTimeSeriesMotifs = list()
    for (k in 1:K) {
      tempTimeSeriesMotifs[[k]] = matrix(0, nrow = shapeFeatureCount, ncol = inputMotifLength)
      for (i in 1:peakCount) {
        #prior[k,i]
        currentMotif = inputData[[i]][motifLocations[k, i]:(motifLocations[k, i] +
                                                                inputMotifLength - 1)]
        tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] + currentMotif *
          prior[k, i]
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
        if (shifttedLocations[i] >= 1 &
            (shifttedLocations[i] + inputMotifLength - 1) <= length(inputData[[i]])) {
          newTimeSeriesMotif = newTimeSeriesMotif + prior[k, i] * inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] +
                                                                                                           inputMotifLength - 1)]
          actualCount = actualCount + prior[k, i]
        }
      }
      newTimeSeriesMotif = newTimeSeriesMotif / actualCount

      l_totalProb = 0
      for (i in 1:peakCount) {
        allBG = sum(log10(dnorm(inputData[[i]], BGmeans, BGsd)))
        if (shifttedLocations[i] >= 1 &
            (shifttedLocations[i] + inputMotifLength - 1) <= length(inputData[[i]])) {
          l_totalProb = l_totalProb +  prior[k, i] * 10 ^ (sum(log10(
            dnorm(inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength -
                                                           1)], mean = newTimeSeriesMotif, sd = BGsd)
          )) - sum(log10(
            dnorm(inputData[[i]][shifttedLocations[i]:(shifttedLocations[i] + inputMotifLength -
                                                           1)], mean = BGmeans, sd = BGsd)
          )))
        } else{
          l_totalProb = l_totalProb + prior[k, i]
        }
      }
      if (l_totalProb > maxLProb) {
        finalLocations[k, ] = shifttedLocations
        finalTimeSeriesMotifs[[k]] = newTimeSeriesMotif
        maxLProb = l_totalProb
      }

      # averagePD = 0
      # actualCount = 0
      # for(i in 1:peakCount){
      #   if(shifttedLocations[i]>=1 & (shifttedLocations[i]+inputMotifLength-1)<=length(inputData[[i]])){
      #     averagePD = averagePD +  prior[k,i]*sum(abs(inputData[[i]][,shifttedLocations[i]:(shifttedLocations[i]+inputMotifLength-1)]-newTimeSeriesMotif)/abs(newTimeSeriesMotif))/prod(dim(newTimeSeriesMotif))
      #     actualCount = actualCount + prior[k,i]
      #   }
      # }
      # if(actualCount>=1){
      #   averagePD = averagePD / actualCount
      # }else{
      #   averagePD = Inf
      # }
      # if(averagePD<minPD){
      #   finalLocations[k,] = shifttedLocations
      #   finalTimeSeriesMotifs[[k]] = newTimeSeriesMotif
      #   minPD = averagePD
      # }


    }
  }

  #print(finalTimeSeriesMotifs)
  #print(timeSeriesMotifs)
  return(
    list(
      "finalTimeSeriesMotifs" = finalTimeSeriesMotifs,
      "finalLocations" = finalLocations
    )
  )
}
