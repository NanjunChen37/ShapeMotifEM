Gibbs <- function(inputData, K, inputMotifLength, tol){

  numOfFeatures = nrow(inputData[[1]])
  shapeFeatureCount = numOfFeatures
  pseudoPosteriorCount = 0

  # Initialization
  initialTimeSeriesMotifs = list()
  for(k in 1:K){
    randomSeqIndex = sample(1:length(inputData),1)
    randomSeqPosition = sample(1:(ncol(inputData[[randomSeqIndex]])-inputMotifLength+1),1)
    initialTimeSeriesMotifs[[k]] = inputData[[randomSeqIndex]][,randomSeqPosition:(randomSeqPosition+inputMotifLength-1)]
  }

  # Loop
  maxIterations = 10*ncol(inputData[[1]])
  currentTol = 1
  iter = 1
  currentTimeSeriesMotifs = initialTimeSeriesMotifs

  pi = rep(1/K,K)
  prior = matrix(0,nrow=K,ncol=peakCount)
  posterior = list()
  for(k in 1:K){
    posterior[[k]] = list()
    for(i in 1:peakCount){
      peakLen = ncol(inputData[[i]])
      posterior[[k]][[i]] = rep(1/peakLen, peakLen)
    }
  }
  tempData = c()
  for(i in 1:peakCount){
    tempData = cbind(tempData,inputData[[i]])
  }
  BGmeans = rowMeans(tempData)
  BGsd = apply(tempData,1,sd)


  while( currentTol>tol && iter<maxIterations){

    print(sprintf("Iteration = %d  | currentTol = %0.10f",iter, currentTol))

    #### E-step ######
    motifLocations = matrix(0,nrow=K, ncol=peakCount)
    for(i in 1:peakCount){
      for(k in 1:K){
        prior[k,i] = 0
        probb = rep(0,(ncol(inputData[[i]])-inputMotifLength+1))
        for(j in 1:(ncol(inputData[[i]])-inputMotifLength+1)){
          tempMotif = inputData[[i]][,j:(j+inputMotifLength-1)]
          probb[j] = prod(dnorm(tempMotif, mean=currentTimeSeriesMotifs[[k]], sd=BGsd)) / prod(dnorm(tempMotif,BGmeans,BGsd))
        }
        prior[k,i] = sum(probb)
        prior[k,i] = prior[k,i]*pi[k] # *dnorm(inputData[[i]],BGmeans,BGsd)
        if(sum(probb)==0 || is.na(sum(probb))){
          probb = rep(1/length(probb),length(probb))  # Numerical Underflow Handling
        }
        probb = probb / sum(probb)
        motifLocations[k,i] = sample(1:length(probb), 1, prob=probb)
      }
    }
    colSumMask = matrix(rep(colSums(prior),K),nrow=K,byrow=TRUE)
    prior = prior / colSumMask
    prior[is.na(colSumMask) | colSumMask==0] = 1/K  # Numerical Underflow Handling

    newMotifLocations = list()
    for(k in 1:K){
      newMotifLocations[[k]] = list()
    }
    for(i in 1:peakCount){
      k = which.max(prior[,i])
      newMotifLocations[[k]] = rbind(newMotifLocations[[k]], list(i,motifLocations[k,i]))
    }

    #### M-step ######
    pi= rowSums(prior)/peakCount
    tempTimeSeriesMotifs = list()
    for(k in 1:K){
      tempTimeSeriesMotifs[[k]] = matrix(0, nrow=shapeFeatureCount, ncol=inputMotifLength)

      summ = 0

      # Reset the motif without any support
      if(length(newMotifLocations[[k]])==0){
        randomSeqIndex = sample(1:length(inputData),1)
        randomSeqPosition = sample(1:(ncol(inputData[[randomSeqIndex]])-inputMotifLength+1),1)
        tempTimeSeriesMotifs[[k]] = inputData[[randomSeqIndex]][,randomSeqPosition:(randomSeqPosition+inputMotifLength-1)]
        next
      }
      for(xx in 1:nrow(newMotifLocations[[k]])){
        bb = newMotifLocations[[k]][xx,]
        i = bb[[1]]
        loc = bb[[2]]
        currentMotif = inputData[[i]][,loc:(loc+inputMotifLength-1)]
        tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] + currentMotif*prior[k,i]
        summ = summ + prior[k,i]
      }
      tempTimeSeriesMotifs[[k]] = tempTimeSeriesMotifs[[k]] / summ
    }

    currentTol = 0
    for(k in 1:K){
      currentTol = currentTol + sum(abs(tempTimeSeriesMotifs[[k]]-currentTimeSeriesMotifs[[k]])) / sum(abs(currentTimeSeriesMotifs[[k]]))
    }
    currentTol = currentTol / K
    currentTimeSeriesMotifs = tempTimeSeriesMotifs
    iter = iter + 1

    if(is.na(currentTol)){
      return()
    }
  }

  # The Phase Shift Problem
  finalLocationsInitial = matrix(0,nrow=K, ncol=peakCount)
  for(k in 1:K){
    for(i in 1:peakCount){
      finalLocationsInitial[k,i] = motifLocations[k,i]
    }
  }
  finalLocations = matrix(0,nrow=K, ncol=peakCount)
  finalTimeSeriesMotifs = list()
  for(k in 1:K){
    minPD = Inf
    maxLProb = -Inf
    finalTimeSeriesMotifs[[k]] = matrix(0, nrow=shapeFeatureCount, ncol=inputMotifLength)
    for(shift in (-inputMotifLength+1):(inputMotifLength-1)){
      shifttedLocations = finalLocationsInitial[k,] + shift

      newTimeSeriesMotif = matrix(0, nrow=shapeFeatureCount, ncol=inputMotifLength)
      actualCount = 0
      for(i in 1:peakCount){
        if(shifttedLocations[i]>=1 & (shifttedLocations[i]+inputMotifLength-1)<=ncol(inputData[[i]])){
          newTimeSeriesMotif = newTimeSeriesMotif + prior[k,i]*inputData[[i]][,shifttedLocations[i]:(shifttedLocations[i]+inputMotifLength-1)]
          actualCount = actualCount + prior[k,i]
        }
      }
      if(actualCount==0){
        next
      }
      newTimeSeriesMotif = newTimeSeriesMotif / actualCount

      l_totalProb = 0
      for(i in 1:peakCount){
        allBG = sum(log10(dnorm(inputData[[i]],BGmeans,BGsd)))
        if(shifttedLocations[i]>=1 & (shifttedLocations[i]+inputMotifLength-1)<=ncol(inputData[[i]])){
          l_totalProb = l_totalProb +  prior[k,i]*10^( sum(log10(dnorm(inputData[[i]][,shifttedLocations[i]:(shifttedLocations[i]+inputMotifLength-1)], mean=newTimeSeriesMotif, sd=BGsd))) - sum(log10(dnorm(inputData[[i]][,shifttedLocations[i]:(shifttedLocations[i]+inputMotifLength-1)], mean=BGmeans, sd=BGsd))) )
        }else{
          l_totalProb = l_totalProb + prior[k,i]
        }
      }
      if(is.na(l_totalProb)){
        next
      }
      if(l_totalProb>maxLProb){
        finalLocations[k,] = shifttedLocations
        finalTimeSeriesMotifs[[k]] = newTimeSeriesMotif
        maxLProb = l_totalProb
      }
    }
  }

  #print(finalTimeSeriesMotifs)
  #print(timeSeriesMotifs)
  return(list("finalTimeSeriesMotifs"=finalTimeSeriesMotifs,"finalLocations"=finalLocations))
}
