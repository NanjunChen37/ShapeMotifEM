
evaluateMatrices <- function(timeSeriesMotifs, finalTimeSeriesMotifs) {

    # evaluateMatrices function

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

    return(list(
      "EDists" = EDists,
      "Pearsons" = Pearsons,
      "MAPEs" = MAPEs
    ))
  }
