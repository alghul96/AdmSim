### MCPOPSIM ###
MCPopSim = setClass(
  "MCPopSim",
  slots = c(posteriorDosages = "list",
            nSims = "numeric")
)

setGeneric(
  name = "MCSim",
  def = function(popList, nGens, nSims = 99,
                 resolution = 100, mc.cores = 4) {
    standardGeneric("MCSim")
  }
)

setMethod(
  f = "MCSim",
  signature = c("list", "numeric"),
  definition = function(popList, nGens, nSims = 99,
                        resolution = 100, mc.cores = 4){

    # Obtaining the populations genomes
    for (pop in popList) {
      genefreq = GetGeneFreq(pop)
      if (nrow(genefreq@chromosome1) != 1 || nrow(genefreq@chromosome2) != 1){
        stop('The list did not contain isolated populations at generation 0 only.')
      }
    }

    # Obtain the populations frequencies at time 0 and the number of children
    nChild = 0
    popSizes = NULL
    for (pop in popList) {
      popSize = length(pop@members)
      popSizes = c(popSizes, popSize)
      nChild = nChild + popSize
    }
    popSizes = t(as.matrix(popSizes))

    cat(nSims, "simulations for", nGens, "generation of", nChild, "individuals \n")

    # first post admixture generation
    ProgressBar(1, nGens)
    simuList = rep(0, nSims)
    simuList = lapply(simuList, MatePopulations.par,
                      popList = popList,
                      nChildren = nChild,
                      cpc = 1,
                      mc.cores = mc.cores)

    # other post admixture generations
    for(j in 2:nGens) {
      ProgressBar(j, nGens)
      simuList = lapply(simuList, c)
      simuList = mclapply(simuList, function(x) {
        MatePopulations.par(
          popList = x,
          nChildren = nChild,
          cpc = 1,
          mc.cores = 1)},
        mc.cores = mc.cores)
    }

    # obtaining the number of individuals for pop. at every locus
    freqs = lapply(simuList, GetGeneFreq, resolution = resolution)
    chrom1Nindivid = lapply(freqs, function(f)
      f@chromosome1 * f@nIndividuals)
    chrom2Nindivid = lapply(freqs, function(f)
      f@chromosome2 * f@nIndividuals)

    popNames = unique(unlist(lapply(chrom1Nindivid, rownames)))

    colnames(popSizes) = popNames

    output = list()
    for (popName in popNames) {

      popZChrom1 = rowMeans(sapply(1:nSims, function(i)
        chrom1Nindivid[[i]][popName,]))
      popZChrom2 = rowMeans(sapply(1:nSims, function(i)
        chrom1Nindivid[[i]][popName,]))

      # getting the beta and alpha parameters of the prior
      mu = popSizes[, popNames] / nChild
      var = nGens * (1 / (popSizes[, popNames]))
      alphaPrior <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      betaPrior <- alphaPrior * (1 / mu - 1)

      # getting the posterior Beta parameters
      alpha = (popZChrom1 + popZChrom2) / 2 + alphaPrior
      beta = (nChild) - (popZChrom1 + popZChrom2) / 2 + betaPrior

      output[[popName]] = data.frame(alpha = alpha, beta = beta)
    }
    return(MCPopSim(posteriorDosages = output, nSims = nSims))
  }
)


setGeneric(
  name = "plot.MCPopSim",
  def = function(popSim, testPop = NULL, alpha = 0.05){
    standardGeneric("plot.MCPopSim")
  }
)

setMethod(
  f = "plot.MCPopSim",
  signature = "MCPopSim",
  definition = function(popSim, testPop, alpha = 0.05) {
    if (is.null(testPop)) {
      stop("Please provvide a population to test!")
    }

    resolution = nrow(popSim@posteriorDosages[[1]])
    popFrequency = GetGeneFreq(testPop, resolution = resolution)
    # unpack the elements in GeneFreq
    chrom1Freq = popFrequency@chromosome1
    chrom2Freq = popFrequency@chromosome2
    resolution = popFrequency@resolution

    # in case of pure individuals
    if (is.null(rownames(chrom1Freq))) {
      nPop <- 1
      rownames(chrom1Freq) = "population"
    } else nPop = nrow(chrom1Freq)

    # Plotting of the frequency between the populations
    totalFreq = chrom1Freq + chrom2Freq

    plot(x = 1:resolution, totalFreq[1, ],
         ylab = "Ancestry Dosages",
         xlab = "Chromosome locus",
         ylim = c(0, 2),
         type = "l", col = 1,
         panel.first = grid(),
         lwd = 2)

    p = 1
    for (mcpopSim in popSim@posteriorDosages) {
      power = 1 - (alpha / 2)
      maxFreq = apply(mcpopSim, 1, function(x)
        qbeta(power, x["alpha"], x["beta"])) * 2
      minFreq = apply(mcpopSim, 1, function(x)
        qbeta(1 - power, x["alpha"], x["beta"])) * 2
      polygon(c(rev(1:resolution), 1:resolution), c(rev(minFreq), maxFreq),
              col = adjustcolor(p, alpha.f = .1), border = NA)
      p = p + 1
    }

    for(p in 1:nPop)
      lines(1:resolution, totalFreq[p, ], col = 1 + p, lwd = 2)
    legend("topright", rownames(chrom1Freq), lwd = 2, lty = 1, col = 2:(nPop + 2))

  }
)
