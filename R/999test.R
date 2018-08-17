setGeneric(
  name = "selTest",
  def = function(pop, window, h0 = 1/2, resolution = 100) {
    standardGeneric("selTest")
  }
)


setMethod(
  f = "selTest",
  signature = "Population",
  definition = function(pop, window, h0 = 1/2, resolution = 100){
    
    freqs = GetGeneFreq(pop, resolution = 100)
    chrom1Freq = freqs@chromosome1
    chrom2Freq = freqs@chromosome2
    
    avgFreq = (chrom1Freq + chrom2Freq) / 2
    
    if (length(h0) == 1) {
      p = h0
    } else if (is.vector(h0)) {
      p = apply(avgFreq[ , h0], 1, mean)
    } else {
      stop("Please provvide a probability or an interval for the null hypthesis...")
    }
    
    genFreq = table(sapply(pop@members, function(x) x@generation))
    maxGen = names(genFreq)[which.max(genFreq)]

    mu = p
    var = as.numeric(maxGen) * p
    alphaPrior <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    betaPrior <- alphaPrior * (1 / mu - 1)
    
    z = round(apply(avgFreq[, window], 1, max) * freqs@nIndividuals)
    
    return(alphaPrior)
    for (anc in names(z)){
      zM0 = dbinom(z[anc], freqs@nIndividuals, p)
      posterior(q) = function(q) {
        dbeta(q, shape1 = (alphaPrior[anc] + z), shape2 = (betaPrior[anc] + freqs@nIndividuals - z))
      }
      zM1 = integrate(posterior, lower = 0, upper = 1)
      cat("For ", anc, "the BF is of", zM0/zM1, "in favour of the H0", "\n")
    }
  }
)

selTest(pop, window = 15:22, h0 = 40:50)
