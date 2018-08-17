#### POPULATION #####
Population <- setClass(
  "Population",
  slots = c(members = "list",
            alleles = "listOrNULL"),
  prototype = list(members = lapply(1:100, function(i)
    Individual()),
    alleles = NULL),
  validity = function(object) {
    if (length(object@members) == 0L)
      return("No members in this population!")
  }
)

# GENE FREQUENCY
GeneFreq <- setClass(
  "GeneFreq",
  #  slots regarding a summary of the two chromosomes in the popoulation
  slots = c(chromosome1 = "matrix",
            chromosome2 = "matrix",
            resolution = "numeric",
            nIndividuals = "numeric")
)


# return the population in a readable format
setMethod(
  f = "show",
  signature = "Population",
  definition = function(object) {
    individuals = object@members
    if (length(individuals) == 0){
      cat("Population has no individuals!")
      return(NULL)
    }
    popLen = paste("Population of", length(individuals), "individuals")
    genFreq = table(sapply(individuals, function(x) x@generation))
    maxGen = names(genFreq)[which.max(genFreq)]
    cat(popLen, "\nMost frequent generation number: ", maxGen)
    for (allele in object@alleles) {
      if (allele@s != 0)
        cat("\nSelection of ", allele@s,
          " for allele type ", allele@type,
          " at distance ", allele@location,
          " with freq. ", round(GetAlleleFreq(object, allele), 3),
          sep = "")
    }
    cat("\n")
  }
)


# Perform selection of the individuals that will reach reproduction
setGeneric(
  name = "SelectIndividuals",
  def = function(population){
    standardGeneric("SelectIndividuals")
  }
)


setMethod(
  f = "SelectIndividuals",
  signature = "Population",
  definition = function(population) {
    individuals = population@members
    survProb = rep(1, length(individuals))
    chrom1List = sapply(individuals, function(x) x@chromosome1)
    chrom2List = sapply(individuals, function(x) x@chromosome2)
    for(allele in population@alleles){
      viability = ComputeViability(allele)
      chrom1 = sapply(chrom1List, checkGenotype, site = allele)
      chrom2 = sapply(chrom2List, checkGenotype, site = allele)
      # this may be either 2, 1 or 0
      genotype = chrom1 + chrom2
      prob = ifelse(genotype == 2,
                    viability$vAA,
                    ifelse(genotype == 1, viability$vAa, viability$vaa))
      # assuming independence, update the total probability
      survProb = survProb * prob
    }
    # creating an index to select individual that will survive
    thres = runif(length(individuals))
    survingIndex = survProb > thres
    output = Population(members = individuals[survingIndex],
                        alleles = population@alleles)
    validObject(output)
    return(output)
  }
)


# Obtaining the dominant allele frequency on a given locus
setGeneric(
  name = "GetAlleleFreq",
  def = function(population, allele) {
    standardGeneric("GetAlleleFreq")
  }
)

setMethod(
  f = "GetAlleleFreq",
  signature = c("Population", "Allele"),
  definition = function(population, allele) {
    individuals = population@members
    chrom1List = sapply(individuals, function(x) x@chromosome1)
    chrom2List = sapply(individuals, function(x) x@chromosome2)
    chromList = unlist(c(chrom1List, chrom2List))
    chrom = sapply(chromList, checkGenotype,
                    site = allele)
    output = sum(chrom)/length(chrom)
    return(output)
  }
)


# Obtaining the frequency of populations per gene
setGeneric(
  name = "GetGeneFreq",
  def = function(population, resolution = 100) {
    standardGeneric("GetGeneFreq")
  }
)


setMethod(
  f = "GetGeneFreq",
  signature = "Population",
  definition = function(population, resolution = 100) {
    individuals = population@members
    popLen = length(individuals)
    chromosomes1 = sapply(individuals, function(x) x@chromosome1)
    chromosomes2 = sapply(individuals, function(x) x@chromosome2)
    # make the chromosomes in vector form
    chromosomes1 = mclapply(chromosomes1,
                            as.vector.Chromosome, resolution = resolution)
    chromosomes2 = mclapply(chromosomes2,
                            as.vector.Chromosome, resolution = resolution)
    popNames = levels(as.factor(unlist(c(chromosomes1, chromosomes2))))
    # make a dataframe with the vectorized chromosomes
    chromosomes1 = do.call(rbind.data.frame, chromosomes1)
    chromosomes2 = do.call(rbind.data.frame, chromosomes2)
    colnames(chromosomes1) = 1:(resolution)
    colnames(chromosomes2) = 1:(resolution)
    # ensure that you at least have a chromosome per haplotype
    toAppend = as.data.frame(t(sapply(popNames, rep, resolution)))
    names(toAppend) = names(chromosomes1)
    chromosomes1 = rbind(toAppend, chromosomes1)
    chromosomes2 = rbind(toAppend, chromosomes2)
    # obtain the frequency of each sequence
    chrom1Freq = (apply(chromosomes1, 2, table) - 1) / popLen
    chrom2Freq = (apply(chromosomes2, 2, table) - 1) / popLen
    # in case we have 1 population on a gene, treat it as a matrix
    if (is.null(nrow(chrom1Freq))) chrom1Freq = t(matrix(chrom1Freq))
    if (is.null(nrow(chrom2Freq))) chrom2Freq = t(matrix(chrom2Freq))
    output = GeneFreq(chromosome1 = chrom1Freq,
      chromosome2 = chrom2Freq,
      resolution = resolution,
      nIndividuals = popLen)
    return(output)
  }
)


# Obtaining the frequency of populations per gene
setGeneric(
  name = "plot.GeneFreq",
  def = function(popFrequency, ...) {
    standardGeneric("plot.GeneFreq")
  }
)

setMethod(
  f = "plot.GeneFreq",
  signature = "GeneFreq",
  definition = function(popFrequency, ...) {
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
      xlab = "Chromosome Locus",
      ylim = c(0, 2),
      type = "l", col = 1,
      panel.first = grid(),
      lwd = 2,
      ...)
    for(p in 1:nPop)
      lines(1:resolution, totalFreq[p, ], col = 1 + p, lwd = 2)
    legend("topright", rownames(chrom1Freq), lty = 1, col = 2:(nPop + 2))
  }
)


# A wrapper around the previous two functions to make it nicer to call
setGeneric(
  name = "plot.Population",
  def = function(population, resolution = 100, ...) {
    standardGeneric("plot.Population")
  }
)

setMethod(
  f = "plot.Population",
  signature = "Population",
  definition = function(population, resolution = 100, ...) {
    popFreq = GetGeneFreq(population, resolution)
    plot.GeneFreq(popFreq, ...)
  }
)
