####    INDIVIDUAL    ####
Individual <- setClass(
  "Individual",
  slots = c(
    generation = "numeric",
    chromosome1 = "Chromosome",
    chromosome2 = "Chromosome"
  ),
  prototype = list(
    generation = 0,
    chromosome1 = Chromosome(),
    chromosome2 = Chromosome()
  ),
  validity = function(object) {
    if (object@generation < 0)
      return("A generation should always be greater then 0")
    return(TRUE)
  }
)


# return the individual in a readable format
setMethod(
  f = "show",
  signature = "Individual",
  definition = function(object) {
    gen = paste("Individual of generation:", object@generation)
    output = cat(gen, "\nChromosomes:\n",
      as.char.chromosome(object@chromosome1), "\n",
      as.char.chromosome(object@chromosome2), "\n")
  }
)


# Perform the crossover and return a combined chromosome
setGeneric(
  name = "DoCrossover",
  def = function(individual, lambda = 2.0) {
    standardGeneric("DoCrossover")
  }
)

setMethod(
  f = "DoCrossover",
  signature = "Individual",
  definition = function(individual, lambda = 2.0) {
    nSplits = rpois(1, lambda) # poisson for regulating the exhanges
    if (nSplits <= 1) {
      splitPoints = NULL
      nSplits <- 1
    }
    else{
      splitPoints = sort(runif(nSplits - 1)) # where to split
    }
    # chromosomes and alleles splits
    splits1pop = splitChromosome(individual@chromosome1, splitPoints)
    splits2pop = splitChromosome(individual@chromosome2, splitPoints)
    alleles1 = splitChromosomeAlleles(individual@chromosome1, splitPoints)
    alleles2 = splitChromosomeAlleles(individual@chromosome2, splitPoints)
    index = sample(1:2, nSplits, replace = T)
    outputComponents = sapply(1:nSplits,
                              function(x) {
                                if (index[x] == 1) splits1pop[[x]]
                                else splits2pop[[x]]
                              })
    outputAlleles = sapply(1:nSplits,
                           function(x) {
                             if (index[x] == 1) alleles1[[x]]
                             else alleles2[[x]]
                           })
    outputChromosome = Chromosome(components = unlist(outputComponents),
                                  alleles = unlist(outputAlleles))
    return(compactChromosome(outputChromosome))
  }
)


# Combine two individuals to generate a child
setGeneric(
  name = "MateIndividuals",
  def = function(parent1, parent2, lambda = 2.0) {
    standardGeneric("MateIndividuals")
  }
)

setMethod(
  f = "MateIndividuals",
  signature = "Individual",
  definition = function(parent1, parent2, lambda = 2.0) {
    # checking that the parents are of the same generation
    if (parent1@generation != parent2@generation)
      warning("\nParents not of the same generation!")
    # performing the crossovers on two individuals
    chromParent1 = DoCrossover(parent1, lambda)
    chromParent2 = DoCrossover(parent2, lambda)
    # merging the two chromosomes into a child
    child = Individual(
      generation = max(parent1@generation, parent2@generation) + 1,
      chromosome1 = chromParent1,
      chromosome2 = chromParent2
    )
    return(child)
  }
)


# Plotting the individual
setGeneric(
  name = "plot.Individual",
  def = function(indiv, resolution = 100, ...) {
    standardGeneric("plot.Individual")
  }
)

setMethod(
  f = "plot.Individual",
  signature = "Individual",
  definition = function(indiv, resolution = 100, ...) {
    chrom1 = as.vector.Chromosome(indiv@chromosome1, resolution)
    chrom2 = as.vector.Chromosome(indiv@chromosome2, resolution)
    # creating an object to store two chromosomes and a white space
    binded = cbind(chrom1, chrom1,
                   rep(NA, resolution),
                   chrom2, chrom2)
    binded = as.numeric(as.factor(binded))
    # creating the matrix to plot
    toPlot <- matrix(data = binded,
                     nc = 5)
    # plotting as an image
    image(toPlot,
          yaxt = "n",
          frame.plot = F,
          ...)
  }
)
