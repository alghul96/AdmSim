#### CHROMOSOME ####
setClassUnion("listOrNULL",members=c("list", "NULL"))
Chromosome <- setClass(
  "Chromosome",
  #  slots containing gene ranges
  slots = c(components = "list",
            alleles = "listOrNULL"),
  # the prototype is a chromosome "Red to 1" with no alleles
  prototype = list(components = c(geneSeq()),
                   alleles = NULL),
  validity = function(object) {
    sequences <- object@components
    nSeq <- length(sequences)
    lastSequence <- sequences[[nSeq]]
    ranges <- sapply(1:nSeq, function(i)
      sequences[[i]]@Rng)
    if (lastSequence@Rng != 1) {
      # checking is last sequence range is a 1
      return(
        paste(
          "A gene sequence in a chromosome should
          always end with a Range of 1, but it's",
          lastSequence@Rng
        ))
    } else if (is.unsorted(ranges)) {
      # checking if the sequences are sorted
      return("Ranges in gene sequences should be sorted!")
    }
    return(TRUE)
  }
)

# return the chromosome in a readable format
setGeneric(
  name = "as.char.chromosome",
  def = function(chrom) {
    standardGeneric("as.char.chromosome")
  }
)

setMethod(
  f = "as.char.chromosome",
  signature = "Chromosome",
  definition = function(chrom) {
    sequences <- chrom@components
    alleles <- chrom@alleles
    nSeq <- length(sequences)
    nAlleles <- length(alleles)
    # obtaining the genSeq in a readable format
    output = sapply(1:nSeq, function(i)
      as.char.geneSeq(sequences[[i]]))
    output = c(output, "  ")
    # adding the alleles
    for (allele in alleles) {
      output = c(output, paste("Allele", allele@type,
      "at distance", allele@location))
    }
    return(paste(output, collapse = " | "))
  }
)

# split a chromosome in smaller gene sequences
setGeneric(
  name = "splitChromosome",
  def = function(chrom, splitPoints) {
    standardGeneric("splitChromosome")
  }
)

setMethod(
  f = "splitChromosome",
  signature = "Chromosome",
  definition = function(chrom, splitPoints) {
    splitPoints <- c(0, splitPoints, 1)
    sequences <- chrom@components
    nSeq <- length(sequences)
    nSplits <- length(splitPoints) + 1
    # obtaining the populations and ranges
    pops <- sapply(1:nSeq, function(i)
      sequences[[i]]@pop)
    rngs <- sapply(1:nSeq, function(i)
      sequences[[i]]@Rng)
    # this function performs a single split!
    SplitSingle <- function(splitPoint) {
      prevSplit = splitPoints[sum(splitPoint > splitPoints)]
      toIncludeBot = sum(rngs < prevSplit) + 1
      toIncludeTop = sum(rngs < splitPoint) + 1
      splitPop = pops[toIncludeBot:toIncludeTop]
      splitRng = sort(c(rngs, splitPoint))[toIncludeBot:toIncludeTop]
      lapply(1:length(splitPop), function(i)
        geneSeq(pop = splitPop[i], Rng = splitRng[i]))
    }
    return(lapply(splitPoints[2:length(splitPoints)], SplitSingle))
  }
)


# The following function is needed to store compactly informations
# e.g.: "red to 0.1 | red to 0.2 | blue to 1" to "red to 0.2 | blue to 1"
setGeneric(
  name = "compactChromosome",
  def = function(chrom) {
    standardGeneric("compactChromosome")
  }
)

setMethod(
  f = "compactChromosome",
  signature = "Chromosome",
  definition = function(chrom) {
    seq <- chrom@components
    # obtaining the populations and ranges
    pops <- sapply(1:length(seq), function(x) seq[[x]]@pop)
    repetitions <- sum(pops[1:length(pops) - 1] == pops[2:length(pops)])
    while (repetitions > 0) {
      i = 1
      seqNew = list()
      # check all the elements in sequences minus the last one
      while (i < length(seq)) {
        if (seq[[i]]@pop == seq[[i + 1]]@pop) {
          seqNew <- c(seqNew, seq[[i + 1]])
          i = i + 2
        } else {
          seqNew <- c(seqNew, seq[[i]])
          i = i + 1
        }
      }
      # check the last element and ensure that it is 1
      if (seqNew[[length(seqNew)]]@pop != seq[[length(seq)]]@pop) {
        seqNew = c(seqNew, seq[[length(seq)]])
      } else {
        seqNew = c(seqNew[1:length(seqNew) - 1], seq[[length(seq)]])
      }
      # update the sequences
      seq <- seqNew
      pops <- sapply(1:length(seq), function(x) seq[[x]]@pop)
      repetitions <- sum(pops[1:length(pops) - 1] == pops[2:length(pops)])
    }
    output <- Chromosome(components = seq, alleles = chrom@alleles)
    validObject(output)
    return(output)
  }
)


# This function is needed for returning which allele is
# on each split of a chromosome
setGeneric(
  name = "splitChromosomeAlleles",
  def = function(chrom, splitPoints) {
    standardGeneric("splitChromosomeAlleles")
  }
)

setMethod(
  f = "splitChromosomeAlleles",
  signature = "Chromosome",
  definition = function(chrom, splitPoints) {
    alleles = chrom@alleles
    nSplits = length(splitPoints) + 1
    alleleList = sapply(1:nSplits, function(i) list(NULL))
    # putting each allele in the correct spot in the list
    for (allele in alleles) {
      index = sum(allele@location >= splitPoints) + 1
      alleleList[[index]] = c(alleleList[[index]], allele)
    }
    return(alleleList)
  }
)


# return the chromosome in a vectorized format
setGeneric(
  name = "as.vector.Chromosome",
  def = function(chrom, resolution = 100) {
    standardGeneric("as.vector.Chromosome")
  }
)

setMethod(
  f = "as.vector.Chromosome",
  signature = c("Chromosome"),
  definition = function(chrom, resolution = 100) {
    sequences <- chrom@components
    nSeq <- length(sequences)
    pops <- sapply(1:nSeq, function(i) sequences[[i]]@pop)
    rngs <- sapply(1:nSeq, function(i) sequences[[i]]@Rng)
    output <- vector(mode = "character", length = resolution)
    rngsScaled = ceiling(rngs * resolution)
    # creating the new object (filling from top to bottom)
    for (i in nSeq:1) {
      output[1:rngsScaled[i]] <- pops[i]
    }
    return(output)
  }
)


# this function return true whether an allele is of a given type
setGeneric(
  name = "checkGenotype",
  def = function(chromosome, site) {
    standardGeneric("checkGenotype")
  }
)

setMethod(
  f = "checkGenotype",
  signature = c("Chromosome", "Allele"),
  definition = function(chromosome, site) {
    alleles <- chromosome@alleles
    locations = sapply(alleles, function(x) x@location)
    types = sapply(alleles, function(x) x@type)
    selection = locations == site@location
    if (sum(selection) == 0) {
      warning("There is no such allele at this locus on the chromosome!")
      return(FALSE)
    }
    return(types[selection] == site@type)
  }
)
