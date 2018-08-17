# find the proportion of individuals with genotype AA, Aa, or aa based on
# Hardy-Weinberg proportions
GenotypeProportions = function(freqA) {

  if (freqA > 1 || freqA < 0)
    return("The allele frequency should be between 0 and 1.")
  freqa = 1 - freqA

  freqAA = freqA^2
  freqAa = 2 * freqA * freqa
  freqaa = freqa^2
  return(list(fAA = freqAA, fAa = freqAa, faa = freqaa))
}

# generate a base population with individuals of a certain genotype in selection
# a list containing the alleles that are favourable for selection, while the
# alleleFreq the frequency on which the dominant allele is present
CreatePopulation = function(popName = "Red", nIndividuals = 100, generation = 0) {

  chrom = Chromosome(components = c(geneSeq(pop = popName)), alleles = NULL)

  indList = lapply(1:nIndividuals, function(i) Individual(generation = generation,
    chromosome1 = chrom, chromosome2 = chrom))

  return(Population(members = indList, alleles = NULL))
}


AddAllele = function(population, selectionAllele, alleleFreq) {

  # add the new allele to the already existing ones
  alleleList = c(population@alleles, selectionAllele)

  # unpack the individuals
  indList = population@members
  nIndivid = length(indList)

  # finding the frequencies of the alleles
  freq = GenotypeProportions(alleleFreq)

  alleleA = Allele(location = selectionAllele@location, type = toupper(selectionAllele@type))
  allelea = Allele(location = selectionAllele@location, type = tolower(selectionAllele@type))

  # getting which genotype each individual will be

  indPerPop = round(c(nIndivid * freq$fAA, nIndivid * freq$fAa, nIndivid * freq$faa))

  if (sum(indPerPop) != nIndivid) {
    warning("Not possible to plug an allele with these frequencies.
      Rounding up the aa genotype individuals.")
    indPerPop[3] = nIndivid - (indPerPop[1] + indPerPop[2])
  }

  genot = rep(c("AA", "Aa", "aa"), times = indPerPop)
  genot = sample(genot)

  for (i in 1:nIndivid) {
    if (genot[i] == "AA") {
      indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles, alleleA)
      indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles, alleleA)
    } else if (genot[i] == "Aa") {
      indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles, alleleA)
      indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles, allelea)
    } else if (genot[i] == "aa") {
      indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles, allelea)
      indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles, allelea)
    } else {
      return("There is an unrecognized genotype")
    }
  }
  output = Population(members = indList, alleles = unlist(alleleList))
  validObject(output)
  return(output)
}


#### PROGRESS BAR #### just a progress bar to keep track of the simulation
ProgressBar = function(current, total, msg = NULL) {
  prop = current/total
  percentage = paste0(round(prop, 3) * 100, "% ")
  text = paste0("(", current, "/", total, ") ")
  bar = c(rep("=", floor(40 * prop)), ">", rep(" ", ceiling(40 * (1 - prop))))
  cat("|", bar, "| ", percentage, text, msg, sep = "", "\r")
  flush.console()
}


#### MATE POPULATION #### combine two or more population toghether to generate a new
#### one popList: 1 or more populations that are expected to reproduce nChildren:
#### expected number of children obtained from this generation cpc: expected number
#### of children per couple

library(parallel)
MatePopulations.par = function(popList, nChildren = 100, cpc = 2, lambda = 2, mc.cores = 1) {
  # merging all the populations individuals at initial generation
  individuals = sapply(popList, function(x) x@members)
  individuals = unlist(individuals)

  # obtaining all the populations allele loci
  alleles = sapply(popList, function(x) x@alleles)
  alleles = unlist(alleles)
  alleles = unique(alleles)

  # creating a vector of expected children for each couple
  nCouples = round(nChildren/cpc)  # see how many couples should reproduce
  childrenPerCouple = rpois(nCouples, cpc)
  childrenPerCouple = childrenPerCouple[childrenPerCouple != 0]

  # function for generating children in a couple
  generateChildren = function(i) {
    parents = sample(individuals, 2)
    children = sapply(1:i, function(x) MateIndividuals(parent1 = parents[[1]],
      parent2 = parents[[2]],
      lambda = lambda))

    return(children)
  }

  # merging all the children
  childrenList = mclapply(childrenPerCouple, generateChildren, mc.cores = mc.cores)

  output = Population(members = unlist(childrenList), alleles = alleles)
  validObject(output)
  return(SelectIndividuals(output))
}
