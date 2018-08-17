#### ALLELE ####
# An allele can be a piece of information sitting on a chromosome of either type A or a, or
# any other type.

Allele <- setClass(
  "Allele",
  
  slots = c(location = "numeric",
            type = "character"),
  
  # the prototype is "Allele of type A located at 0"
  prototype = list(location = 0,
                   type = "A"),
  
  validity = function(object) {
    if ((object@location > 1) || (object@location < 0)) {
      return("The location should be between 0 and 1")
    }
    return(TRUE)
  }
)


#### GENE SITE ####
# A gene site is an allele on a chromosome that can have an advantageous (or disadvantageous) allele.
# Hence, the type field refears to a type of an allele (usually A dominant, a recessive).
# The s is the selection coefficient, while l is the locus of the allele on a chromosome.
# A selection coefficient of 0 will mean that no advantage
# is carried by the presence/absence of any population genome on that locus.

geneSite <- setClass(
  "geneSite",
  
  slots = c(s = "numeric",
            selectionCase = "character"),
  contains = "Allele",
  
  # the prototype is "No advantage, it is located at 0"
  prototype = list(s = 0,
                   selectionCase = "additive"),
  
  validity = function(object) {
    if ((object@location > 1) || (object@location < 0)) {
      return("The location should be between 0 and 1")
    }
    if ((object@s > 1) || (object@s < 0)) {
      return("The selection coefficient should be between 0 and 1")
    }
    return(TRUE)
  }
)


# This function returns the viability corresponding to the three different genotypes
# that can occur from the mating of two individuals
setGeneric(
  name = "ComputeViability",
  def = function(genesite) {
    standardGeneric("ComputeViability")
  }
)

setMethod(
  f = "ComputeViability",
  signature = "geneSite",
  definition = function(genesite) {
    
    switch (genesite@selectionCase,
      additive = return(list(
        vAA = 1,
        vAa = 1 - genesite@s,
        vaa = 1 - 2 * genesite@s
      )),
      dominant = return(list(
        vAA = 1,
        vAa = 1,
        vaa = 1 - genesite@s
      )),
      recessive = return(list(
        vAA = 1,
        vAa = 1 - genesite@s,
        vaa = 1 - genesite@s
      ))
    )

  }
)
