#### GENE SEQUENCE ####

# a gene sequence is needed for coding conseguent
# populations haplotypes
geneSeq <- setClass(
  "geneSeq",
  
  slots = c(pop = "character",
            Rng = "numeric"),
  
  # the prototype is "Red to 1"
  prototype = list(pop = "red",
                   Rng = c(1)),
  
  validity = function(object) {
    if ((object@Rng > 1) || (object@Rng) < 0) {
      return("The Rng should be between 0 and 1")
    } else if (length(object@pop) > 1) {
      return("A sequence should only refer to one population")
    }
    return(TRUE)
  }
)


# METHODS ASSOCIATED WITH GENESEQ
# return the sequence in a readable format
setGeneric(
  name = "as.char.geneSeq",
  def = function(geneseq) {
    standardGeneric("as.char.geneSeq")
  }
)

setMethod(
  f = "as.char.geneSeq",
  signature = "geneSeq",
  definition = function(geneseq) {
    validObject(geneseq)
    output = paste(geneseq@pop, "to", round(geneseq@Rng, 3))
    return(output)
  }
)
