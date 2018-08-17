#### ALLELE ####
Allele <- setClass(
  "Allele",
  slots = c(location = "numeric",
            type = "character"),
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
geneSite <- setClass(
  "geneSite",
  slots = c(s = "numeric",
            selectionCase = "character"),
  contains = "Allele",
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


# compute viabilities according to directional selection
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
