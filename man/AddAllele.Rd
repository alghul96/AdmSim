\name{AddAllele}
\alias{AddAllele}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AddAllele
}
\description{
Add an allele the members of a population. It is possible to add dominant and recessive alleles to each member of a population based on Hardy-Weinberg proportions.
Moreover, it add that specific allele was added to the population in the form a \code{geneSite} object.
}
\usage{
AddAllele(population, selectionAllele, alleleFreq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{population}{
  An object of type \code{population}.
}
  \item{selectionAllele}{
  An object of type \code{geneSite}. Notice, here a selection coefficient should be specified.
}
  \item{alleleFreq}{
  The frequency of the dominant allele in the population.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (population, selectionAllele, alleleFreq)
{
    alleleList = c(population@alleles, selectionAllele)
    indList = population@members
    nIndivid = length(indList)
    freq = GenotypeProportions(alleleFreq)
    alleleA = Allele(location = selectionAllele@location, type = toupper(selectionAllele@type))
    allelea = Allele(location = selectionAllele@location, type = tolower(selectionAllele@type))
    indPerPop = round(c(nIndivid * freq$fAA, nIndivid * freq$fAa,
        nIndivid * freq$faa))
    if (sum(indPerPop) != nIndivid) {
        warning("Not possible to plug an allele with these frequencies.\n      Rounding up the aa genotype individuals.")
        indPerPop[3] = nIndivid - (indPerPop[1] + indPerPop[2])
    }
    genot = rep(c("AA", "Aa", "aa"), times = indPerPop)
    genot = sample(genot)
    for (i in 1:nIndivid) {
        if (genot[i] == "AA") {
            indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles,
                alleleA)
            indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles,
                alleleA)
        }
        else if (genot[i] == "Aa") {
            indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles,
                alleleA)
            indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles,
                allelea)
        }
        else if (genot[i] == "aa") {
            indList[[i]]@chromosome1@alleles = c(indList[[i]]@chromosome1@alleles,
                allelea)
            indList[[i]]@chromosome2@alleles = c(indList[[i]]@chromosome2@alleles,
                allelea)
        }
        else {
            return("There is an unrecognized genotype")
        }
    }
    output = Population(members = indList, alleles = unlist(alleleList))
    validObject(output)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
