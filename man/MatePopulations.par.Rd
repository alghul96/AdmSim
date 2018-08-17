\name{MatePopulations.par}
\alias{MatePopulations.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a new generation out of 1 or more populations. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MatePopulations.par(popList, nChildren = 100, cpc = 2, lambda = 2, mc.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popList}{
%%     ~~Describe \code{popList} here~~
}
  \item{nChildren}{
%%     ~~Describe \code{nChildren} here~~
}
  \item{cpc}{
%%     ~~Describe \code{cpc} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
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
function (popList, nChildren = 100, cpc = 2, lambda = 2, mc.cores = 1)
{
    individuals = sapply(popList, function(x) x@members)
    individuals = unlist(individuals)
    alleles = sapply(popList, function(x) x@alleles)
    alleles = unlist(alleles)
    alleles = unique(alleles)
    nCouples = round(nChildren/cpc)
    childrenPerCouple = rpois(nCouples, cpc)
    childrenPerCouple = childrenPerCouple[childrenPerCouple !=
        0]
    generateChildren = function(i) {
        parents = sample(individuals, 2)
        children = sapply(1:i, function(x) MateIndividuals(parent1 = parents[[1]],
            parent2 = parents[[2]], lambda = lambda))
        return(children)
    }
    childrenList = mclapply(childrenPerCouple, generateChildren,
        mc.cores = mc.cores)
    output = Population(members = unlist(childrenList), alleles = alleles)
    return(SelectIndividuals(output))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
