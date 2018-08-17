\name{GetAlleleFreq}
\alias{GetAlleleFreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the frequency of an allele in a population.
}
\description{
Takes a population as input an the allele of interest, and return its frequency.
}
\usage{
GetAlleleFreq(population, allele)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{population}{
%%     ~~Describe \code{population} here~~
}
  \item{allele}{
%%     ~~Describe \code{allele} here~~
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
structure(function (population, allele)
{
    standardGeneric("GetAlleleFreq")
  }, generic = structure("GetAlleleFreq", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = c("population",
"allele"), default = `\001NULL\001`, skeleton = (function (population,
    allele)
stop("invalid call in method dispatch to 'GetAlleleFreq' (no default method)",
    domain = NA))(population, allele), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
