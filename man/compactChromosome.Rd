\name{compactChromosome}
\alias{compactChromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduces the dimention of a chromosome.
}
\description{
The following function is needed to store compactly informations since it
merge toghether conseguent gene sequences of the same popoulation.
For instance: "red to 0.1 | red to 0.2 | blue to 1" becomes "red to 0.2 | blue to 1"
}
\usage{
compactChromosome(chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{
%%     ~~Describe \code{chrom} here~~
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
structure(function (chrom)
{
    standardGeneric("compactChromosome")
  }, generic = structure("compactChromosome", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = "chrom", default = `\001NULL\001`, skeleton = (function (chrom)
stop("invalid call in method dispatch to 'compactChromosome' (no default method)",
    domain = NA))(chrom), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
