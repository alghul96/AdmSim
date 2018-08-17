\name{MateIndividuals}
\alias{MateIndividuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for mating two individuals in a population.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MateIndividuals(parent1, parent2, lambda = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parent1}{
%%     ~~Describe \code{parent1} here~~
}
  \item{parent2}{
%%     ~~Describe \code{parent2} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
structure(function (parent1, parent2, lambda = 2)
{
    standardGeneric("MateIndividuals")
  }, generic = structure("MateIndividuals", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = c("parent1",
"parent2", "lambda"), default = `\001NULL\001`, skeleton = (function (parent1,
    parent2, lambda = 2)
stop("invalid call in method dispatch to 'MateIndividuals' (no default method)",
    domain = NA))(parent1, parent2, lambda), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
