\name{plot.Population}
\alias{plot.Population}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the populations frequencies inside a population. A wrapper around plot.GeneFreq.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.Population(population, resolution = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{population}{
%%     ~~Describe \code{population} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
structure(function (population, resolution = 100, ...)
{
    standardGeneric("plot.Population")
  }, generic = structure("plot.Population", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = c("population",
"resolution"), default = `\001NULL\001`, skeleton = (function (population,
    resolution = 100, ...)
stop("invalid call in method dispatch to 'plot.Population' (no default method)",
    domain = NA))(population, resolution, ...), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
