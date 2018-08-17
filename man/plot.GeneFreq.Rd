\name{plot.GeneFreq}
\alias{plot.GeneFreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the frequency of population genes inside a population.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.GeneFreq(popFrequency, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popFrequency}{
%%     ~~Describe \code{popFrequency} here~~
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
structure(function (popFrequency, ...)
{
    standardGeneric("plot.GeneFreq")
  }, generic = structure("plot.GeneFreq", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = "popFrequency", default = `\001NULL\001`, skeleton = (function (popFrequency,
    ...)
stop("invalid call in method dispatch to 'plot.GeneFreq' (no default method)",
    domain = NA))(popFrequency, ...), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
