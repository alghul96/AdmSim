\name{DoCrossover}
\alias{DoCrossover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs crossing over for an Individual and return a chromosome.
}
\description{
This function splits and combines two chromosome of each individual, returning a resulting
chromosome that is then eventually passed to its child.
}
\usage{
DoCrossover(individual, lambda = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individual}{
The individual on which crossing over should be performed.
}
  \item{lambda}{
The expected number of pieces in which a chromosome will split.
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
structure(function (individual, lambda = 2)
{
    standardGeneric("DoCrossover")
  }, generic = structure("DoCrossover", package = "AdmSim"), package = "AdmSim", group = list(), valueClass = character(0), signature = c("individual",
"lambda"), default = `\001NULL\001`, skeleton = (function (individual,
    lambda = 2)
stop("invalid call in method dispatch to 'DoCrossover' (no default method)",
    domain = NA))(individual, lambda), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
