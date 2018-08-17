\name{Individual-class}
\Rdversion{1.1}
\docType{class}
\alias{Individual-class}
\alias{Individual}
\alias{DoCrossover,Individual-method}
\alias{MateIndividuals,Individual-method}
\alias{plot.Individual,Individual-method}
\alias{show,Individual-method}

\title{Class \code{"Individual"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Individual(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{generation}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{chromosome1}:}{Object of class \code{"Chromosome"} ~~ }
    \item{\code{chromosome2}:}{Object of class \code{"Chromosome"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{DoCrossover}{\code{signature(individual = "Individual")}: ... }
    \item{MateIndividuals}{\code{signature(parent1 = "Individual")}: ... }
    \item{plot.Individual}{\code{signature(indiv = "Individual")}: ... }
    \item{show}{\code{signature(object = "Individual")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Individual")
}
\keyword{classes}
\usage{
Individual(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
