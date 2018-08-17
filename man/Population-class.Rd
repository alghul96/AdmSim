\name{Population-class}
\Rdversion{1.1}
\docType{class}
\alias{Population-class}
\alias{Population}
\alias{AveragePairwiseDistances,Population-method}
\alias{GetAlleleFreq,Population,Allele-method}
\alias{GetGeneFreq,Population-method}
\alias{plot.Population,Population-method}
\alias{SelectIndividuals,Population-method}
\alias{show,Population-method}

\title{Class \code{"Population"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Population(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{members}:}{Object of class \code{"list"} ~~ }
    \item{\code{alleles}:}{Object of class \code{"listOrNULL"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{AveragePairwiseDistances}{\code{signature(population = "Population")}: ... }
    \item{GetAlleleFreq}{\code{signature(population = "Population", allele = "Allele")}: ... }
    \item{GetGeneFreq}{\code{signature(population = "Population")}: ... }
    \item{plot.Population}{\code{signature(population = "Population")}: ... }
    \item{SelectIndividuals}{\code{signature(population = "Population")}: ... }
    \item{show}{\code{signature(object = "Population")}: ... }
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
showClass("Population")
}
\keyword{classes}
\usage{
Population(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
