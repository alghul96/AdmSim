\name{Allele-class}
\Rdversion{1.1}
\docType{class}
\alias{Allele-class}
\alias{Allele}
\alias{checkGenotype,Chromosome,Allele-method}
\alias{GetAlleleFreq,Population,Allele-method}

\title{Class \code{"Allele"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Allele(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{location}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{checkGenotype}{\code{signature(chromosome = "Chromosome", site = "Allele")}: ... }
    \item{GetAlleleFreq}{\code{signature(population = "Population", allele = "Allele")}: ... }
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
showClass("Allele")
}
\keyword{classes}
\usage{
Allele(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
