\name{Chromosome-class}
\Rdversion{1.1}
\docType{class}
\alias{Chromosome-class}
\alias{Chromosome}
\alias{as.char.chromosome,Chromosome-method}
\alias{as.vector.Chromosome,Chromosome-method}
\alias{checkGenotype,Chromosome,Allele-method}
\alias{compactChromosome,Chromosome-method}
\alias{splitChromosome,Chromosome-method}
\alias{splitChromosomeAlleles,Chromosome-method}

\title{Class \code{"Chromosome"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Chromosome(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{components}:}{Object of class \code{"list"} ~~ }
    \item{\code{alleles}:}{Object of class \code{"listOrNULL"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.char.chromosome}{\code{signature(chrom = "Chromosome")}: ... }
    \item{as.vector.Chromosome}{\code{signature(chrom = "Chromosome")}: ... }
    \item{checkGenotype}{\code{signature(chromosome = "Chromosome", site = "Allele")}: ... }
    \item{compactChromosome}{\code{signature(chrom = "Chromosome")}: ... }
    \item{splitChromosome}{\code{signature(chrom = "Chromosome")}: ... }
    \item{splitChromosomeAlleles}{\code{signature(chrom = "Chromosome")}: ... }
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
showClass("Chromosome")
}
\keyword{classes}
\usage{
Chromosome(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
