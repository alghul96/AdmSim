\name{CreatePopulation}
\alias{CreatePopulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function intialize a new population without any advantageous allele.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CreatePopulation(popName = "Red", nIndividuals = 100, generation = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popName}{
The population name. For instance "Red" or "Blue". Needed to keep track of the
source of genome.
}
  \item{nIndividuals}{
The numbers of individuals.
}
  \item{generation}{
The generation of each one of the individuals.
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
function (popName = "Red", nIndividuals = 100, generation = 0)
{
    chrom = Chromosome(components = c(geneSeq(pop = popName)),
        alleles = NULL)
    indList = lapply(1:nIndividuals, function(i) Individual(generation = generation,
        chromosome1 = chrom, chromosome2 = chrom))
    return(Population(members = indList, alleles = NULL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
