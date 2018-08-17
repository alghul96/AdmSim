\name{GenotypeProportions}
\alias{GenotypeProportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the proportion of individuals with genotype AA, Aa, or aa based on
Hardy-Weinberg proportions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenotypeProportions(freqA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freqA}{
The frequency of a dominant allele.
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
function (freqA)
{
    if (freqA > 1 || freqA < 0)
        return("The allele frequency should be between 0 and 1.")
    freqa = 1 - freqA
    freqAA = freqA^2
    freqAa = 2 * freqA * freqa
    freqaa = freqa^2
    return(list(fAA = freqAA, fAa = freqAa, faa = freqaa))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
