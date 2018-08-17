\name{ProgressBar}
\alias{ProgressBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A simple percentage bar to keep track of a simulation status.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProgressBar(current, total, msg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{current}{
%%     ~~Describe \code{current} here~~
}
  \item{total}{
%%     ~~Describe \code{total} here~~
}
  \item{msg}{
%%     ~~Describe \code{msg} here~~
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
function (current, total, msg = NULL)
{
    prop = current/total
    percentage = paste0(round(prop, 3) * 100, "\% ")
    text = paste0("(", current, "/", total, ") ")
    bar = c(rep("=", floor(40 * prop)), ">", rep(" ", ceiling(40 *
        (1 - prop))))
    cat("|", bar, "| ", percentage, text, msg, sep = "", "\r")
    flush.console()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
