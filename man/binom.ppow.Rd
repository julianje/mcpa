% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomial.R
\name{binom.ppow}
\alias{binom.ppow}
\title{Compute the power of a binomial experiment using pilot data.}
\usage{
binom.ppow(n, pilotdata, r = 10000, alternative = c("two.sided", "less",
  "greater"), alpha = 0.05, nullp = 0.5, conf.level = 0.95)
}
\arguments{
\item{n}{sample size.}

\item{pilotdata}{a list with numerical or categorical pilot data. Largest number of second category is treated as success.}

\item{r}{number of simulations to compute power.}

\item{alternative}{type of alternative hypothesis in binomial test. Must be "\code{two.sided}" (default), "\code{greater}", or "\code{less}".}

\item{alpha}{significance threshhold.}

\item{nullp}{probability of success in null hypothesis.}

\item{conf.level}{size of confidence intervals.}
}
\value{
The probability of finding \eqn{p < \alpha} with the experiment description and a 95% Clopper-Pearson confidence interval.
}
\description{
\code{binom.ppow} computes (through simulation) the power of a binomial experiment with a specified sample size, using pilot data.
}
\examples{
binom.ppow(n=16, pilotdata = (1, 1, 0, 0, 1, 1), p=0.8)
binom.ppow(n=16, pilotdata = ("a", "b", "b", "b"), p=0.8, alternative="greater")
binom.ppow(n=16, pilotdata = ("a", "b", "b", "b"), r=5000, nullp=0.25, alternative="greater")
}
\seealso{
\code{\link{binom.power}}, \code{\link{binom.ppow}}, \code{\link{binom.explore}}, and \code{\link{binom.pexplore}}.
}
