% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomial.R
\name{binom.pow}
\alias{binom.pow}
\title{Compute the power of a binomial experiment.}
\usage{
binom.pow(n, p, r = 10000, alternative = c("two.sided", "less", "greater"),
  alpha = 0.05, nullp = 0.5, conf.level = 0.95)
}
\arguments{
\item{n}{sample size.}

\item{p}{predicted probability of success.}

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
\code{binom.pow} computes (via simulation) the power of a binomial experiment with a specified sample size and probability of success.
}
\examples{
binom.pow(n=16, p=0.8)
binom.pow(n=16, p=0.8, alternative="greater")
binom.pow(n=32, p=0.6, r=5000, nullp=0.25)
}
\seealso{
\code{\link{binom.power}}, \code{\link{binom.ppow}}, \code{\link{binom.explore}}, and \code{\link{binom.pexplore}}.
}
