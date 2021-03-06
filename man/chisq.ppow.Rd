% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chisq.R
\name{chisq.ppow}
\alias{chisq.ppow}
\title{Compute the power using a chi-squared test by bootstrapping pilot data.}
\usage{
chisq.ppow(x, y = NULL, n, r = 10000, alpha = 0.05, conf.level = 0.95)
}
\arguments{
\item{x}{a list containing count data of pilot observations.}

\item{y}{an optional list with the count observations of a second group.}

\item{n}{sample size.}

\item{r}{number of simulations to compute power.}

\item{alpha}{significance threshhold.}

\item{conf.level}{size of confidence intervals.}
}
\value{
The probability of finding \eqn{p < \alpha} with the experiment description.
}
\description{
\code{chisq.ppow} computes (via simulation) the power of an experiment that will be analyzed using a chi-square test by bootstrapping pilot data.
This function directly calls \code{chisq.pow} but is included to maintaing consistency with other power functions in the package.
}
\examples{
chisq.ppow(x=c(3, 5, 2), n=32) # Chi-squared distribution against chance with 32 observations.
chisq.ppow(x=c(2, 2, 2), y=c(8, 14, 1)) # Compare both group with 16 observations in each group.
}
\seealso{
\code{\link{chisq.pow}}, \code{\link{chisq.ppow}}, \code{\link{chisq.explore}}, and \code{\link{chisq.pexplore}}.
}
