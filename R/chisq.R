chisq.pow <- function(x, y = NULL, n, r = 10000, alpha = 0.05, conf.level = 0.95){
  #' Compute the power using a chi-squared test.
  #'
  #' \code{chisq.pow} computes (via simulation) the power of an experiment that will be analyzed using a chi-square test.
  #'
  #' @param x a list containing probabilities of each category.
  #' @param y an optional list with the probabilities of the second group.
  #' @param n sample size.
  #' @param r number of simulations to compute power.
  #' @param alpha significance threshhold.
  #' @param mu mean value according to null hypothesis (default = \code{0}). Only used in one sample t-tests.
  #' @param conf.level size of confidence intervals.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' chisq.pow(x=c(0.2, 0.5, 0.3), n=32) # Chi-squared distribution against chance with 32 observations.
  #' chisq.pow(x=c(0.2, 0.5, 0.3), y=c(0.8, 0.1, 0.1)) # Compare both group with 16 observations in each group.
  #' @seealso \code{\link{chisq.pow}}, \code{\link{chisq.ppow}}, \code{\link{chisq.explore}}, and \code{\link{chisq.pexplore}}.
  if (any(is.na(n) | (n < 0))) stop("'n' must be a nonnegative integer")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  if (any(is.na(conf.level) | (conf.level < 0) | (conf.level > 1))) stop("'alpha' must be between 0 and 1")
  x <- x/sum(x)
  if (!is.null(y)){
    y <- y/sum(y)
    if (length(y) != length(x)) stop("'x' and 'y' must have the same length")
  } else{
    y <- rep(1/length(x),length(x))
  }
  xsamples <- rmultinom(r, n, x)
  ysamples <- rmultinom(r, n, y)
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r, clear = FALSE, width = 60
  )
  power <- sum(sapply(1:r,function(sim){
    pb$tick()
    return(chisq.test(data.frame(xsamples[,sim],ysamples[,sim]))$p.value)}) < alpha)/r
  CI <- exactci(round(power*r), r, conf.level = conf.level)$conf.int
  structure(list(power = power, conf.int = CI))
}

chisq.explore <- function(x, y = NULL, lown, topn, r = 10000, alpha = 0.05, conf.level = 0.95, plotit = TRUE){
  #' Compute the power using a chi-squared test by bootstrapping pilot data.
  #'
  #' \code{chisq.pexplore} computes (via simulation) the power of an experiment that will be analyzed using a chi-square test using a range of potential sample sizes.
  #'
  #' @inherit chisq.pow
  #' @param lown smallest sample size to test.
  #' @param topn largest sample size to test.
  #' @examples
  #' chisq.explore(c(0.5,0.25,0.25),lown=18,topn=20,r=500)
  if (any(is.na(lown) | (lown < 0))) stop("'lown' must be a nonnegative integer")
  if (any(is.na(topn) | (topn < 0))) stop("'topn' must be a nonnegative integer")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  if (any(is.na(conf.level) | (conf.level < 0) | (conf.level > 1))) stop("'alpha' must be between 0 and 1")
  x <- x/sum(x)
  if (!is.null(y)){
    y <- y/sum(y)
    if (length(y) != length(x)) stop("'x' and 'y' must have the same length")
  } else {
    y <- rep(1/length(x),length(x))
  }
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r*length(lown:topn), clear = FALSE, width = 60
  )
  GetPower <- function(samplesize){
    xsamples <- rmultinom(r, samplesize, x)
    ysamples <- rmultinom(r, samplesize, y)
    return(sum(sapply(1:r,function(sim){
      return(chisq.test(data.frame(xsamples[,sim],ysamples[,sim]))$p.value)}) < alpha)/r)
  }
  powerlist <- sapply(lown:topn, function(samplesize){
    pb$tick()
    return(GetPower(samplesize))
  })
  CIs <- data.frame(sapply(powerlist, function(x){return(exactci(round(x*r), r, conf.level = conf.level)$conf.int)}))
  results <- data.frame(
    samplesize = lown:topn,
    power = powerlist,
    LowCI = as.numeric(CIs[1,]),
    TopCI = as.numeric(CIs[2,])
  )
  if (plotit){
    return(ggplot(results, aes(x = samplesize, y = power))+geom_line()+geom_point()+theme_bw()+
             scale_x_continuous("Sample size", breaks = c(lown:topn))+scale_y_continuous("Power")+geom_errorbar(aes(ymin=LowCI,ymax=TopCI),width=0.2))
  } else {
    return(results)
  }
}

chisq.ppow <- function(x, y = NULL, n, r = 10000, alpha = 0.05, conf.level = 0.95){
  #' Compute the power using a chi-squared test by bootstrapping pilot data.
  #'
  #' \code{chisq.ppow} computes (via simulation) the power of an experiment that will be analyzed using a chi-square test by bootstrapping pilot data.
  #' This function directly calls \code{chisq.pow} but is included to maintaing consistency with other power functions in the package.
  #'
  #' @inherit chisq.pow
  #' @param x a list containing count data of pilot observations.
  #' @param y an optional list with the count observations of a second group.
  #' @examples
  #' chisq.ppow(x=c(3, 5, 2), n=32) # Chi-squared distribution against chance with 32 observations.
  #' chisq.ppow(x=c(2, 2, 2), y=c(8, 14, 1)) # Compare both group with 16 observations in each group.
  return(chisq.pow(x, y, n, r, alpha, conf.level))
}

chisq.pexplore <- function(x, y = NULL, lown, topn, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, plotit = TRUE){
  #' Explore power as a function of sample sie using a one- or two-sample unpaired t-test by bootstrapping pilot data.
  #'
  #' \code{ttest.pexplore} computes (via simulation) the power of an experiment that will be analyzed using a t-test for a range of sample sizes.
  #' Rather than taking a theoretical distribution, this function takes empirical data and bootstraps them to calculate the power.
  #' For an equivalent function that does not rely on pilot data see \link{ttestpower}.
  #'
  #' @inherit ttest.pow
  #' @param lown smallest sample size to explore.
  #' @param topn largest sample size to explore.
  #' @param plotit logical (default=\code{TRUE}) value. Function generates a plot when \code{TRUE} and returns a data frame otherwise.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' ttest.pexplore(x=c(0, 5, 10), lown=16, topn=24) # Power for a one-sample t-test with n in 16-24. Pilot data consists of three data points.
  #' ttest.pexplore(x=c(0, 5, 10), lown=16, topn=24,mu = -5) # Same as above, changing the avarege under the null to -5.
  #' ttest.pexplore(x=c(0, 5, 10), lown=16, topn=24, y=c(9, 3, 2, 1)) # Power for a two-sample t-test with n=16-24 (per condition) using unbalanced pilot data.
  if (any(sapply(x, class)=="factor")) stop("Data frame must not have categorical data")
  if (any(is.na(lown) | (lown < 0))) stop("'lown' must be nonnegative and integer")
  if (any(is.na(topn) | (topn < 0))) stop("'topn' must be nonnegative and integer")
  if (is.data.frame(x)){
    if (ncol(x) != 2) stop("When x is a data frame it must have exactly two columns")
    xvals = x[,1]
    yvals = y[,2]
    twotail = TRUE
  }
  else {
    xvals = x
    if (is.null(y)){
      twotail = FALSE
    } else{
      yvals = y
      twotail = TRUE
    }
  }
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r*length(lown:topn), clear = FALSE, width = 60
  )
  if (twotail){
    getpower <- function(samplesize){
      return(sum(sapply(1:r,function(z){
        pb$tick()
        return(t.test(sample(xvals, size = samplesize, replace = TRUE), sample(yvals, size = samplesize, replace = TRUE), alternative = alternative)$p.value)
      }) < alpha)/r)
    }
    powerlist = sapply(lown:topn, getpower)
  }
  else{
    if (is.null(mu)) {
      message("setting 'mu' to 0 (value according to null hypothesis). Change using 'mu = x'")
      mu = 0
    }
    getpower <- function(samplesize){
      return(sum(sapply(1:r,function(x){
      pb$tick()
      return(t.test(sample(xvals, size = samplesize, replace = TRUE), alternative = alternative, mu = mu)$p.value)
      }) < alpha)/r)}
    powerlist = sapply(lown:topn, getpower)
  }
  CIs <- data.frame(sapply(powerlist, function(x){return(exactci(round(x*r), r, conf.level = 0.95)$conf.int)}))
  results <- data.frame(
    samplesize = lown:topn,
    power = powerlist,
    LowCI = as.numeric(CIs[1,]),
    TopCI = as.numeric(CIs[2,])
  )
  if (plotit){
    return(ggplot(results, aes(x = samplesize, y = power))+geom_line()+geom_point()+theme_bw()+
             scale_x_continuous("Sample size", breaks = c(lown:topn))+scale_y_continuous("Power")+geom_errorbar(aes(ymin=LowCI, ymax=TopCI), width=0.2))
  } else {
    return(results)
  }
}
