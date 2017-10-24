ttest.pow <- function(means, var, n, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, conf.level = 0.95){
  #' Compute the power using a one- or two-sample unpaired t-test.
  #'
  #' \code{ttest.pow} computes (via simulation) the power of an experiment that will be analyzed using a t-test. When two means are provided,
  #' function assumes a two-sample unpaired t-test, and \code{n} is interpreted as the sample size of each group (for a total sample size or \code{2n}).
  #'
  #' @param means either a list with two average values (computes a two-sample t-test) or a single value (computes a one-sample t-test).
  #' @param var expected variance in each group.
  #' @param n sample size.
  #' @param r number of simulations to compute power.
  #' @param alternative type of alternative hypothesis in binomial test. Must be "\code{two.sided}" (default), "\code{greater}", or "\code{less}".
  #' @param alpha significance threshhold.
  #' @param mu mean value according to null hypothesis (default = \code{0}). Only used in one sample t-tests.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' ttest.pow(means=c(5, 10), var=10, n=16) # two-sample t-test. n=16 refers to each condition, for a total of 32.
  #' ttest.pow(means=20, var=10, n=16) # one-sample t-test. Comparing if average is different from 0. Because there is only condition, the total sample isze is 16.
  #' ttest.pow(means=20, var=10, n=16, mu=10, alternative="higher") # one-sample t-test. Comparing if average is higher than 10.
  #' @seealso \code{\link{ttest.pow}}, \code{\link{ttest.ppow}}, \code{\link{ttest.explore}}, and \code{\link{ttest.pexplore}}.
  if (any(is.na(n) | (n < 0))) stop("'n' must be a nonnegative integer")
  if (any(is.na(var) | (var < 0))) stop("'var' must be positive")
    if (length(means) != 2) {
        if (length(means) != 1 || is.null(mu)){
          message("setting 'mu' to 0 (value according to null hypothesis). Change using 'mu = x'")
          mu = 0
        }
  }
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  alternative <- match.arg(alternative)
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r, clear = FALSE, width = 60
  )
  if (length(means) == 2){
    power <- (sum(sapply(1:r,function(x){
      pb$tick()
      return(t.test(rnorm(n, means[1], sqrt(var)), rnorm(n, means[2], sqrt(var)), alternative=alternative)$p.value)
    }) <= alpha)/r)
  }
  else{
    power <- (sum(sapply(1:r,function(x){return(t.test(rnorm(n, means, sqrt(var)), mu = mu, alternative = alternative)$p.value)}) <= alpha)/r)
  }
  CI <- exactci(round(power*r), r, conf.level = conf.level)$conf.int
  structure(list(power = power, conf.int = CI))
}

ttest.explore <- function(lown, topn, means, var, r = 10000, alternative = c("two.sided", "less", "greater"), mu = 0, alpha = 0.05, conf.level = 0.95, plotit = TRUE){
  #' Explore the power using a t-test under different sample sizes.
  #'
  #' \code{ttest.explore} computes (through simulation) the power of an experiment that will be analyzed using a t-test for a set of potential sample sizes. When two means are provided,
  #' function assumes a two-sample unpaired t-test, and \code{n} is interpreted as the sample size of each group (for a total sample size or \code{2n}).
  #'
  #' @inherit ttest.pow
  #' @param lown smallest sample size to explore.
  #' @param topn largest sample size to explore.
  #' @param plotit logical (default=\code{TRUE}) value. Function generates a plot when \code{TRUE} and returns a data frame otherwise.
  #' @examples
  #' ttest.explore(lown=10, topn=15, means=c(5, 10), var=10) # two-sample t-test. Effective sample sizes are 20 to 30 (10 to 15 per group)
  #' ttest.explore(lown=10, topn=15, means=20, var=10) # one-sample t-test. Comparing if average is different from 0.
  #' ttest.explore(lown=10, topn=15, means=20, var=10, mu=10, alternative="higher") # one-sample t-test. Comparing if average is higher than 10.
  if (any(is.na(lown) | (lown < 0))) stop("'lown' must be nonnegative and integer")
  if (any(is.na(topn) | (topn < 0))) stop("'topn' must be nonnegative and integer")
  if (topn < lown) stop("'topn' must be higher than 'lown'")
  if (any(is.na(var) | (var < 0))) stop("'var' must be positive")
  if (length(means) != 2) {
    if (length(means) != 1 || is.null(mu)){
      stop("'means' should have two values, or 'means' should have one value and 'mu' should be provided")
    }
  }
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  alternative <- match.arg(alternative)
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r*length(lown:topn), clear = FALSE, width = 60
  )
  if (length(means) == 2){
    getpower <- function(n, means, var, alternative, alpha, r){
      return(sum(sapply(1:r, function(x){
        pb$tick()
        return(t.test(rnorm(n, means[1], sqrt(var)), rnorm(n, means[2], sqrt(var)), alternative=alternative)$p.value)
      }) <= alpha)/r)
    }
    powerlist <- sapply(lown:topn, function(samplesize){
      return(getpower(samplesize, means, var, alternative, alpha, r))})
  }
  else{
    getpower <- function(n, meanval, var, mu, alternative, alpha, r){
      simulations <- sapply(1:r, function(x){
        pb$tick()
        return(t.test(rnorm(n, meanval, sqrt(var)), mu = mu, alternative = alternative)$p.value)
      })
      return(sum(simulations < alpha)/r)
    }
    powerlist <- sapply(lown:topn, function(samplesize){
      return(getpower(samplesize, means[1], var, mu, alternative, alpha, r))
    })
  }
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

ttest.ppow <- function(x, y = NULL, n, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, conf.level = 0.95){
  #' Compute the power using a one- or two-sample unpaired t-test by bootstrapping pilot data.
  #'
  #' \code{ttest.ppow} computes (via simulation) the power of an experiment that will be analyzed using a t-test.
  #' Rather than taking a theoretical distribution, this function takes empirical data and bootstraps them to calculate the power.
  #' For an equivalent function that does not rely on pilot data see \link{ttest.pow}.
  #'
  #' @inherit ttest.pow
  #' @param x a data frame with two columns, or a list with pilot data.
  #' @param y a list with pilot data. When x is a list and y is not provided, a one-tailed t-test is used.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' ttest.ppow(x=c(0, 5, 10), n=16) # Power for a one-sample t-test with n=16. Pilot data consists of three data points.
  #' ttest.ppow(x=c(0, 5, 10), n=16, mu = -5) # Same as above, changing the avarege under the null to -5.
  #' ttest.ppow(x=c(0, 5, 10), y=c(9, 3, 2, 1), n=30) # Power for a two-sample t-test with n=30 (per condition) using unbalanced pilot data.
  if (any(sapply(x, class)=="factor")) stop("Data frame must not have categorical data")
  if (any(is.na(n) | (n < 0))) stop("'n' must be a nonnegative integer")
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
    total = r, clear = FALSE, width = 60
  )
  if (twotail){
    power <- (sum(sapply(1:r,function(z){
      pb$tick()
      return(t.test(sample(xvals, size = n, replace = TRUE), sample(yvals, size = n, replace = TRUE), alternative = alternative)$p.value)
      }) < alpha)/r)
  }
  else{
    if (is.null(mu)) {
      message("setting 'mu' to 0 (value according to null hypothesis). Change using 'mu = x'")
      mu = 0
    }
    power <- (sum(sapply(1:r,function(x){
      pb$tick()
      return(t.test(sample(xvals, size = n, replace = TRUE), alternative = alternative, mu = mu)$p.value)
      }) < alpha)/r)
  }
  CI <- exactci(power*r, r, conf.level = conf.level)$conf.int
  structure(list(power = power, conf.int = CI))
}

ttest.pexplore <- function(x, y = NULL, lown, topn, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, conf.level = 0.95, plotit = TRUE){
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
  CIs <- data.frame(sapply(powerlist, function(x){return(exactci(round(x*r), r, conf.level = conf.level)$conf.int)}))
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
