utest.ppow <- function(x, y = NULL, n, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, conf.level = 0.95){
  #' Compute the power using a one- or two-sample unpaired u test by bootstrapping pilot data.
  #'
  #' \code{utest.ppow} computes (via simulation) the power of an experiment that will be analyzed using a u test (also called Wilcoxon test).
  #' Rather than taking a theoretical distribution, this function takes empirical data and bootstraps them to calculate the power.
  #'
  #' @inherit utest.pow
  #' @param x a data frame with two columns, or a list with pilot data.
  #' @param y a list with pilot data. When x is a list and y is not provided, a one-tailed t-test is used.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' utest.ppow(x=c(0, 5, 10), n=16) # Power for a one-sample u test with n=16. Pilot data consists of three data points.
  #' utest.ppow(x=c(0, 5, 10), n=16, mu = -5) # Same as above, changing the avarege under the null to -5.
  #' utest.ppow(x=c(0, 5, 10), y=c(9, 3, 2, 1), n=30) # Power for a two-sample u test with n=30 (per condition) using unbalanced pilot data.
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

utest.pexplore <- function(x, y = NULL, lown, topn, r = 10000, alternative = c("two.sided", "less", "greater"), mu = NULL, alpha = 0.05, plotit = TRUE, conf.level = 0.95){
  #' Explore power as a function of sample sie using a one- or two-sample unpaired u test by bootstrapping pilot data.
  #'
  #' \code{ttest.pexplore} computes (via simulation) the power of an experiment that will be analyzed using a u test  (also called Wilcoxon test) for a range of sample sizes.
  #' Rather than taking a theoretical distribution, this function takes empirical data and bootstraps them to calculate the power.
  #'
  #' @inherit utest.pow
  #' @param lown smallest sample size to explore.
  #' @param topn largest sample size to explore.
  #' @param plotit logical (default=\code{TRUE}) value. Function generates a plot when \code{TRUE} and returns a data frame otherwise.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description.
  #' @examples
  #' utest.pexplore(x=c(0, 5, 10), lown=16, topn=24) # Power for a one-sample u test with n in 16-24. Pilot data consists of three data points.
  #' utest.pexplore(x=c(0, 5, 10), lown=16, topn=24,mu = -5) # Same as above, changing the avarege under the null to -5.
  #' utest.pexplore(x=c(0, 5, 10), lown=16, topn=24, y=c(9, 3, 2, 1)) # Power for a two-sample u test with n=16-24 (per condition) using unbalanced pilot data.
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
