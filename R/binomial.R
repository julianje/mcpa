binom.pow <- function(n, p, r = 10000, alternative = c("two.sided", "less", "greater"), alpha = 0.05, nullp=0.5, conf.level = 0.95){
  #' Compute the power of a binomial experiment.
  #'
  #' \code{binom.pow} computes (via simulation) the power of a binomial experiment with a specified sample size and probability of success.
  #'
  #' @param n sample size.
  #' @param p predicted probability of success.
  #' @param r number of simulations to compute power.
  #' @param alternative type of alternative hypothesis in binomial test. Must be "\code{two.sided}" (default), "\code{greater}", or "\code{less}".
  #' @param alpha significance threshhold.
  #' @param nullp probability of success in null hypothesis.
  #' @param conf.level size of confidence intervals.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description and a 95% Clopper-Pearson confidence interval.
  #' @examples
  #' binom.pow(n=16, p=0.8)
  #' binom.pow(n=16, p=0.8, alternative="greater")
  #' binom.pow(n=32, p=0.6, r=5000, nullp=0.25)
  #' @seealso \code{\link{binom.power}}, \code{\link{binom.ppow}}, \code{\link{binom.explore}}, and \code{\link{binom.pexplore}}.
  if (any(is.na(n) | (n < 0))) stop("'n' must be a nonnegative integer")
  if (any(is.na(p) | (p < 0) | (p > 1))) stop("'p' must be between 0 and 1")
  if (any(is.na(nullp) | (nullp < 0) | (nullp > 1))) stop("'nullp' must be between 0 and 1")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  alternative <- match.arg(alternative)
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  samples <- rbinom(r, n, p)
  power <- switch(alternative,
         "less" = {sum(samples < qbinom(alpha, n, nullp))/r},
         "greater" = {sum(samples > qbinom(alpha, n, nullp, lower.tail = FALSE))/r},
         "two.sided" = {(sum(samples < qbinom(alpha/2, n, nullp))+sum(samples > qbinom(alpha/2, n, nullp, lower.tail = FALSE)))/r})
  CI <- exactci(round(power*r), r, conf.level = conf.level)$conf.int
  message(paste("Power estimate using",r,"samples."))
  structure(list(power = power, conf.int = CI))
}

binom.ppow <- function(n, pilotdata, r = 10000, alternative = c("two.sided", "less", "greater"), alpha = 0.05, nullp=0.5, conf.level = 0.95){
  #' Compute the power of a binomial experiment using pilot data.
  #'
  #' \code{binom.ppow} computes (through simulation) the power of a binomial experiment with a specified sample size, using pilot data.
  #'
  #' @param pilotdata a list with numerical or categorical pilot data. Largest number of second category is treated as success.
  #' @inherit binom.pow
  #' @examples
  #' binom.ppow(n=16, pilotdata = (1, 1, 0, 0, 1, 1), p=0.8)
  #' binom.ppow(n=16, pilotdata = ("a", "b", "b", "b"), p=0.8, alternative="greater")
  #' binom.ppow(n=16, pilotdata = ("a", "b", "b", "b"), r=5000, nullp=0.25, alternative="greater")
  if (any(is.na(pilotdata) | length(unique(pilotdata)) != 2)) stop("'pilotdata' must be a list with two types of values (cannot simulate power if pilot data has no variance).")
  probs = prop.table(table(pilotdata))
  if (names(probs)[0]!="0" || names(probs)[1]!="1"){
    message(paste("Interpreting '",names(probs)[2],"' as success and '",names(probs)[1],"' as failure.",sep=""))
  }
  binom.pow(n, probs[2], r, alternative, alpha, nullp, conf.level)
}

binom.explore <- function(lown, topn, p, r = 10000, alternative = c("two.sided", "less", "greater"), alpha = 0.05, nullp=0.5, conf.level = 0.95, plotit=TRUE){
  #' Explore the power of a binomial experiment under different sample sizes.
  #'
  #' \code{binom.explore} computes (through simulation) the power of a binomial experiment under different sample sizes.
  #'
  #' @inherit binom.pow
  #' @param lown smallest sample size to explore.
  #' @param topn largest sample size to explore.
  #' @param plotit logical (default=\code{TRUE}) value. Function generates a plot when \code{TRUE} and returns a data frame otherwise.
  #' @examples
  #' binom.explore(lown=16, topn=24, p=0.8)
  #' binom.explore(lown=16, topn=24, p=0.8, alternative="greater")
  #' binom.explore(lown=16, topn=24, p=0.6, r=5000, nullp=0.25)
  #'
  if (any(is.na(lown) | (lown < 0))) stop("'lown' must be nonnegative and integer")
  if (any(is.na(topn) | (topn < 0))) stop("'topn' must be nonnegative and integer")
  if (topn < lown) stop("'topn' must be higher than 'lown'")
  if (any(is.na(p) | (p < 0) | (p > 1))) stop("'p' must be between 0 and 1")
  if (any(is.na(nullp) | (nullp < 0) | (nullp > 1))) stop("'nullp' must be between 0 and 1")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  alternative <- match.arg(alternative)
  if (any(is.na(alpha) | (alpha < 0) | (alpha > 1))) stop("'alpha' must be between 0 and 1")
  getpower <- function(n, r, p, nullp, alternative){
    samples <- rbinom(r, n, p)
    power <- switch(alternative,
                    "less" = {sum(samples < qbinom(alpha, n, nullp))/r},
                    "greater" = {sum(samples > qbinom(alpha, n, nullp, lower.tail = FALSE))/r},
                    "two.sided" = {(sum(samples < qbinom(alpha/2, n, nullp))+sum(samples > qbinom(alpha/2, n, nullp, lower.tail = FALSE)))/r})
    return(power)
  }
  powervalues <- sapply(lown:topn,function(x){return(getpower(x, r, p, nullp, alternative))})
  CIs <- data.frame(sapply(powervalues, function(x){return(exactci(round(x*r), r, conf.level = conf.level)$conf.int)}))
  results <- data.frame(samplesize = lown:topn,
             power = powervalues,
             LowCI = as.numeric(CIs[1,]),
             TopCI = as.numeric(CIs[2,]))
  message(paste("Power estimate using",r,"samples per group."))
  if (plotit){
    return(ggplot(results, aes(x = samplesize, y = power))+
             geom_line()+geom_point()+theme_bw()+
      scale_x_continuous("Sample size", breaks=c(lown:topn))+scale_y_continuous("Power")+geom_errorbar(aes(ymin=LowCI,ymax=TopCI),width=0.2))
  } else {

    return(results)
  }
}

binom.pexplore <- function(lown, topn, pilotdata, r = 10000, alternative = c("two.sided", "less", "greater"), alpha = 0.05, nullp=0.5, conf.level = 0.95, plotit=TRUE){
  #' Explore the power of a binomial experiment under different sample sizes using pilot data.
  #'
  #' \code{binom.pexplore} computes (through simulation) the power of a binomial experiment under different sample sizes.
  #' Rather than taking a probability of success (like \code{binom.explore}), \code{binom.pexplore} takes a vector of pilot data.
  #'
  #' @inherit binom.pow
  #' @param lown smallest sample size to explore.
  #' @param topn largest sample size to explore.
  #' @param plotit logical (default=\code{TRUE}) value. Function generates a plot when \code{TRUE} and returns a data frame otherwise.
  #' @examples
  #' binom.pexplore(lown=16, topn=24, pilotdata = (1, 1, 0, 0, 1, 1), p=0.8)
  #' binom.pexplore(lown=16, topn=24, pilotdata = ("a", "b", "b", "b"), p=0.8, alternative="greater")
  #' binom.pexplore(lown=16, topn=24, pilotdata = ("a", "b", "b", "b"), r=5000, nullp=0.25, alternative="greater")
  if (any(is.na(pilotdata) | length(unique(pilotdata)) != 2)) stop("'pilotdata' must be a list with two types of values (cannot simulate power if pilot data has no variance).")
  probs = prop.table(table(pilotdata))
  if (names(probs)[0]!="0" || names(probs)[1]!="1"){
    message(paste("Interpreting '",names(probs)[2],"' as success and '",names(probs)[1],"' as failure.",sep=""))
  }
  return(binom.explore(lown, topn, probs[2], r, alternative, alpha, nullp, conf.level, plotit))
}
