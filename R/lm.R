lm.pow <- function(lmodel, n, r = 10000, alpha = 0.05, conf.level = 0.95){
  #' Compute the power of a linear model.
  #'
  #' \code{lm.pow} computes (via simulation) the power of a linear model trained with pilot data.
  #'
  #' @param n sample size.
  #' @param lmodel linear model trained through lm object.
  #' @param r number of simulations to compute power.
  #' @param alpha significance threshhold.
  #' @param conf.level size of confidence intervals.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description and a 95\% Clopper-Pearson confidence interval for each item in the equation.
  #' @examples
  #' library(datasets)
  #' Model <- lm(mpg ~ hp + vs, data = mtcars)
  #' lm.pow(n=16, Model)
  if (any(is.na(n) | (n < 0))) stop("'n' must be a nonnegative integer")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  PilotData <-lmodel$model
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r, clear = FALSE, width = 60
  )
  simoutcomes <- replicate(r, {pb$tick(); as.numeric(summary(update(lmodel, data = dplyr::sample_n(PilotData, size = n, replace = TRUE)))$coefficients[,4]) < alpha})
  powerlist <- sapply(1:nrow(simoutcomes), function(x){return(sum(simoutcomes[x,]))})/r
  CIs <- data.frame(sapply(powerlist, function(x){return(exactci(round(x*r), r, conf.level = conf.level)$conf.int)}))
  results <- data.frame(
    power = powerlist,
    LowCI = as.numeric(CIs[1,]),
    TopCI = as.numeric(CIs[2,])
  )
  rownames(results) = names(lmodel$coefficients)
  message(paste("Power estimate using",r,"samples."))
  return(results)
}

lm.explore <- function(lmodel, lown, topn, r = 10000, alpha = 0.05, conf.level = 0.95, plotit = TRUE){
  #' Compute the power of a linear model for a range of sample sizes.
  #'
  #' \code{lm.explore} computes (via simulation) the power of a linear model trained with pilot data.
  #'
  #' @inherit lm.pow
  #' @param lown lowest sample size.
  #' @param topn highest sample size.
  #' @return The probability of finding \eqn{p < \alpha} with the experiment description and a 95\% Clopper-Pearson confidence interval for each item in the equation.
  #' @examples
  #' library(datasets)
  #' Model <- lm(mpg ~ hp + vs, data = mtcars)
  #' lm.pow(n=16, Model)
  if (any(is.na(lown) | (lown < 0))) stop("'lown' must be a nonnegative integer")
  if (any(is.na(topn) | (topn < 0))) stop("'lown' must be a nonnegative integer")
  if (topn <= lown) stop("topn must be greater than lown")
  if (any(is.na(r) | (r < 0))) stop("'r' must be nonnegative and integer")
  PilotData <-lmodel$model
  pb <- progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = r*length(lown:topn), clear = FALSE, width = 60
  )
  modelrownames = names(lmodel$coefficients)
  GetLmPower <- function(n){
    simoutcomes <- replicate(r, {pb$tick(); as.numeric(summary(update(lmodel, data = dplyr::sample_n(PilotData, size = n, replace = TRUE)))$coefficients[,4]) < alpha})
    powerlist <- sapply(1:nrow(simoutcomes), function(x){return(sum(simoutcomes[x,]))})/r
    CIs <- data.frame(sapply(powerlist, function(x){return(exactci(round(x*r), r, conf.level = conf.level)$conf.int)}))
    return(results <- data.frame(
      variable = modelrownames,
      samplesize = rep(n,length(powerlist)),
      power = powerlist,
      LowCI = as.numeric(CIs[1,]),
      TopCI = as.numeric(CIs[2,])
    ))
  }
  results <- do.call(rbind,lapply(lown:topn,function(x){return(GetLmPower(x))}))
  message(paste("Power estimate using",r,"samples."))
  if (plotit){
    return(ggplot(results, aes(x=samplesize, y=power,color=variable,group=variable))+geom_point()+geom_line()+geom_errorbar(aes(ymin=LowCI,ymax=TopCI))+theme_bw())+
    scale_x_discrete("sample size")
  } else{
    return(results)
  }
}
