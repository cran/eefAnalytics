% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eefAnalytics_14_09_2020.R
\name{crtBayes}
\alias{crtBayes}
\title{Bayesian analysis of cluster randomised education trials using Vague Priors.}
\usage{
crtBayes(
  formula,
  random,
  intervention,
  adaptD = NULL,
  nsim = 2000,
  data,
  threshold = 1:10/10,
  ...
)
}
\arguments{
\item{formula}{the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.}

\item{random}{a string variable specifying the "clustering variable" as contained in the data. See example below.}

\item{intervention}{a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.}

\item{adaptD}{As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.}

\item{nsim}{number of MCMC iterations per chain. Default is 2000.}

\item{data}{data frame containing the data to be analysed.}

\item{threshold}{a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).}

\item{...}{additional arguments of \code{\link[rstanarm]{stan_lmer}} to be passed to the function.}
}
\value{
S3 object; a list consisting of
\itemize{
\item \code{Beta}: Estimates and credible intervals for variables specified in the model.
\item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
\item \code{covParm}: A vector of variance decomposition into between cluster variance (Schools) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
\item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
\item \code{ProbES}: A matrix of Bayesian Posterior Probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
\item \code{Unconditional}: A list of unconditional effect sizes, covParm and ProbES obtained based on between and within cluster variances from the unconditional model (model with only the intercept as a fixed effect).
}
}
\description{
\code{crtBayes} performs analysis of cluster randomised education trials using a multilevel model under a Bayesian setting,
assuming vague priors.
}
\examples{
if(interactive()){

  data(crtData)

  ########################################################
  ## Bayesian analysis of cluster randomised trials     ##
  ########################################################

  output <- crtBayes(Posttest~ Intervention+Prettest,random="School",
                     intervention="Intervention",nsim=2000,data=crtData)

  ### Fixed effects
  beta <- output$Beta
  beta

  ### Effect size
  ES1 <- output$ES
  ES1

  ## Covariance matrix
  covParm <- output$covParm
  covParm

  ### plot random effects for schools

  plot(output)

  ### plot posterior probability of an effect size to be bigger than a pre-specified threshold

  plot(output,group=1)
}
}