% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crtFREQ.R
\name{crtFREQ}
\alias{crtFREQ}
\title{Analysis of Cluster Randomised Education Trials using Multilevel Model under a Frequentist Setting.}
\usage{
crtFREQ(
  formula,
  random,
  intervention,
  baseln,
  nPerm,
  nBoot,
  type,
  ci,
  seed,
  data
)
}
\arguments{
\item{formula}{the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.}

\item{random}{a string variable specifying the "clustering variable" as contained in the data. See example below.}

\item{intervention}{a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.}

\item{baseln}{A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.}

\item{nPerm}{number of permutations required to generate a permutated p-value.}

\item{nBoot}{number of bootstraps required to generate bootstrap confidence intervals.}

\item{type}{method of bootstrapping including case re-sampling at student level "case(1)", case re-sampling at school level "case(2)", case re-sampling at both levels "case(1,2)" and residual bootstrapping using "residual". If not provided, default will be case re-sampling at student level.}

\item{ci}{method for bootstrap confidence interval calculations; options are the Basic (Hall's) confidence interval "basic" or the simple percentile confidence interval "percentile". If not provided default will be percentile.}

\item{seed}{seed required for bootstrapping and permutation procedure, if not provided default seed will be used.}

\item{data}{data frame containing the data to be analysed.}
}
\value{
S3 object; a list consisting of
\itemize{
\item \code{Beta}: Estimates and confidence intervals for variables specified in the model.
\item \code{ES}: Conditional Hedges' g effect size and its 95\% confidence intervals. If nBoot is not specified, 95\% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
\item \code{covParm}: A vector of variance decomposition into between cluster variance (Schools) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
\item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
\item \code{Perm}: A "nPerm x 2w" matrix containing permutated effect sizes using residual variance and total variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only when \code{nPerm} is specified.
\item \code{Bootstrap}: A "nBoot x 2w" matrix containing the bootstrapped effect sizes using residual variance (Within) and total variance (Total). "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is only produced when \code{nBoot} is specified.
\item \code{Unconditional}: A list of unconditional effect sizes, covParm, Perm and Bootstrap obtained based on variances from the unconditional model (model with only the intercept as a fixed effect).
}
}
\description{
\code{crtFREQ} performs analysis of cluster randomised education trials using a multilevel model under a frequentist setting.
}
\examples{
if(interactive()){

data(crtData)

########################################################
## MLM analysis of cluster randomised trials + 1.96SE ##
########################################################

output1 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=crtData)


### Fixed effects
beta <- output1$Beta
beta

### Effect size
ES1 <- output1$ES
ES1

## Covariance matrix
covParm <- output1$covParm
covParm

### plot random effects for schools

plot(output1)

##################################################
## MLM analysis of cluster randomised trials    ##
## with residual bootstrap confidence intervals ##
##################################################

output2 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,type="residual",data=crtData)


### Effect size

ES2 <- output2$ES
ES2

### plot bootstrapped values

plot(output2, group=1)

#######################################################################
## MLM analysis of cluster randomised trials with permutation p-value##
#######################################################################

output3 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=crtData)

### Effect size

ES3 <- output3$ES
ES3


### plot permutated values

plot(output3, group=1)
}
}
