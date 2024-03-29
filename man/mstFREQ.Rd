% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mstFREQ.R
\name{mstFREQ}
\alias{mstFREQ}
\title{Analysis of Multisite Randomised Education Trials using Multilevel Model under a Frequentist Setting.}
\usage{
mstFREQ(
  formula,
  random,
  intervention,
  baseln,
  nPerm,
  data,
  type,
  ci,
  seed,
  nBoot
)
}
\arguments{
\item{formula}{the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.}

\item{random}{a string variable specifying the "clustering variable" as contained in the data. See example below.}

\item{intervention}{a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.}

\item{baseln}{A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.}

\item{nPerm}{number of permutations required to generate permutated p-value.}

\item{data}{data frame containing the data to be analysed.}

\item{type}{method of bootstrapping including case re-sampling at student level "case(1)", case re-sampling at school level "case(2)", case re-sampling at both levels "case(1,2)" and residual bootstrapping using "residual". If not provided, default will be case re-sampling at student level.}

\item{ci}{method for bootstrap confidence interval calculations; options are the Basic (Hall's) confidence interval "basic" or the simple percentile confidence interval "percentile". If not provided default will be percentile.}

\item{seed}{seed required for bootstrapping and permutation procedure, if not provided default seed will be used.}

\item{nBoot}{number of bootstraps required to generate bootstrap confidence intervals.}
}
\value{
S3 object; a list consisting of
\itemize{
\item \code{Beta}: Estimates and confidence intervals for variables specified in the model.
\item \code{ES}: Conditional Hedge's g effect size (ES) and its 95\% confidence intervals. If nBoot is not specified, 95\% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
\item \code{covParm}: A list of variance decomposition into between cluster variance-covariance matrix (schools and school by intervention) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
\item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept and intervention slope.
\item \code{Perm}: A "nPerm x 2w" matrix containing permutated effect sizes using residual variance and total variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only when \code{nPerm} is specified.
\item \code{Bootstrap}: A "nBoot x 2w" matrix containing the bootstrapped effect sizes using residual variance (Within) and total variance (Total). "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is only produced when \code{nBoot} is specified.
\item \code{Unconditional}: A list of unconditional effect sizes, covParm, Perm and Bootstrap obtained based on variances from the unconditional model (model with only the intercept as a fixed effect).
}
}
\description{
\code{mstFREQ} performs analysis of multisite randomised education trials using a multilevel model under a frequentist setting.
}
\examples{
if(interactive()){

data(mstData)

###############################################
## MLM analysis of multisite trials + 1.96SE ##
###############################################

output1 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=mstData)


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
## MLM analysis of multisite trials             ##
## with residual bootstrap confidence intervals ##
##################################################

output2 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,type="residual",data=mstData)

tp <- output2$Bootstrap
### Effect size

ES2 <- output2$ES
ES2

### plot bootstrapped values

plot(output2, group=1)

#######################################################################
## MLM analysis of mutltisite trials with permutation p-value##
#######################################################################

output3 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=mstData)

ES3 <- output3$ES
ES3

#### plot permutated values

plot(output3, group=1)
}
}
