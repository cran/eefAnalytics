
#' Multisite Trial Data.
#'
#' A multisite trial dataset containing 54 schools. This data contains a random sample of test data of pupils and not actual trial data.
#'
#' \itemize{
#'   \item Posttest: posttest scores
#'   \item Prettest: prettest scores
#'   \item Intervention: the indicator for the intervention groups in a two arm trial,
#' coded as 1 for intervention group and 0 for control group. The intervention variable should always be numeric with the control group as the first value.
#'   \item Intervention2: a simulated indicator for intervention groups in a three arm trial. The intervention variable should always be numeric with the control group as the first value.
#'   \item School: numeric school identifier
#' }
#'
#' @format A data frame with 210 rows and 5 variables
#' @name mstData
NULL

###iwq

#' Cluster Randomised Trial Data.
#'
#' A cluster randomised trial dataset containing 22 schools. The data contains a random sample of test data of pupils and not actual trial data.
#'
#' \itemize{
#'   \item Posttest: posttest scores
#'   \item Prettest: prettest scores
#'   \item Intervention: the indicator for intervention groups in a two arm trial,
#'coded as 1 for intervention group and 0 for control group. The intervention variable should always be numeric with the control group as the first value.
#'   \item Intervention2: a simulated indicator for intervention groups in a three arm trial. The intervention variable should always be numeric with the control group as the first value.

#'   \item School: numeric school identifier
#' }
#'
#' @format A data frame with 265 rows and 5 variables
#' @name crtData
NULL


## IMPORTS ##
#' @importFrom lme4 lmer ranef VarCorr
#' @importFrom mvtnorm rmvnorm
#' @importFrom metafor forest
#' @importFrom rstanarm stan_glm stan_lmer stan_glmer
#' @importFrom graphics abline barplot hist legend lines par plot points text mtext title
#' @importFrom stats confint lm model.frame model.matrix model.response quantile rnorm na.omit update na.omit as.formula
NULL


#############################################################################
############# SRT main functions ################################################


#' Analysis of Simple Randomised Education Trial using Linear Regression Model.
#'
#' \code{srtFREQ} performs analysis of educational trials under the assumption of independent errors among pupils.
#' This can also be used with schools as fixed effects.
#'
#' @export
#' @param formula the model to be analysed is of the form y~x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals. Default is NULL.
#' @param nPerm number of permutations required to generate permutated p-value.  Default is NULL.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and confidence intervals for the variables specified in the model.
#' \item \code{ES}: Conditional Hedges'g effect size and its 95% confidence intervals. If nBoot is not specified, 95% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{sigma2}: Residual variance.
#' \item \code{Perm}: A "nPerm x w" matrix containing permutated effect sizes using residual variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only if \code{nPerm} is specified.
#' \item \code{Bootstrap}: A "nBoot x w" matrix containing the bootstrapped effect sizes using residual variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only if \code{nBoot} is specified.
#' \item \code{Unconditional}: A list of unconditional effect size, sigma2, Perm and Bootstrap obtained based on variances from the unconditional model (model with only intercept as fixed effect).
#' }
#' @example inst/examples/srtExample.R
srtFREQ <- function(formula,intervention,nBoot=NULL,nPerm=NULL,data)UseMethod("srtFREQ")

#' @export
srtFREQ.default <- function(formula,intervention,nBoot=NULL,nPerm=NULL,data){stop("No correct formula input given.")}

#' @export
srtFREQ.formula <- function(formula,intervention,nBoot=NULL,nPerm=NULL,data=data){

  if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(is.null(nPerm)){nPerm <-0}
  if(is.null(nBoot)){nBoot <-0}
  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])
  mf <- model.frame(formula=formula, data=data)
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]

  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  output <- srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)
  output$ES <- sapply(output$ES, function(x) round(x,2), simplify = F)

  if(nPerm>0){

    if(nPerm<1000){stop("nPerm must be greater than 1000")}

    permES1<- matrix(NA,nPerm,(length(unique(trt))-1))
    permES2<- matrix(NA,nPerm,(length(unique(trt))-1))
    set.seed(1020252)
    for (i in 1:nPerm){

      data[,which(colnames(data)==intervention)]<- sample(trt)
      fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))

      p2CRTFREQ <-srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)

      permES1[i,]  <-  p2CRTFREQ$ES$Conditional[,1]
      permES2[i,]  <-  p2CRTFREQ$ES$Unconditional[,1]
    }
    permES1=data.frame(permES1)
    permES2=data.frame(permES2)
    names(permES1) <- row.names(output$ES$Conditional)
    names(permES2) <- row.names(output$ES$Conditional)

    Perm<- data.frame(cbind(permES_cond=permES1, permES_uncond=permES2))
    perm.names <- paste0( rep(c("cond", "uncond"),each=dim(permES1)[2]),"_" ,rep(row.names(output$ES$Conditional),dim(permES1)[2]))
    names(Perm) <- perm.names
    output$Perm<-Perm
  }

  if(nBoot > 0){

    if(nBoot<1000){stop("nBoot must be greater than 1000")}


    tid <- c(1:nrow(fixedDesignMatrix))

    set.seed(1020252)
    bootSamples <- sapply(c(1:nBoot),function(x)sample(tid,replace=TRUE))

    bootResults <- sapply(1:dim(bootSamples)[2],function(bt)srt.srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,bt=bootSamples[,bt]), simplify = F)
    bootResults2 <- data.frame(do.call(rbind, bootResults))

    bootES <- apply(bootResults2,2,function(x)quantile(x,prob=c(0.025,0.975),na.rm=TRUE))
    all.ES <- as.data.frame(do.call(rbind, output$ES))[,1]
    bootES2.1 <- round(data.frame(t(rbind(all.ES,bootES))),2)
    names(bootES2.1)=colnames(output$ES$Conditional)

    condName <- grep("uncond", names(bootResults2), invert = T)
    UncondName <- grep("uncond", names(bootResults2), invert = F)
    bootES2  <- list(Conditional=bootES2.1[condName,], Unconditional=bootES2.1[UncondName,])
    bootES2 <- lapply(bootES2, function(x){ row.names(x)<-row.names(output$ES$Conditional); x})

    bootR1 <- data.frame(bootResults2[,condName])
    bootR2  <- data.frame(bootResults2[,UncondName])
    names(bootR1)<- row.names(output$ES$Conditional)
    names(bootR2)<- row.names(output$ES$Conditional)

    output$ES <-bootES2
    output$Bootstrap <- bootResults2
  }


  output1 <- list()
  output1$Beta   <- output$Beta
  output1$ES     <- output$ES$Conditional
  output1$sigma2 <- output$sigma2["Conditional"]
  if(nPerm > 0){output1$permES    <- permES1}
  if(nBoot > 0){output1$Bootstrap <- round(bootR1,2)}

  output1$Unconditional$ES     <- output$ES$Unconditional
  output1$Unconditional$sigma2 <- output$sigma2["Unconditional"]
  if(nPerm > 0){output1$Unconditional$permES <- permES2}
  if(nBoot > 0){output1$Unconditional$Bootstrap <- round(bootR2,2)}


  output1$Method <- "LM"
  class(output1) <- "eefAnalytics"
  return(output1)
}


#############################################################################
############# CRT main functions ################################################
#' Analysis of Cluster Randomised Education Trials using Multilevel Model under a Frequentist Setting.
#'
#' \code{crtFREQ} performs analysis of cluster randomised education trials using a multilevel model under a frequentist setting.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals. Default is NULL.
#' @param nPerm number of permutations required to generate a permutated p-value.  Default is NULL.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and confidence intervals for variables specified in the model.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% confidence intervals. If nBoot is not specified, 95% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{covParm}: A vector of variance decomposition into between cluster variance (Schools) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
#' \item \code{Perm}: A "nPerm x 2w" matrix containing permutated effect sizes using residual variance and total variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only when \code{nPerm} is specified.
#' \item \code{Bootstrap}: A "nBoot x 2w" matrix containing the bootstrapped effect sizes using residual variance (Within) and total variance (Total). "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is only produced when \code{nBoot} is specified.
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm, Perm and Bootstrap obtained based on variances from the unconditional model (model with only the intercept as a fixed effect).
#'  }
#' @example inst/examples/crtExample.R
crtFREQ<- function(formula,random,intervention,nPerm=NULL,nBoot=NULL,data)UseMethod("crtFREQ")
#' @export
crtFREQ.default<- function(formula,random,intervention,nPerm=NULL,nBoot=NULL,data){stop("No correct formula input given.")}


#' @export
crtFREQ.formula <- function(formula,random,intervention,nPerm=NULL,nBoot=NULL,data){

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  trt<-as.factor(trt)
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk >0){stop("This is not a CRT design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  stp2 <- as.numeric(apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))

  if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(is.null(nPerm)){nPerm <-0}
  if(is.null(nBoot)){nBoot <-0}
  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])

  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)


  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


  output <- crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)


  if(nPerm>0){
    if(nPerm<999){stop("nPerm must be greater than 1000")}
    output$Perm<- crt.perm(formula,data,stp,stp2,intervention,cluster,nPerm,random)
    output$Condtional$Perm <- round(data.frame(output$Perm$condtional),2)
    output$Unconditional$Perm <- round(data.frame(output$Perm$unconditional),2)
  }


  if(nBoot >0){

     if(nBoot<1000){stop("nBoot must be greater than 1000")}

    tid <- c(1:nrow(fixedDesignMatrix))

    set.seed(1020252)
    bootSamples <- NULL

    for(ii in 1:length(unique(cluster))){
      selID <- tid[cluster==unique(cluster)[ii]]
      if(length(selID)>0){
        selID2<- sapply(c(1:nBoot),function(x)sample(selID,length(selID),replace=TRUE))
        bootSamples <- rbind(bootSamples ,selID2)
      }

    }


    bootResults <- apply(bootSamples ,2,function(bt)crt.crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster,bt=bt))
    bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults,intervention=intervention)
    output$ES <- bootES

    output$Bootstrap <- bootResults
    output$Condtional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Conditional))))
    output$Unconditional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Unconditional))))
    Bnames <- gsub("Estimate1","Within", gsub("Estimate2", "Total",names(output$Condtional$Bootstrap )))
    names(output$Condtional$Bootstrap) <- Bnames
    names(output$Unconditional$Bootstrap) <- Bnames
  }


  output1 <- list()
  output1$Beta   <- output$Beta
  output1$covParm <- output$covParm["Conditional",]
  output1$ES     <- output$ES$Conditional
  output1$SchEffects     <- output$SchEffects
  if(nPerm > 0){output1$permES    <- output$Condtional$Perm}
  if(nBoot > 0){output1$Bootstrap <- output$Condtional$Bootstrap}

  output1$Unconditional$ES     <- output$ES$Unconditional
  output1$Unconditional$covParm <- output$covParm["Unconditional",]
  if(nPerm > 0){output1$Unconditional$permES <- output$Unconditional$Perm}
  if(nBoot > 0){output1$Unconditional$Bootstrap <-  output$Unconditional$Bootstrap}


  output1$Method <- "MLM"

  class(output1) <- "eefAnalytics"
  return(output1)
}



#############################################################################
############# MST main functions ################################################

#############################################################################
############# MST main functions ################################################
#' Analysis of Multisite Randomised Education Trials using Multilevel Model under a Frequentist Setting.
#'
#' \code{mstFREQ} performs analysis of multisite randomised education trials using a multilevel model under a frequentist setting.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals. Default is NULL.
#' @param nPerm number of permutations required to generate permutated p-value.  Default is NULL.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and confidence intervals for variables specified in the model.
#' \item \code{ES}: Conditional Hedge's g effect size (ES) and its 95% confidence intervals. If nBoot is not specified, 95% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{covParm}: A list of variance decomposition into between cluster variance-covariance matrix (schools and school by intervention) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept and intervention slope.
#' \item \code{Perm}: A "nPerm x 2w" matrix containing permutated effect sizes using residual variance and total variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only when \code{nPerm} is specified.
#' \item \code{Bootstrap}: A "nBoot x 2w" matrix containing the bootstrapped effect sizes using residual variance (Within) and total variance (Total). "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is only prduced when \code{nBoot} is specified.
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm, Perm and Bootstrap obtained based on variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/mstExample.R
mstFREQ<- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL)UseMethod("mstFREQ")

#'@export
mstFREQ.default<- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL){stop("No correct formula input given.")}


#' @export
mstFREQ.formula <- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL){

  data <- na.omit(data[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data.frame(data)[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
  trt <- data[,which(colnames(data)==intervention)]
  trt <- as.factor(trt)
  tmp2 <- which(colnames(data)==random)
  cluster2 <- data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk ==0){stop("This is not a MST design")}

  if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(is.null(nPerm)){nPerm <-0}
  if(is.null(nBoot)){nBoot <-0}


  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])
  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention


  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


  output <- rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster)

  if(nPerm > 0){

    if(nPerm<999){stop("nPerm must be greater than 1000")}
    output$Perm <- mst.perm(formula,data,trt,intervention,nPerm,random,cluster)
    output$Conditional$Perm <- round(data.frame(output$Perm$Conditional),2)
    output$Unconditional$Perm <- round(data.frame(output$Perm$Unconditional),2)
  }

  if(nBoot>0){


    tid <- c(1:nrow(fixedDesignMatrix))

    set.seed(1020252)
    bootSamples <- NULL

    for(ii in 1:length(unique(cluster))){
      selID <- tid[cluster==unique(cluster)[ii]]
      if(length(selID)>0){
        selID2<- sapply(c(1:nBoot),function(x)sample(selID,length(selID),replace=TRUE))
        bootSamples <- rbind(bootSamples ,selID2)
      }

    }

    bootResults <- apply(bootSamples ,2,function(bt)rbd.rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,bt=bt))

    bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults,intervention=intervention)
    output$ES <- bootES
    output$Bootstrap <- bootResults
    output$Conditional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Conditional))))
    output$Unconditional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Unconditional))))
    Bnames <- gsub("Estimate1","Within", gsub("Estimate2", "Total",names(output$Conditional$Bootstrap)))
    names(output$Conditional$Bootstrap) <- Bnames
    names(output$Unconditional$Bootstrap) <- Bnames
  }
  output1 <- list()
  output1$Beta   <- output$Beta
  output1$covParm <- output$covParm$Conditional
  output1$ES     <- output$ES$Conditional
  output1$SchEffects     <- output$SchEffects
  if(nPerm > 0){output1$permES    <- output$Conditional$Perm}
  if(nBoot > 0){output1$Bootstrap <- output$Conditional$Bootstrap}

  output1$Unconditional$ES     <- output$ES$Unconditional
  output1$Unconditional$covParm <- output$covParm$Unconditional
  if(nPerm > 0){output1$Unconditional$permES <- output$Unconditional$Perm}
  if(nBoot > 0){output1$Unconditional$Bootstrap <-  output$Unconditional$Bootstrap}


  output1$Method <- "MLM"

  class(output1) <- "eefAnalytics"
  return(output1)
}


#############################################################################
############# SRT Bayesian main functions ################################################


#' Analysis of Simple Randomised Education Trials using Bayesian Linear Regression Model with Vague Priors.
#'
#' \code{srtBayes} performs analysis of educational trials under the assumption of independent errors among pupils using Bayesian framework with Stan.
#' This can also be used with schools as fixed effects.
#'
#' @export
#' @param formula the model to be analysed is of the form y~x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param adaptD As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param ... additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed to the function.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for the variables specified in the model.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{sigma2}: Residual variance.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Unconditional}: A list of unconditional effect sizes, sigma2 and ProbES obtained based on residual variance from the unconditional model (model with only the intercept as a fixed effect).
#'  }
#' @example inst/examples/srtBExample.R
srtBayes <- function(formula,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...)UseMethod("srtBayes")

#' @export
srtBayes.default <- function(formula,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
srtBayes.formula <- function(formula,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){

  if(nsim<2000){stop("nsim must be at least 2000")}

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])
  mf <- model.frame(formula=formula, data=data)
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]

  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  new.formula0 <- as.formula(post~1)
  new.formula  <- as.formula(paste0(LHS.formula,"~",RHS.formula ))
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1])

  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  output <- srtB(formula=new.formula, unc.formula=new.formula0,intervention=intervention, adaptD=adaptD,nsim=nsim,data1=new.data,threshold=threshold,...)
  output$Method <- "MLM"
  class(output) <- "eefAnalytics"
  return(output)
}




###################################################################################################
############# Bayesian Multilevel Analysis of Cluster Randomised Education Trials
##################################################################################################

#' Bayesian analysis of cluster randomised education trials using Vague Priors.
#'
#' \code{crtBayes} performs analysis of cluster randomised education trials using a multilevel model under a Bayesian setting,
#' assuming vague priors.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param adaptD As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param ... additional arguments of \code{\link[rstanarm]{stan_lmer}} to be passed to the function.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for variables specified in the model.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{covParm}: A vector of variance decomposition into between cluster variance (Schools) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
#' \item \code{ProbES}: A matrix of Bayesian Posterior Probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm and ProbES obtained based on between and within cluster variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/crtBExample.R
crtBayes <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...)UseMethod("crtBayes")

#' @export
crtBayes.default <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
crtBayes.formula <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){

  if(nsim<2000){stop("nsim must be at least 2000")}

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk >0){stop("This is not a CRT design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  stp2 <- as.numeric(apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))


  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])

  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1],cluster=cluster)

  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  RHS.cluster <-"+(1|cluster)"
  new.formula <-  as.formula(paste0(LHS.formula,"~",RHS.formula,RHS.cluster))
  new.formula0 <- as.formula(paste0(LHS.formula,"~",RHS.cluster))


  nsim=nsim
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  #if(nsim < 2000){stop("nsim >= 10000 is recommended")}

  output <-  erantBAYES(formula=new.formula, unc.formula=new.formula0,intervention=intervention,data1=new.data,cluster=random, adaptD=adaptD,nsim=nsim,threshold=threshold,...)
  output$Method <- "MLM"
  class(output) <- "eefAnalytics"
  return(output)
}





###################################################################################################
############# Bayesian Multilevel Analysis of Multisite Randomised Education Trials
##################################################################################################

#' Bayesian analysis of Multisite Randomised Education Trials using Vague Priors.
#'
#' \code{mstBayes} performs analysis of multisite randomised education trials using a multilevel model under a Bayesian setting
#' assuming vague priors.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param adaptD As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability that the observed effect size is greater than or equal to the threshold(s).
#' @param ... additional arguments of \code{\link[rstanarm]{stan_lmer}} to be passed to the function.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for variables specified in the model.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{covParm}: A list of variance decomposition into between cluster variance-covariance matrix (schools and school by intervention) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept and intervention slope.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm and ProbES obtained based on between and within cluster variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/mstBExample.R
mstBayes <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...)UseMethod("mstBayes")

#' @export
mstBayes.default <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
mstBayes.formula <- function(formula,random,intervention,adaptD=NULL,nsim=2000,data,threshold=1:10/10,...){
  if(nsim<2000){stop("nsim must be at least 2000")}

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk ==0){stop("This is not a MST design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  #stp2 <- as.numeric(apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))


  tmp3 <- which(colnames(data)==intervention)
  data[,tmp3] <- as.factor(data[,tmp3])
  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  btp  <- which(substring(tmp,1,nchar(intervention))==intervention & nchar(tmp)==(nchar(intervention)+1))
  btp1 <- tmp[btp]
  btp2 <- gsub(intervention,"trt", btp1)
  btp3 <- data.frame(fixedDesignMatrix[,btp1])
  colnames(btp3) <- btp2
  posttest <- model.response(mf)
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1],btp3,cluster=cluster)

  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  RHS.cluster0 <-"+(1|cluster)"
  RHS.cluster <-paste0("+(1+", paste0(btp2, collapse = "+"),"|cluster)")
  new.formula <-  as.formula(paste0(LHS.formula,"~",RHS.formula,RHS.cluster))
  new.formula0 <- as.formula(paste0(LHS.formula,"~",RHS.cluster0))


  nsim=nsim
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  #if(nsim < 2000){stop("nsim >= 10000 is recommended")}

  output <-  erantMstBAYES(formula=new.formula, unc.formula=new.formula0,intervention=intervention,data1=new.data,cluster=random, adaptD=adaptD,nsim=nsim,threshold=threshold,btp=btp,btp2=btp2,...)
  output$Method <- "MLM"
  class(output) <- "eefAnalytics"
  return(output)
}




## - internal SRT functions
srt <- function(posttest,fixedDesignMatrix,intervention,trt){

  freqFit <- lm(posttest~ fixedDesignMatrix-1)
  cit <- confint(freqFit)
  citt <- rowSums(is.na(cit))
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[which(citt==0),1],cit[which(citt==0),]))
  row.names(betaB)<- colnames(fixedDesignMatrix)[which(citt==0)]
  colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB

  btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention&
                 nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))
  btp2 <- which(substring(colnames(fixedDesignMatrix),1,nchar(intervention))==intervention &
                  nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

  tmpTRT <- table(trt)
  sigma1 <- c(Conditional=summary(freqFit)$sigma,
              Unconditional=summary(lm(posttest~1))$sigma)

  output2.0 <- matrix(NA,length(btp),3 )
  colnames(output2.0)<- c("Estimate","95% LB","95% UB")
  row.names(output2.0) <- row.names(betaB)[btp]

  output2 <- list()
  for( j in names(sigma1)){
    sd.pool <- sigma1[j]
    for( i in 1:length(btp)){
      beta <- betaB[btp[i],1]
      cd <- (beta/sd.pool)
      trt2 <- unique(fixedDesignMatrix[,btp2[i]])
      n.c <- tmpTRT[names(tmpTRT)==trt2[2]]
      n.t <- tmpTRT[names(tmpTRT)==trt2[1]]
      var.cd <- ((n.t+n.c)/(n.t*n.c)+cd^2/(2*(n.t+n.c)))
      se.cd <- sqrt(var.cd)
      cd.lb <- (cd - 1.96*se.cd)
      cd.ub <- (cd + 1.96*se.cd)
      j.df <- (1 - (3/(4*(n.t+n.c-2)-1)))
      g <- (j.df*cd)
      var.g <- (j.df^2 * var.cd)
      se.g <- sqrt(var.g)
      g.lb <- (g - 1.96*se.g)
      g.ub <- (g + 1.96*se.g)

      output2.0[i,] <- c(g, g.lb, g.ub)
    }
    output2[[j]] <- output2.0
  }

  sigma2=round(sigma1^2,2)
  output <- list(Beta=round(betaB,2),ES=output2,sigma2=sigma2)
  return(output)

}


## - internal
srt.srt<- function(posttest,fixedDesignMatrix,intervention,bt,variances=variances){

  posttest2 <- posttest[bt]
  fixedDesignMatrix2 <- fixedDesignMatrix[bt,]



  freqFit <- try(lm(posttest2~ fixedDesignMatrix2-1),silent=TRUE)
  if(attr(freqFit,"class")!="try-error"){

    betaB <- data.frame(summary(freqFit)$coefficients[,1])
    ntpp <- as.character(sapply(as.character(rownames(betaB)),function(x)substring(x,(nchar("fixedDesignMatrix2")+1),nchar(x))))
    row.names(betaB)<- ntpp
    betaB <- betaB

    sigma1 <- c(Conditional=summary(freqFit)$sigma,
                Unconditional=summary(lm(posttest~1))$sigma)

    btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention &
                   nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

    output1<- sapply(names(sigma1), function(x) sapply(1:length(btp ), function (i) betaB[btp[i],1]/sigma1[x] ))
  }

  boot.names <- paste0( rep(c("cond", "uncond"),each=length(btp)),"_" ,rep(row.names(betaB)[btp],length(btp)))
  output1<- matrix(output1,1,(length(btp)*2))
  colnames(output1)<- boot.names
  return(output1)
}


#srt.srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,bt=bt, conditional=T)










## random intercept model - internal
crt <- function(posttest,fixedDesignMatrix,intervention,cluster){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+ (1|cluster))

  np<- row.names(summary(freqFit)$coef)
  cit <- confint(freqFit,np)
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
  row.names(betaB)<- colnames(fixedDesignMatrix)
  colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B<- as.numeric(summary(freqFit)$varcor)
  var.W<- summary(freqFit)$sigma^2
  var.tt <- var.W+var.B
  ICC1 <- var.B/var.tt
  sigmaBE1 <- round(c(var.B,var.W,var.B+var.W,(var.B/(var.B+var.W))),2)
  names(sigmaBE1)<- c("Schools","Pupils","Total","ICC")
  var.B1<- as.numeric(summary(lmer(posttest~ 1+ (1|cluster)))$varcor)
  var.W1<- summary(lmer(posttest~ 1+(1|cluster)))$sigma^2
  var.tt1 <- var.W1+var.B1
  ICC <-(c(Conditional=ICC1,Unconditional=var.B1/var.tt1))
  sigmaBE2 <- round(c(var.B1,var.W1,var.B1+var.W1,(var.B1/(var.B1+var.W1))),2)
  names(sigmaBE2)<- c("Schools","Pupils","Total","ICC")
  sigmaBE <- data.frame(rbind(Conditional=sigmaBE1,
                              Unconditional=sigmaBE2))
  sigmaBE <- sigmaBE
  sigma.W<-c(Conditional=var.W,Unconditional=var.W1)
  sigma.tt<-c(Conditional=var.tt,Unconditional=var.tt1)
  schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
  names(schRand)<- c("Schools","Estimate")
  btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention &
                 nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

  output2 <- list()
  for( j in names(sigma.tt)){
    var.w <-sigma.W[j]
    var.tt<-sigma.tt[j]
    icc<-ICC[j]
    output.1<-list()
    for( i in 1:length(btp )){
      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      esWithin <- g.within(var.w=var.w, beta=beta, icc=icc, group=group, schoolID=cluster)
      esTotal <- g.total(var.tt=var.tt, beta=beta, icc=icc, group=group, schoolID=cluster)

      output1 <- data.frame(rbind(esWithin,esTotal))
      colnames(output1) <- c("Estimate","95% LB","95% UB")
      rownames(output1) <- c("Within","Total")
      output.1[[i]] <- round(output1,2)
    }
    names(output.1) <- row.names(betaB)[btp]
    output2[[j]] <- output.1
  }
  output <- list(Beta=round(betaB,2),covParm=sigmaBE,ES=output2,SchEffects=round(schRand,2))

  return(output)

}


## internal
crtP <- function(posttest,fixedDesignMatrix,intervention,cluster){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+ (1|cluster))

  betaB <- data.frame(summary(freqFit)$coefficients[,1])
  row.names(betaB)<- colnames(fixedDesignMatrix)
  #colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B<- as.numeric(summary(freqFit)$varcor)
  var.W<- summary(freqFit)$sigma^2
  var.tt <- var.W+var.B
  ICC1 <- var.B/var.tt
  sigmaBE1 <- c(var.B,var.W)

  freqFit1 <- lmer(posttest~ 1+ (1|cluster))
  var.B1<- as.numeric(summary(freqFit1)$varcor)
  var.W1<- summary(freqFit1)$sigma^2
  var.tt1 <- var.W1+var.B1
  ICC <-(c(Conditional=ICC1,Unconditional=var.B1/var.tt1))
  sigmaBE2 <- c(var.B1,var.W1)
  sigmaBE <- data.frame(rbind(Conditional=sigmaBE1,Unconditional=sigmaBE2))
  sigmaBE <- sigmaBE
  sigma.W<-c(Conditional=var.W,Unconditional=var.W1)
  sigma.tt<-c(Conditional=var.tt,Unconditional=var.tt1)
  btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention &
                 nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

  output2 <- list()
  for( j in names(sigma.tt)){
    var.w <-sigma.W[j]
    var.tt<-sigma.tt[j]
    icc<-ICC[j]
    output.1<-list()
    for( i in 1:length(btp )){
      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      esWithin <- g.within(var.w=var.w, beta=beta, icc=icc, group=group, schoolID=cluster)
      esTotal <- g.total(var.tt=var.tt, beta=beta, icc=icc, group=group, schoolID=cluster)
      output1 <- data.frame(rbind(esWithin,esTotal))
      colnames(output1) <- c("Estimate","95% LB","95% UB")
      rownames(output1) <- c("Within","Total")
      output.1[[i]] <- round(output1,2)
    }
    names(output.1) <- row.names(betaB)[btp]
    output2[[j]] <- output.1
  }

  output <- list(ES=output2)

  return(output)

}


## - internal
crt.perm <- function(formula,data,stp,stp2,intervention,cluster,nPerm,random){

  data2 <- data[,-which(colnames(data)==intervention)]

  g <- matrix(NA,nPerm,2*(length(unique(stp2))-1))
  g.unc <- matrix(NA,nPerm,2*(length(unique(stp2))-1))

  for(i in 1:nPerm){
    set.seed(12890*i+1)
    tp3 <- data.frame(stp,sample(stp2))
    names(tp3) <- c(paste(random),paste(intervention))
    data.tp4 <- merge(data2,tp3,by=random)
    data.tp4 <- data.tp4[order(data.tp4[,which(colnames(data.tp4)==random)]),]
    cluster = data.tp4[,which(colnames(data.tp4)==random)]
    tmp34 <- which(colnames(data.tp4)==intervention)
    data.tp4[,tmp34] <- as.factor(data.tp4[,tmp34])
    mf <- model.frame(formula=formula, data=data.tp4)
    fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data.tp4)))
    tmp <- colnames(fixedDesignMatrix )
    tmp[1]  <- "Intercept"
    colnames(fixedDesignMatrix)<- tmp
    posttest <- model.response(mf)
    intervention <- intervention

    p2CRTFREQ <-crtP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)
    chkppp <- data.frame(cond=unlist(p2CRTFREQ$ES$Conditional),uncond=unlist(p2CRTFREQ$ES$Unconditional))
    chkppp2 <- c(seq(1,6*(length(unique(tp3[,2]))-1),6),seq(2,6*(length(unique(tp3[,2]))-1),6))
    chkppp3 <- chkppp2[order(chkppp2)]
    g[i,]  <-  chkppp[chkppp3,"cond"]
    g.unc[i,]  <-  chkppp[chkppp3,"uncond"]

  }
  ntpp <- rep(names(p2CRTFREQ$ES$Conditional),2)
  ntpp <- ntpp[order(ntpp )]
  wt <- rep(c("Within","Total"),length(names(p2CRTFREQ$ES$Conditional)))
  colnames(g) <- paste(ntpp ,wt,sep="")
  colnames(g.unc) <- paste(ntpp ,wt,sep="")
  g1 <- list(condtional=g,unconditional=g.unc)
  return(g1)
}


## - internal
crt.crt<- function(posttest,fixedDesignMatrix,intervention,cluster,bt){

  posttest2 <- posttest[bt]
  fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
  cluster2 <- cluster[bt]


  freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1|cluster2)),silent=TRUE)
  output2 <- NULL
  if(attr(freqFit,"class")!="try-error"){

    betaB <- data.frame(summary(freqFit)$coefficients[,1])
    row.names(betaB)<- colnames(fixedDesignMatrix)
    betaB <- betaB
    var.B2<- as.matrix(summary(freqFit)$varcor)
    var.B3 <- c(matrix(attr(var.B2[[1]],"stddev")))
    var.B <- var.B3^2
    var.W<- summary(freqFit)$sigma^2
    var.tt <- var.W+var.B
    ICC1 <- var.B/var.tt
    sigmaBE1 <- c(var.B,var.W)
    names(sigmaBE1)<- c("Schools","Pupils")
    var.B21<- as.matrix(summary(lmer(posttest2~ 1+ (1|cluster2)))$varcor)
    var.B31 <- c(matrix(attr(var.B21[[1]],"stddev")))
    var.B1 <- var.B31^2
    var.W1<- summary(lmer(posttest2~ 1+(1|cluster2)))$sigma^2
    var.tt1 <- var.W1+var.B1
    ICC <-(c(Conditional=ICC1, Unconditional=var.B1/var.tt1))
    sigmaBE2 <- c(var.B1,var.W1)
    names(sigmaBE2)<- c("Schools","Pupils")
    sigmaBE <- data.frame(rbind(Conditional=sigmaBE1,
                                Unconditional=sigmaBE2))
    sigmaBE <- sigmaBE
    sigma.W<-c(Conditional=var.W,Unconditional=var.W1)
    sigma.tt<-c(Conditional=var.tt,Unconditional=var.tt1)
    btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention&
                   nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

    output2 <- list()
    for( j in names(sigma.tt)){
      var.w <-sigma.W[j]
      var.tt<-sigma.tt[j]
      output.1<-list()
      for( i in 1:length(btp )){

        beta <- betaB[btp[i],1]
        group <- fixedDesignMatrix[,btp[i]]

        esWithin <- beta/sqrt(var.w)
        esTotal <- beta/sqrt(var.tt)

        output1 <- data.frame(rbind(esWithin,esTotal))
        names(output1) <- c("Estimate")
        rownames(output1) <- c("Within","Total")
        output.1[[i]] <- round(output1,2)
      }
      names(output.1) <- row.names(betaB)[btp]
      output2[[j]] <- output.1
    }

  }
  return(output2)
}


## - internal

g.within <- function(var.w, beta, icc, group, schoolID){
  t <- group; id <- schoolID
  d.w <- (beta/sqrt(var.w))
  n.it <- table(id[t==1]); n.ic <- table(id[t==0])
  m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
  M <- (m.t + m.c)
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  n.sim.1 <- ((N.c * sum(n.it^2))/(as.numeric(N.t)*as.numeric(N)))
  n.sim.2 <- ((N.t * sum(n.ic^2))/(as.numeric(N.c)*as.numeric(N)))
  n.sim <- (n.sim.1 + n.sim.2)
  vterm1 <- ((N.t+N.c)/(as.numeric(N.t)*as.numeric(N.c)))
  vterm2 <- (((1+(n.sim-1)*icc))/(1-icc))
  vterm3 <- ((d.w^2)/(2*(N-M)))
  se <- sqrt(vterm1*vterm2+vterm3)
  LB <- (d.w-1.96*se); UB <- (d.w+1.96*se)
  output <- data.frame(d.w, LB, UB)
  names(output) <- c("g", "LB", "UB")
  return(output)
}

## - internal

g.total <- function(var.tt, beta, icc, group, schoolID){
  t <- group; id <- schoolID
  n.it <- table(id[t==1]); n.ic <- table(id[t==0])
  m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
  M <- (m.t + m.c)
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  n.ut <- ((N.t^2-sum(n.it^2))/(as.numeric(N.t)*as.numeric(m.t-1)))
  n.uc <- ((N.c^2-sum(n.ic^2))/(as.numeric(N.c)*as.numeric(m.c-1)))
  dt.1 <- (beta/sqrt(var.tt))
  dt.2 <- sqrt(1-icc*(((N-n.ut*m.t-n.uc*m.c)+n.ut+n.uc-2)/(N-2)))
  d.t <- (dt.1*dt.2)

  n.sim.1 <- ((as.numeric(N.c) * sum(n.it^2))/(as.numeric(N.t)*as.numeric(N)))
  n.sim.2 <- ((as.numeric(N.t) * sum(n.ic^2))/(as.numeric(N.c)*as.numeric(N)))
  n.sim <- (n.sim.1 + n.sim.2)
  B <- (n.ut*(m.t-1)+n.uc*(m.c-1))
  A.t <- ((as.numeric(N.t)^2*sum(n.it^2)+(sum(n.it^2))^2-2*as.numeric(N.t)*sum(n.it^3))/as.numeric(N.t)^2)
  A.c <- ((as.numeric(N.c)^2*sum(n.ic^2)+(sum(n.ic^2))^2-2*as.numeric(N.c)*sum(n.ic^3))/as.numeric(N.c)^2)
  A <- (A.t + A.c)

  vterm1 <- (((N.t+N.c)/(as.numeric(N.t)*as.numeric(N.c)))*(1+(n.sim-1)*icc))
  vterm2 <- (((N-2)*(1-icc)^2+A*icc^2+2*B*icc*(1-icc))*d.t^2)
  vterm3 <- (2*(N-2)*((N-2)-icc*(N-2-B)))
  se <- sqrt(vterm1+vterm2/vterm3)
  LB <- (d.t-1.96*se); UB <- (d.t+1.96*se)
  output <- data.frame(d.t, LB, UB)
  names(output)<- c("g", "LB", "UB")
  return(output)
}


## compile bootstrap results - internal
bootCompile <- function(output,trt,bootResults,intervention){
  Con_withinBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
  Con_totalBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
  Un_withinBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
  Un_totalBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
  for(k in 1:length(bootResults)){
    tmp <- bootResults[[k]]$Conditional
    tmpR <- NULL
    for(j in 1:length(tmp)){

      tmpR  <- c(tmpR,tmp[[j]][,1])
    }
    Con_withinBoot[k,] <- tmpR[seq(1,2*length(tmp),2)]
    Con_totalBoot[k,] <- tmpR[seq(2,2*length(tmp),2)]
  }
  for(k in 1:length(bootResults)){
    tmp1 <- bootResults[[k]]$Unconditional
    tmpR1 <- NULL
    for(j in 1:length(tmp1)){

      tmpR1  <- c(tmpR1,tmp1[[j]][,1])
    }
    Un_withinBoot[k,] <- tmpR1[seq(1,2*length(tmp1),2)]
    Un_totalBoot[k,] <- tmpR1[seq(2,2*length(tmp1),2)]
  }

  withinCI_Conditional <- apply(Con_withinBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
  TotalCI_Conditional <- apply(Con_totalBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
  withinCI_Unconditional <- apply(Un_withinBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
  TotalCI_Unconditional <- apply(Un_totalBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
  withinCI<-list(withinCI_Conditional,withinCI_Unconditional)
  TotalCI<-list(TotalCI_Conditional,TotalCI_Unconditional)
  ES<-list(output$ES$Conditional,output$ES$Unconditional)
  names(ES)<-c("Conditional","Unconditional")
  names(withinCI)<-c("Conditional","Unconditional")
  names(TotalCI)<-c("Conditional","Unconditional")
  btp <- which(substring(row.names(output$Beta),1,nchar(intervention))==intervention)
  tmpES <- list()
  for( j in names(ES)){
    withinCI1<-withinCI[[j]]
    TotalCI1<-TotalCI[[j]]
    ES1<-ES[[j]]
    temp<-list()
    for(kk in 1:length(ES1)){
      tmp1 <- round(rbind(withinCI1[,kk], TotalCI1[,kk]),2)
      tmp2 <- cbind(ES1[[kk]][,1],tmp1)
      colnames(tmp2)<- c("Estimate","95% LB","95% UB")
      row.names(tmp2)	<- c("Within","Total")
      temp[[kk]] <- tmp2
    }
    names(temp) <- names(ES1)
    tmpES[[j]] <- temp
    #names(tmpES) <- names(output$ES)

  }
  return(tmpES)
}



## - internal

rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))
  np<- row.names(summary(freqFit)$coef)
  cit <- confint(freqFit,np)
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
  row.names(betaB)<- colnames(fixedDesignMatrix)
  colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B2<- as.matrix(VarCorr(freqFit)$cluster)
  var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
  var.B2_2<-as.data.frame(VarCorr(freqFit))
  var.B3 <- diag(var.B2)[-1]
  var.sch <- var.B2[1,1]
  vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
  var.E <- var.B3
  var.W<- summary(freqFit)$sigma^2
  N <- as.numeric(length(summary(freqFit)$res))
  N.t <- as.matrix(summary(trt))[-1];
  var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
  ICC1 <- sum(var.B2)/var.tt
  sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
  names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
  freqFit1<- lmer(posttest~1+(1|cluster))
  var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
  var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
  var.sch1 <- var.B21[1,1]
  var.W1<- summary(freqFit1)$sigma^2
  var.tt1 <- var.sch1+var.W1
  ICC2 <- sum(var.B21)/var.tt1
  sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
  names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
  ICC <-(c(conditional=ICC1,unconditional=ICC2))
  sigmaBE <- list(sigmaBE1,sigmaBE2)
  names(sigmaBE)<-c('Conditional','Unconditional')
  sigmaBE <- sigmaBE
  schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
  names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
  btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)

  output2 <- list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    var.e<- var.E[i]
    esWithin <- g.within.mst(var.w=var.W, var.e=var.e, beta=beta, group=group, schoolID=cluster)
    esTotal <- g.total.mst(var.w=var.W, var.e=var.e, var.tt=var.tt, beta=beta, group=group, schoolID=cluster)
    outputc <- data.frame(rbind(esWithin,esTotal))
    colnames(outputc) <- c("Estimate","95% LB","95% UB")
    rownames(outputc) <- c("Within","Total")
    output2[[i]] <- round(outputc,2)
  }
  names(output2) <- row.names(betaB)[btp]
  output3<-list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    esWithin1 <- g.within(var.w=var.W1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    esTotal1 <- g.total(var.tt=var.tt1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    outputu <- data.frame(rbind(esWithin1,esTotal1))
    colnames(outputu) <- c("Estimate","95% LB","95% UB")
    rownames(outputu) <- c("Within","Total")
    output3[[i]] <- round(outputu,2)
  }
  names(output3) <- row.names(betaB)[btp]
  output2.3<-list(Conditional=output2,Unconditional=output3)
  output <- list(Beta=round(betaB,2),covParm=sigmaBE,ES=output2.3,SchEffects=round(schRand,2))

  return(output)
}


## - internal

rbdP <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))
  np<- row.names(summary(freqFit)$coef)
  cit <- confint(freqFit,np)
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
  row.names(betaB)<- colnames(fixedDesignMatrix)
  #colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B2<- as.matrix(VarCorr(freqFit)$cluster)
  var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
  var.B2_2<-as.data.frame(VarCorr(freqFit))
  var.B3 <- diag(var.B2)[-1]
  var.sch <- var.B2[1,1]
  vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
  var.E <- var.B3
  var.W<- summary(freqFit)$sigma^2
  N <- as.numeric(length(summary(freqFit)$res))
  N.t <- as.matrix(summary(trt))[-1];
  var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
  ICC1 <- sum(var.B2)/var.tt
  sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
  names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
  freqFit1<- lmer(posttest~1+(1|cluster))
  var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
  var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
  var.sch1 <- var.B21[1,1]
  var.W1<- summary(freqFit1)$sigma^2
  var.tt1 <- var.sch1+var.W1
  ICC2 <- sum(var.B21)/var.tt1
  sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
  names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
  ICC <-(c(conditional=ICC1,unconditional=ICC2))
  sigmaBE <- list(sigmaBE1,sigmaBE2)
  names(sigmaBE)<-c('Conditional','Unconditional')
  sigmaBE <- sigmaBE
  schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
  names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
  btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)

  output2 <- list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    var.e<- var.E[i]
    esWithin <- g.within.mst(var.w=var.W, var.e=var.e, beta=beta, group=group, schoolID=cluster)
    esTotal <- g.total.mst(var.w=var.W, var.e=var.e, var.tt=var.tt, beta=beta, group=group, schoolID=cluster)
    outputc <- round(data.frame(rbind(esWithin,esTotal)),2)
    colnames(outputc) <- c("Estimate","95% LB","95% UB")
    rownames(outputc) <- c("Within","Total")
    output2[[i]] <- round(outputc,2)
  }
  names(output2) <- row.names(betaB)[btp]
  output3<-list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    esWithin1 <- g.within(var.w=var.W1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    esTotal1 <- g.total(var.tt=var.tt1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    outputu <- data.frame(rbind(esWithin1,esTotal1))
    colnames(outputu) <- c("Estimate","95% LB","95% UB")
    rownames(outputu) <- c("Within","Total")
    output3[[i]] <- round(outputu,2)
  }
  names(output3) <- row.names(betaB)[btp]
  output2.3<-list(Conditional=output2,Unconditional=output3)
  output <- list(ES=output2.3)

  return(output)
}

## - internal

rbd.rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,bt){

  posttest2 <- posttest[bt]
  fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
  trt2 <- trt[bt]
  cluster2 <- cluster[bt]


  freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1+trt2|cluster2)),silent=TRUE)
  output2 <- NULL
  if(attr(freqFit,"class")!="try-error"){

    betaB <- data.frame(summary(freqFit)$coefficients[,1])
    row.names(betaB)<- colnames(fixedDesignMatrix)
    betaB <- betaB
    var.B2<- as.matrix(VarCorr(freqFit)$cluster)
    var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
    var.B2_2<-as.data.frame(VarCorr(freqFit))
    var.B3 <- diag(var.B2)[-1]
    var.sch <- var.B2[1,1]
    vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
    var.E <- var.B3
    var.W<- summary(freqFit)$sigma^2
    N <- as.numeric(length(summary(freqFit)$res))
    N.t <- as.matrix(summary(trt))[-1];
    var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
    ICC1 <- sum(var.B2)/var.tt
    sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
    names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
    freqFit1<- lmer(posttest~1+(1|cluster))
    var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
    var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
    var.sch1 <- var.B21[1,1]
    var.W1<- summary(freqFit1)$sigma^2
    var.tt1 <- var.sch1+var.W1
    ICC2 <- sum(var.B21)/var.tt1
    sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
    names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
    ICC <-(c(conditional=ICC1,unconditional=ICC2))
    sigmaBE <- list(sigmaBE1,sigmaBE2)
    names(sigmaBE)<-c('Conditional','Unconditional')
    sigmaBE <- sigmaBE
    schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
    names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
    btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)


    output3 <- list()
    for( i in 1:length(btp)){

      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      var.e<- var.E[i]
      esWithin <- beta/sqrt(var.W)
      esTotal <- beta/sqrt(var.tt)

      outputc <- data.frame(rbind(esWithin,esTotal))
      names(outputc) <- c("Estimate")
      rownames(outputc) <- c("Within","Total")
      output3[[i]] <- round(outputc,2)
    }
    names(output3) <- row.names(betaB)[btp]
    output4<-list()
    for( i in 1:length(btp)){

      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      esWithin1 <- beta/sqrt(var.W1)
      esTotal1 <- beta/sqrt(var.tt1)

      outputu <- data.frame(rbind(esWithin1,esTotal1))
      names(outputu) <- c("Estimate")
      rownames(outputu) <- c("Within","Total")
      output4[[i]] <- round(outputu,2)
    }
    names(output4) <- row.names(betaB)[btp]
    output2<-list(Conditional=output3,Unconditional=output4)

  }
  return(output2)
}



## - internal

g.within.mst <- function(var.w, var.e, beta, group, schoolID){
  t <- group; id <- schoolID
  d.w <- (beta/sqrt(var.w))
  n.itc <- as.data.frame.matrix(table(id,t))
  n.it <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "1"]
  n.ic <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "0"]
  M <- length(unique(id))
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  nin<-(n.it*n.ic)/(n.it+n.ic)
  ni <- n.it+n.ic
  vterm1 <- 1/(sum(var.w/(var.e+var.w/nin)))
  vterm2 <- ((d.w^2)/((2*N-4*M)))
  se <- sqrt(vterm1+vterm2)
  LB <- (d.w-1.96*se); UB <- (d.w+1.96*se)
  output <- data.frame(d.w, LB, UB)
  names(output) <- c("g", "LB", "UB")
  return(output)
}

## - internal

g.total.mst <- function(var.w, var.e, var.tt, beta, group, schoolID){
  t <- group; id <- schoolID
  n.itc <- as.data.frame.matrix(table(id,t))
  n.it <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "1"]
  n.ic <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "0"]
  M <- length(unique(id))
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  nin<-(n.it*n.ic)/(n.it+n.ic)
  ni <- n.it+n.ic
  d.t <- (beta/sqrt(var.tt))
  vterm1 <- 1/(sum(var.tt/(var.e+var.w/nin)))
  vterm2 <- ((d.t^2)/((2*N-4*M)))
  se <- sqrt(vterm1+vterm2)
  LB <- (d.t-1.96*se); UB <- (d.t+1.96*se)
  output <- data.frame(d.t, LB, UB)
  names(output)<- c("g", "LB", "UB")
  return(output)
}


## - internal

mst.perm <- function(formula,data,trt,intervention,nPerm,random,cluster){

  data2 <- data
  g <- matrix(NA,nPerm,2*(length(unique(trt))-1))
  g.unc <- matrix(NA,nPerm,2*(length(unique(trt))-1))

  for(i in 1:nPerm){
    set.seed(12890*i+1)
    tryCatch({data2[,which(colnames(data)==intervention)]<-unlist(tapply(trt,cluster,function(x)sample(x)))
    data3 <- data2[order(data2[,which(colnames(data2)==random)],data2[,which(colnames(data2)==intervention)]),]
    cluster = data3[,which(colnames(data3)==random)]
    mf <- model.frame(formula=formula, data=data3)
    fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data3)))
    tmp <- colnames(fixedDesignMatrix )
    tmp[1]  <- "Intercept"
    colnames(fixedDesignMatrix)<- tmp
    posttest <- model.response(mf)
    intervention <- intervention
    trt2 <- data3[,which(colnames(data3)==intervention)]
    p2CRTFREQ <-rbdP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt2,cluster=cluster)

    chkppp <- data.frame(cond=unlist(p2CRTFREQ$ES$Conditional),uncond=unlist(p2CRTFREQ$ES$Unconditional))
    chkppp2 <- c(seq(1,6*(length(unique(trt))-1),6),seq(2,6*(length(unique(trt))-1),6))
    chkppp3 <- chkppp2[order(chkppp2)]
    g[i,]  <-  chkppp[chkppp3,"cond"]
    g.unc[i,]  <-  chkppp[chkppp3,"uncond"]}, error=function(e){})
  }
  ntpp <- rep(names(p2CRTFREQ$ES$Conditional),2)
  ntpp <- ntpp[order(ntpp)]
  wt <- rep(c("Within","Total"),length(names(p2CRTFREQ$ES$Conditional)))
  colnames(g) <- paste(ntpp ,wt,sep="")
  colnames(g.unc) <- paste(ntpp ,wt,sep="")
  g1 <- list(Conditional=g,Unconditional=g.unc)
  return(g1)
}



## - internal
srtB <- function(formula,unc.formula,intervention,tmp, adaptD=adaptD,nsim=nsim,data1,threshold,...){
  freqFit <- stan_glm(formula=formula, data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,... )
  freqFit0 <-stan_glm(formula=unc.formula,data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,... )

  Sim_Beta <- data.frame(as.matrix(freqFit))
  btp <- which(substring(names(Sim_Beta),1,nchar(intervention))==intervention &
                 nchar(names(Sim_Beta))==(nchar(intervention)+1))

  output1 <- NULL
  output2 <- NULL
  ProbES  <- list()
  for( i in 1:length(btp )){

    Output0=Bsrt.Summary(freqFit=freqFit,freqFit0=freqFit0,Sim_Beta= Sim_Beta, btp= btp[i])
    Beta <-  Output0$Beta
    Sim_Beta_t <-  Output0$Sim_Beta_t
    sigma2.1 <-  Output0$sigma2.1
    sigma2.2 <-  Output0$sigma2.2
    Sim_sigma2.1 <-  Output0$Sim_sigma2.1
    Sim_sigma2.2 <-  Output0$Sim_sigma2.2
    sim_ES1 <-  Output0$sim_ES1
    sim_ES2 <-  Output0$sim_ES2

    ProbES0 <- data.frame(cbind(Cond=sapply(threshold, function(j) mean(sim_ES1 > j)),
                                Uncond=sapply(threshold, function(j) mean(sim_ES2 > j))))

    row.names(ProbES0) <- paste0("P(ES>",threshold,")")
    ProbES[[i]] <- round(ProbES0,3)


    gtmp1 <- round(c("ES"=mean(sim_ES1), quantile(sim_ES1,probs=c(0.025,0.975))),2)
    gtmp2 <- round(c("ES"=mean(sim_ES2), quantile(sim_ES2,probs=c(0.025,0.975))),2)
    output1 <- rbind(output1,gtmp1)
    output2 <- rbind(output2,gtmp2)
  }

  row.names(output1) <- names(Sim_Beta)[btp]
  row.names(output2) <- names(Sim_Beta)[btp]
  output1<- data.frame(output1)
  output2<- data.frame(output2)
  colnames(output1)<- c("Estimate","95% LB","95% UB")
  colnames(output2)<- c("Estimate","95% LB","95% UB")
  output <- list(Beta=round(Beta,2),ES=round(output1,2),sigma2=sigma2.1, ProbES=lapply(ProbES,function(x) x[, "Cond"]),
                 unconditional= list(ES=round(output2,2),sigma2=sigma2.2,ProbES=lapply(ProbES,function(x) x[, "Uncond"])))
  return(output)
}

## summarise beta parameters - internal
Bsrt.Summary <- function(freqFit,freqFit0, Sim_Beta, btp){

  Sim_Beta_t <- as.matrix(freqFit,pars=names(Sim_Beta[btp])) #treatment effect
  Beta <- data.frame( summary(freqFit,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
  names(Beta)<- c("Estimate","95% LB","95% UB")

  Sim_sigma2.1 <- as.matrix(freqFit,pars="sigma")  #sigma(pupil): sqrt of residual
  Sim_sigma2.2 <- as.matrix(freqFit0,pars="sigma")
  sigma2.1 <- mean(Sim_sigma2.1^2)
  sigma2.2 <- mean(Sim_sigma2.2^2)
  sim_ES1 <- Sim_Beta_t/Sim_sigma2.1
  sim_ES2 <- Sim_Beta_t/Sim_sigma2.2


  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t, sigma2.1=sigma2.1, Sim_sigma2.1=Sim_sigma2.1, sigma2.2=sigma2.2, Sim_sigma2.2=Sim_sigma2.2, sim_ES1=sim_ES1,sim_ES2=sim_ES2)
  return(ModelSmry)
}





## - internal

erantBAYES<- function(formula,unc.formula,intervention,data1,cluster,  adaptD=adaptD,nsim=nsim,threshold,...){

  freqFit <- stan_lmer(formula=formula, data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,...)
  freqFit0 <-stan_lmer(formula=unc.formula,data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,...)

  Sim_Beta <- data.frame(as.matrix(freqFit))
  btp <- which(substring(names(Sim_Beta),1,nchar(intervention))==intervention &
                 nchar(names(Sim_Beta))==(nchar(intervention)+1))


  output1 <- list()
  output2 <- list()
  ProbES <- list()
  for( i in 1:length(btp )){

    Output0=errantSummary(freqFit=freqFit,freqFit0=freqFit0,Sim_Beta= Sim_Beta, btp= btp[i])
    Beta <-  Output0$Beta
    ICC1 <- Output0$ICC1
    Sim_Beta_t <-  Output0$Sim_Beta_t
    var.B1 <- Output0$var.B1
    var.W1 <- Output0$var.W1
    var.tt1 <- Output0$var.tt1

    sim_ES.tt1 <- Sim_Beta_t/sqrt(var.tt1)
    sim_ES.W1 <- Sim_Beta_t/sqrt(var.W1)

    ICC2 <- Output0$ICC2
    var.W2 <- Output0$var.W2
    var.tt2 <- Output0$var.tt2
    sim_ES.tt2 <- Sim_Beta_t/sqrt(var.tt2)
    sim_ES.W2 <- Sim_Beta_t/sqrt(var.W2)


    ProbES0 <- data.frame(cbind(Total1=sapply(threshold, function(j) mean(sim_ES.tt1 > j)),
                                Within1 =sapply(threshold, function(j) mean(sim_ES.W1 > j)),
                                Total2=sapply(threshold, function(j) mean(sim_ES.tt2 > j)),
                                Within2 =sapply(threshold, function(j) mean(sim_ES.W2 > j))))

    row.names(ProbES0) <- paste0("P(ES>",threshold,")")
    ProbES[[i]] <- round(ProbES0,3)


    gtmp.W1 <- round(c("ES"=mean(sim_ES.W1), quantile(sim_ES.W1,probs=c(0.025,0.975))),2)
    gtmp.W2 <- round(c("ES"=mean(sim_ES.W2), quantile(sim_ES.W2,probs=c(0.025,0.975))),2)
    gtmp.tt1 <- round(c("ES"=mean(sim_ES.tt1), quantile(sim_ES.tt1,probs=c(0.025,0.975))),2)
    gtmp.tt2 <- round(c("ES"=mean(sim_ES.tt2), quantile(sim_ES.tt2,probs=c(0.025,0.975))),2)
    output.1 <- rbind(Within=gtmp.W1,Total=gtmp.tt1)
    output.2 <- rbind(Within=gtmp.W2,Total=gtmp.tt2)
    colnames(output.1) <- c("Estimate","95% LB","95% UB")
    colnames(output.2) <- c("Estimate","95% LB","95% UB")
    output1[[i]] <- round(output.1,2)
    output2[[i]] <- round(output.2,2)
  }

  names(output1) <- names(Sim_Beta)[btp]
  names(output2) <- names(Sim_Beta)[btp]
  Output <- list(Beta=round(Beta,2), ES=data.frame(output1), covParm=Output0$sigmaBE1, SchEffects=Output0$SchEffects, ProbES=lapply(ProbES, function(x) x[, c("Within1","Total1")]),
                 unconditional= list(ES=data.frame(output2), covParm=Output0$sigmaBE2, ProbES=lapply(ProbES, function(x) x[, c("Within2","Total2")])))
  return(Output)
}



## - internal

errantSummary <- function(freqFit,freqFit0,Sim_Beta, btp){

  Sim_Beta_t <- as.matrix(freqFit,pars=names(Sim_Beta[btp])) #treatment effect
  Beta <- data.frame( summary(freqFit,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
  names(Beta)<- c("Estimate","95% LB","95% UB")

  #conditional ouptut
  Sim_var.W1 <- as.matrix(freqFit,pars="sigma")^2  #sigma(pupil)
  Sim_var.B1 <- as.matrix(freqFit,pars="Sigma[cluster:(Intercept),(Intercept)]")
  Sim_var.tt1 <- Sim_var.B1+Sim_var.W1
  Sim_ICC1 <- Sim_var.B1/Sim_var.tt1
  sigmaBE1 <- c(mean(Sim_var.B1),mean(Sim_var.W1),mean(Sim_var.tt1),mean(Sim_ICC1))
  names(sigmaBE1)<- c("Schools","Pupils","Total","ICC")

  #conditional ouptut
  Sim_var.W2 <- as.matrix(freqFit0,pars="sigma")^2  #sigma(pupil)
  Sim_var.B2 <- as.matrix(freqFit0,pars="Sigma[cluster:(Intercept),(Intercept)]")
  Sim_var.tt2 <- Sim_var.B2+Sim_var.W2
  Sim_ICC2 <- Sim_var.B2/Sim_var.tt2
  sigmaBE2 <- c(mean(Sim_var.B2),mean(Sim_var.W2),mean(Sim_var.tt2),mean(Sim_ICC2))
  names(sigmaBE2)<- c("Schools","Pupils","Total","ICC")

  #random effects
  SchEffects0 <- summary(freqFit, regex_pars = "b\\[\\(Intercept\\) cluster\\:", probs = c(0.025, 0.975),digits = 2)
  SchEffects1 <- data.frame(SchEffects0)[,"mean"]
  schl_ID <-  as.numeric(as.character(substr(rownames(SchEffects0), nchar("b[(Intercept) cluster\\:"), nchar(rownames(SchEffects0))-1)))
  SchEffects <- round(data.frame(cbind(Schools=schl_ID,Intercept=SchEffects1)),2)



  #all output
  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t,
                    var.B1=Sim_var.B1, var.W1=Sim_var.W1, var.tt1=Sim_var.tt1, ICC1=Sim_ICC1, sigmaBE1=round(sigmaBE1,2),
                    var.B2=Sim_var.B2, var.W2=Sim_var.W2, var.tt2=Sim_var.tt2, ICC2=Sim_ICC2, sigmaBE2=round(sigmaBE2,2),
                    SchEffects =SchEffects)
  return(ModelSmry)
}





## - internal

erantMstBAYES<- function(formula,unc.formula,intervention,data1,cluster,  adaptD=adaptD,nsim=nsim,threshold,btp,btp2,...){

  freqFit <- stan_lmer(formula=formula, data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,...)
  freqFit0 <-stan_lmer(formula=unc.formula,data=data1,adapt_delta =adaptD, seed = 1234,iter=nsim,...)

  Sim_Beta <- data.frame(as.matrix(freqFit))

  output1 <- list()
  output2 <- list()
  ProbES <- list()
  for( i in 1:length(btp )){

    Output0=errantMStSummary(freqFits=freqFit,Cond=F,Sim_Beta= Sim_Beta, btp= btp[i],btp1=btp,btp2=btp2 )
    Beta <-  Output0$Beta
    Sim_Beta_t <-  Output0$Sim_Beta_t
    var.B1 <- Output0$var.B1
    var.W1 <- Output0$var.W1
    var.tt1 <- Output0$var.tt1
    sim_ES.tt1 <- Sim_Beta_t/sqrt(var.tt1)
    sim_ES.W1 <- Sim_Beta_t/sqrt(var.W1)

    Output01=errantMStSummary(freqFits=freqFit0,Cond=T,Sim_Beta= Sim_Beta, btp= btp[i],btp1=btp,btp2=btp2 )

    var.W2 <- Output01$var.W1
    var.tt2 <- Output01$var.tt1
    sim_ES.tt2 <- Sim_Beta_t/sqrt(var.tt2)
    sim_ES.W2 <- Sim_Beta_t/sqrt(var.W2)


    ProbES0 <- data.frame(cbind(Total1=sapply(threshold, function(j) mean(sim_ES.tt1 > j)),
                                Within1 =sapply(threshold, function(j) mean(sim_ES.W1 > j)),
                                Total2=sapply(threshold, function(j) mean(sim_ES.tt2 > j)),
                                Within2 =sapply(threshold, function(j) mean(sim_ES.W2 > j))))

    row.names(ProbES0) <- paste0("P(ES>",threshold,")")
    ProbES[[i]] <- round(ProbES0,3)

    gtmp.W1 <- round(c("ES"=mean(sim_ES.W1), quantile(sim_ES.W1,probs=c(0.025,0.975))),2)
    gtmp.W2 <- round(c("ES"=mean(sim_ES.W2), quantile(sim_ES.W2,probs=c(0.025,0.975))),2)
    gtmp.tt1 <- round(c("ES"=mean(sim_ES.tt1), quantile(sim_ES.tt1,probs=c(0.025,0.975))),2)
    gtmp.tt2 <- round(c("ES"=mean(sim_ES.tt2), quantile(sim_ES.tt2,probs=c(0.025,0.975))),2)
    output.1 <- rbind(Within=gtmp.W1,Total=gtmp.tt1)
    output.2 <- rbind(Within=gtmp.W2,Total=gtmp.tt2)
    colnames(output.1) <- c("Estimate","95% LB","95% UB")
    colnames(output.2) <- c("Estimate","95% LB","95% UB")
    output1[[i]] <- round(output.1,2)
    output2[[i]] <- round(output.2,2)
  }

  names(ProbES) <- names(Sim_Beta)[btp]
  names(output1) <- names(Sim_Beta)[btp]
  names(output2) <- names(Sim_Beta)[btp]
  Output <- list(Beta=round(Beta,2), ES=data.frame(output1), covParm=Output0$sigmaBE1, SchEffects=Output0$SchEffects, ProbES=lapply(ProbES, function(x) x[, c("Within1","Total1")]),
                 unconditional= list(ES=data.frame(output2), covParm=Output01$sigmaBE1, ProbES=lapply(ProbES, function(x) x[, c("Within2","Total2")])))
  return(Output)
}



## - internal


errantMStSummary <- function(freqFits, Cond=F,Sim_Beta, btp, btp1,btp2){

  data1 <- freqFits$data
  N <- length(freqFits$residuals)
  N.t <- sapply(btp1, function(i) sum(data1[,i])) # need to check how to add this

  if(Cond==T){
    # Unconditional results
    Beta <- NULL
    Sim_Beta_t <- NULL
    SchEffects <- NULL
    Sim_var.B1 <- as.matrix(freqFits,pars="Sigma[cluster:(Intercept),(Intercept)]")
    Sim_var.W1 <- as.matrix(freqFits,pars="sigma")^2  #sigma(pupil)
    Sim_var.tt1 <- Sim_var.B1+Sim_var.W1
    Sim_ICC1 <- Sim_var.B1/Sim_var.tt1
    sigmaBE1 <- round(c(mean(Sim_var.B1),mean(Sim_var.W1),mean(Sim_var.tt1),mean(Sim_ICC1)),2)
    names(sigmaBE1)<- c("Schools","Pupils","Total","ICC")
  }




  if(Cond==F){
    #conditional ouptut
    Sim_Beta_t <- as.matrix(freqFits,pars=names(Sim_Beta[btp])) #treatment effect
    Beta <- data.frame( summary(freqFits,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
    names(Beta)<- c("Estimate","95% LB","95% UB")

    Sim_var.B   <- as.matrix(freqFits, regex_pars ="Sigma\\[\\cluster\\:trt")

    name.varW <- "sigma"
    name.Vsch <- "Sigma[cluster:(Intercept),(Intercept)]"
    name.diag <- grep(paste(paste0(btp2,",",btp2),collapse = "|"), colnames(Sim_var.B), value = T)
    name.diag0<- grep(paste(paste0(btp2,",",btp2),collapse = "|"), colnames(Sim_var.B), value = T, invert = T)
    name.cov  <- grep(paste(paste0(btp2,",","\\(Intercept\\)"),collapse = "|"), colnames(Sim_var.B), value = T)


    Sim_var.W1 <- as.matrix(freqFits, pars=name.varW)^2  #sigma(pupil)
    Sim_var.sch1 <- as.matrix(freqFits, pars=name.Vsch)
    Sim_var.B1  <- Sim_var.B[, name.diag]
    Sim_cov.B1  <- Sim_var.B[, name.cov]
    Sim_Term1.1 <-  Sim_var.W1+Sim_var.sch1
    Sim_Term2.1 <- 2*Sim_cov.B1+Sim_var.B1
    sim_Nt_B1  <- rowSums(matrix(Sim_Term2.1, ncol = length(btp1))%*%matrix(N.t))
    Sim_var.tt1 <- (N*Sim_Term1.1 + sim_Nt_B1)/N
    Sim_ICC1 <- rowSums(matrix(Sim_var.B1, ncol = length(btp1)))/Sim_var.tt1


    ## variance covariance matrix
    varAll <- data.frame(summary(freqFits, regex_pars = "Sigma\\[cluster\\:", probs = 0.025,digits = 2))
    diag <- varAll[c(name.Vsch,name.diag),"mean"]
    offdiag <- varAll[name.diag0,"mean"]
    Vcov1 <- matrix(NA, ncol = length(diag), nrow = length(diag))
    rownames(Vcov1) <- c("Intercept",btp2)
    colnames(Vcov1) <- c("Intercept",btp2)
    Vcov1[lower.tri(Vcov1)] <- offdiag
    Vcov1[upper.tri(Vcov1)] <- t(Vcov1)[upper.tri(t(Vcov1))]
    diag(Vcov1) <- diag
    sigmaBE1 <- list(round(Vcov1,2),round(mean(Sim_var.W1),2),round(mean(Sim_var.tt1),2),round(mean(Sim_ICC1),2))
    names(sigmaBE1)<- c("School:trt","Pupils","Total","ICC")

    #random effects
    SchEffects00 <- data.frame(summary(freqFits, regex_pars = "b\\[\\(Intercept\\) cluster\\:", probs = c(0.025, 0.975),digits = 2))
    SchEffects0 <- SchEffects00[,"mean"]
    SchEffects1 <- sapply(btp2,function(x) data.frame(summary(freqFits, regex_pars = paste0("b\\[",x, " cluster\\:"),digits = 2))[,"mean"])

    schl_ID <-  as.numeric(as.character(substr(rownames(SchEffects00), nchar("b[(Intercept) cluster\\:"), nchar(rownames(SchEffects00))-1)))
    SchEffects <- round(data.frame(cbind(Schools=schl_ID,Intercept=SchEffects0,SchEffects1)),2)
  }


  #all output
  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t,
                    var.B1=Sim_var.B1, var.W1=Sim_var.W1, var.tt1=Sim_var.tt1, ICC1=Sim_ICC1, sigmaBE1=sigmaBE1,
                    SchEffects =SchEffects)
  return(ModelSmry)
}
