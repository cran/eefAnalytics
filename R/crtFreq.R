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
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals.
#' @param nPerm number of permutations required to generate a permutated p-value.
#' @param seed seed required for bootstrapping and permutation procedure, if not provided default seed will be used.
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
crtFREQ<- function(formula,random,intervention,baseln,nPerm,nBoot,seed,data)UseMethod("crtFREQ")
#' @export
crtFREQ.default<- function(formula,random,intervention,baseln,nPerm,nBoot,seed,data){stop("No correct formula input given.")}


#' @export
crtFREQ.formula <- function(formula,random,intervention,baseln,nPerm,nBoot,seed,data){

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  #trt<-as.factor(trt)
  if( missing(baseln)){trt <- as.factor(trt)}
  if(!missing(baseln)){trt <- relevel(as.factor(trt),baseln)}
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk >0){stop("This is not a CRT design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  control<-colnames(table(cluster2,trt))[1]
  stp2 <- (apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))

  if(!missing(nPerm) & !missing(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(missing(nPerm)){nPerm <-0}
  if(missing(nBoot)){nBoot <-0}
  tmp3 <- which(colnames(data)==intervention)
  if( missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),baseln)}
  #data[,tmp3] <- as.factor(data[,tmp3])

  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)

  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))

  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


  output <- crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster, btp=btp)


  if(nPerm>0){
    if(nPerm<999){stop("nPerm must be greater than 1000")}
    output$Perm<- crt.perm(formula,data,stp,stp2,intervention,cluster,nPerm,random, btp,seed,baseln)
    output$Condtional$Perm <- round(data.frame(output$Perm$condtional),2)
    output$Unconditional$Perm <- round(data.frame(output$Perm$unconditional),2)
  }


  if(nBoot >0){

    #if(nBoot<1000){stop("nBoot must be greater than 1000")}

    tid <- c(1:nrow(fixedDesignMatrix))

    #set.seed(1020252)
    if(!missing(seed)){set.seed(seed)}
    bootSamples <- NULL

    for(ii in 1:length(unique(cluster))){
      selID <- tid[cluster==unique(cluster)[ii]]
      if(length(selID)>0){
        selID2<- sapply(c(1:nBoot),function(x) selID [sample(1:length(selID), length(selID),replace=TRUE)])
        bootSamples <- rbind(bootSamples ,selID2)
      }

    }


    bootResults <- apply(bootSamples ,2,function(bt)crt.crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster,bt=bt, btp=btp))
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
  output1$Function <- "srtFREQ"
  class(output1) <- "eefAnalytics"
  return(output1)
}



#######################################################################################

## random intercept model - internal
crt <- function(posttest,fixedDesignMatrix,intervention,cluster,btp){

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
  #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention &
  #               nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

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
crtP <- function(posttest,fixedDesignMatrix,intervention,cluster,btp){

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
  #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention &
  #               nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

  output2 <- list()
  for( j in names(sigma.tt)){
    var.w <-sigma.W[j]
    var.tt<-sigma.tt[j]
    icc<-ICC[j]
    output.1<-list()
    for( i in 1:length(btp)){
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
crt.perm <- function(formula,data,stp,stp2,intervention,cluster,nPerm,random,btp,seed,baseln){

  data2 <- data[,-which(colnames(data)==intervention)]

  g <- matrix(NA,nPerm,2*(length(unique(stp2))-1))
  g.unc <- matrix(NA,nPerm,2*(length(unique(stp2))-1))

  for(i in 1:nPerm){
    #set.seed(12890*i+1)
    if(!missing(seed)){set.seed(seed*i+1)}
    tp3 <- data.frame(stp,sample(stp2))
    names(tp3) <- c(paste(random),paste(intervention))
    data.tp4 <- merge(data2,tp3,by=random)
    data.tp4 <- data.tp4[order(data.tp4[,which(colnames(data.tp4)==random)]),]
    cluster = data.tp4[,which(colnames(data.tp4)==random)]
    tmp34 <- which(colnames(data.tp4)==intervention)
    #data.tp4[,tmp34] <- as.factor(data.tp4[,tmp34])
    if( missing(baseln)){data.tp4[,tmp34] <- as.factor(data.tp4[,tmp34])}
    if(!missing(baseln)){data.tp4[,tmp34] <- relevel(as.factor(data.tp4[,tmp34]),baseln)}
    mf <- model.frame(formula=formula, data=data.tp4)
    fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data.tp4)))
    tmp <- colnames(fixedDesignMatrix )
    tmp[1]  <- "Intercept"
    colnames(fixedDesignMatrix)<- tmp
    posttest <- model.response(mf)
    intervention <- intervention

    p2CRTFREQ <-crtP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster, btp=btp)
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
crt.crt<- function(posttest,fixedDesignMatrix,intervention,cluster,bt,btp){

  posttest2 <- posttest[bt]
  fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
  cluster2 <- cluster[bt]


  freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1|cluster2)),silent=TRUE)
  output2 <- NULL
  if(!is(freqFit, "try-error")){

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
    #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention&
    #               nchar(colnames(fixedDesignMatrix))==(nchar(intervention)+1))

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
bootCompile <- function(output,trt,bootResults,intervention,btp){
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
  #btp <- which(substring(row.names(output$Beta),1,nchar(intervention))==intervention)
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


