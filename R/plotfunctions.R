#' A plot method for an eefAnalytics S3 object obtained from the eefAnalytics package.
#'
#' Plots different figures based on output from eefAnalytics package.
#'
#' @export
#' @param x an output object from the eefAnalytics package.
#' @param group a scalar value indicating which intervention to plot.
#' This must not be greater than the number of intervention groups excluding the control group.
#' For a two arm trial, the maximum value is 1 and a maximum value of 2 for three arm trial.
#' @param Conditional a logical value to indicate whether to plot the conditional effect size.
#' The default is Conditional=TRUE, otherwise Conditional=FALSE should be specified for plot based on the unconditional effect size.
#'  Conditional variance is total or residual variance from a multilevel model with fixed effects, whilst unconditional variance is total variance or residual variance from a multilevel model with only intercept as fixed effect.
#' @param ES_Total A logical value indicating whether to plot the effect size based on total variance or within school variance.
#' The default is ES_Total=TRUE, to plot the effect size using total variance.
#' ES_Total=FALSE should be specified for the effect size based on within school or residuals variance.
#' @param slop A logical value indicating whether to return the plot of random intercept (default is slop=FALSE).
#' return other school-by-intervention interaction random slope (s) is slop=TRUE.
#' This argument is suitable only for mstBayes and mstFREQ functions.
#' @param ... arguments passed to \code{\link[graphics]{plot.default}}
#' @details Plot produces a graphical visualisation depending on which model is fitted:
#' \itemize{
#' \item For \code{srtFREQ()}, plot can only be used when \code{nBoot} or \code{nPerm} is specified to visualise the distribution of bootstrapped or permutated values.
#' \item For \code{crtFREQ()} or \code{mstFREQ()}, plot shows the distribution of random intercepts when \code{group=NULL}.
#' It produces histogram of permutated or bootstrapped values when \code{group} is specified and either \code{nBoot} or \code{nPerm} is also specified.
#' }
#' @return Returns relevant plots for each model.
#' @example inst/examples/plotExample.R
plot.eefAnalytics <- function(x,group=NULL, Conditional=TRUE,ES_Total=TRUE,slop=FALSE,...){
    plotObject(analyticObject=x,group=group, Conditional=Conditional,ES_Total=ES_Total, compare=FALSE,modelNames=FALSE,...)

}



#' A plot function to compare different eefAnalytics S3 objects from the eefAnalytics package.
#'
#' @description It generates bar plot that compares the effect size from eefAnalytics' methods.
#'
#' @export
#' @param eefAnalyticsList A list of eefAnalytics S3 objects from eefAnalytics package.
#' @param group a scalar value indicating which intervention to plot.
#' This must not be greater than the number of intervention groups excluding the control group.
#' For a two arm trial, the maximum value is 1 and a maximum value of 2 for three arm trial.
#' @param Conditional  a logical value to indicate whether to plot conditional effect size.
#' The default is Conditional=TRUE, otherwise Conditional=FALSE should be specified for plot based on unconditional effect size.
#'  Conditional variance is total or residual variance a multilevel model with fixed effects, whilst unconditional variance is total variance or residual variance from a multilevel model with only intercept as fixed effect.
#' @param ES_Total A logical value indicating whether to plot the effect size based on total variance or within school variance.
#' The default is ES_Total=TRUE, to plot effect size using total variance.
#' ES_Total=FALSE should be specified for effect size based on within school or residuals variance.
#' @param modelNames a string factor containing the names of model to compare. See examples below.
#' @details \code{ComparePlot} produces a bar plot which compares the effect sizes and the associated confidence intervals from the different models.
#' For a multilevel model, it shows the effect size based on residual variance and total variance.
#'
#' @return Returns a bar plot to compare the different methods.
#' @example inst/examples/compareExample.R
ComparePlot <- function(eefAnalyticsList,group=NULL, Conditional=TRUE,ES_Total=TRUE,modelNames=NULL){
  if(class(eefAnalyticsList)!="list"){stop("eefAnalyticsList is not a list.")}

  if(!all(unlist(lapply(eefAnalyticsList,FUN=class))=="eefAnalytics")){stop("Not all list objects are a eefAnalytics class object.")}

  if(is.null(modelNames)){stop("modelNames must be specified.")}
  if(is.null(group)){stop("group must be specified.")}

  plotObject(analyticObject=eefAnalyticsList,group=group, Conditional=Conditional,ES_Total=ES_Total, compare=TRUE,modelNames=modelNames)

}





# Internal plot function

plotObject <- function(analyticObject,group=NULL, Conditional=TRUE,ES_Total=TRUE,slop=FALSE, compare=FALSE,modelNames=FALSE,...){

  if(Conditional ==TRUE){analyticObject2=analyticObject; Condname="Conditional"}
  if(Conditional ==FALSE){analyticObject2=analyticObject$Unconditional; Condname="Unconditional"}
  if(ES_Total ==TRUE){ES_TW<-"Total"}
  if(ES_Total ==FALSE){ES_TW<-"Within"}

  if(compare==TRUE & !is.null(names(analyticObject))){stop("Specify the list of objects to compare")}

  if(sum(analyticObject$Method=="LM")==1){
    if(sum(names(analyticObject)=="Bootstrap"|names(analyticObject)=="permES")==0){stop("Only relevant for bootstrapped or permutated values")}
    if(sum(names(analyticObject)=="Bootstrap")==1){
      ntp <- nrow(as.matrix(analyticObject$ES))

      if(is.null(group)){stop("Group number must be defined")}
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      obs.est <- analyticObject2$ES[group,1]
      tmp2 <- as.numeric(analyticObject2$Bootstrap[,group])
      xlabs=paste0(Condname," Bootstrap estimates")
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }

    if(sum(names(analyticObject)=="permES")==1){
      ntp <- nrow(as.matrix(analyticObject$ES))
      if(is.null(group)){stop("Group number must be defined")}
      if(group > ntp){stop("Group number must be less than the number of intervention")}


      Perm.names <-names(analyticObject2$permES)

      obs.est <- analyticObject2$ES[group,1]
      tmp2 <- as.numeric(analyticObject2$permES[,group])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Permutation values (PermES) based on ",Condname, " ES")

      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
    }


  }


  if(sum(analyticObject$Method=="MLM")==1){

    if(is.null(group)){

      tmp000  <- data.frame(analyticObject$SchEffects)
      if(slop==FALSE){
        tmp00 <- tmp000[,grep("Schools|Intercept|Estimate",names(tmp000))]
        mar11 <-  c(5, 4, 4, 2) + 0.1
      }
      if(slop==TRUE & dim(tmp000)[2] ==2){stop("x must be mstFREQ or mstBAyes object")}
      if(slop==TRUE & dim(tmp000)[2] >2){
        tmp00 <- tmp000[,!(names(tmp000) %in% "Intercept")]
        if(dim(tmp00)[2]==2){mar11 <- c(5, 4, 4, 2) + 0.1}
        if(dim(tmp00)[2]==3){mar11 <- c(5, 2, 4, 0) + 1.0}
        if(dim(tmp00)[2] >3){mar11 <- c(3, 2, 0, 0) + 1.0}}

      op <- par(mfrow = c(floor(dim(tmp00)[2]/2),round(dim(tmp00)[2]/2)),
                mar = mar11)


      for(i in 2:dim(tmp00)[2]){
        tmp <- data.frame(y=tmp00[,i],x=c(1:length(tmp00[,i])))
        tmp2 <- tmp[order(tmp$y),]
        ylabs=gsub("trt", "Intervention ",gsub("Estimate","Intercept", names(tmp00)[i]))
        barplot(tmp2$y,names.arg=tmp2$x,las=2,col="cornflowerblue",border="cornflowerblue")
        if(dim(tmp00)[2]<=2){mtext(ylabs, side = 2.5, line = 2)}
        if(dim(tmp00)[2] >2){mtext(ylabs, side = 2, line = 1.7, cex = 0.8)}
      }
      lines1=-2.5
      if(dim(tmp00)[2] >3){lines1=-1}
      title(xlab="School labels", outer = TRUE, line = lines1,cex.lab = 1.2)
      par(op)
    }

    if( !is.null(group) & sum(names(analyticObject)=="Bootstrap")>0){
      ntp <- length(analyticObject$ES)
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      Boot.names <-names(analyticObject2$Bootstrap)
      obs.est <- analyticObject2$ES[[group]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$Bootstrap[,grep(ES_TW, Boot.names, ignore.case = T)[group]])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Bootstrap estimates for ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }


    if( !is.null(group) & sum(names(analyticObject)=="permES")>0){
      ntp <- ifelse(is.list(analyticObject$ES),length(analyticObject$ES),1)
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      Perm.names <-names(analyticObject2$permES)
      obs.est <- analyticObject2$ES[[group]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$permES[,grep(ES_TW, Perm.names, ignore.case = T)[group]])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Permutation values(PermES) based on ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
    }


    if( !is.null(group) &sum(names(analyticObject)=="ProbES")> 0 ){
      ntp <- length(analyticObject$ES)
      nmax <- group*3
      nmin <- nmax -2
      tmp <- analyticObject$ProbES[,-1]
      nntmp <- ncol(tmp)/3
      if(group > nntmp){stop("Group number must be less than the number of intervention")}
      es <- analyticObject$ProbES[,1]
      tmp2 <- tmp[,c(nmin:nmax)]

      par_original <- par()[c("mar","xpd")]

      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(es,tmp2[,1],ylim=c(0,max(tmp2)),ylab="Probability",cex.lab=1,cex.axis=1,type="n", xlab=expression("Effect size" >= "x"),cex=1,...)
      lines(es,tmp2[,1],col="chartreuse3",cex=1.5,lwd=1.5,lty=2)
      lines(es,tmp2[,2],col="violetred",cex=1.5,lwd=1.5,lty=3)
      lines(es,tmp2[,3],col="cornflowerblue",cex=1.5,lwd=1.5,lty=1)
      points(es,tmp2[,1],col="chartreuse3",cex=1.5,lwd=1.5,pch=7)
      points(es,tmp2[,2],col="violetred",cex=1.5,lwd=1.5,pch=1)
      points(es,tmp2[,3],col="cornflowerblue",cex=1.5,lwd=1.5,pch=12)
      legend("topright",inset=c(-.2,0) ,legend=c("Within ","Between ","Total "),bty="n",lty=c(2,3,1), cex=.8,
             pch=c(7,1,12),col=c("chartreuse3","violetred","cornflowerblue"),title="Variance")

      par(par_original)
    }

  }

  if(is.null(names(analyticObject))){

    ltp <- names(analyticObject)
    if(!is.null(ltp)){stop("Specify list of eefAnalytics objects for comparison")}
    ntp <- length(analyticObject)
    if(length(modelNames)!= ntp){stop("Names must be equal to the number of eefAnalytics objects")}
    es.mean <- es.lower <- es.upper <- p.name <- var.name <- NULL
    for(k in 1:ntp){
      tmp <- analyticObject[[k]]

      if(tmp$Method=="LM"){

        tmp2 <- as.matrix(tmp$ES)
        if(is.null(group)){stop("Group number must be defined")}
        if(group > nrow(tmp2 )){stop("Group number must be less than the number of intervention")}
        es.mean1 <-tmp2[group,1]
        es.lower1 <-tmp2[group,2]
        es.upper1 <-tmp2[group,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep("Within",length(es.mean1))

      }

      if(tmp$Method=="MLM"){
        nlist <- 1
        if(is.null(group)){stop("Group number must be defined")}
        if(group > length(tmp$ES)){stop("Group number must be less than the number of intervention")}
        tmp2 <- unlist(tmp$ES[[group]])
        estmp <-c(sapply(c(1:nlist),function(x)c(6*(x-1)+1,6*(x-1)+2)))
        lbtmp <-c(sapply(c(1:nlist),function(x)c(6*(x-1)+3,6*(x-1)+4)))
        ubtmp <-c(sapply(c(1:nlist),function(x)c(6*(x-1)+5,6*(x-1)+6)))
        es.mean1 <- tmp2[estmp]
        es.lower1 <-tmp2[lbtmp]
        es.upper1 <-tmp2[ubtmp]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep(c("Within","Total"),nlist)

      }


      if(tmp$Method=="MLM"&sum(names(tmp)=="ProbES")==1){
        if(is.null(group)){stop("Group number must be defined")}
        if(group > length(tmp$ES)){stop("Group number must be less than the number of intervention")}
        tmp0 <- unlist(tmp$ES)
        if(length(tmp0)>9){tmp2 <- as.matrix(tmp$ES[[group]][-2,])}
        if(!length(tmp0)>9){tmp2 <- as.matrix(data.frame(tmp$ES)[-2,])}
        es.mean1 <-tmp2[,1]
        es.lower1 <-tmp2[,2]
        es.upper1 <-tmp2[,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep(c("Within","Total"),1)

      }

      es.mean <- c(es.mean,es.mean1)
      es.lower <- c(es.lower,es.lower1)
      es.upper <- c(es.upper,es.upper1)
      p.name <- c(p.name,p.name1)
      var.name <- c(var.name,var.name1)

    }


    xtp <- max(range(c(es.mean,es.lower,es.upper)))-1.2
    xub <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, 2+xtp ,2)
    xlb <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, -2+xtp,-2)
    tlb <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, -2+xtp ,-2)
    vlb <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, -1.25+xtp ,-1)
    tub <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, 1.25+xtp ,1.5)
    alb <- ifelse(max(range(c(es.mean,es.lower,es.upper))) > 1.3, vlb ,-1)


    alb <- ifelse(min(range(c(es.mean,es.lower,es.upper)))<=(alb+0.1),alb-(min(range(c(es.mean,es.lower,es.upper)))- alb ),alb)

    par_original <- par()[c("cex","font")]

    op <- par(cex=1,font=1)
    forest(x=es.mean,ci.lb=es.lower,ci.ub=es.upper,xlab="Hedge's g", ylim=c(0,(length(es.upper)+3)),
           slab=p.name,xlim=c(xlb,xub),alim = c(-1,1.5),ilab=var.name,ilab.xpos=alb,
           pch=as.numeric(as.factor(p.name)))



    par(font=2)
    text(xlb, (length(es.upper)+1.5),"Model", pos=4)
    text((alb-0.25), (length(es.upper)+1.5),"Variance", pos=4)
    text(tub, (length(es.upper)+1.5),"95% CI", pos=4)

    par(par_original)
  }

}
