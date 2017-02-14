#' A plot method for the an eefAnalytics S3 object obtained  from the eefAnalytic package.
#' 
#' plots different figures based on output from eefAnalytics package. 
#' 
#' @export
#' @param x an output object from the eefAnalytics package.
#' @param group a value indicating which intervention to plot.   
#' This must not be greater than the number of intervention excluding the control group.
#' For a two arm trial, the maximum value is 1 and a maximum value of 2 for three arm trial.
#' @param ... arguments passed to \code{\link[graphics]{plot.default}}
#' @details Plot produces graphical visualisation depending on which model is fitted:
#' \itemize{
#' \item For \code{srtFREQ()}, plot can only be used when \code{nBoot} or \code{nPerm} is specified to visualise the distribution of boostrapped or permutated values. 
#' \item For \code{crtFREQ()} or \code{mstFREQ}, plot shows the distribution of random intercepts when \code{group=NULL}. 
#' It produces histogram of permutated or bootstrapped values when \code{group} is specified and either \code{nBoot} or \code{nPerm} is also specified. 
#' \item For \code{mlmBayes()}, plot produces the distrbution of random intercepts when \code{group = NULL}. 
#' It produces the probability of effect size to be greater than a pre-specified threshold when group is specified. 
#' \item Lastly, plot produces forest plots to compare CACE estimated for different level of compliance when \code{caceSRTBoot()} or 
#' \code{caceCRTBoot()} or \code{caceMSTBoot()} is used.
#' } 
#' @return Returns relevant plots for each model.
#' @example inst/examples/plotExample.R
plot.eefAnalytics <- function(x,group=NULL,...){
    plotObject(analyticObject=x,group=group, compare=FALSE,modelNames=FALSE,...)
  
}



#' A plot function to compare diferent eefAnalytics S3 objects from the eefAnalytics package.
#' 
#' @description A forest plot comparing the different eefAnalytics results.
#' 
#' @export
#' @param eefAnalyticsList A list of eefAnalytics S3 objects from eefAnalytics package.
#' @param group a value indicating which intervention to plot.   
#' This must not be greater than the number of intervention excluding the control group.
#' For a two arm trial, the maximum value is 1 and a maximum value of 2 for three arm trial.
#' @param modelNames a string factor containing the names of model to compare
#' @details \code{ComparePlot} produces a forest plot which compares the effect size and the associated confidence interval from the different model. 
#' For a multilevel model, it shows effect size based on residual variance and total variance.
#' 
#' @return Returns a forest plot to compare the different models
#' @example inst/examples/compareExample.R
ComparePlot <- function(eefAnalyticsList,group=NULL,modelNames=NULL){
  if(class(eefAnalyticsList)!="list"){stop("eefAnalyticsList is not a list.")}
  
  if(!all(unlist(lapply(eefAnalyticsList,FUN=class))=="eefAnalytics")){stop("Not all list objects are a eefAnalytics class object.")}
  
  if(is.null(modelNames)){stop("modelNames must be specified.")}
  if(is.null(group)){stop("group must be specified.")}
  
  plotObject(analyticObject=eefAnalyticsList,group=group, compare=TRUE,modelNames=modelNames)
    
}





# Internal plot function

plotObject <- function(analyticObject,group=NULL, compare=FALSE,modelNames=FALSE,...){
  
  if(compare==TRUE & !is.null(names(analyticObject))){stop("Specify the list of objects to compare")}
  
  if(sum(analyticObject$Method=="LM")==1 & sum(names(analyticObject)=="CACE")!=1){
    if(sum(names(analyticObject)=="Bootstrap"|names(analyticObject)=="Perm")==0){stop("Only relevant for bootstrapped or permutated values")}
    if(sum(names(analyticObject)=="Bootstrap")==1){
      ntp <- nrow(as.matrix(analyticObject$ES))
      
      if(is.null(group)){stop("Group number must be defined")}
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      obs.est <- analyticObject$ES[group,1]
      tmp2 <- as.numeric(analyticObject$Bootstrap[,group])
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab="Bootstrap estimates",main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
      
    }
    
    if(sum(names(analyticObject)=="Perm")==1){
      ntp <- nrow(as.matrix(analyticObject$ES))
      if(is.null(group)){stop("Group number must be defined")}
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      obs.est <- analyticObject$ES[group,1]
      tmp2 <- as.numeric(analyticObject$Perm[,group])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab="Permutation values",main=paste("P(X|NULL)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
      
    }
    
    
  }
  
  if(sum(names(analyticObject)=="CACE")==1){
    tmp <- analyticObject$CACE
    par_original <- par()[c("cex","font")]
    
    op <- par(cex=1,font=1)
    xtp <- max(range(c(tmp[,-1])))-1.2
    xub <- ifelse(max(range(c(tmp[,-1]))) > 1.5, 2+xtp ,2)
    xlb <- ifelse(max(range(c(tmp[,-1]))) > 1.5, -1-xtp,-1)
    tlb <- ifelse(max(range(c(tmp[,-1]))) > 1.5, -1-xtp ,-1)
    tub <- ifelse(max(range(c(tmp[,-1]))) > 1.5, 1.25+xtp ,1.5)
    forest(x=tmp[,2],ci.lb=tmp[,3],ci.ub=tmp[,4],xlab="Hedge's g", ylim=c(0,(length(tmp[,3])+3)),
           slab=tmp[,1],xlim=c(xlb,xub),alim = c(-0.5,1.5))
    par(font=2)
    text(tlb, (length(tmp[,3])+1.5),"Compliance", pos=4)
    text(tub, (length(tmp[,3])+1.5),"95% CI", pos=4)
    
    par(par_original)
  }
  
  if(sum(analyticObject$Method=="MLM")==1 & sum(names(analyticObject)=="CACE")!=1){
    if(is.null(group)){
      tmp <- data.frame(y=analyticObject$SchEffects$Estimate,x=c(1:length(analyticObject$SchEffects$Estimate)))
      tmp2 <- tmp[order(tmp$y),]
      barplot(tmp2$y,ylab="Random intercepts for schools",names.arg=tmp2$x,las=2,col="cornflowerblue",border="cornflowerblue",...)
      
    }
    
    if( !is.null(group) & sum(names(analyticObject)=="Bootstrap")>0){
      ntp <- length(analyticObject$ES)
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      obs.est <- analyticObject$ES[[group]][2,1]
      tmp <- as.numeric(analyticObject$Bootstrap[2,])
      dim(tmp) <- c(ntp,floor(length(tmp)/ntp))
      tmp <- t(tmp)
      tmp2 <- tmp[,group]
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab="Bootstrap estimates",main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }
    
    
    if( !is.null(group) & sum(names(analyticObject)=="Perm")>0){
      ntp <- ifelse(is.list(analyticObject$ES),length(analyticObject$ES),1)
      if(group > ntp){stop("Group number must be less than the number of intervention")}
      obs.est <- analyticObject$ES[[group]][2,1]
      tmp2 <- as.numeric(analyticObject$Perm[,group])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab="Permutation values",main=paste("P(X|NULL)=",pvalue,sep=""))
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
        if(!length(tmp0)>9){tmp2 <- as.matrix(tmp$ES[-2,])}				
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