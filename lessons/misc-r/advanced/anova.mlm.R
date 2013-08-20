#######################################################################
anova.mlm=function(object,...,force.int=FALSE) {
  n=num=length(obj <- list(object, ...))
  if(num>1) {
# order objects by increasing rank
    ord=object$rank
    for (i in 2:num) ord=c(ord,obj[[i]]$rank)
    obj=obj[order(ord)]
# check all have the same response variables
    y=model.response(model.frame(obj[[1]]))
    for (i in 2:num) {
      y1=model.response(model.frame(obj[[i]]))
      if (any(dim(y)!=dim(y1)))
        stop("Response Variables different\n")
      if(any(y!=y1))
        stop("Response Variables different\n")
      }
# force the intercept 
    for (i in 1:num) if (attributes(terms(obj[[i]]))$intercept==0) {
      m0=update(obj[[i]],y[,1]~.)
      m1=update(obj[[i]],y[,1]~.+1)
      if ((m1$rank==m0$rank) | force.int)
        obj[[i]]=update(obj[[i]],~.+1)
    }
# find the Base model: smallest common model
    y1=y[,1]
    z1=model.matrix(obj[[1]])
    m1=lm(y1~0+z1)
    z1=z1[,!is.na(coef(m1)),drop=FALSE] 
    for (i in 2:num) {
      z2=model.matrix(obj[[i]])
      m1=lm(y1~0+z2)
      z2=z2[,!is.na(coef(m1)),drop=FALSE] 
      z1=cbind(z2,z1)
      m1=lm(y1~0+z1)
      z1=z1[,is.na(coef(m1)),drop=FALSE]
    }
    if(dim(z1)[2]>0) {
      m1=lm(y~0+z1)
      SS0=deviance(m1)
      DF0=m1$df.residual
      fm="y ~ 1"
    }
    else {
      SS0=colSums(y^2)
      DF0=dim(z1)[1]
      fm="y ~ 0"
    }
# compute sum of squares and degrees of freedom
    SSE=SSR=matrix(NA,length(SS0),num)
    DFE=DFR=numeric(num)
    for (i in 1:num) {
      SSE[,i]=deviance(obj[[i]])
      SSR[,i]=SS0-SSE[,i]
      DFE[i]=obj[[i]]$df.residual
      DFR[i]=DF0-DFE[i]
    }
# check if all models are Base models
    if (all(DFR==0)) num=1
  }  
  if (num>1) {  
# check if any models are the Base model
    if (any(DFR==0)) {
      m1=obj[[1]]
      keep=(DFR>0)
      SSE=SSE[,keep,drop=FALSE]
      SSR=SSR[,keep,drop=FALSE]
      DFE=DFE[keep]
      DFR=DFR[keep]
      n=sum(keep)
      obj=obj[keep]
      fm=deparse(formula(m1))
    }
# if none are Base determine the form of the Base model
    if (n==num) {
      Terms=attributes(terms(obj[[1]]))$term.labels
      for (i in 2:num) Terms=c(Terms,attributes(terms(obj[[i]]))$term.labels)
      Terms=table(Terms)==num
      Terms=names(Terms[Terms])
      if (length(Terms)==0) {
        z0=model.matrix(as.formula(fm))
        if (any(add <- !(colnames(z1)%in%colnames(z0))))
          fm=paste("y ~",paste(colnames(z1)[add],collapse=" + "))
        }
      else {
        fm=paste("y ~",paste(Terms,collapse=" + "))
        z0=model.matrix(as.formula(fm),data=model.frame(obj[[1]])[-1])
        if (any(add <- !(colnames(z1)%in%colnames(z0))))
          fm=paste(fm,paste(colnames(z1)[add],collapse=" + "),sep=" + ")
        }
    }
# compute F statistics, P values and finish it up
    nm1=paste("M",1:n,sep="")
    Fstat=t(t(SSR)*DFE/t(SSE)/DFR)
    Pval=t(pf(t(Fstat),DFR,DFE,lower.tail=FALSE))
    colnames(Fstat)=paste(nm1,"Fstat",sep=":")
    colnames(Pval)=paste(nm1,"Pval",sep=":")
    DF=cbind(Model=DFR,Residuals=DFE)
    rownames(DF)=nm1
    name="\nF Test compared to Base Model\n\n"
    betaparms=cbind(t(apply(Pval,2,est.beta)),mu=colMeans(Pval))
    rownames(betaparms)=nm1
    z=cbind(Fstat,Pval)[,order(c(2*(1:n)-1,2*(1:n)))]
    rownames(z)=names(SS0)
    form=character(n+1)
    form[1] <- fm
    for (i in 1:n) form[i+1]=deparse(formula(obj[[i]]))
    nm1=format(c("Base",paste("Model ",1:n,sep="")),justify="right")
    form=paste(nm1,": ",form,"\n",sep="")
    SS0=list(SSR,SSE)
  }
# If only 1 model supplied.
  if(num==1) {
    y=model.response(model.frame(object))
    x=model.matrix(object)
    ord=attributes(x)$assign
# force the intercept if required
    if (force.int & min(ord)>0) {
      m1=lm(y~x,singular.ok=TRUE)
      if (diff(range(deviance(object)/deviance(m1)))<1e-6)
      x=cbind("(Intercept)"=1,x[,!is.na(coef(m1)[-1,1])])
      ord=c(0,ord[!is.na(coef(m1)[-1,1])])
      }
# compute sum of squares and degrees of freedom for each term
    if(min(ord)==1) {SS0=colSums(y^2)} 
    else {SS0=colSums(scale(y,scale=FALSE)^2)} 
    DF=object$rank+object$df.residual-1+min(ord)
    if ((n0 <- max(ord))>1)
    for (i in 2:n0-1) {
      m1=lm(y~x[,ord<=i])
      SS0=cbind(SS0,deviance(m1))
      DF=c(DF,m1$df.residual)
    }
    SS0=cbind(SS0,SSE<-deviance(object))
    DF=c(DF,dfE<-df.residual(object))
    for (i in 1:n0) {
      SS0[,i]=SS0[,i]-SS0[,i+1]
      DF[i]=DF[i]-DF[i+1]
    }
    nm1=attributes(terms(object))$term.labels[DF[1:n0]>0]
    n0=n0-sum(DF==0)
    SS0=SS0[,DF>0]
    DF=DF[DF>0]
    colnames(SS0)=names(DF)=c(nm1,"Residuals")
# Compute F stats, p values and finish up
    Fstat=t(t(SS0[,1:n0])*dfE/DF[1:n0])/SSE
    Pval=t(pf(t(Fstat),DF[1:n0],dfE,lower.tail=FALSE))
    colnames(Fstat)=paste(nm1,"Fstat",sep=":")
    colnames(Pval)=paste(nm1,"Pval",sep=":")
    name="\nF Test for Term Significance\nTerms added sequentially\n\nModel: "
    form=paste(deparse(formula(object)),"\n")
    betaparms=cbind(t(apply(Pval,2,est.beta)),mu=colMeans(Pval))
    rownames(betaparms)=nm1
    z=cbind(Fstat,Pval)[,order(c(2*(1:n0)-1,2*(1:n0)))]
  }
# return the object
  structure(z,SS=SS0,DF=DF,class=c("anova.mlm","matrix"),formula=form,betaparms=betaparms,heading=name)
}
########################################################################
plot.anova.mlm=function(object) {
# if(dev.cur()==1) X11()
# save the old par values and setup the graphsheet for output
  old.par=par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(oma=c(0,0,5,0),mfcol=c(2,2))
  n=dim(object)[2]%/%2
  parms=attr(object,"betaparms")
  DF=attr(object,"DF")
# for each term or model
  for(i in 1:n) {
# plot histogram of p-values and impose beta density
    z=hist(object[,2*i],plot=FALSE)
    x1=seq(min(object[,2*i]),max(object[,2*i]),length=10000)
    y1=dbeta(x1,parms[i],parms[n+i])
    hist(object[,2*i],freq=FALSE,main="Histogram of P values",xlab="P Values")
    lines(x1,y1,col="red")
    box()
# make quantile plot for beta density
    lab1=paste("Quantiles from a Beta(",round(parms[i],2),",",
      round(parms[n+i],2),") Distribution",sep="")
    x1=qbeta((1:dim(object)[1])/(dim(object)[1]+1),parms[i],parms[n+i])
    qqplot(x1,object[,2*i],xlab=lab1,ylab="Empirical Quantiles",type="l",lwd=2)
    title("Quantile plot of P values vs Beta Dist.")
# make histogram on G log scale for F statistic include Null distribution
    z=hist(asinh(object[,2*i-1]),plot=FALSE)
    x1=seq(min(object[,2*i-1]),max(object[,2*i-1]),length=10000)
    if (length(DF)==n+1) y1=df(x1,DF[i],DF[n+1])
    if (length(DF)==2*n) y1=df(x1,DF[i],DF[n+i])
    tp=ifelse(DF[i]==1,max(1,z$density),max(c(y1,z$density))) 
    hist(asinh(object[,2*i-1]),freq=FALSE,axes=FALSE,ylim=c(0,tp),
      main="Histogram of F Statistic",xlab="F Statistic: G. Log Scale")
    lines(asinh(x1),y1,col="red")
    axis(side=2); box()
    axis(side=1,at=quantile(asinh(object[,2*i-1]),na.rm=TRUE),
      label=round(sinh(quantile(asinh(object[,2*i-1]))),1))
# make quantile plot of observed F versus NULL F
    lab1=paste("Quantiles from a F(",DF[i],",",DF[n+1],") Distribution",sep="")
    if (length(DF)==n+1) {
      x1=qf((1:dim(object)[1])/(dim(object)[1]+1),DF[i],DF[n+1])
      lab1=paste("Quantiles from a F(",DF[i],",",DF[n+1],") Distribution",
        sep="") 
    }
    if (length(DF)==n*2) {
      x1=qf((1:dim(object)[1])/(dim(object)[1]+1),DF[i],DF[n+i])
      lab1=paste("Quantiles from a F(",DF[i],",",DF[n+i],") Distribution",
        sep="")
    }
    qqplot(x1,object[,2*i-1],xlab=lab1,ylab="Empirical Quantiles",type="l",
      lwd=2)
    title("Quantile plot of F Stat vs F Dist.")
# label sheet before moving on
    mtext(paste("Graphics for",rownames(parms)[i]),outer=TRUE,
      font=par()$font.main)
    }
}
###############################################################################
print.anova.mlm=function(object,show=1:2,digits=max(3,getOption("digits")-2),signif.stars = getOption("show.signif.stars")) {
# determine which model to display
  jj=class(show)
  ii=1:dim(object)[1]
  if (jj=="logical") show=ii[show]
  if (jj=="character") if (show=="All") show=ii
    else show=ii[match(show,rownames(object[[2]]),nomatch=0)]
  SS=attributes(object)$SS
# if a single model show anova tables for seperate terms
  if (class(SS)=="matrix") {
    DF=attributes(object)$DF
    n0=length(DF)-1
    cat("Analysis of Variance Table\n\n----------------------------\n")
    for (i in show) {
      cat(paste("Response[",i,"]:",sep=""),rownames(object)[i],"\n")
      mat=cbind.data.frame(Df=DF,"Sum Sq"=SS[i,],"Mean Sq"=SS[i,]/DF)
      F=c(mat[1:n0,3]/mat[n0+1,3],NA)
      P=pf(F,DF,DF[n0+1],lower.tail=FALSE)
      mat=cbind.data.frame(mat,"F value"=F,"Pr(>F)"=P)
      print.anova(mat,digits=digits)
      cat("----------------------------\n")
      }
    }
# if more than one model compare to smallest common model
  else if (class(SS)=="list") {
    dm=dim(object)[2]/2
    fm=attributes(object)$form
    DF=attributes(object)$DF
    cat("Analysis of Variance Table\n\n")
    for (j in 1:length(fm)) cat(fm[j])
    cat("----------------------------\n")
    for (i in show) {
      cat(paste("Response[",i,"]:",sep=""),rownames(object)[i],"\n\n")
      Res.Df=c(sum(DF[1,]),DF[,2])
      names(Res.Df)=substring(fm,1,8)
      RSS=c(SS[[1]][i,1]+SS[[2]][i,1],SS[[2]][i,])
      Df=c(NA,DF[,1])
      SSR=c(NA,SS[[1]][i,])
      F=c(NA,object[i,(1:dm)*2-1])
      P=c(NA,object[i,(1:dm)*2])
      mat=cbind.data.frame(Res.Df,RSS,Df,"Sum of Sq"=SSR,F,"Pr(>F)"=P)
      print.anova(mat,digits=digits,signif.stars=signif.stars)
      cat("----------------------------\n")
      }
    }
  if (dim(object)[1]>length(show)) 
    cat(dim(object)[1]-length(show),"models not displayed\n\n")
  invisible()
}
######################################################################
summary.mlm=function(object,correlation=FALSE) {
  Beta=t(coef(object))
  fit=t(object$fitted.values)
  dm=c(dim(Beta),4)
  x1=array(NA,dm)
  x1[,,1]=Beta 
  dimnames(x1)=list(rownames(Beta),colnames(Beta),
    c("Estimate","Std. Error","t value","Pr(>|t|)"))
  x=model.matrix(object)[,!is.na(Beta[1,])]
  des=solve(crossprod(x))
  SSE=deviance(object)
  dfE=object$df.residual
  x1[,!is.na(Beta[1,]),2]=sqrt(outer(SSE/dfE,diag(des)))
  x1[,,3]=x1[,,1]/x1[,,2]
  x1[,,4]=2*pt(abs(x1[,,3]),dfE,lower.tail=FALSE)
  SSM=rowSums((fit-rowMeans(fit))^2)
  dfM=object$rank-attr(object$terms,"intercept")
  Fstat=SSM*dfE/SSE/dfM
  Pval=pf(Fstat,dfM,dfE,lower.tail=FALSE)
  R2=SSM/(SSE+SSM)
  aR2=1-(1-R2)*(dim(fit)[2]-attr(object$terms,"intercept"))/dfE
  x2=cbind("Res. SE"=sqrt(SSE/dfE),Fstat,Pval,"R-Squared"=R2,
    "Adj R-Sqr"=aR2)
  z=list(Coefficients=x1,"FullModelStats"=x2,Residuals=t(object$residuals))
  if (correlation) {
    se=sqrt(diag(des))
    z[[4]]=des/outer(se,se)
    names(z)[4]="Correlations"
    
  }
  structure(z,DF=c(dfM,dfE),class="summary.mlm",call=object$call)
}
######################################################################
print.summary.mlm=function(object,show=1:2,digits=max(3,getOption("digits")-3),signif.stars = getOption("show.signif.stars")) {
  jj=class(show)
  ii=1:dim(object[[2]])[1]
  if (jj=="logical") show=ii[show]
  if (jj=="character") if (show=="All") show=ii
    else show=ii[match(show,rownames(object[[2]]),nomatch=0)]
  cat("Call:\n",deparse(attributes(object)$call))
  cat("\n\n----------------------------\n")
  DF=attributes(object)$DF
  for (i in show) {
    cat(paste("Response[",i,"]:",sep=""),dimnames(object[[1]])[[1]][i],"\n")
    z=quantile(object[[3]][i,])
    names(z)=c("Min","Q1","Median","Q3","Max")
    cat("\nResiduals:\n")
    print(z,digits=digits)
    cat("\nCoefficients:\n")
    printCoefmat(object[[1]][i,,],digits=digits,signif.stars=signif.stars)
    z=object[[2]][i,]
    cat("\nResidual standard error:",format(signif(z[1],digits)))
    cat(" on",DF[2],"degrees of freedom\n")
    cat("Multiple R-Squared:",formatC(z[4],digits=digits))
    cat(",\tAdjusted R-squared:",formatC(z[5],digits=digits))
    cat("\nF-statistic:",formatC(z[2],digits=digits),"on",DF[1])
    cat(" and",DF[2],"DF,  p-value:",format.pval(z[3],digits=digits))
    cat("\n----------------------------\n")
  }
  if (dim(object[[1]])[1]>length(show)) 
    cat(dim(object[[1]])[1]-length(show),"models not displayed\n\n")
  if (length(object)==4) {
    cat("Correlation of Coefficients:\n")
    correl <- format(round(object[[4]],2),nsmall=2,digits=digits)
    correl[!lower.tri(correl)] <- ""
    print(correl[-1, -NCOL(correl),drop=FALSE],quote=FALSE)
    }
  invisible()
}
######################################################################
est.beta<-function(y,method=c("moments","mle"),Var=FALSE) {
# compute moment estimates
  method<-match.arg(method)
  mu<-mean(y)
  sig2<-var(y)
  x<-c(alpha=mu,beta=1-mu)*((mu-mu^2)/sig2-1)
  if (method=="moments") return(x)
# mle function
  f1<-function(x,y) {
    z<-(x[1]-1)*log(y)+(x[2]-1)*log(1-y)-lbeta(x[1],x[2])
    return(-sum(z))
  }
# mle gradient
  g1<-function(x,y) {
    z1<-log(y)+digamma(sum(x))-digamma(x[1])
    z2<-log(1-y)+digamma(sum(x))-digamma(x[2])
    return(-c(sum(z1),sum(z2)))
  }
# if variance required estimate parameters, compute hessian and return
  if (Var) {
    z<-optim(x,f1,g1,method="L-BFGS-B",lower=rep(.Machine$double.eps,2),
      hessian=TRUE,y=y)
    cov.par<-solve(z$hessian)
    par<-cbind(z$par,sqrt(diag(cov.par)))
    colnames(par)<-c("Est","SE")
    return(list(Estimates=par,Correlation=cov.par[1,2]/
      sqrt(cov.par[1,1]*cov.par[2,2])))
  }
# otherwise just find solution and return
  z<-optim(x,f1,g1,method="L-BFGS-B",lower=rep(.Machine$double.eps,2),y=y)
  return(z$par)
}

