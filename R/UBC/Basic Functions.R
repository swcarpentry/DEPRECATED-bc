#writing functions in R

#function.name <- function(arguments of the form: tag or tag = value, may be many arguments) {
#             function body is a set of commands 
#  }

#help("function")

# methods of return

x<-rnorm(5)
y<-rnorm(5)

mp1<-function(x,y) {
  sum<-x+y
  diff<-x-y
  list(sum=sum,diff=diff)
}

mp2<-function(x,y) {
  sum<-x+y
  diff<-x-y
  invisible(list(sum=sum,diff=diff))
}

mp3<-function(x,y) {
  sum<-x+y
  diff<-x-y
  print(list(sum=sum,diff=diff))
}

mp4<-function(x,y) {
  sum<-x+y
  diff<-x-y
  print(list(sum=sum,diff=diff))
  invisible()
}

# run this code with each example above

mp1(x,y)
z=mp1(x,y)
z

mp2(x,y)
z=mp2(x,y)
z

mp3(x,y)
z=mp3(x,y)
z

mp4(x,y)
z=mp4(x,y)
z

# mp5 acts like which previous version of the function?

mp5<-function(x,y) {
  sum<-x+y
  diff<-x-y
  z<-list(sum=sum,diff=diff)
}


# behaviour to be aware of local versus global variables

fn1<-function() 
  return(x)

fn2<-function(x) 
  return(x)

fn3<-function(x=z) 
  return(x)

# What will be retuned in each case

rm(x)
fn1()
fn2()
fn3()

x=1
fn1()
fn2(x)
fn3(x)

x=1
z=2
fn3()

# arguments passed by value not address
#"what happens in Vegas stays in Vegas"

fn <- function() {
  x=x+100
  print(x)
}

x=1
fn()
x

# Default parameter values
# variable initialized first time it is used inside function

sx1<-function(x=z) {
  z=3
  return(x)
}

sx2<-function(x=z) {
  y=x
  z=3
  return(x)
}

z=2
x=1
sx1()
sx2()


# specifying x,y,z as arguments makes them local variables global versions are ignored
sx<-function(x=z,y=x,z=y) {
  return(x+y+z)
}

x=1
y=2
z=3
sx()
sx(1)
sx(1,2)
sx(2,1)
sx(y=2,x=1)

# careful when using <- in function call not the same as = in this case
# see help("=")

rm(x,y,z)
x
y
sx(2,1)
sx(1,2)
sx(y<-2,x<-1)
x
y

# careful of <<-  notice what happens locally and globally

sx<-function() {
  x<-1
  x<<-2
  return(x)
}

rm(x)
x
sx()
x

#--------------------------------------------
# Advanced Functions 
# notice ... and match.fun

porder<-function(q,i=1,n=1,dist="norm",...) {
  FUN<-match.fun(paste("p",dist,sep=""))
  return(1-pbinom(i-1,n,FUN(q,...))) 
}

# get with mode="function" also works for match.fun
qorder<-function(p,i=1,n=1,dist="norm",...) {
print(c(NULL,...))
  FUN<-get(paste("q",dist,sep=""),mode="function")
  p1<-1-qbeta(1-p,n-i+1,i)
  return(FUN(p1,...)) 
}


qorder(.5,,10)
qorder(.5,n=10)
qorder(.5,n=10,i=1)
qorder(.5,mean=3,sd=10)
qorder(.5,dist="pois",lambda=10)

#----------------------------------------------
# Advanced Functions #2 
# notice match.arg for method to set finite set of default values
# notice function has 3 return points, first one hit leaves function
# notice functions defined inside other function is okay

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
    z<-optim(x,f1,g1,method="L-BFGS-B",lower=rep(.Machine$double.eps,2),hessian=TRUE,y=y)
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

x<-rbeta(1000,1.2,.8)
est.beta(x)
est.beta(x,"mle")
est.beta(x,"mle",T)

est.beta(x,"bob")

#--------------------------------------------------------
# Special Function forms

# only 2 arguments are possible
"%+%" <- function(x,y)
  return(paste(x,y,sep=" + "))

x <- 1:10
y <- 10:1
x%+%y

# value is the required argument name
"inc<-" = function(x,value)
  return(x+value)

inc(x) <- 3
x
