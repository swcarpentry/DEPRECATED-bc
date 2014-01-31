# modes: numeric, character, logical, complex, list, function (raw and complex)
# basic classes: numeric, integer, character, logical, list, function (raw and complex)
x = numeric(10)
x = character(10)
x = logical(10)
x = vector("list",10) 
x = function() {}

# What is an attribute?
# basic attributes: class, dim, dimnames, names, row.names, comment

x = vector("list",10)
names(x) = 101:110
comment(x)="This is a silly list"

for (i in 1:10) x[[i]]=rnorm(10)

#-----------------------------------------------------------------------
# create a matrix of lists

dim(x)=c(5,2)
dimnames(x)=list(letters[1:5],LETTERS[1:2])

x[1,1]
x[,2]
x[1,]

# look at
# help(attributes)
# help(attr)
# help(structure)

# create your own array without using the array command

x=1:100
dim(x)<-c(2,5,10)
class(x)
mode(x)
attributes(x)


# use gl to create a factor
x=gl(5,1,20,labels=LETTERS[1:5])
# examine the object
class(x)
mode(x)
attributes(x)

# create your own factor without using gl or factor
x1=rep(1:5,4)
levels(x1)=LETTERS[1:5]
class(x1)="factor"
x1
unclass(x1)

# lets look at a data frame
x=vector("list",5)
names(x)=LETTERS[1:5]
for (i in 1:5) x[[i]]=rnorm(20)
x1=data.frame(x,row.names=1:20)
mode(x1)
attributes(x1)

# turn x into a dataframe without using the data.frame command

class(x)="data.frame"
x
row.names(x)=1:20
x
mode(x)
class(x)

# Look at more involved classes

m1=lm(A~B+C,data=x1)
s1=summary(m1)
a1=anova(m1)

class(m1)
class(s1)
class(a1)

mode(m1)
mode(s1)
mode(a1)

# S3 object oriented nature of R

methods(print)
methods(class=lm)

print.summary.lm
# use getAnywhere(print.summary.lm) to view

#UseMethod, NextMethod
print
UseMethod

y=matrix(rnorm(100000),ncol=100,dimnames=list(NULL,outer(LETTERS[1:25],1:4,paste,sep="")))
A=gl(2,500)
B=gl(2,250,1000)
m1=lm(y~A)
m2=lm(y~B)
m3=lm(y~A+B)

class(m1)
class(m1) <- c("bob",class(m1))
class(m1)

anova(m1)
anova.bob <- function(...) {
  "HaHaHa"
}
anova(m1)
class(m1) <- class(m1)[-1]


s1=summary(m1)
class(s1)
names(s1)
s1[[1]]
class(s1[[1]])

# look at my version of anova.mlm
source("anova.mlm.R")

s2=summary(m1)
class(s2)
names(s2)
class(s2[[1]])
dimnames(s2$Coefficients)
s2$Coefficients[1:5,2,]
class(s2[[2]])
dimnames(s2[[2]])
print(s2,show=sample(100,2))
s2$FullModelStats[1:5,]

a0=stats::anova.mlm(m1)
a1=anova(m1)
print(a1,show=sample(100,2))

anova(m1,m3)
anova(m1,m2,m3)