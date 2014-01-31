y <- matrix(rnorm(100000),ncol=100,dimnames=list(NULL,outer(LETTERS[1:25],1:4,paste,sep="")))
A <- gl(2,500)
B <- gl(2,250,1000)
m1 <- lm(y~A)
m2 <- lm(y~B)
m3 <- lm(y~A+B)

anova(m1)
summary(m1)

library(anova.mlm)

summary(m1)

a1=anova(m1)

stats::anova.mlm(m1)
