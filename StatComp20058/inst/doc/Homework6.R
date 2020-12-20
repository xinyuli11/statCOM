## ----warning=FALSE------------------------------------------------------------
set.seed(1234)
library(bootstrap)
b.cor <- function(x,i) cor(x[i,1],x[i,2])
b<-data.matrix(law)
theta.hat <- b.cor(b,1:15)
theta.jack <- numeric(15)
for(i in 1:15){
    theta.jack[i] <- b.cor(b,(1:15)[-i])
}
bias.jack <- (15-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((15-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
        se.jack=se.jack),3)

## -----------------------------------------------------------------------------
library(boot)

## -----------------------------------------------------------------------------
n=1e1
R<-rexp(n,rate=1)
data(R, package = "boot") 
boot.obj <- boot(aircondit, R = 2000, statistic = function(x, i){mean(x[i,1])})
print(boot.ci(boot.obj, type=c("basic","norm","perc","bca")))

## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
mu<-1;n<-1e1;m<-1e4;library(boot);set.seed(1)
boot.mean <- function(x,i) 1/mean(x[i])
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m){
  R<-rexp(n,rate=1)
  de <- boot(data=R,statistic=boot.mean, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
  ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5] }
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),
    'BCa =',mean(ci.bca[,1]<=mu & ci.bca[,2]>=mu))

## -----------------------------------------------------------------------------
n<-nrow(scor)
jack<-numeric(n) 
for (i in 1:n) {
  scor<-scor[-i,]
  cov.e<-eigen(cov(scor))
  lameda<-cov.e$values
  jack[i]<-lameda[1]/sum(lameda)
}
theta.hat<-eigen(cov(scor))$values[1]/sum(eigen(cov(scor))$values)
bias.jack<-(n-1)*(mean(jack)-theta.hat)
se.jack<-sqrt(n-1)*mean((jack - mean(jack))^2)
print(c(bias.jack<-(n-1)*(mean(jack)-theta.hat),se.jack<-sqrt(n-1)*mean((jack - mean(jack))^2)))

## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
library(lattice)
library(DAAG)
attach(ironslag)
n <- length(magnetic)-1  #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  y<-y[-k]
  x <- chemical[-k]
  x<-x[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  yhat12 <- J1$coef[1] + J1$coef[2] * chemical[k+1]
  e1[k] <- (magnetic[k] - yhat1)^2+(magnetic[k+1] - yhat12)^2
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
  yhat22 <- J2$coef[1] + J2$coef[2] * chemical[k+1] +J2$coef[3] * chemical[k+1]^2
  e2[k] <- (magnetic[k] - yhat2)^2+(magnetic[k+1] - yhat22)^2
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[k+1]
  yhat32 <- exp(logyhat32)
  e3[k] <- (magnetic[k] - yhat3)^2+(magnetic[k+1] - yhat32)^2
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
  yhat4 <- exp(logyhat4)
  logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[k+1])
  yhat42 <- exp(logyhat42)
  e4[k] <- (magnetic[k] - yhat4)^2+(magnetic[k+1] - yhat42)^2
}

## -----------------------------------------------------------------------------
c(mean(e1), mean(e2), mean(e3), mean(e4))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(DAAG); attach(ironslag)
L1 <- lm(magnetic ~ chemical)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L3 <- lm(log(magnetic) ~ chemical)
L4 <- lm(log(magnetic) ~ log(chemical))

## ----eval=FALSE, include=FALSE------------------------------------------------
#  par(mfrow = c(1, 2))  #layout for graphs
#  plot(L2$fit, L2$res)  #residuals vs fitted values
#  abline(0, 0) #reference line
#  qqnorm(L2$res) #normal probability plot
#  qqline(L2$res) #reference line

## ----message=FALSE, warning=FALSE---------------------------------------------
library(DAAG)
attach(ironslag)
a <- seq(10, 40, .1)
n <- length(magnetic)
c<-combn(n,2)
m<-ncol(c)
e1 <- e2 <- e3 <- e4 <- numeric(m)
for(i in 1:m){
  x<-magnetic[-c[,i]]
  y<-chemical[-c[,i]]
  J1 <- lm(y ~ x)
  yhat11 <- J1$coef[1] + J1$coef[2] * chemical[c[1,i]] 
  yhat12 <- J1$coef[1] + J1$coef[2] * chemical[c[2,i]] 
  e1[i] <- (magnetic[c[1,i]]-yhat11)^2+(magnetic[c[2,i]] - yhat12)^2
  
  J2 <- lm(y ~ x + I(x^2))
  yhat21 <- J2$coef[1] + J2$coef[2] * chemical[c[1,i]] + J2$coef[3] * chemical[c[1,i]]^2
  yhat22 <- J2$coef[1] + J2$coef[2] * chemical[c[2,i]] + J2$coef[3] * chemical[c[2,i]]^2
  e2[i] <- (magnetic[c[1,i]]-yhat21)^2+(magnetic[c[2,i]] - yhat22)^2
  
  J3 <- lm(log(y) ~ x)
  logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[c[1,i]] 
  logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[c[2,i]] 
  yhat31<-exp(logyhat31)
  yhat32<-exp(logyhat32)
  e3[i] <- (magnetic[c[1,i]]-yhat31)^2+(magnetic[c[2,i]] - yhat32)^2
  
  J4 <- lm(log(y) ~ log(x))
  logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[c[1,i]])
  logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[c[2,i]]) 
  yhat41 <- exp(logyhat41) 
  yhat42 <- exp(logyhat42) 
  e4[i] <- (magnetic[c[1,i]] - yhat41)^2+(magnetic[c[2,i]] - yhat42)^2

}
c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
summary(L3)

