## ----warning=FALSE------------------------------------------------------------
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 1000
nofbeta<- c(seq(0, 100, 1))
N <- length(nofbeta)
pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  e <- nofbeta[j]
  sktests <- numeric(m)
  for (i in 1:m) { 
    x <- rbeta(n,e,e)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(nofbeta, pwr, ylim = c(0,1))
abline(h = .1, lty = 1)
se <- sqrt(pwr * (1-pwr) / m)
lines(nofbeta, pwr+se, lty = 3)
lines(nofbeta, pwr-se, lty = 3)

## ----warning=FALSE------------------------------------------------------------
alpha <- .1
n <- 30
m <- 1000
noft<- c(seq(0, 100, 1))
N <- length(noft)
pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  e <- noft[j]
  sktests <- numeric(m)
  for (i in 1:m) { 
    x <- rt(n,e)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(noft, pwr, ylim = c(0,1))
abline(h = .1, lty = 1)
se <- sqrt(pwr * (1-pwr) / m)
lines(noft, pwr+se, lty = 3)
lines(noft, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
power1<-numeric(2)
set.seed(2334)
m=1000
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5.5))
}
sigma1 <- 1
sigma2 <- 1.5
power <-(replicate(m, expr={
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  count5test(x, y)
}))
power1[1]<-mean(power)

sigma1 <- 1
sigma2 <- 1.5
pvalues <- (replicate(m, expr={
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  Ftest = var.test(y, x, ratio = 1)
  Ftest$p.value
}))
power1[2] = mean(pvalues <= .055)
rbind(power1)


## ----warning=FALSE------------------------------------------------------------
n <- c(10, 20, 30, 40, 50, 60,70,80,90,100,150)
m=1000
sigma1 <- 1
sigma2 <- 1.5
power<-power2 <- numeric(length(n)) 
for (i in 1:length(n)) {
  sktests <- numeric(m)  
  for (j in 1:m) {
    x <- rnorm(n[i],0,sigma1)
    y<-rnorm(n[i],0,sigma2)
    sktests[j]= count5test(x, y)
  }
  power[i] <- mean(sktests) 
}
for (i in 1:length(n)) {
  sktests <- numeric(m)  
  for (j in 1:m) {
    x <- rnorm(n[i],0,sigma1)
    y<-rnorm(n[i],0,sigma2)
    Ftest = var.test(y, x, ratio = 1)
    sktests[j]= Ftest$p.value
  }
  power2[i] <- mean(sktests<0.055) 
}
rbind(n,power,power2)
library(ggplot2)
ggplot(ylab="POWER")+geom_line(data=data.frame(n,power),aes(n,power),col = 'red')+geom_line(data=data.frame(n,power2),aes(n,power2),col = 'black')

## ----warning=FALSE------------------------------------------------------------
library(MVN)
library(MASS)
n <- c(10, 50, 100, 150, 200, 500)
p.values <-Statistic<- numeric(6)
for (i in 1:length(n)) { 
  Sigma <- matrix(c(1,0.2,0.2,100),2,2)
  U <- mvrnorm(n[i], rep(0, 2), Sigma)
  sk = mvn(U, mvnTest = c("mardia"))
  sk1 = sk$multivariateNormality[3]
  sk2 = sk$multivariateNormality[2]
  a = data.frame(sk1[1])
  b = data.frame(sk2[1])
  pvalues1 =as.numeric(as.character(a$p.value[1]))
  pvalues2 =as.numeric(as.character(b$Statistic[1]))
  p.values[i]<- pvalues1
  Statistic[i]<- pvalues2
}
rbind(n,Statistic,p.values)

## ----warning=FALSE------------------------------------------------------------
library(MASS)
library(MVN)
set.seed(1234)
alpha <-.05
n <-30
m <- 2500
epsilon <- c(seq(0,1,.1))
N <- length(epsilon)
pwr <- numeric(N)
cv <- qchisq(1-alpha/2,(2*(2+1)*(2+2)/6))
for (j in 1:N) { 
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) {  
    Sigma <- matrix(c(1,e,e,100),2,2)
    U <- mvrnorm(n, rep(0, 2), Sigma)
    sk = mvn(U, mvnTest = c("mardia"))
    sk2 = sk$multivariateNormality[2]
    b = data.frame(sk2[1])
    pvalues2 =as.numeric(as.character(b$Statistic[1]))
    sktests[i] <- as.integer(abs(n*pvalues2/6) >= cv)
  }
  pwr[j] <- mean(sktests)
}
n<-c(seq(0,1,.1))
rbind(n,pwr)

## -----------------------------------------------------------------------------
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

