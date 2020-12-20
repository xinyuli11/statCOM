## -----------------------------------------------------------------------------
library(boot)
library(energy)
library(MASS)
set.seed(1234)
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy))) 
}
n1 <-20
n2 <-30
mu1 <- 1
mu2 <- 4
sigma1 <-1 
sigma2 <- 6
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
d0<-maxout(x, y)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:50
reps <- numeric(R) #storage for replicates
for (i in 1:R) {
  k <- sample(K, size = 10, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  reps[i] <- maxout(x1, y1)
}
p <- mean(c(d0, reps) >= d0)
p

## -----------------------------------------------------------------------------
hist(reps, breaks = "scott",main = "Histogram",xlab = "a",col="black")
points(d0, 0, cex = 1, pch = 16,col="red")

## -----------------------------------------------------------------------------
library(Ball)
library(RANN)

## -----------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2;set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n2*p,0,2),ncol=p);
  y <- cbind(rnorm(n2,0,2),rnorm(n2,0,1));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

## -----------------------------------------------------------------------------
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
matrix(pow, nrow=1,dimnames=list("power", c("eqdist.nn","eqdist.etest","bd.test")))

## -----------------------------------------------------------------------------
library(Ball)
library(RANN)

## -----------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2;set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n2*p),ncol=p);
  y <- cbind(rnorm(n2,0,1),rnorm(n2,1,1));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

## -----------------------------------------------------------------------------
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
matrix(pow, nrow=1,dimnames=list("power", c("eqdist.nn","eqdist.etest","bd.test")))

## -----------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2;set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
Sigma <- matrix(c(1,0.1,0.1,1),2,2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rt(n2*p,1),ncol=p);
  y <- mvrnorm(n=n2, rep(0,2), Sigma);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

## -----------------------------------------------------------------------------
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
matrix(pow, nrow=1,dimnames=list("power", c("eqdist.nn","eqdist.etest","bd.test")))

## -----------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2;set.seed(12345)
n1 <- 20
n2 <- 200; R<-999; n <- n1+n2; N = c(n1,n2)
Sigma <- matrix(c(1,0.1,0.1,1),2,2)
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- matrix(rnorm(n2*p),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

## -----------------------------------------------------------------------------
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
matrix(pow, nrow=1,dimnames=list("power", c("eqdist.nn","eqdist.etest","bd.test")))

