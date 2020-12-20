## -----------------------------------------------------------------------------
library(actuar)
set.seed(1234)
x<-runif(100)
a=2
b=2
U=runif(100)
X<-2/((1-U)^(1/2))
hist(X,prob=T,col=gray(0.9),main="pareto(2,2) from uniform")
curve(dpareto(x,2,2),add=T,col="red",lwd=2)

## -----------------------------------------------------------------------------
set.seed(12345)
x1mn<-runif(100000,-1,1)
x2mn<-runif(100000,-1,1)
x3mn<-runif(100000,-1,1)
data<-data.frame(x1=x1mn,x2=x2mn,x3=x3mn)
head(data,5)
result3=""
m<-1
for(i in 1:length(x1mn)){
  if(abs(x3mn[i])>=abs(x2mn[i]) & abs(x3mn[i])>=abs(x1mn[i])){
    result3[m]<-x2mn[i]
  }else{
    result3[m]<-x3mn[i]
  }
  m<-m+1
}
head(result3,5)

## -----------------------------------------------------------------------------
result3<-as.numeric(result3)
hist(result3,prob=T,col=gray(0.9),main="Epanechnikov kernel from uniform")
lines(density(result3),lwd=2,col="red")

## -----------------------------------------------------------------------------
set.seed(1234)
x<-runif(100000,-1,1)
hist(result3,prob=T,col=gray(0.9),main="pareto(2,2,4) from uniform")
curve(3*(1-x^2)/4,add=T,col="red",lwd=2)

## -----------------------------------------------------------------------------
set.seed(1234)
a=2
b=4
U=runif(1000)
X1<-(a/((1-U)^(1/b)))-a
hist(X1,prob=T,col=gray(0.9),main="the density histogram of the sample")
lines(density(X1,kernel="epanechnikov"),col="red",lwd=2)

## ----eval=FALSE, message=FALSE, include=FALSE---------------------------------
#  library(qcc)
#  pareto.chart(X1, ylab = "Frequency",ylab2 ="per%", main='curve',lwd=1)

