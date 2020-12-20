## -----------------------------------------------------------------------------
library(StatComp20058)
x0<-25
sigma<-c(0.05,1,2,16)
N<-3000
M1=MetropolisC(sigma[1],x0,N)
M2=MetropolisC(sigma[2],x0,N)
M3=MetropolisC(sigma[3],x0,N)
M4=MetropolisC(sigma[4],x0,N)
acceptanceC<-1-c(M1[[2]][1],M2[[2]][1],M3[[2]][1],M4[[2]][1])/N
index<-1:2000
My1<-M1[[1]][index]
My2<-M2[[1]][index]
My3<-M3[[1]][index]
My4<-M4[[1]][index]
opar<-par(no.readonly = T)
par(mar = c(4,3,1,1)+0.1)
plot(index,My1,type = "l",xlab = "sigma=0.05")
plot(index,My2,type = "l",xlab = "sigma=1")
plot(index,My3,type = "l",xlab = "sigma=2")
plot(index,My4,type = "l",xlab = "sigam=16")

## -----------------------------------------------------------------------------
library(StatComp20058)
MM1=MetropolisR(sigma[1],x0,N)
MM2=MetropolisR(sigma[2],x0,N)
MM3=MetropolisR(sigma[3],x0,N)
MM4=MetropolisR(sigma[4],x0,N)
acceptanceR<-1-c(M1$k,M2$k,M3$k,M4$k)/N
list(acceptanceC=acceptanceC,acceptanceR=acceptanceR)
opar<-par(no.readonly = T)
par(mar = c(4,3,1,1)+0.1)
qqplot(M1[[1]],MM1$x)
qqplot(M2[[1]],MM2$x)
qqplot(M3[[1]],MM3$x)
qqplot(M4[[1]],MM4$x)

## -----------------------------------------------------------------------------
library(microbenchmark)
library(StatComp20058)
tS1<-microbenchmark(MC=MetropolisC(sigma[1],x0,N),MR=MetropolisR(sigma[1],x0,N))
tS2<-microbenchmark(MC=MetropolisC(sigma[2],x0,N),MR=MetropolisR(sigma[2],x0,N))
tS3<-microbenchmark(MC=MetropolisC(sigma[3],x0,N),MR=MetropolisR(sigma[3],x0,N))
tS4<-microbenchmark(MC=MetropolisC(sigma[4],x0,N),MR=MetropolisR(sigma[4],x0,N))
TI1<-summary(tS1)
TI2<-summary(tS2)
TI3<-summary(tS3)
TI4<-summary(tS4)

## -----------------------------------------------------------------------------
print(TI1)

## -----------------------------------------------------------------------------
print(TI2)

## -----------------------------------------------------------------------------
print(TI3)

## -----------------------------------------------------------------------------
print(TI4)

