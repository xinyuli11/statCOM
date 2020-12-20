## -----------------------------------------------------------------------------
Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=(exp(-abs(y))/exp(-abs(x[i-1]))))
      x[i]<-y else{
        x[i]<-x[i-1]
        k<-k+1
    }
  }
  return(list(x=x,k=k))
}
N<-3000
sigma<-c(0.05,1,5,16)
x0<-25
r1<-Metropolis(sigma[1],x0,N)
r2<-Metropolis(sigma[2],x0,N)
r3<-Metropolis(sigma[3],x0,N)
r4<-Metropolis(sigma[4],x0,N)
print(1-c(r1$k,r2$k,r3$k,r4$k)/N)
index<-1:3000
X1<-Metropolis(sigma[1],x0,N)$x[index]
X2<-Metropolis(sigma[2],x0,N)$x[index]
X3<-Metropolis(sigma[3],x0,N)$x[index]
X4<-Metropolis(sigma[4],x0,N)$x[index]
opar<-par(no.readonly = T)
plot(index,X1,type = "l",xlab = "sigma=0.05")
plot(index,X2,type = "l",xlab = "sigma=1")
plot(index,X3,type = "l",xlab = "sigma=5")
plot(index,X4,type = "l",xlab = "sigam=16")

## -----------------------------------------------------------------------------
R<-function(psi){
 psi<-as.matrix(psi)
 n<-ncol(psi)
 k<-nrow(psi)
 psi.means<-rowMeans(psi)
 B<-n*var(psi.means)
 psi.w<-apply(psi,1,"var")
 W<-mean(psi.w)
 v.hat<-W*(n-1)/n+(B/n)
 r.hat<-v.hat/W
 return(r.hat)
}
Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=(exp(-abs(y))/exp(-abs(x[i-1]))))
      x[i]<-y else{
        x[i]<-x[i-1]
        k<-k+1
      }
  }
  return(x)
}



## -----------------------------------------------------------------------------
sigma1<-0.05
N<-10000
b<-2000
K<-4
x0<-c(-8,-4,4,8)
X<-matrix(0,K,N)
for(i in 1:K)
  X[i,]<-Metropolis(sigma1,x0[i],N)
psi1<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi1))
  psi1[i,]<-psi1[i,]/(1:ncol(psi1))
par(mar = c(4,3,1,1)+0.1)
for(i in 1:K)
  plot(psi1[i,(b+1):N],type="l",xlab = i,ylab = bquote(psi1))

sigma2<-1
for(i in 1:K)
  X[i,]<-Metropolis(sigma2,x0[i],N)
psi2<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi2))
  psi2[i,]<-psi2[i,]/(1:ncol(psi2))
print(R(psi2))
par(mar = c(4,3,1,1)+0.1)
for(i in 1:K)
  plot(psi2[i,(b+1):N],type="l",xlab = i,ylab = bquote(psi2))

sigma3<-5
for(i in 1:K)
  X[i,]<-Metropolis(sigma3,x0[i],N)
psi3<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi3))
  psi3[i,]<-psi3[i,]/(1:ncol(psi3))
print(R(psi3))
par(mar = c(4,3,1,1)+0.1)
for(i in 1:K)
  plot(psi3[i,(b+1):N],type="l",xlab = i,ylab = bquote(psi))

sigma4<-16
for(i in 1:K)
  X[i,]<-Metropolis(sigma4,x0[i],N)
psi4<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi4))
  psi4[i,]<-psi4[i,]/(1:ncol(psi4))
print(R(psi4))
par(mar = c(4,3,1,1)+0.1)
for(i in 1:K)
  plot(psi4[i,(b+1):N],type="l",xlab = i,ylab = bquote(psi4))

## -----------------------------------------------------------------------------
cn<-c("sigma=0.05","sigma=1","sigma=5","sigma=16")
rn<-c("R")
sigma1<-0.05
for(i in 1:K)
  X[i,]<-Metropolis(sigma1,x0[i],N)
psi1<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi1))
  psi1[i,]<-psi1[i,]/(1:ncol(psi1))

sigma2<-1
for(i in 1:K)
  X[i,]<-Metropolis(sigma2,x0[i],N)
psi2<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi2))
  psi2[i,]<-psi2[i,]/(1:ncol(psi2))

sigma3<-5
for(i in 1:K)
  X[i,]<-Metropolis(sigma3,x0[i],N)
psi3<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi3))
  psi3[i,]<-psi3[i,]/(1:ncol(psi3))

sigma4<-16
for(i in 1:K)
  X[i,]<-Metropolis(sigma4,x0[i],N)
psi4<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi4))
  psi4[i,]<-psi4[i,]/(1:ncol(psi4))

data<-c(R(psi1),R(psi2),R(psi3),R(psi4))
matrix(data,1,4,dimnames =list(rn,cn))
rh1<-rep(0,N)
for(j in (b+1):N)
  rh1[j]<-R(psi1[,1:j])
rh2<-rep(0,N)
for(j in (b+1):N)
  rh2[j]<-R(psi2[,1:j])
rh3<-rep(0,N)
for(j in (b+1):N)
  rh3[j]<-R(psi3[,1:j])
rh4<-rep(0,N)
for(j in (b+1):N)
  rh4[j]<-R(psi4[,1:j])

## -----------------------------------------------------------------------------
plot(rh1[(b+1):N],type = "l",xlab = "sigma=0.05",ylab = "R")
abline(h=1.2,lty=2)
plot(rh2[(b+1):N],type = "l",xlab = "sigma=1",ylab = "R")
abline(h=1.2,lty=2)
plot(rh3[(b+1):N],type = "l",xlab = "sigma=5",ylab = "R")
abline(h=1.2,lty=2)
plot(rh4[(b+1):N],type = "l",xlab = "sigma=16",ylab = "R")
abline(h=1.2,lty=2)

