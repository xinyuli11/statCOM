## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20058)

## ----message=FALSE------------------------------------------------------------
function(x, conf) {
  se <- sd(x) / sqrt(length(x))
  alpha <- 1 - conf
  mean(x) + se * qnorm(c(alpha / 2, 1 - alpha / 2))
}
data("data")
attach(data)
myst(age)

## ----message=FALSE, warning=FALSE---------------------------------------------
function(x){
  m<-mean(x)
  n<-length(x)
  s<-sd(x)
  skew<-sum((x-m^3/s^3))/n
  kurt<-sum((x-m^4/s^4))/n-3
  return(c(n=n,mean=m,stdev=s,skew=skew,kurtosis=kurt))
}
data("data")
attach(data)
mean_ci(age,0.95)

## -----------------------------------------------------------------------------
function(x,min,q25,q50,q75,max){
  if (x%%2 == 0 & (x/2)%%2 == 0) {
    tmp1 = runif(x/4-2,min,q25-1)
    tmp2 = runif(x/4-2,q25+1,q50-1)
    tmp3 = runif(x/4-2,q50+1,q75-1)
    tmp4 = runif(x/4-2,q75+1,max)
    y = c(min,tmp1,q25-1,q25+1,tmp2,q50-1,q50+1,tmp3,q75-1,q75+1,tmp4,max)
  } else if(x%%2 == 0 & (x/2)%%2 != 0){
    tmp1 = runif(ceiling(x/4-2),min,q25)
    tmp2 = runif(ceiling(x/4-2),q25,q50-1)
    tmp3 = runif(ceiling(x/4-2),q50+1,q75)
    tmp4 = runif(ceiling(x/4-2),q75,max)
    y = c(min,tmp1,q25,tmp2,q50-1,q50+1,tmp3,q75,tmp4,max)
  } else if(x%%2 != 0 & (x-1)%%4 == 0){
    tmp1 = runif((x-1)/4-1,min,q25)
    tmp2 = runif((x-1)/4-1,q25,q50)
    tmp3 = runif((x-1)/4-1,q50,q75)
    tmp4 = runif((x-1)/4-1,q75,max)
    y = c(min,tmp1,q25,tmp2,q50,tmp3,q75,tmp4,max)
  } else if(x%%2 != 0 & (x-1)%%4 != 0){
    tmp1 = runif(ceiling((x-1)/4-2),min,q25-1)
    tmp2 = runif(ceiling((x-1)/4-2),q25+1,q50-1)
    tmp3 = runif(ceiling((x-1)/4-2),q50+1,q75-1)
    tmp4 = runif(ceiling((x-1)/4-2),q75+1,max)
    y = c(min,tmp1,q25-1,q25+1,tmp2,q50,tmp3,q75-1,q75+1,tmp4,max)
  }
  return(y)
}
z = getRandomNum(x = 5,min = 2230,q25 = 18123,q50 = 52213,q75 = 78312,max = 234234123)
quantile(z)
fivenum(z)

## ----message=FALSE------------------------------------------------------------
function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  tot <- sum(!is.na(var_name))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "\n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / tot*100, 1), "\n")
  cat("Mean of the outliers:", round(mo, 2), "\n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "\n")
  cat("Mean if we remove outliers:", round(m2, 2), "\n")
}
data("data")
attach(data)
par(mfrow=c(4,1),mar=rep(2,4))
outlierKD(data,age)

## -----------------------------------------------------------------------------
function(x,y,error,maxiter,stepmethod,step,alpha,beta)
{
  m<-nrow(x)
  x<-cbind(matrix(1,m,1),x)
  n<-ncol(x)
  theta<-matrix(rep(0,n),n,1)  
  iter<-1
  newerror<-1
  while((newerror>error)|(iter<maxiter)){
    iter<-iter+1
    h<-x%*%theta  
    des<-t(t(h-y)%*%x)  
    if(stepmethod==T){
      sstep=1
      new_theta<-theta-sstep*des
      new_h<-x%*%new_theta
      costfunction<-t(h-y)%*%(h-y)  
      new_costfunction<-t(new_h-y)%*%(new_h-y)
      while(new_costfunction>costfunction-alpha*sstep*sum(des*des)){
        sstep<-sstep*beta
        new_theta<-theta-sstep*des
        new_h<-x%*%new_theta
        new_costfunction<-t(new_h-y)%*%(new_h-y)  
      }
      newerror<-t(theta-new_theta)%*%(theta-new_theta)       
      theta<-new_theta     
    }
    
    if(stepmethod==F){        
      new_theta<-theta-step*des
      new_h<-x%*%new_theta
      newerror<-t(theta-new_theta)%*%(theta-new_theta)
      theta<-new_theta 
    }
    
  }
  costfunction<-t(x%*%theta-y)%*%(x%*%theta-y)
  result<-list(theta,iter,costfunction)
  names(result)<-c('Coefficient','iteration number','error')
  result
}
data("irisl")
x<-matrix(iris[1:50,1],50,1)
y<-matrix(iris[1:50,2],50,1)
l<-lm(y~x)
summary(l)

## -----------------------------------------------------------------------------
GradientDescent(x,y,1e-14,1000,stepmethod=T,step=0.001,alpha=0.25,beta=0.8)

