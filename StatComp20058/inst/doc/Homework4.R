## -----------------------------------------------------------------------------
set.seed(123342)
theta.hat <- se <- numeric(3)

u <- runif(100000) 
x <- 1-log(1-u)
g<-function(x){
  ((1/(2*pi)^(1/2))*(x^2)*exp(-(x^2)/2))*(x>1)
}
fg <- g(x) / (exp(-x+1))
theta.hat[1]=mean(fg)
se[1]=sd(fg)

M<-runif(10000)
c<-(-2*log(1-M))^(1/2)
g<-function(x){
  ((1/(2*pi)^(1/2))*((x+1)^2)*exp(-((x+1)^2)/2))
}
fg <- g(c) / (c*exp(-(c^2)/2))
theta.hat[2]=mean(fg)
se[2]=sd(fg)

m<-10000
g<-function(x){
  (x^2/(2*pi)^0.5)*exp(-0.5*x^2)*(x>1)
}
u<-runif(m)
x<-tan(pi*u/2)+1
fg<-g(x)/((2/pi)*(1+(x-1)^2)^(-1))
theta.hat[3]<-mean(fg)
se[3]<-sd(fg)

rbind(theta.hat, se)


## ----paged.print=TRUE---------------------------------------------------------
c1=d1=numeric(5)
u <- runif(2000)
x1 <- - log((1 - u * (1- exp(-1/5))))
g1<-function(x){
  exp(-x - log(1+x^2))*(x > 0) * (x < 1/5)
}
fg1 <- g1(x1)/(((exp(-x1)))/(1 - exp(-1/5)))
c1[1] = mean(fg1)
d1[1] = var(fg1)


x2 <- - log((exp(-1/5) - u * (exp(-1/5)- exp(-2/5))))
g2<-function(x){
  exp(-x - log(1+x^2))*(x > 1/5) * (x < 2/5)
}
fg2 <- g2(x2)/(((exp(-x2)))/(exp(-1/5) - exp(-2/5)))
c1[2] = mean(fg2)
d1[2] = var(fg2)


x3 <- - log((exp(-2/5) - u * (exp(-2/5)- exp(-3/5))))
g3<-function(x){
  exp(-x - log(1+x^2))*(x > 2/5) * (x < 3/5)
}
fg3 <- g3(x3)/(((exp(-x3)))/(exp(-2/5)- exp(-3/5)))
c1[3] = mean(fg3)
d1[3] = var(fg3)


x4 <- - log((exp(-3/5) - u * (exp(-3/5)- exp(-4/5))))
g4<-function(x){
  exp(-x - log(1+x^2))*(x > 3/5) * (x < 4/5)
}
fg4 <- g4(x4)/(((exp(-x4)))/(exp(-3/5)- exp(-4/5)))
c1[4] = mean(fg4)
d1[4] = var(fg4)


x5 <- - log(exp(-4/5) - u * (exp(-4/5)- exp(-1)))
g5<-function(x){
  exp(-x - log(1+x^2))*(x >4/5) * (x < 1)
}
fg5 <- g5(x5)/(((exp(-x5)))/(exp(-4/5)- exp(-1)))
c1[5] = mean(fg5)
d1[5] = var(fg5)

theta.hat <- se <- numeric(2)
theta.hat[2]=round(sum(c1),5)
se[2]=round(sum(d1),5)
theta.hat[1]=round(0.5257801,5)
se[1]=round(0.0970,5)
rbind(theta.hat, se)

## -----------------------------------------------------------------------------
set.seed=123
UCL <- replicate(10000, expr = {
  s = rlnorm(20, meanlog = 0, sdlog = 1)
  a = log(s)
  sqrt(20)*abs(mean(a))/((20/19)*var(a))^(1/2)
} )
print(c(sum(UCL<qt(0.975,19)),mean(UCL<qt(0.975,19))))


## -----------------------------------------------------------------------------
n<-20
m<-1000
set.seed(12)
#problem6.5
CL1=matrix(0,m,2)
for(i in 1:m){
  x<-rchisq(n,2)
  CL1[i,2]=mean(x)+((((20/19)*var(x))^(1/2))*qt(0.975,n-1))/(n^0.5)
  CL1[i,1]=mean(x)-((((20/19)*var(x))^(1/2))*qt(0.975,n-1))/(n^0.5)
}
a1=mean(CL1[,1]<=2&CL1[,2]>=2)
a2=sd(CL1[,1]<=2&CL1[,2]>=2)
print(c(a1,a2))

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rchisq(n, df = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
print(c(mean(UCL>4),sd(UCL)))

