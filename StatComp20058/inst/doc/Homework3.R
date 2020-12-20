## -----------------------------------------------------------------------------
m <- 1e6; x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * (pi/3)
print(c(theta.hat,-cos(pi/3) + cos(0)))

## -----------------------------------------------------------------------------
#简单蒙特卡罗模拟积分
m <- 1e6; x <- runif(m, min=0, max=1)
theta.hat <- mean(exp(x)) 
print(c(theta.hat,exp(1)-exp(0)))

## -----------------------------------------------------------------------------
#antithetic variate approach
seed=1234
m <- 1e8
u <- runif(m/2,0,1)
v <- 1 - u
print(c(theta.hat,(sum(exp(u))+sum(exp(v)))/m,exp(1)-exp(0)))

## -----------------------------------------------------------------------------
#计算方差
sd1<-sd(exp(u)+exp(v))
sd2<-sd(exp(x))
sd3<-sd1/sd2
print(c(sd1,sd2,sd3))

