#' @title A illustration dataset
#' @name data
#' @description A dataset used to illustrate the performance of \code{vaccR} and \code{vaccC}.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' tm <- microbenchmark::microbenchmark(
#'   vR = vaccR(age,female,ily),
#'   vC = vaccC(age,female,ily)
#' )
#' print(summary(tm)[,c(1,3,5,6)])
#' }
NULL


#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{gibbsR} and \code{vaccR}) and Cpp functions (\code{gibbsC} and \code{vaccC}).
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma qnorm runif sd
#' @importFrom grDevices boxplot.stats
#' @importFrom graphics boxplot hist par title
#' @useDynLib StatComp20058
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = gibbsR(100,10),
#'   rnC = gibbsC(100,10)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' 
#' tm2 <- microbenchmark::microbenchmark(
#'   vR = vaccR(age,female,ily),
#'   vC = vaccC(age,female,ily)
#' )
#' print(summary(tm2)[,c(1,3,5,6)])
#' }
NULL

#' @title Use three inputs to predict response using R.
#' @description The prediction model is described in http://www.babelgraph.org/wp/?p=358.
#' @param age the first predictor (numeric)
#' @param female the second predictor (logical)
#' @param ily the third predictor (logical)
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' res <- vaccR(age,female,ily)
#' }
#' @export
vaccR <- function(age, female, ily) {
  p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
  p <- p * ifelse(female, 1.25, 0.75)
  p <- pmax(0, p)
  p <- pmin(1, p)
  p
}

#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsR <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rgamma(1, 3, y * y + 4)
      y <- rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))
    }
    mat[i, ] <- c(x, y)
  }
  mat
}


#' @title A Metropolis sampler using R
#' @description A Metropolis sampler using R
#' @param sigma the sigma of samples
#' @param x0 a parameter
#' @param N sample of size
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#' x0<-25
#' sigma<-c(0.05,1,2,16)
#' N<-3000
#' M1=MetropolisR(sigma[1],x0,N)
#' M2=MetropolisR(sigma[2],x0,N)
#' M3=MetropolisR(sigma[3],x0,N)
#' M4=MetropolisR(sigma[4],x0,N)
#' acceptanceC<-1-c(M1[[2]][1],M2[[2]][1],M3[[2]][1],M4[[2]][1])/N
#' index<-1:2000
#' My1<-M1[[1]][index]
#' My2<-M2[[1]][index]
#' My3<-M3[[1]][index]
#' My4<-M4[[1]][index]
#'opar<-par(no.readonly = T)
#' par(mfrow=c(2,2))
#' plot(index,My1,type = "l",xlab = "sigma=0.05")
#' plot(index,My2,type = "l",xlab = "sigma=1")
#' plot(index,My3,type = "l",xlab = "sigma=2")
#' plot(index,My4,type = "l",xlab = "sigam=16")
#' }
#' @export
MetropolisR<-function(sigma,x0,N){
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

#' @title The confidence interval calculate
#' @description The confidence interval between the two ends of the mean is calculated using an approximate normal distribution
#' @param x   a input list confidence 
#' @param conf  a parameter of confidence
#' @return A confidence interval
#' @examples
#' \dontrun{
#' conf=0.95
#' x <- runif(100)
#' mean_ci(x,conf)
#' mean_ci(x, conf = 0.99)
#' }
#' @export
mean_ci <- function(x, conf) {
  se <- sd(x) / sqrt(length(x))
  alpha <- 1 - conf
  mean(x) + se * qnorm(c(alpha / 2, 1 - alpha / 2))
}


#' @title Basic statistic
#' @description Mean, variance, sample size, kurtosis, skewness
#' @param x  a input list
#' @return A list about some basic statistics
#' @examples
#' \dontrun{
#' x=c(1,2,3,4)
#' mystat(x)
#' }
#' @export
myst<-function(x){
  m<-mean(x)
  n<-length(x)
  s<-sd(x)
  skew<-sum((x-m^3/s^3))/n
  kurt<-sum((x-m^4/s^4))/n-3
  return(c(n=n,mean=m,stdev=s,skew=skew,kurtosis=kurt))
}



#' @title Random number generation based on order statistics
#' @description Random number generation based on order statistics
#' @param x   Sample size
#' @param min  minimum
#' @param q25 quantile
#' @param q50 quantile
#' @param q75 quantile
#' @param max maximum
#' @return A Random sequence
#' @examples
#' \dontrun{
#' z = getRandomNum(x = 5,min = 2230,q25 = 18123,q50 = 52213,q75 = 78312,max = 234234123)
#' quantile(z)
#' fivenum(z)
#' }
#' @export
getRandomNum <-function(x,min,q25,q50,q75,max){
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


#' @title Outlier analysis
#' @description Outliers analysis  for a sequence
#' @param dt  data Frame 
#' @param var  variable
#' @return Outlier analysis result and some plots
#' @examples
#' \dontrun{
#' data("data")
#' attach(data)
#' outlierKD(data,age)
#' }
#' @export
outlierKD <- function(dt, var) {
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
