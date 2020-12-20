## -----------------------------------------------------------------------------
library(lattice)
n <- seq(20, 50, 10)
x <- rnorm(sum(n))
y <- factor(rep(n, n), labels=paste("n =", n))
densityplot(~ x | y,
panel = function(x, ...) {
panel.densityplot(x, col="red", ...)
panel.mathdensity(dmath=dnorm,
args=list(mean=mean(x), sd=sd(x)),
col="darkblue")
})

## -----------------------------------------------------------------------------
x1<-c(0.10, 0.11, 0.12, 0.13, 0.14, 0.15,
0.16, 0.17, 0.18, 0.20)
x2<-c(1,2,3,4,5,5,4,3,2,1)
y<-c(42.0, 43.5, 49.0, 45.5, 45.0, 47.5,
49.0, 53.0, 53.0, 55.0)
lm.sol<-lm(y ~ x1+x2)
knitr::kable(summary(lm.sol)$coef)

