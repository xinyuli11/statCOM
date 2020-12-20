## -----------------------------------------------------------------------------
dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
                 Frequency=c('p^2','q^2','r^2','2pr','2qr','2pq','1'),
                 Count=c('AA','BB','OO','AO','BO','AB','n'))
dat

## -----------------------------------------------------------------------------
est <- function(ep = 1e-8){
  M = 0.2; N= 0.1
  k = 1
  repeat{
    k = k+1
    M[k]=(444*M[k-1]/(2-M[k-1]-2*N[k-1])+507)/2000
    N[k]=(132*N[k-1]/(2-N[k-1]-2*M[k-1])+195)/2000
    if(abs(M[k]-M[k-1])<ep | abs(N[k]-N[k-1])<ep)break
  }
  list(M=M[k], N=N[k], iter=k)
}
est()

## ----message=FALSE, warning=FALSE---------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
attach(mtcars) 
formulas <- as.character(formulas)
for (i in 1:length(formulas)){
  lm <- lm(formulas[i])
  print(lm)
}
lapply(formulas, function(x) lm(x))

## -----------------------------------------------------------------------------
t <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(t, function(f) f$p.value)
sapply(t, "[[", 3)

