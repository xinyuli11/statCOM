## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20058)

## ----message=FALSE------------------------------------------------------------
data("data")
attach(data)
myst(age)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("data")
attach(data)
mean_ci(age,0.95)

## -----------------------------------------------------------------------------
z = getRandomNum(x = 5,min = 2230,q25 = 18123,q50 = 52213,q75 = 78312,max = 234234123)
quantile(z)
fivenum(z)

## ----fig.height=30, fig.width=30, include=FALSE-------------------------------
data("data")
attach(data)
par(mfrow=c(4,1),mar=rep(0,4))
outlierKD(data,age)

