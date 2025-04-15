getwd()
setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")

# Neyman repeated sampling

rm(list=ls())
library(Matching)
data(lalonde)
z = lalonde$treat
y = lalonde$re78

## Neymanian inference

# calculate point estimator and standard error based on the formulas in Theorem 4.1
n1 = sum(z)
n0 = length(z) - n1
tauhat = mean(y[z==1]) - mean(y[z==0])
vhat = var(y[z==1])/n1 + var(y[z==0])/n0
sehat = sqrt(vhat)

tauhat
# [1] 1794.343
sehat
# [1] 670.9967

## OLS
# Using OLS to estimate the ATE which also gives a standard error.
# However, this se seems too small compared to the one based on Theorem 4.1

olsfit = lm(y ~ z)
summary(olsfit)$coef[2, 1:2]
# Estimate  Std. Error 
# 1794.3431   632.8536 

## Eicker-Huber-White se
# EHW robust standard error.

library(car)
sqrt(hccm(olsfit)[2, 2])
# [1] 672.6823

sqrt(hccm(olsfit, type="hc0")[2, 2])
# [1] 669.3155
sqrt(hccm(olsfit, type="hc2")[2, 2])
# [1] 670.9967

