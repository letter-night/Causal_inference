rm(list=ls())
library(car)

## function for SRE
## estimation
Neyman_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K = length(xlevels)
  PiK = rep(0, K)
  TauK = rep(0, K)
  varK = rep(0, K)
  for (k in 1:K)
  {
    xk = xlevels[k]
    zk = z[x == xk]
    yk = y[x == xk]
    PiK[k] = length(zk) / length(z)
    TauK[k] = mean(yk[zk==1]) - mean(yk[zk==0])
    varK[k] = var(yk[zk==1])/sum(zk) + var(yk[zk==0])/sum(1-zk)
  }
  return (c(sum(PiK*TauK), sqrt(sum(PiK^2*varK))))
}

## Data "nhanes_bmi" from the "ATE" package in R
# treatment: School_meal - indicator for participation in the school mean plan
# outcome: BMI 
nhanes_bmi = read.csv("nhanes_bmi.csv")[, -1]
z = nhanes_bmi$School_meal 
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

## simple regression analysis
# Without adjusting for any covariates / Fisher's / Lin's estimators
DiM = lm(y ~ z)
Fisher = lm(y ~ z + x)
Lin = lm(y ~ z + x + z*x)

res.regression = c(coef(DiM)[2], hccm(DiM)[2, 2]^0.5,
                   coef(Fisher)[2], hccm(Fisher)[2, 2]^0.5,
                   coef(Lin)[2], hccm(Lin)[2, 2]^0.5)
res.regression = matrix(res.regression, nrow=2, ncol=3)
rownames(res.regression) = c("est", "se")
colnames(res.regression) = c("naive", "fisher", "lin")
round(res.regression, 3)
#     naive fisher    lin
# est 0.534  0.061 -0.017
# se  0.225  0.227  0.226

# The naive difference in means differs greatly from the other methods
# due to the large imbalance in covariates.
# Fisher's estimator still gives a positive point estimator although it is not significant.
# Lin's estimator and the propensity score stratification estimators give qualitatively the same
# results.


## Propensity score stratification
# Based on propensity score stratification, we can calculate the point estimators &
# standard errors for different choice of K in {5,10,20,50,80}


## discretized by the quantiles of the pscore
pscore = glm(z ~ x, family=binomial)$fitted.values
n.strata = c(5, 10, 20, 50, 80)
strat.res = sapply(n.strata, FUN=function(nn){
  q.pscore = quantile(pscore, (1:(nn-1))/nn)
  ps.strata = cut(pscore, breaks = c(0, q.pscore, 1),
                  labels = 1:nn)
  Neyman_SRE(z, y, ps.strata)})

rownames(strat.res) = c("est", "se")
colnames(strat.res) = n.strata
round(strat.res, 3)
#          5     10     20     50     80
# est -0.116 -0.178 -0.200 -0.265 -0.204
# se   0.283  0.282  0.279  0.272     NA

# Increasing K from 5 to 50 reduces the standard error.
# However, we cannot go as extreme as K = 80 because the se is not well-defined in some strata
# with only one treated or one control unit.
# Estimators show a negative but insignificant effect of the meal program on BMI.

# The propensity score stratification estimators are stable across different choices of K.

pdf("pscore_stratified_hist.pdf", height=4, width=8)

pscore = glm(z ~ x, family=binomial)$fitted.values
par(mfrow = c(1, 3), mar = c(2, 0.1, 1, 0.1))
hist(pscore[z==1], freq=FALSE, col="grey", border=NA, xlab="", ylab="",
     yaxt="n", breaks=5, main="breaks=5", xlim=c(0, 1), ylim=c(0, 4.5))
hist(pscore[z==0], freq=FALSE, add=TRUE, breaks = 5)

hist(pscore[z==1], freq=FALSE, col="grey", border=NA, xlab="", ylab="",
     yaxt="n", breaks=10, main="breaks=10", xlim=c(0, 1), ylim=c(0, 4.5))
hist(pscore[z==0], freq=FALSE, add=TRUE, breaks = 10)

hist(pscore[z==1], freq=FALSE, col="grey", border=NA, xlab="", ylab="",
     yaxt="n", breaks=30, main="breaks=30", xlim=c(0, 1), ylim=c(0, 4.5))
hist(pscore[z==0], freq=FALSE, add=TRUE, breaks = 30)
dev.off()







