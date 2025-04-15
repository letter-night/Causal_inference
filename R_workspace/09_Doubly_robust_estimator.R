setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")
getwd()

rm(list=ls())

# Outcome regression, HT, Hajek, Doubly robust estimators for Tau.
# Default choice for the propensity score model: logistic model
# Default choice for the outcome model: linear model 

OS_est = function(z, y, x, out.family=gaussian,
                  truncps=c(0,1))
{
  ## fitted propensity score
  pscore = glm(z ~ x, family = binomial)$fitted.values
  pscore = pmax(truncps[1], pmin(truncps[2], pscore))
  
  ## fitted potential outcomes
  outcome1 = glm(y ~ x, weights = z, 
                 family = out.family)$fitted.values
  outcome0 = glm(y ~ x, weights = (1 - z),
                 family = out.family)$fitted.values
  
  ## outcome regression estimator
  ace.reg = mean(outcome1 - outcome0)
  
  ## IPW estimators
  y.treat = mean(z*y/pscore)
  y.control = mean((1-z)*y / (1-pscore))
  one.treat = mean(z/pscore)
  one.control = mean((1-z)/(1-pscore))
  ace.ipw0 = y.treat - y.control
  ace.ipw = y.treat/one.treat - y.control/one.control
  
  ## doubly robust estimator
  res1 = y - outcome1
  res0 = y - outcome0
  r.treat = mean(z*res1/pscore)
  r.control = mean((1-z)*res0 / (1-pscore))
  ace.dr = ace.reg + r.treat - r.control 
  
  return(c(ace.reg, ace.ipw0, ace.ipw, ace.dr))
}

# Bootstrap approximation to the variances based on resampling from {Zi, Xi, Yi}.
# Funtion returning point estimators & bootstrap standard errors.
OS_ATE = function(z, y, x, n.boot = 2*10^2,
                  out.family = gaussian, truncps = c(0, 1))
{
  point.est = OS_est(z, y, x, out.family, truncps)
  
  ## nonparametric bootstrap
  n = length(z)
  x = as.matrix(x)
  boot.est = replicate(n.boot, {
    id.boot = sample(1:n, n, replace=TRUE)
    OS_est(z[id.boot], y[id.boot], x[id.boot, ],
           out.family, truncps)
  })
  boot.se = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg", "HT", "Hajek", "DR")
  
  return(res)
}

library(parallel)
numCores = detectCores()

# set sample size n = 500 & generate 500 independent data sets according to the DGP
n.sim = 500
n = 500

## ------------------ Simulation ---------------- ##
# simulation to evaluate the finite-sample properties of the estimators
# 1. both the propensity score & outcome models are correct
# 2. the propensity score model is wrong but the outcome model is correct
# 3. the propensity score model is correct but the outcome model is wrong
# 4. both the propensity score & outcome models are wrong

# average bias, true standard error, average estimated standard error of estimators 

## data generating process
# 1. simulation with correct models
simu.11 = function(mc)
{
  x = matrix(rnorm(n*2), n, 2)
  x1 = cbind(1, x)
  beta.z = c(0, 1, 1)
  pscore = 1 / (1 + exp(- as.vector(x1%*%beta.z)))
  z = rbinom(n, 1, pscore)
  beta.y1 = c(1, 2, 1)
  beta.y0 = c(1, 2, 1)
  y1 = rnorm(n, x1%*%beta.y1)
  y0 = rnorm(n, x1%*%beta.y0)
  y = z*y1 + (1-z)*y0
  ce = OS_ATE(z, y, x)
  
  ## true ACE, est, se
  c(mean(y1-y0), ce[1, ], ce[2,])
}

# 2. simulation with an incorrect propensity score model
# modify propensity score model to be nonlinear
simu.01 = function(mc)
{
  x = matrix(rnorm(n*2), n, 2)
  x1 = cbind(1, x, exp(x))
  beta.z = c(-1, 0, 0, 1, -1)
  pscore = 1 / (1 + exp(- as.vector(x1%*%beta.z)))
  z = rbinom(n, 1, pscore)
  beta.y1 = c(1, 2, 1, 0, 0)
  beta.y0 = c(1, 1, 1, 0, 0)
  y1 = rnorm(n, x1%*%beta.y1)
  y0 = rnorm(n, x1%*%beta.y0)
  y = z*y1 + (1-z)*y0
  ce = OS_ATE(z, y, x)
  
  ## true ACE, est, se
  c(mean(y1-y0), ce[1, ], ce[2, ])
}

# 3. simulation with an incorrect outcome model
# modify outcome model to be nonlinear
simu.10 = function(mc)
{
  x = matrix(rnorm(n*2), n, 2)
  x1 = cbind(1, x, exp(x))
  beta.z = c(0, 1, 1, 0, 0)
  pscore = 1 / (1 + exp(- as.vector(x1%*%beta.z)))
  z = rbinom(n, 1, pscore)
  beta.y1 = c(1, 0, 0, 0.2, -0.1)
  beta.y0 = c(1, 0, 0, -0.2, 0.1)
  y1 = rnorm(n, x1%*%beta.y1)
  y0 = rnorm(n, x1%*%beta.y0)
  y = z*y1 + (1-z)*y0
  ce = OS_ATE(z, y, x)
  
  ## true ACE, est, se
  c(mean(y1-y0), ce[1, ], ce[2, ])
}

# 4. simulation without correct models
# modify both the propensity score and the outcome model
simu.00 = function(mc)
{
  x = matrix(rnorm(n*2), n, 2)
  x1 = cbind(1, x, exp(x))
  beta.z = c(-1, 0, 0, 1, -1)
  pscore = 1 / (1 + exp(- as.vector(x1%*%beta.z)))
  z = rbinom(n, 1, pscore)
  beta.y1 = c(1, 0, 0, 0.2, -0.1)
  beta.y0 = c(1, 0, 0, -0.2, 0.1)
  y1 = rnorm(n, x1%*%beta.y1)
  y0 = rnorm(n, x1%*%beta.y0)
  y = z*y1 + (1-z)*y0
  ce = OS_ATE(z, y, x)
  
  ## true ACE, est, se
  c(mean(y1-y0), ce[1, ], ce[2, ])
}

## simulations
#res = mclapply(1:n.sim, simu.11, mc.cores = numCores)
res = mclapply(1:n.sim, simu.11)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)


res = mclapply(1:n.sim, simu.01)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)


res = mclapply(1:n.sim, simu.10)
res = simplify2array(res)

bias = res[2:5, ] - 0.2*exp(1/2)
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)


res = mclapply(1:n.sim, simu.00)
res = simplify2array(res)

bias = res[2:5, ] - 0.2*exp(1/2)
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

# case 1: All estimators are nearly unbiased.
# The two weighting esetimators have larger variances.
#            reg   HT Hajek    DR
# ave.bias -0.01 0.01  0.03 -0.01
# true.se   0.11 0.30  0.27  0.13
# est.se    0.10 0.26  0.23  0.12

# case2: the two weighting estimators are severely biased due to the misspecification of
# the propensity score model.
# The outcome regression & doubly robust estimators are nearly unbiased.
#            reg    HT Hajek   DR
# ave.bias -0.01 -0.80 -0.76 0.00
# true.se   0.13  0.78  0.53 0.34
# est.se    0.13  0.55  0.39 0.20

# case3: the outcome regression estimator has a larger bias than the other three estimators
# dut eo the misspecification of the outcome model.
# The weighting & doubly robust estimators are nearly unbiased.
#            reg   HT Hajek   DR
# ave.bias -0.05 0.00  0.01 0.00
# true.se   0.11 0.16  0.15 0.14
# est.se    0.11 0.15  0.14 0.14

# case4: all estimators are biased because both the propensity score & outcome models
# are wrong. 
# The HT & doubly robust estimator has the largest bias.
# When both models are wrong, the doubly robust estimator appears to be doubly fragile.
#            reg   HT Hajek   DR
# ave.bias -0.07 0.10 -0.06 0.14
# true.se   0.13 0.26  0.19 0.28
# est.se    0.13 0.23  0.16 0.23

# The bootstrap standard errors are close to the true ones when the estimators are
# nearly unbiased for the true average causal effect.



rm(list=ls())
library(car)

## Data "nhanes_bmi" from the "ATE" package in R
# treatment: School_meal - indicator for participation in the school mean plan
# outcome: BMI 
nhanes_bmi = read.csv("nhanes_bmi.csv")[, -1]
z = nhanes_bmi$School_meal 
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

causaleffects = OS_ATE(z, y, x, n.boot=10^3)
round(causaleffects, 3)
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

#        reg     HT  Hajek     DR
# est -0.017 -1.516 -0.156 -0.019
# se   0.230  0.533  0.262  0.239

#            reg         HT      Hajek         DR
#[1,] -0.4673063 -2.5611292 -0.6687146 -0.4868155
#[2,]  0.4333985 -0.4714384  0.3573769  0.4482284

# The two weighting estimators are much larger than the other two estimators.

## checking the data
pscore = glm(z ~ x, family = binomial)$fitted.values
hist(pscore[z==1], col="grey", border=NA, freq=FALSE, 
     ylim=c(0, 4.5), breaks=30, main="",
     xlab=expression(hat(e)(X)), ylab="")

## truncated propensity score
# Truncating the estimated propensity score at [0.1, 0.9], we obtain the following
# estimators and bootstrap standard errors.
causaleffects = OS_ATE(z, y, x, n.boot=10^3,
                       truncps=c(0.1, 0.9))
round(causaleffects, 3)
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

#        reg     HT  Hajek     DR
# est -0.017 -0.713 -0.054 -0.043
# se   0.226  0.414  0.235  0.230

#             reg          HT      Hajek         DR
# [1,] -0.4606230 -1.52506527 -0.5139619 -0.4951546
# [2,]  0.4267152  0.09826943  0.4068225  0.4083924

# The Hajek estimator becomes much closer to the outcome regression & dr estimators,
# while Horvitz-Thompson estimator is still an outlier.

