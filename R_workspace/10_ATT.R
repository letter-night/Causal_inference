
# Two outcome regression estimators, two IPW estimators, and the dr estimator for Tau_T


ATT.est = function(z, y, x, out.family = gaussian, Utruncps=1)
{
  ## sample size
  nn = length(z)
  nn1 = sum(z)
  
  ## fitted propensity score
  pscore = glm(z ~ x, family = binomial)$fitted.values
  pscore = pmin(Utruncps, pscore)
  odds = pscore/(1-pscore)
  
  ## fitted potential outcomes
  outcome0 = glm(y ~ x, weights = (1 - z),
                 family = out.family)$fitted.values
  
  ## outcome regression estimator
  ace.reg0 = lm(y ~ z + x)$coef[2]
  ace.reg = mean(y[z==1]) - mean(outcome0[z==1])
  
  ## propensity score weighting estimator
  ace.ipw0 = mean(y[z==1]) - mean(odds*(1-z)*y)*nn/nn1
  ace.ipw = mean(y[z==1]) - mean(odds*(1-z)*y) / mean(odds*(1-z))
  
  ## doubly robust estimator
  res0 = y - outcome0
  ace.dr = ace.reg - mean(odds*(1-z)*res0)*nn/nn1
  
  return(c(ace.reg0, ace.reg, ace.ipw0, ace.ipw, ace.dr))
}

# function implementing bootstrap variance estimators
OS_ATT = function(z, y, x, n.boot=10^2, out.family=gaussian, Utruncps=1)
{
  point.est = ATT.est(z, y, x, out.family, Utruncps)
  
  ## nonparametric bootstrap
  n = length(z)
  x = as.matrix(x)
  boot.est = replicate(n.boot, {
    id.boot = sample(1:n, n, replace=TRUE)
    ATT.est(z[id.boot], y[id.boot], x[id.boot, ], out.family, Utruncps)
  })
  
  boot.se = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg0", "reg", "HT", "Hajak", "DR")
  
  return(res)
}


## an application: the NHANES BMI dataset
## Data "nhanes_bmi" from the "ATE" package in R
# Estimate Tau_T

nhanes_bmi = read.csv("nhanes_bmi.csv")[, -1]
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

ATTs = OS_ATT(z, y, x, n.boot=10^3)
round(ATTs, 3)

#      reg0    reg     HT  Hajak     DR
# est 0.061 -0.351 -1.992 -0.351 -0.187
# se  0.217  0.256  0.711  0.325  0.276

## truncated propensity score
# truncating the estimated propensity scores from the above at 0.9
ATTs = OS_ATT(z, y, x, n.boot=10^3, Utruncps=0.9)
round(ATTs, 3)

#      reg0    reg     HT  Hajak     DR
# est 0.061 -0.351 -0.597 -0.192 -0.230
# se  0.227  0.254  0.572  0.304  0.275

# The HT estimator is sensitive to the truncation.
# The regression estimator in exp13.1 is quite different from other estimators.
# It imposes an unnecessary assumption that the regression functions in the treatment & control
# group share the same coefficient of X.
# The regression estimator in exp13.2 is much closer to the Hajek and dr estimators.
# The estimates are slightly different from those in Section 12.3.3, suggesting 
# the existence of treatment effect heterogeneity across Tau_T and Tau.


## Simulation
## ------------------ Simulation ---------------- ##
# simulation to evaluate the finite-sample properties of the estimators
# 1. both the propensity score & outcome models are correct
# 2. the propensity score model is wrong but the outcome model is correct
# 3. the propensity score model is correct but the outcome model is wrong
# 4. both the propensity score & outcome models are wrong

# average bias, true standard error, average estimated standard error of estimators 

# set sample size n = 500 & generate 500 independent data sets according to the DGP
n.sim = 500
n = 500

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
  ce = OS_ATT(z, y, x)
  
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
  ce = OS_ATT(z, y, x)
  
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
  ce = OS_ATT(z, y, x)
  
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
  ce = OS_ATT(z, y, x)
  
  ## true ACE, est, se
  c(mean(y1-y0), ce[1, ], ce[2, ])
}

## simulations
res = mclapply(1:n.sim, simu.11)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

##
res = mclapply(1:n.sim, simu.01)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

##
res = mclapply(1:n.sim, simu.10)
res = simplify2array(res)

bias = res[2:5, ] - 0.2*exp(1/2)
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

##
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
#           reg0  reg    HT Hajak
# ave.bias -0.01 0.00 -0.01  0.02
# true.se   0.10 0.11  0.54  0.35
# est.se    0.00 0.10  0.12  0.44

# case2: the two weighting estimators are severely biased due to the misspecification of
# the propensity score model.
# The outcome regression & doubly robust estimators are nearly unbiased.
#          reg0  reg   HT Hajak
# ave.bias 0.22 0.60 1.10  1.04
# true.se  0.12 0.14 0.21  0.21
# est.se   0.60 0.13 0.14  0.21

# case3: the outcome regression estimator has a larger bias than the other three estimators
# dut eo the misspecification of the outcome model.
# The weighting & doubly robust estimators are nearly unbiased.
#           reg0  reg   HT Hajak
# ave.bias -0.04 0.08 0.12  0.11
# true.se   0.11 0.13 0.18  0.18
# est.se    0.45 0.11 0.13  0.18

# case4: all estimators are biased because both the propensity score & outcome models
# are wrong. 
# The HT & doubly robust estimator has the largest bias.
# When both models are wrong, the doubly robust estimator appears to be doubly fragile.
#          reg0  reg   HT Hajak
# ave.bias 0.15 0.44 0.46  0.40
# true.se  0.12 0.14 0.15  0.15
# est.se   0.75 0.12 0.13  0.15

# The bootstrap standard errors are close to the true ones when the estimators are
# nearly unbiased for the true average causal effect.
