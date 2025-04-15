
# Problem 5 #######################################################################################

## Simulation studies for Tau_T with either correct or incorrect propensity score or outcome models.


# Estimators for Tau_T: Outcome regression, HT, Hajek, Doubly robust estimators. ------------ #

# Default choice for the propensity score model: logistic model
# Default choice for the outcome model: linear model

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
  outcome0 = glm(y ~ x, weights = (1 - z), family = out.family)$fitted.values
  
  ## outcome regression estimator
  ace.reg = mean(y[z==1]) - mean(outcome0[z==1])
  
  ## propensity score weighting estimator
  ace.ipw0 = mean(y[z==1]) - mean(odds*(1-z)*y)*nn/nn1
  ace.ipw = mean(y[z==1]) - mean(odds*(1-z)*y) / mean(odds*(1-z))
  
  ## doubly robust estimator
  res0 = y - outcome0
  ace.dr = ace.reg - mean(odds*(1-z)*res0)*nn/nn1
  
  return (c(ace.reg, ace.ipw0, ace.ipw, ace.dr))
}


# Bootstrap approximation to the variances based on resampling from {Zi, Xi, Yi}. --------- #
# function returning point estimators & bootstrap standard errors.

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
  colnames(res) = c("reg", "HT", "Hajek", "DR")
  
  return(res)
}

# Simulation ------------------------------------------------------------------------------ #
# simulation to evaluate the finite-sample properties of the estimators
# 1. both the propensity score & outcome models are correct
# 2. the propensity score model is wrong but the outcome model is correct
# 3. the propensity score model is correct but the outcome model is wrong
# 4. both the propensity score & outcome models are wrong

# Report average bias, true standard error, average estimated standard error of estimators.

# set sample size n = 500 & generate 500 independent datasets according to the DGP.
n.sim = 500
n = 500

# data generating process ----------------------------------- ##
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
  
  ## true ACE, est, s.e.
  c(mean(y1-y0), ce[1, ], ce[2, ])
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
  
  ## true ACE, est, s.e.
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
  
  ## true ACE, est, s.e.
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
  
  ## true ACE, est, s.e.
  c(mean(y1-y0), ce[1, ], ce[2, ])
}


# simulation results ---------------------------------------- ##
library(parallel)
# 1. simulation with correct models
res = mclapply(1:n.sim, simu.11)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

# 2. simulation with an incorrect propensity score model
res = mclapply(1:n.sim, simu.01)
res = simplify2array(res)

bias = res[2:5, ] - 0
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

# 3. simulation with an incorrect outcome model
res = mclapply(1:n.sim, simu.10)
res = simplify2array(res)

bias = res[2:5, ] - 0.2*exp(1/2)
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

# 4. simulation without correct models
res = mclapply(1:n.sim, simu.00)
res = simplify2array(res)

bias = res[2:5, ] - 0.2*exp(1/2)
res = rbind(apply(bias, 1, mean),
            apply(bias, 1, sd),
            apply(res[6:9, ], 1, mean))
row.names(res) = c("ave.bias", "true.se", "est.se")
round(res, 2)

# case 1 ---------------------------------- ##
#           reg   HT Hajek   DR
# ave.bias 0.00 0.04  0.05 0.00
# true.se  0.12 0.54  0.35 0.14
# est.se   0.12 0.41  0.27 0.14

# All estimators are nearly unbiased. The two weighting estimators have larger variances.

# case2 ---------------------------------- ##
#           reg   HT Hajek   DR
# ave.bias 0.59 1.10  1.04 0.59
# true.se  0.14 0.21  0.21 0.16
# est.se   0.14 0.20  0.20 0.16

# The two weighting estimators are severely biased 
# due to the misspecification of the propensity score model.

# Even though outcome regression model is correct, DR estimator has a larger variance
# compared to the direct regression estimator. 

# case3 ---------------------------------- ##
#           reg   HT Hajek   DR
# ave.bias 0.09 0.12  0.12 0.12
# true.se  0.14 0.20  0.19 0.18
# est.se   0.13 0.19  0.17 0.17
# The outcome regression estimator has a large bias due to the misspecification of the outcome model.

# case4 ---------------------------------- ##
#           reg   HT Hajek   DR
# ave.bias 0.42 0.45  0.39 0.41
# true.se  0.14 0.15  0.16 0.16
# est.se   0.13 0.15  0.15 0.15

# All estimators are biased because both the propensity score & outcome models are wrong.
# When both models are wrong, the doubly robust estimator appears to have large bias.

# The bootstrap standard errors are close to the true ones when the estimators are
# nearly unbiased for the true average causal effect.



# The end of the Problem 5. solutions. #############################################################

