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

# Compute the IPW estimators & bootstrap standard errors.

## functions for IPW estimation
ipw.est = function(z, y, x, truncps = c(0, 1))
{
  ## fitted propensity score
  pscore = glm(z ~ x, family = binomial)$fitted.values
  pscore = pmax(truncps[1], pmin(truncps[2], pscore))
  
  ace.ipw0 = mean(z*y/pscore - (1-z)*y/(1-pscore))
  ace.ipw = mean(z*y/pscore) / mean(z/pscore) - mean((1-z)*y/(1-pscore)) / mean((1-z)/(1-pscore))
  
  return(c(ace.ipw0, ace.ipw))
}

ipw.boot = function(z, y, x, n.boot = 500, truncps = c(0, 1))
{
  point.est = ipw.est(z, y, x, truncps)
  
  ## nonparametric bootstrap
  n.sample = length(z)
  x = as.matrix(x)
  boot.est = replicate(n.boot, {
    id.boot = sample(1:n.sample, n.sample, replace=TRUE)
    ipw.est(z[id.boot], y[id.boot], x[id.boot,], truncps)
  })
  boot.se = apply(boot.est, 1, sd)
  
  res = cbind(point.est, boot.se)
  colnames(res) = c("est", "se")
  rownames(res) = c("HT", "Hajek")
  
  return(res)
}


# IPW estimators based on different truncations of the estimated propensity scores.
trunc.list = list(trunc0 = c(0, 1),
                  trunc.01 = c(0.01, 0.99),
                  trunc.05 = c(0.05, 0.95),
                  trunc.1 = c(0.1, 0.9))
trunc.est = lapply(trunc.list,
                   function(t){
                     est = ipw.boot(z, y, x, truncps = t)
                     round(est, 3)
                   })

trunc.est
# $trunc0
#          est    se
# HT    -1.516 0.479
# Hajek -0.156 0.248

# $trunc.01
#          est    se
# HT    -1.516 0.488
# Hajek -0.156 0.256

# $trunc.05
#          est    se
# HT    -1.499 0.455
# Hajek -0.152 0.252

# $trunc.1
#          est    se
# HT    -0.713 0.424
# Hajek -0.054 0.233

# HT estimator gives results far away from all other estimators we discussed so far.
# The point estimates seem too large and they are negatively significant unless 
# we truncate the estimated propensity scores at (0.1, 0.9).
# This shows the instability of HT estimator.

