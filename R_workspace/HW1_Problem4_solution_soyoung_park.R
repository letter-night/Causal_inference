
# Problem 4 (a), (b), (c) ###########################################################################

## Estimating the treatment effect using the outcome regression / IPW / doubly-robust estimator.

# cps1re74: observational study examining the effect of a job training program.

# z: binary treatment vector indicating whether a unit was randomly assigned to the job program or not.
# y: outcome vector representing a unit's earnings in 1978.
dat = read.table("./data/cps1re74.csv", header=TRUE)

# Binary treatment
z = dat$treat

# Outcome
y = dat$re78

# Covariate
dat$u74 <- as.numeric(dat$re74==0) # unemployed
dat$u75 <- as.numeric(dat$re75==0) # unemployed
x = as.matrix(dat[, c("age", "educ", "black", "hispan", "married",
                      "nodegree", "re74", "re75", "u74", "u75")]) 

x = scale(x)

# Estimators for Tau: Outcome regression, HT, Hajek, Doubly robust estimators. ------------ #

# Default choice for the propensity score model: logistic model
# Default choice for the outcome model: linear model

OS_est = function(z, y, x, out.family=gaussian, truncps=c(0, 1))
{
  ## fitted propensity score
  pscore = glm(z ~ x, family = binomial)$fitted.values
  pscore = pmax(truncps[1], pmin(truncps[2], pscore))
  
  ## fitted potential outcomes
  outcome1 = glm(y ~ x, weights = z, family = out.family)$fitted.values
  outcome0 = glm(y ~ x, weights = (1 - z), family = out.family)$fitted.values
  
  ## outcome regression estimator
  ace.reg = mean(outcome1 - outcome0)
  
  ## IPW estimators
  y.treat = mean(z*y/pscore)
  y.control = mean((1-z)*y/(1-pscore))
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

# Bootstrap approximation to the variances based on resampling from {Zi, Xi, Yi}. --------- #
# function returning point estimators & bootstrap standard errors.

OS_ATE = function(z, y, x, n.boot=2*10^2, out.family=gaussian, truncps=c(0, 1))
{
  point.est = OS_est(z, y, x, out.family, truncps)
  
  ## nonparametric bootstrap
  n = length(z)
  x = as.matrix(x)
  boot.est = replicate(n.boot, {
    id.boot = sample(1:n, n, replace=TRUE)
    OS_est(z[id.boot], y[id.boot], x[id.boot, ], out.family, truncps)
  })
  boot.se = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg", "HT", "Hajek", "DR")
  
  return (res)
}

# Estimate the treatment effects using outcome regression / IPW / doubly-robust estimator -- #
causaleffects = OS_ATE(z, y, x, n.boot=10^3)
round(causaleffects, 3)
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

#           reg         HT   Hajek        DR
# est -4265.801 -10557.510 -6950.4 -4207.565
# se   3509.803   1682.104  1709.6  3568.536

#             reg         HT      Hajek         DR
# [1,] -11145.014 -13854.434 -10301.217 -11201.897
# [2,]   2613.413  -7260.586  -3599.583   2786.766

# The two weighting estimators are much larger than the other two estimators.
# Hajek shows more stable estimates compared to the HT, due to the normalization.


# truncated propensity score - truncating the estimated propensity score at [0.1, 0.9]
causaleffects = OS_ATE(z, y, x, n.boot=10^3, truncps=c(0.01, 0.99))
round(causaleffects, 3)
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

#           reg         HT     Hajek        DR
# est -4265.801 -15972.091 -8141.602 -4258.749
# se   3422.331     98.796   643.089  3424.690

#             reg        HT     Hajek         DR
# [1,] -10973.569 -16165.73 -9402.057 -10971.141
# [2,]   2441.968 -15778.45 -6881.147   2453.643


# The end of the Problem 4. (a), (b), (c) solutions. ###############################################

