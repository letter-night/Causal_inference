setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")
getwd()

## for parallel computing
library("car")
library(parallel)
numCores = detectCores()

## IV point estimator ------------------------------------------------------------------------- #
IV_Wald = function(Z, D, Y)
{
  tau_D = mean(D[Z==1]) - mean(D[Z==0])
  tau_Y = mean(Y[Z==1]) - mean(Y[Z==0])
  CACE = tau_Y / tau_D
  
  c(tau_D, tau_Y, CACE)
}

## IV se via the delta method ----------------------------------------------------------------- #
IV_Wald_delta = function(Z, D, Y)
{
  est = IV_Wald(Z, D, Y)
  AdjustedY = Y - D*est[3]
  VarAdj = var(AdjustedY[Z==1])/sum(Z) + var(AdjustedY[Z==0])/sum(1-Z)
  
  c(est[3], sqrt(VarAdj)/abs(est[1]))
}

## IV se via the bootstrap -------------------------------------------------------------------- #
IV_Wald_bootstrap = function(Z, D, Y, n.boot=200)
{
  est = IV_Wald(Z, D, Y)
  
  CACEBoot = replicate(n.boot, {
    id.boot = sample(1:length(Z), replace=TRUE)
    IV_Wald(Z[id.boot], D[id.boot], Y[id.boot])[3]
  })
  
  c(est[3], sd(CACEBoot))
}

## covariate adjustment in IV analysis -------------------------------------------------------- #
IV_Lin = function(Z, D, Y, X)
{
  X = scale(as.matrix(X))
  tau_D = lm(D ~ Z + X + Z*X)$coef[2]
  tau_Y = lm(Y ~ Z + X + Z*X)$coef[2]
  CACE = tau_Y / tau_D
  
  c(tau_D, tau_Y, CACE)
}

## IV_adj se via the bootsrap ----------------------------------------------------------------- #
IV_Lin_bootstrap = function(Z, D, Y, X, n.boot=200)
{
  X = scale(as.matrix(X))
  est = IV_Lin(Z, D, Y, X)
  CACEboot = replicate(n.boot, {
    id.boot = sample(1:length(Z), replace=TRUE)
    IV_Lin(Z[id.boot], D[id.boot], Y[id.boot], X[id.boot])[3]
  })
  
  c(est[3], sd(CACEboot))
}

#################################################################################################
## Examine whether it is better for patients to be treated in a large or small volume hospital.

# 158 cancer patients - 79 diagnosed at large hospital / 79 diagnosed at small hospital.

# (IV) z: indicator of whether a patient was diagnosed at a large volume hospital 
# (treatment) d: indicator of whether a patient received treatment at a large volume hospital.
# (outcome) y: whether the patient survived longer than 1 year after the diagnosis.
# (covariates) x: information about age, rural area, male.

# assuming the IV is randomly assigned conditional on observed covariates.

karolinska = read.table("./data/karolinska.txt", header = TRUE)
z = karolinska$hvdiag
d = karolinska$hvtreat
y = 1 - (karolinska$year.survival == 1)
x = as.matrix(karolinska[, c(3, 4, 5)])

getX = lm(hvtreat ~ age + rural + male, data = karolinska)
X = model.matrix(getX)[, -1]

## Estimate Tau_c_hat
## without covariates
res = rbind(IV_Wald_delta(z, d, y), IV_Wald_bootstrap(z, d, y, n.boot = 10^3))

## with covariates
res = rbind(res, IV_Lin_bootstrap(z, d, y, X, n.boot = 10^3))
res = cbind(res, res[, 1] - 1.96*res[, 2], res[, 1] + 1.96*res[, 2])

row.names(res) = c("delta", "bootstrap", "with covariates")
colnames(res) = c("est", "se", "lower CI", "upper CI")
round(res, 3)

#                    est    se lower CI upper CI
# delta           -0.023 0.139   -0.295    0.249
# bootstrap       -0.023 0.145   -0.306    0.261
# with covariates  0.065 0.156   -0.242    0.371

# Receiving treatment at a large volume hospital positively affected the outcome of whether
# the patient survived longer than 1 year after the diagnosis.

#################################################################################################

