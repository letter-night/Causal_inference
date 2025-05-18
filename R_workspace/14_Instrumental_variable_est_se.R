setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")
getwd()

rm(list=ls())

## for parallel computing
library("car")
library(parallel)
numCores = detectCores()

## IV point estimator ------------------------------------------------------------------ #
IV_Wald = function(Z, D, Y)
{
  tau_D = mean(D[Z==1]) - mean(D[Z==0])
  tau_Y = mean(Y[Z==1]) - mean(Y[Z==0])
  CACE = tau_Y / tau_D
  
  c(tau_D, tau_Y, CACE)
}

## IV se via the delta method
IV_Wald_delta = function(Z, D, Y)
{
  est = IV_Wald(Z, D, Y)
  AdjustedY = Y - D*est[3]
  VarAdj = var(AdjustedY[Z==1])/sum(Z) + var(AdjustedY[Z==0])/sum(1-Z)
  
  c(est[3], sqrt(VarAdj)/abs(est[1]))
}

## IV se via the bootstrap
IV_Wald_bootstrap = function(Z, D, Y, n.boot=200)
{
  est = IV_Wald(Z, D, Y)
  
  CACEBoot = replicate(n.boot, {
    id.boot = sample(1:length(Z), replace=TRUE)
    IV_Wald(Z[id.boot], D[id.boot], Y[id.boot])[3]
  })
  
  c(est[3], sd(CACEBoot))
}


## covariate adjustment in IV analysis ------------------------------------------------- #

IV_Lin = function(Z, D, Y, X)
{
  X = scale(as.matrix(X))
  tau_D = lm(D ~ Z + X + Z*X)$coef[2]
  tau_Y = lm(Y ~ Z + X + Z*X)$coef[2]
  CACE = tau_Y / tau_D
  
  c(tau_D, tau_Y, CACE)
}

## IV_adj se via the bootstrap
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


## One sided noncompliance: D(0) = 0 --------------------------------------------------- #

# Efficacy of a job training intervention on unemployed workers 
# randomized experiment
# treat: indicator for whether a participant was randomly selected for training program.
# comply: indicator for whether a participant actually participated in the program.
# job_seek: outcome for measuring the level of job-search self-efficacy (1 ~ 5)
# Covariates: sex, age, marital, nonwhite, educ, income

jobsdata = read.csv("./data/jobsdata.csv")
Z = jobsdata$treat
D = jobsdata$comply
Y = jobsdata$job_seek

getX = lm(treat ~ sex + age + marital + nonwhite + educ + income,
          data = jobsdata)

X = model.matrix(getX)[, -1]

## Estimate Tau_c_hat 
## without covariates
res = rbind(IV_Wald_delta(Z, D, Y), IV_Wald_bootstrap(Z, D, Y, n.boot=10^3))

## with covariates
res = rbind(res, IV_Lin_bootstrap(Z, D, Y, X, n.boot = 10^3))
res = cbind(res, res[, 1] - 1.96*res[, 2], 
            res[, 1] + 1.96*res[, 2])

row.names(res) = c("delta", "bootstrap", "with covariates")
colnames(res) = c("est", "se", "lower CI", "upper CI")
round(res, 3)
#                   est    se lower CI upper CI
# delta           0.109 0.081   -0.050    0.268
# bootstrap       0.109 0.082   -0.052    0.269
# with covariates 0.118 0.082   -0.044    0.279


## FAR CI --------------------------------------------------------------------------------- #
## Fieller-Anderson-Rubin confidence interval
## without covariates
FARci = function(Z, D, Y, Lower, Upper, grid)
{
  CIrange = seq(Lower, Upper, grid)
  Pvalue  = sapply(CIrange, function(t){
    Y_t       = Y - t*D
    Tauadj    = mean(Y_t[Z==1]) - mean(Y_t[Z==0])
    VarAdj    = var(Y_t[Z==1])/sum(Z) + 
      var(Y_t[Z==0])/sum(1 - Z)
    Tstat     = Tauadj/sqrt(VarAdj)
    (1 - pnorm(abs(Tstat)))*2
  })
  
  return(list(CIrange = CIrange, Pvalue  = Pvalue))
}

linestimator = function(Z, Y, X){
  ## standardize X
  X      = scale(X)
  n      = dim(X)[1]
  p      = dim(X)[2]
  
  ## fully interacted OLS
  linreg = lm(Y ~ Z*X)
  est    = coef(linreg)[2]
  vehw   = hccm(linreg)[2, 2]
  
  ## super population correction
  inter  = coef(linreg)[(p+3):(2*p+2)]
  vsuper = vehw + sum(inter*(cov(X)%*%inter))/n
  
  c(est, sqrt(vehw), sqrt(vsuper))
}

FARciX = function(Z, D, Y, X, Lower, Upper, grid)
{
  CIrange = seq(Lower, Upper, grid)
  X       = scale(X)
  Pvalue  = sapply(CIrange, function(t){
    Y_t       = Y - t*D
    linest    = linestimator(Z, Y_t, X)
    Tstat     = linest[1]/linest[3]
    (1 - pnorm(abs(Tstat)))*2
  })
  
  return(list(CIrange = CIrange, Pvalue  = Pvalue))
}


pdf("farci_jobs.pdf", height = 8.5, width = 8.5)
par(mfrow = c(2, 1), mai = c(0.8, 0.8, 0.3, 0.3))
FARplot = FARci(Z, D, Y, Lower = -0.2, Upper = 0.4, grid = 0.001)
plot(FARplot$Pvalue ~ FARplot$CIrange,
     type = "l", 
     xlab = expression(tau[c]), 
     ylab = "p-value",
     main = "without covariates")
abline(h = 0.05, lty = 2)  
ARCI = range(FARplot$CIrange[FARplot$Pvalue >= 0.05])

FARplot = FARciX(Z, D, Y, X, Lower = -0.2, Upper = 0.4, grid = 0.001)
plot(FARplot$Pvalue ~ FARplot$CIrange,
     type = "l", 
     xlab = expression(tau[c]), 
     ylab = "p-value",
     main = "with covariates")
abline(h = 0.05, lty = 2)  
ARCI.x = range(FARplot$CIrange[FARplot$Pvalue >= 0.05])

ARCIs = rbind(ARCI, ARCI.x)
row.names(ARCIs) = c("without covariates", "with covariates")
colnames(ARCIs) = c("lower CI", "upper CI")
round(ARCIs, 3)
dev.off()

#                    lower CI upper CI
# without covariates   -0.050    0.267
# with covariates      -0.047    0.282









