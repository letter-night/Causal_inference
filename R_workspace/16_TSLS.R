rm(list=ls())

library("car")

## Card Data
card.data = read.csv("./data/card1995.csv")

Y = card.data[, "lwage"]
D = card.data[, "educ"]
Z = card.data[, "nearc4"]
X = card.data[, c("exper", "expersq", "black", "south",
                  "smsa", "reg661", "reg662", "reg663",
                  "reg664", "reg665", "reg666",
                  "reg667", "reg668", "smsa66")]
X = as.matrix(X)

## TSLS ------------------------------------------------------------------------------- #
# Based on TSLS, obtain the point estimator and 95% confidence interval.
Dhat = lm(D ~ Z + X)$fitted.values
tslsreg = lm(Y ~ Dhat + X)
tslsest = coef(tslsreg)[2]

## correct se by changing the residuals
res.correct = Y - cbind(1, D, X)%*%coef(tslsreg)
tslsreg$residuals = as.vector(res.correct)
tslsse = sqrt(hccm(tslsreg, type = "hc0")[2, 2])

res = c(tslsest, tslsest - 1.96*tslsse, tslsest + 1.96*tslsse)

names(res) = c("TSLS", "lower CI", "upper CI")
round(res, 3)
#  TSLS lower CI upper CI 
# 0.132    0.026    0.237 


## FAR confidence set ------------------------------------------------------------------- #
# Using FAR confidence set, compute p-value as a function of b.
BetaAR   = seq(-0.1, 0.4, 0.001)
PvalueAR = sapply(BetaAR, function(b){
  Y_b   = Y - b*D
  ARreg = lm(Y_b ~ Z + X)
  coefZ = coef(ARreg)[2]
  seZ   = sqrt(hccm(ARreg)[2, 2])
  Tstat = coefZ/seZ
  (1 - pnorm(abs(Tstat)))*2
})

# Figure shows the p-values for a sequence of tests for the coefficient of D 
pdf("FAR_Carddata.pdf", height = 5, width = 8.5)
plot(PvalueAR ~ BetaAR, type = "l",
     xlab = "coefficient of D",
     ylab = "p-value",
     main = "Fieller-Anderson-Rubin interval based on Card's data")
point.est = BetaAR[which.max(PvalueAR)]
abline(h = 0.05, lty = 2, col = "grey")
abline(v = point.est, lty = 2, col = "grey")
ARCI = range(BetaAR[PvalueAR >= 0.05])
abline(v = ARCI[1], lty = 2, col = "grey")
abline(v = ARCI[2], lty = 2, col = "grey")

## FAR results
# point estimate as the value of b with the largest p-value .
# confidence interval as the region of b with p-values larger than 0.05.
FARres = c(point.est, ARCI)
names(FARres) = c("FAR est", "lower CI", "upper CI")
round(FARres, 3)
dev.off()
# FAR est lower CI upper CI 
#   0.132    0.028    0.282 

# Comparing the TSLS and FAR methods, the lower confidence limits are very close
# but the upper confidence limits are slightly different
# due to the possibly heavy right tail of the distribution of the TSLS estimator.
# Overall, the TSLS and FAR methods give similar results in this example
# because the IV is not weak.
