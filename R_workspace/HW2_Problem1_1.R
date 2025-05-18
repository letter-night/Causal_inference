rm(list=ls())
library("car")
library("Matching")

## Problem 1. NHANES data to study whether participation in school meal programs led to an increase in BMI for school children.

## Data "nhanes_bmi"
# treatment: School_meal - indicator for participation in the school meal plan.
# outcome: BMI
nhanes_bmi = read.csv("./data/nhanes_bmi.csv")[, -1]
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

## (a) Fit a logistic regression for Pr(Z = 1 | X) using all covariates. & 
##     Plot the estimated propensity scores for treated and control units. 

pdf("pscore_hist.pdf", height=4, width=8)

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


## (b) Nearest-neighbor matching without replacement using the estimated propensity scores.
##     1) k = 1, 2 matches per treated unit.
##     2) assess covariate balance before and after matching & graphical diagnostics. 
##     3) alternative approaches for improving covariate balance.

# matching (k = 1) without replacement using the estimated propensity scores ------------ #
matchest = Match(Y = y, Tr = z, X = pscore, replace = FALSE)
summary(matchest)
# Estimate...  0.50175 
# SE.........  0.2382 
# T-stat.....  2.1064 
# p.val......  0.035167 

# Original number of observations..............  2330 
# Original number of treated obs...............  1284 
# Matched number of observations...............  1046 
# Matched number of observations  (unweighted).  1046 

# balance checking before and after matching -------------------------------------------- #

# Use simple OLS to check covariate balance.
# Before matching, the covariates are highly imbalanced,
# signified by many stars associated with the coefficients.

# However, after matching, the covariates are well-balanced, 
# signified by the absence of stars for all coefficients.

lm.before = lm(z ~ x)
summary(lm.before)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.5510730  0.0088522  62.253  < 2e-16 ***
# xage          0.0384998  0.0094972   4.054 5.21e-05 ***
# xChildSex    -0.0009824  0.0088683  -0.111  0.91181    
# xblack        0.0871794  0.0096867   9.000  < 2e-16 ***
# xmexam        0.0933521  0.0098366   9.490  < 2e-16 ***
# xpir200_plus -0.1498185  0.0103655 -14.454  < 2e-16 ***
# xWIC          0.0169171  0.0096372   1.755  0.07932 .  
# xFood_Stamp   0.0939160  0.0104951   8.949  < 2e-16 ***
# xfsdchbi      0.0287162  0.0093979   3.056  0.00227 ** 
# xAnyIns      -0.0026647  0.0092550  -0.288  0.77343    
# xRefSex       0.0015053  0.0093552   0.161  0.87218    
# xRefAge       0.0018197  0.0096494   0.189  0.85044    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lm.after = lm(z ~ x, subset = c(matchest$index.treated, 
                                matchest$index.control))
summary(lm.after)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.5101639  0.0093286  54.688  < 2e-16 ***
# xage          0.0417804  0.0099332   4.206 2.71e-05 ***
# xChildSex     0.0055282  0.0093447   0.592   0.5542    
# xblack        0.0848452  0.0102128   8.308  < 2e-16 ***
# xmexam        0.0961833  0.0104150   9.235  < 2e-16 ***
# xpir200_plus -0.1501435  0.0109386 -13.726  < 2e-16 ***
# xWIC          0.0205085  0.0102636   1.998   0.0458 *  
# xFood_Stamp   0.1015770  0.0112791   9.006  < 2e-16 ***
# xfsdchbi      0.0237824  0.0100948   2.356   0.0186 *  
# xAnyIns       0.0026992  0.0099108   0.272   0.7854    
# xRefSex      -0.0008433  0.0098706  -0.085   0.9319    
# xRefAge       0.0058418  0.0102744   0.569   0.5697    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# matching (k = 2) without replacement using the estimated propensity scores ------------ #
matchest = Match(Y = y, Tr = z, X = pscore, M = 2, replace = FALSE)
summary(matchest)
# Estimate...  0.43054 
# SE.........  0.32457 
# T-stat.....  1.3265 
# p.val......  0.18467 

# Original number of observations..............  2330 
# Original number of treated obs...............  1284 
# Matched number of observations...............  523 
# Matched number of observations  (unweighted).  1046 

# balance checking before and after matching --------------------------------------------- #
lm.after = lm(z ~ x, subset = c(matchest$index.treated,
                                matchest$index.control))
summary(lm.after)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.507301   0.009273  54.709  < 2e-16 ***
# xage          0.037656   0.009790   3.846 0.000123 ***
# xChildSex     0.004211   0.009282   0.454 0.650138    
# xblack        0.085901   0.010243   8.386  < 2e-16 ***
# xmexam        0.102398   0.010405   9.841  < 2e-16 ***
# xpir200_plus -0.148334   0.010836 -13.689  < 2e-16 ***
# xWIC          0.027055   0.010082   2.684 0.007343 ** 
# xFood_Stamp   0.110747   0.011173   9.912  < 2e-16 ***
# xfsdchbi      0.011571   0.010081   1.148 0.251189    
# xAnyIns      -0.012739   0.009672  -1.317 0.187970    
# xRefSex       0.009751   0.009736   1.002 0.316672    
# xRefAge       0.009755   0.010367   0.941 0.346814    
# ---

matchest = Match(Y = y, Tr = z, X = pscore, M = 3, replace = TRUE, BiasAdjust = TRUE)
summary(matchest)

lm.after = lm(z ~ x, subset = c(matchest$index.treated,
                                matchest$index.control))
summary(lm.after)

############################################################################################
# MatchIt 














