getwd()
setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")

rm(list=ls())
library("car")
library("Matching")

## experimental data ####################################################################
data("lalonde")
head(lalonde)
#   age educ black hisp married nodegr re74 re75     re78 u74 u75 treat
# 1  37   11     1    0       1      1    0    0  9930.05   1   1     1
# 2  22    9     0    1       0      1    0    0  3595.89   1   1     1
# 3  30   12     1    0       0      0    0    0 24909.50   1   1     1
# 4  27   11     1    0       0      1    0    0  7506.15   1   1     1
# 5  33    8     1    0       0      1    0    0   289.79   1   1     1
# 6  22    9     1    0       0      1    0    0  4056.49   1   1     1

y = lalonde$re78
z = lalonde$treat
x = as.matrix(lalonde[, c("age", "educ", "black", "hisp", "married", "nodegr", 
                          "re74", "re75")])

## analysis the randomized experiment -------------------------------------------------- #
neymanols = lm(y ~ z)
fisherols = lm(y ~ z + x)
xc = scale(x)
linols = lm(y ~ z*xc)

resols = c(neymanols$coef[2],
           fisherols$coef[2],
           linols$coef[2],
           sqrt(hccm(neymanols, type = "hc2")[2, 2]),
           sqrt(hccm(fisherols, type = "hc2")[2, 2]),
           sqrt(hccm(linols, type = "hc2")[2, 2]))
resols = matrix(resols, 3, 2)
rownames(resols) = c("neyman", "fisher", "lin")
colnames(resols) = c("est", "se")
resols
#             est       se
# neyman 1794.343 670.9967
# fisher 1676.343 677.0493
# lin    1621.584 694.7217

# All regression estimators show positive significant results on the job training program.
# We can analyze the data as if it is an observational study based on 1-1 matching, 
# yielding the following results:

## analysis as if it is an observational study ----------------------------------------- #
matchest.adj = Match(Y = y, Tr = z, X = x, BiasAdjust = TRUE)
summary(matchest.adj)
# Estimate...  2119.7 
# AI SE......  876.42 
# T-stat.....  2.4185 
# p.val......  0.015583 

# Original number of observations..............  445 
# Original number of treated obs...............  185 
# Matched number of observations...............  185 
# Matched number of observations  (unweighted).  268 

# Both the point estimator and standard error increase, 
# but qualitatively, the conclusion remains the same.


## observational data #####################################################################
dat <- read.table("data/cps1re74.csv", header=TRUE)
# unemployed
dat$u74 <- as.numeric(data$re74==0)
dat$u75 <- as.numeric(data$re75==0)

head(dat)
#   treat age educ black hispan married nodegree re74 re75       re78 u74 u75
# 1     1  37   11     1      0       1        1    0    0  9930.0460   1   1
# 2     1  22    9     0      1       0        1    0    0  3595.8940   1   1
# 3     1  30   12     1      0       0        0    0    0 24909.4500   1   1
# 4     1  27   11     1      0       0        1    0    0  7506.1460   1   1
# 5     1  33    8     1      0       0        1    0    0   289.7899   1   1
# 6     1  22    9     1      0       0        1    0    0  4056.4940   1   1

y = dat$re78
z = dat$treat
x = as.matrix(dat[, c("age", "educ", "black", "hispan", "married", "nodegree", 
                      "re74", "re75", "u74", "u75")])

## analyze as if it is from a randomized experiment ------------------------------------- #
neymanols = lm(y ~ z)
fisherols = lm(y ~ z + x)
xc = scale(x)
linols = lm(y ~ z*xc)

resols = c(neymanols$coef[2],
           fisherols$coef[2],
           linols$coef[2],
           sqrt(hccm(neymanols, type = "hc2")[2, 2]),
           sqrt(hccm(fisherols, type = "hc2")[2, 2]),
           sqrt(hccm(linols, type = "hc2")[2, 2]))
resols = matrix(resols, 3, 2)
rownames(resols) = c("neyman", "fisher", "lin")
colnames(resols) = c("est", "se")
resols
#              est        se
# neyman -8506.495  583.4426
# fisher  1067.546  628.4389
# lin    -4265.801 3211.7718

# If we use simple OLS estimators, the results are far from the experimental benchmark
# and are sensitive to the specification of the regression.

# However, if we use 1-1 matching, 
# the results almost recover those based on the experimental data:

## analyze the observational study ------------------------------------------------------- #
matchest = Match(Y = y, Tr = z, X = x, BiasAdjust = TRUE)
summary(matchest)
# Estimate...  1747.8 
# AI SE......  916.59 
# T-stat.....  1.9068 
# p.val......  0.056543 

# Original number of observations..............  16177 
# Original number of treated obs...............  185 
# Matched number of observations...............  185 
# Matched number of observations  (unweighted).  248 


## matched pairs analysis

# Ignoring the ties in the matched data, we can also use the matched-pairs analysis,
# which again yields results similar to those based on the experimental data:
diff = y[matchest$index.treated] - y[matchest$index.control]
round(summary(lm(diff ~ 1))$coef[1, ], 2)
# Estimate Std. Error    t value   Pr(>|t|) 
#  1581.44     558.55       2.83       0.01 


## balance checking before and after matching ###############################################

# We can use simple OLS to check covariate balance.
# Before matching, the covariates are highly imbalanced, 
# signified by many stars associated with the coefficients.

# However, after matching, the covariates are well-balanced, 
# signified by the absence of stars for all coefficients.

lm.before = lm(z ~ x)
summary(lm.before)
# Call:
# lm(formula = z ~ x)

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.18508 -0.01057  0.00303  0.01018  1.01355 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.404e-03  6.326e-03   0.222   0.8243    
# xage        -4.043e-04  8.512e-05  -4.750 2.05e-06 ***
# xeduc        3.220e-04  4.073e-04   0.790   0.4293    
# xblack       1.070e-01  2.902e-03  36.871  < 2e-16 ***
# xhispan      6.377e-03  3.103e-03   2.055   0.0399 *  
# xmarried    -1.525e-02  2.023e-03  -7.537 5.06e-14 ***
# xnodegree    1.345e-02  2.523e-03   5.331 9.89e-08 ***
# xre74        7.601e-07  1.806e-07   4.208 2.59e-05 ***
# xre75       -1.231e-07  1.829e-07  -0.673   0.5011    
# xu74         4.224e-02  3.271e-03  12.914  < 2e-16 ***
# xu75         2.424e-02  3.399e-03   7.133 1.02e-12 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.09935 on 16166 degrees of freedom
# Multiple R-squared:  0.1274,	Adjusted R-squared:  0.1269 
# F-statistic: 236.1 on 10 and 16166 DF,  p-value: < 2.2e-16


lm.after = lm(z ~ x, subset = c(matchest$index.treated,
                                matchest$index.control))
summary(lm.after)
# Call:
# lm(formula = z ~ x, subset = c(matchest$index.treated, matchest$index.control))

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.66864 -0.49161 -0.03679  0.50378  0.65122 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  6.003e-01  2.427e-01   2.474   0.0137 *
# xage         3.199e-03  3.427e-03   0.933   0.3511  
# xeduc       -1.501e-02  1.634e-02  -0.918   0.3590  
# xblack       6.141e-05  7.408e-02   0.001   0.9993  
# xhispan      1.391e-02  1.208e-01   0.115   0.9084  
# xmarried    -1.328e-02  6.729e-02  -0.197   0.8437  
# xnodegree   -3.023e-02  7.144e-02  -0.423   0.6723  
# xre74        6.754e-06  9.864e-06   0.685   0.4939  
# xre75       -9.848e-06  1.279e-05  -0.770   0.4417  
# xu74         2.179e-02  1.027e-01   0.212   0.8321  
# xu75        -2.642e-02  8.327e-02  -0.317   0.7512  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5043 on 485 degrees of freedom
# Multiple R-squared:  0.005101,	Adjusted R-squared:  -0.01541 
# F-statistic: 0.2487 on 10 and 485 DF,  p-value: 0.9909





















