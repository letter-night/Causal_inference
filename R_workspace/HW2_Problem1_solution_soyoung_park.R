rm(list=ls())
library("car")
library("Matching")

## Problem 1. NHANES data to study whether participation in school meal programs led to an increase in BMI for school children.

## Data "nhanes_bmi"
# treatment: School_meal - indicator for participation in the school meal plan.
# outcome: BMI
nhanes_bmi = read.csv("./data/nhanes_bmi.csv")[, -1]
head(nhanes_bmi)
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

#############################################################################################
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


#############################################################################################
## (b) Nearest-neighbor matching without replacement using the estimated propensity scores.
##     1) k = 1 matches per treated unit, covariate balance check before & after matching
##     2) k = 2 matches per treated unit, covariate balance check before & after matching
##     3) alternative approaches for improving covariate balance.

library("MatchIt")

## Check Initial Imbalance ---------------------------------------------------------------- #
m.out0 <- matchit(School_meal ~ age + ChildSex + black + mexam + pir200_plus
                  + WIC + Food_Stamp + fsdchbi + AnyIns + RefSex + RefAge,
                  data=nhanes_bmi,
                  method=NULL,
                  distance = "glm")

# Checking balance prior to matching
summary(m.out0)
# Summary of Balance for All Data:
#             Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
# distance           0.6722        0.4023          1.2702     0.8596    0.3029   0.4539
# age               10.0732        9.8528          0.0628     0.6264    0.0702   0.1607
# ChildSex           0.5117        0.5249         -0.0264          .    0.0132   0.0132
# black              0.3084        0.1989          0.2372          .    0.1096   0.1096
# mexam              0.3279        0.1778          0.3197          .    0.1501   0.1501
# pir200_plus        0.2469        0.6616         -0.9617          .    0.4147   0.4147
# WIC                0.2555        0.1099          0.3336          .    0.1455   0.1455
# Food_Stamp         0.4408        0.1166          0.6529          .    0.3242   0.3242
# fsdchbi            0.3255        0.1482          0.3785          .    0.1774   0.1774
# AnyIns             0.8380        0.8862         -0.1309          .    0.0482   0.0482
# RefSex             0.3941        0.5029         -0.2226          .    0.1088   0.1088
# RefAge            38.5506       40.2859         -0.1675     1.1467    0.0368   0.1518

# severe imbalance as measured by the standardized mean differences (Std. Mean Diff.),
# variance ratios (Var. Ratio), empirical cumulative distribution function (eCDF) statistics.

# Values of standardized mean differences and eCDF statistics close to zero and
# values of variance ratios close to one indicate good balance.

# Here, many of them are far from their ideal values.

## Matching (k = 1) ------------------------------------------------------------------------ #
# 1:1 nearest neighbor (NN) matching on the propensity score
# which is appropriate for estimating the ATT.
# One by one, each treated unit is paired with an available control unit 
# that has the closest propensity score to it. 
# Any remaining control units are left unmatched and excluded from further analysis.

# Propensity score matching can be an effective way to achieve covariate balance
# in treatment groups.

# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(School_meal ~ age + ChildSex + black + mexam + pir200_plus
                  + WIC + Food_Stamp + fsdchbi + AnyIns + RefSex + RefAge,
                  data=nhanes_bmi,
                  method="nearest",
                  distance = "glm") # logistic regression propensity score

# Checking balance after 1:1 NN matching
summary(m.out1, un = FALSE) # `un = TRUE` display the balance before & after matching
# Summary of Balance for Matched Data:
#         Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance       0.7502        0.4023          1.6372     0.3886    0.3921   0.6291          1.6372
# age           10.2992        9.8528          0.1271     0.6373    0.0693   0.1616          1.4533
# ChildSex       0.5086        0.5249         -0.0325          .    0.0163   0.0163          0.9697
# black          0.3050        0.1989          0.2298          .    0.1061   0.1061          0.8632
# mexam          0.3595        0.1778          0.3869          .    0.1816   0.1816          0.8228
# pir200_plus    0.1042        0.6616         -1.2926          .    0.5574   0.5574          1.3015
# WIC            0.2964        0.1099          0.4275          .    0.1864   0.1864          0.6861
# Food_Stamp     0.5315        0.1166          0.8357          .    0.4149   0.4149          0.8396
# fsdchbi        0.3795        0.1482          0.4937          .    0.2314   0.2314          0.8528
# AnyIns         0.8260        0.8862         -0.1635          .    0.0602   0.0602          0.6876
# RefSex         0.3671        0.5029         -0.2778          .    0.1358   0.1358          1.0369
# RefAge        38.3709       40.2859         -0.1848     1.1513    0.0379   0.1539          1.0822

# Although balance has improved for some covariates, in general balance is still quite poor,
# indicating that nearest neighbor propensity score matching is not sufficient for
# removing confounding in this dataset.

# Distribution of propensity scores of those who were matched 
plot(m.out1, type = "jitter", interactive = FALSE)

# Examine balance on the covariates
plot(m.out1, type = "density", interactive = FALSE, which.xs = ~age + ChildSex + black)
plot(m.out1, type = "density", interactive = FALSE, which.xs = ~mexam + pir200_plus + WIC)
plot(m.out1, type = "density", interactive = FALSE, which.xs = ~Food_Stamp + fsdchbi + AnyIns)
plot(m.out1, type = "density", interactive = FALSE, which.xs = ~RefSex + RefAge)

# Imbalances are represented by the differences between the black(treated) & gray(control)
# distributions.

vignette("matching-methods")

## Matching (k = 2) ------------------------------------------------------------------------ #
# 2:1 NN PS matching w/o replacement
m.out2 <- matchit(School_meal ~ age + ChildSex + black + mexam + pir200_plus
                  + WIC + Food_Stamp + fsdchbi + AnyIns + RefSex + RefAge,
                  data=nhanes_bmi,
                  method="nearest",
                  ratio = 2,
                  distance = "glm") # logistic regression propensity score

# Checking balance after 2:1 NN matching
summary(m.out2, un = FALSE) # `un = TRUE` display the balance before & after matching
# Summary of Balance for Matched Data:
#         Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance       0.7502        0.4023          1.6372     0.3886    0.3921   0.6291          1.6372
# age           10.2992        9.8528          0.1271     0.6373    0.0693   0.1616          1.4533
# ChildSex       0.5086        0.5249         -0.0325          .    0.0163   0.0163          0.9697
# black          0.3050        0.1989          0.2298          .    0.1061   0.1061          0.8632
# mexam          0.3595        0.1778          0.3869          .    0.1816   0.1816          0.8228
# pir200_plus    0.1042        0.6616         -1.2926          .    0.5574   0.5574          1.3015
# WIC            0.2964        0.1099          0.4275          .    0.1864   0.1864          0.6861
# Food_Stamp     0.5315        0.1166          0.8357          .    0.4149   0.4149          0.8396
# fsdchbi        0.3795        0.1482          0.4937          .    0.2314   0.2314          0.8528
# AnyIns         0.8260        0.8862         -0.1635          .    0.0602   0.0602          0.6876
# RefSex         0.3671        0.5029         -0.2778          .    0.1358   0.1358          1.0369
# RefAge        38.3709       40.2859         -0.1848     1.1513    0.0379   0.1539          1.0822

# Although balance has improved for some covariates, in general balance is still quite poor,
# indicating that nearest neighbor propensity score matching is not sufficient for
# removing confounding in this dataset.

# Distribution of propensity scores of those who were matched 
plot(m.out2, type = "jitter", interactive = FALSE)

# Examine balance on the covariates
plot(m.out2, type = "density", interactive = FALSE, which.xs = ~age + ChildSex + black)
plot(m.out2, type = "density", interactive = FALSE, which.xs = ~mexam + pir200_plus + WIC)
plot(m.out2, type = "density", interactive = FALSE, which.xs = ~Food_Stamp + fsdchbi + AnyIns)
plot(m.out2, type = "density", interactive = FALSE, which.xs = ~RefSex + RefAge)

# Imbalances are represented by the differences between the black(treated) & gray(control)
# distributions.

## Trying a Different Matching Specification ----------------------------------------------- #
# Given the poor performance of nearest neighbor matching in this example, 
# try a different matching method or make other changes to the matching algorithm 
# or distance specification.

# Full matching: matches every treated unit to at least one control and every control to 
# at least one treated unit.

# Probit for propensity score model.

# Full matching on a probit PS
m.out3 <- matchit(School_meal ~ age + ChildSex + black + mexam + pir200_plus
                  + WIC + Food_Stamp + fsdchbi + AnyIns + RefSex + RefAge,
                  data=nhanes_bmi,
                  method="full",
                  distance = "glm",
                  link = "probit")
m.out3

# Checking balance after full matching
summary(m.out3, un = FALSE)
# Summary of Balance for Matched Data:
#            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance          0.6694        0.6694          0.0001     0.9965    0.0014   0.0125          0.0044
# age              10.0732       10.0989         -0.0073     0.5872    0.0840   0.1510          0.8170
# ChildSex          0.5117        0.4960          0.0313          .    0.0157   0.0157          0.9383
# black             0.3084        0.3355         -0.0586          .    0.0270   0.0270          0.6401
# mexam             0.3279        0.3061          0.0463          .    0.0217   0.0217          0.7137
# pir200_plus       0.2469        0.2516         -0.0109          .    0.0047   0.0047          0.2546
# WIC               0.2555        0.2432          0.0280          .    0.0122   0.0122          0.6662
# Food_Stamp        0.4408        0.4361          0.0095          .    0.0047   0.0047          0.3402
# fsdchbi           0.3255        0.3495         -0.0511          .    0.0239   0.0239          0.5564
# AnyIns            0.8380        0.8324          0.0152          .    0.0056   0.0056          0.6449
# RefSex            0.3941        0.3957         -0.0033          .    0.0016   0.0016          0.9356
# RefAge           38.5506       38.9914         -0.0425     0.9512    0.0202   0.0920          0.9843

# Balance is far better, as determined by the lower standardized mean differences and eCDF statistics.
# The balance should be reported when publishing the results of a matching analysis.
# This can be done either in a table, using the values resulting from summary(), or in a plot,
# such as a Love plot.
plot(summary(m.out3))

#############################################################################################
## (c) Using matched sample, estimate the average treatment effect on the treated (ATT)
##     1) Simple difference in matched outcomes
##     2) Regression adjustment on the matched sample

# Extract the matched dataset: contains matched units, distance, weights, subclass
m.data <- match_data(m.out3)
head(m.data)

## 1) Simple difference in matched outcomes ------------------------------------------------ #
fit <- lm(BMI ~ School_meal,
          data=m.data,
          weights=weights)

avg_comparisons(fit, variables = "School_meal",
                newdata = subset(School_meal == 1))

# Estimate Std. Error      z Pr(>|z|)   S  2.5 % 97.5 %
#  -0.0875      0.228 -0.384    0.701 0.5 -0.535   0.36

## 2) Regression adjustment on the matched sample ------------------------------------------ #
library("marginaleffects")

fit <- lm(BMI ~ School_meal * (age + ChildSex + black + mexam + pir200_plus + WIC + 
                                 Food_Stamp + fsdchbi + AnyIns + RefSex + RefAge),
          data=m.data,
          weights=weights)

avg_comparisons(fit, variables = "School_meal",
                vcov = ~subclass,
                newdata = subset(School_meal == 1))

# Estimate Std. Error      z Pr(>|z|)   S  2.5 % 97.5 %
#  -0.0481       0.31 -0.155    0.876 0.2 -0.655  0.559

################################################################################################
