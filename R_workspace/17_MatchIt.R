library("MatchIt")
data("lalonde")

head(lalonde)

## Check Initial Imbalance ------------------------------------------------------------------ #
m.out0 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                  data=lalonde,
                  method=NULL,
                  distance = "glm")

# Checking balance prior to matching
summary(m.out0)
# Summary of Balance for All Data:
#            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
# distance          0.5774        0.1822          1.7941     0.9211    0.3774   0.6444
# age              25.8162       28.0303         -0.3094     0.4400    0.0813   0.1577
# educ             10.3459       10.2354          0.0550     0.4959    0.0347   0.1114
# raceblack         0.8432        0.2028          1.7615          .    0.6404   0.6404
# racehispan        0.0595        0.1422         -0.3498          .    0.0827   0.0827
# racewhite         0.0973        0.6550         -1.8819          .    0.5577   0.5577
# married           0.1892        0.5128         -0.8263          .    0.3236   0.3236
# nodegree          0.7081        0.5967          0.2450          .    0.1114   0.1114
# re74           2095.5737     5619.2365         -0.7211     0.5181    0.2248   0.4470
# re75           1532.0553     2466.4844         -0.2903     0.9563    0.1342   0.2876

# severe imbalance as measured by the standardized mean differences (Std. Mean Diff.),
# variance ratios (Var. Ratio), empirical cumulative distribution function (eCDF) statistics.

# Values of standardized mean differences and eCDF statistics close to zero and
# values of variance ratios close to one indicate good balance.

# Here, many of them are far from their ideal values.

## Matching -------------------------------------------------------------------------------- #
# 1:1 nearest neighbor (NN) matching on the propensity score
# which is appropriate for estimating the ATT.
# One by one, each treated unit is paired with an available control unit 
# that has the closest propensity score to it. 
# Any remaining control units are left unmatched and excluded from further analysis.

# Propensity score matching can be an effective way to achieve covariate balance
# in treatment groups.

# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(treat ~ age + educ + race + married +
                    nodegree + re74 + re75,
                  data=lalonde,
                  method = "nearest",
                  distance = "glm") # logistic regression propensity score.
# m.out1
# A `matchit` object
# - method: 1:1 nearest neighbor matching without replacement
# - distance: Propensity score
# - estimated with logistic regression
# - number of obs.: 614 (original), 370 (matched)
# - target estimand: ATT
# - covariates: age, educ, race, married, nodegree, re74, re75

# Checking balance after NN matching
summary(m.out1, un = FALSE) # `un = TRUE` display the balance before & after matching

# Summary of Balance for Matched Data:
#            Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
# distance          0.5774        0.3629          0.9739     0.7566    0.1321   0.4216
# age              25.8162       25.3027          0.0718     0.4568    0.0847   0.2541
# educ             10.3459       10.6054         -0.1290     0.5721    0.0239   0.0757
# raceblack         0.8432        0.4703          1.0259          .    0.3730   0.3730
# racehispan        0.0595        0.2162         -0.6629          .    0.1568   0.1568
# racewhite         0.0973        0.3135         -0.7296          .    0.2162   0.2162
# married           0.1892        0.2108         -0.0552          .    0.0216   0.0216
# nodegree          0.7081        0.6378          0.1546          .    0.0703   0.0703
# re74           2095.5737     2342.1076         -0.0505     1.3289    0.0469   0.2757
# re75           1532.0553     1614.7451         -0.0257     1.4956    0.0452   0.2054
#            Std. Pair Dist.
# distance            0.9740
# age                 1.3938
# educ                1.2474
# raceblack           1.0259
# racehispan          1.0743
# racewhite           0.8390
# married             0.8281
# nodegree            1.0106
# re74                0.7965
# re75                0.7381

# Although balance has improved for some covariates, in general balance is still quite poor,
# indicating that nearest neighbor propensity score matching is not sufficient for
# removing confounding in this dataset.

# Distribution of propensity scores of those who were matched 
plot(m.out1, type = "jitter", interactive = FALSE)

# Examine balance on the covariates
plot(m.out1, type = "density", interactive = FALSE, which.xs = ~age + married + re75)

# Imbalances are represented by the differences between the black(treated) & gray(control)
# distributions.
# Although married and re75 appear to have improved balance after matching, the case is mixed for age.

## Trying a Different Matching Specification ----------------------------------------------- #
# Given the poor performance of nearest neighbor matching in this example, 
# try a different matching method or make other changes to the matching algorithm 
# or distance specification.

# Full matching: matches every treated unit to at least one control and every control to 
# at least one treated unit.

# Probit for propensity score model.

# Full matching on a probit PS
m.out2 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                  data=lalonde,
                  method = "full",
                  distance = "glm",
                  link = "probit")
m.out2
# A `matchit` object
# - method: Optimal full matching
# - distance: Propensity score
# - estimated with probit regression
# - number of obs.: 614 (original), 614 (matched)
# - target estimand: ATT
# - covariates: age, educ, race, married, nodegree, re74, re75

# Checking balance after full matching
summary(m.out2, un = FALSE)
# Summary of Balance for Matched Data:
#          Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance        0.5773        0.5764          0.0045     0.9949    0.0043   0.0486          0.0198
# age            25.8162       25.5347          0.0393     0.4790    0.0787   0.2742          1.2843
# educ           10.3459       10.5381         -0.0956     0.6192    0.0253   0.0730          1.2179
# raceblack       0.8432        0.8389          0.0119          .    0.0043   0.0043          0.0162
# racehispan      0.0595        0.0492          0.0435          .    0.0103   0.0103          0.4412
# racewhite       0.0973        0.1119         -0.0493          .    0.0146   0.0146          0.3454
# married         0.1892        0.1633          0.0660          .    0.0259   0.0259          0.4473
# nodegree        0.7081        0.6577          0.1110          .    0.0504   0.0504          0.9872
# re74         2095.5737     2100.2150         -0.0009     1.3467    0.0314   0.1881          0.8387
# re75         1532.0553     1561.4420         -0.0091     1.5906    0.0536   0.1984          0.8240

# Balance is far better, as determined by the lower standardized mean differences and eCDF statistics.
# The balance should be reported when publishing the results of a matching analysis.
# This can be done either in a table, using the values resulting from summary(), or in a plot,
# such as a Love plot.
plot(summary(m.out2))

## Estimating the Treatment Effect --------------------------------------------------------- #

# Run a regression of the outcome on the treatment and covariates in the matched sample
# and estimate the treatment effect
# Including the covariates used in the matching provide additional robustness to slight imbalances
# remaining after matching and can improve precision.

# 1. Extract the matched dataset - only contains the matched units + distance, weights, subclass
m.data <- match_data(m.out2)
head(m.data)

# 2. Model the outcome in this dataset using the standard regression functions
# Be sure to include the matching weights in the estimation
library("marginaleffects")

fit <- lm(re78 ~ treat * (age + educ + race + married + nodegree + re74 + re75),
          data = m.data,
          weights = weights)

# 3. G-computation to estimate the ATT.
avg_comparisons(fit,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(treat == 1))

# Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#     1977        704 2.81  0.00501 7.6   596   3357


