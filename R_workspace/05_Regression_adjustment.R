rm(list=ls())
library("ggplot2")

library("car")
library("foreign")

# Regression Adjustment



## Angrist 2009 data: Canadian university
# Experiment to evaluate strategies to improve academic performance 
# treatment group offered academic support services and financial incentives.
# outcome: GPA

# Impute missing outcomes with the observed average.
# Two covariates: gender, baseline GPA

angrist = read.dta("star.dta")
angrist2 = subset(angrist, control == 1|sfsp == 1)

## imputing missing outcomes
y = angrist2$GPA_year1

meany = mean(y, na.rm = TRUE)
y = ifelse(is.na(y), meany, y)
z = angrist2$sfsp
x = angrist2[, c("female", 'gpa0')]

## unadjusted estimator
fit_unadj = lm(y ~ z)
ace_unadj = coef(fit_unadj)[2]
se_unadj = sqrt(hccm(fit_unadj, type="hc2")[2, 2])

## regression adjustment
x = scale(x)
fit_adj = lm(y ~ z*x)
ace_adj = coef(fit_adj)[2]
se_adj = sqrt(hccm(fit_adj, type="hc2")[2, 2])

res = c(ace_unadj, ace_adj, se_unadj, se_adj)
dim(res) = c(2, 2)
t.stat = res[, 1] / res[, 2]
p.value = 2*pnorm(abs(t.stat), lower.tail = FALSE)
res = cbind(res, t.stat, p.value)
rownames(res) = c("Neyman",  "Lin")
colnames(res) = c("estimate", "s.e.", "t-stat", "p-value")
round(res, 3)

#        estimate  s.e. t-stat p-value
# Neyman    0.052 0.078  0.669   0.504
# Lin       0.068 0.074  0.925   0.355

# The adjusted estimator has a smaller standard error although it gives the 
# same insignificant result as the unadjusted estimator.
