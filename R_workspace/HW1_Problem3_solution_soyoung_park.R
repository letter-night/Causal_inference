
# Problem 3 (a) ######################################################################################

## Fisher Randomization Test

# LaLonde's experimental data to illustrate the FRT.
# the effect of a job training program using a randomized experiment.

# z: binary treatment vector indicating whether a unit was randomly assigned to the job program or not.
# y: outcome vector representing a unit's earnings in 1978.

library(Matching)
data("lalonde")

# Binary treatment
z = lalonde$treat

# Outcome
y = lalonde$re78


# Computing observed values of the test statistics
tauhat = t.test(y[z==1], y[z==0], var.equal = TRUE)$statistic
tauhat
# t = 2.835321

student = t.test(y[z==1], y[z==0], var.equal = FALSE)$statistic
student
# t = 2.674146

W = wilcox.test(y[z==1], y[z==0])$statistic
W
# W = 27402.5

D = ks.test(y[z==1], y[z==0])$statistic
D
# D = 0.1321206

# Monte Carlo approximation of the randomization distributions of the test statistics
# via randomly permuting the treatment vector.
MC = 10^4
Tauhat = rep(0, MC)
Student = rep(0, MC)
Wilcox = rep(0, MC)
Ks = rep(0, MC)

for (mc in 1:MC)
{
  zperm = sample(z)
  Tauhat[mc] = t.test(y[zperm==1], y[zperm==0], var.equal = TRUE)$statistic
  Student[mc] = t.test(y[zperm==1], y[zperm==0], var.equal = FALSE)$statistic
  Wilcox[mc] = wilcox.test(y[zperm==1], y[zperm==0])$statistic
  Ks[mc] = ks.test(y[zperm==1], y[zperm==0])$statistic
}

# One-sided p-values based on FRT
exact.pv = c(mean(Tauhat >= tauhat),
             mean(Student >= student),
             mean(Wilcox >= W),
             mean(Ks >= D))
round(exact.pv, 3)
# [1] 0.002 0.002 0.004 0.036
# p-values are all smaller than 0.05

# Asymptotic p-values without using Monte Carlo
asym.pv = c(t.test(y[z==1], y[z==0], var.equal=TRUE)$p.value,
            t.test(y[z==1], y[z==0], var.equal=FALSE)$p.value,
            wilcox.test(y[z==1], y[z==0])$p.value,
            ks.test(y[z==1], y[z==0])$p.value)
round(asym.pv, 3)
# [1] 0.005 0.008 0.011 0.046


# Problem 3 (b) ######################################################################################

## Covariate-adjusted FRT
## Pseudo-outcome strategy for covariate-adjusted FRT

# Using test statistics based on residuals from the linear regression.

# Covariates ------------------------------------------------------------------------- #
x = as.matrix(lalonde[, c("age", "educ", "hisp", "black", "married",
                          "nodegr", "re74", "re75", "u74", "u75")])
x = scale(x)

# Pseudo outcomes -------------------------------------------------------------------- #
# Run a linear regression of the outcomes on the covariates, obtains the residuals.
resid = lm(y ~ x)$resid

# Test statistics based on the residuals --------------------------------------------- #
tauhat = t.test(resid[z==1], resid[z==0], var.equal = TRUE)$statistic
tauhat
# t = 2.576048

student = t.test(resid[z==1], resid[z==0], var.equal = FALSE)$statistic
student
# t = 2.438348

W = wilcox.test(resid[z==1], resid[z==0])$statistic
W
# W = 26723

D = ks.test(resid[z==1], resid[z==0])$statistic
D
# D = 0.1390852

# FRT -------------------------------------------------------------------------------- #
# Monte Carlo approximation of the randomization distributions of the test statistics
# via randomly permuting the treatment vector.
MC = 10^4
Tauhat = rep(0, MC)
Student = rep(0, MC)
Wilcox = rep(0, MC)
Ks = rep(0, MC)

for (mc in 1:MC)
{
  zperm = sample(z)
  Tauhat[mc] = t.test(resid[zperm==1], resid[zperm==0], var.equal = TRUE)$statistic
  Student[mc] = t.test(resid[zperm==1], resid[zperm==0], var.equal = FALSE)$statistic
  Wilcox[mc] = wilcox.test(resid[zperm==1], resid[zperm==0])$statistic
  Ks[mc] = ks.test(resid[zperm==1], resid[zperm==0])$statistic
}

# One-sided p-values based on FRT ----------------------------------------------------- #
exact.pv = c(mean(Tauhat >= tauhat),
             mean(Student >= student),
             mean(Wilcox >= W),
             mean(Ks >= D))
round(exact.pv, 3)
# [1] 0.005 0.005 0.022 0.029


# Problem 3 (c) ######################################################################################

## Covariate-adjusted FRT
## Model-outcome strategy for covariate-adjusted FRT

# Using a regression coefficient as a test statistic.
# regress Yi on (Zi, Xi) to obtain the coefficient of Zi as the test statistic.

# Linear regression of the outcomes on the treatment & covariates -------------------- #
fit_adj = lm(y ~ z * x)

# Test statistics based on the coefficient ------------------------------------------- #
stat.obs = coef(fit_adj)[2]
stat.obs
# 1583.468

# FRT -------------------------------------------------------------------------------- #
# Monte Carlo approximation of the randomization distributions of the test statistics
# via randomly permuting the treatment vector.
MC = 10^4
statMC = rep(0, MC)

for (mc in 1:MC)
{
  zperm = sample(z)
  fit_adj = lm(y ~ zperm * x)
  ace_adj = coef(fit_adj)[2]
  statMC[mc] = ace_adj
}

# One-sided p-values based on FRT ----------------------------------------------------- #
exact.pv = mean(statMC >= stat.obs)
round(exact.pv, 3)
# [1] 0.008

# The end of the Problem 3. (a), (b), (c) solutions. ###############################################

