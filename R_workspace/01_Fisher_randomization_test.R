
# Fisher Randomization Test

# LaLonde's experimental data to illustrate the FRT.
library(Matching)
data("lalonde")

z = lalonde$treat
y = lalonde$re78

# z: binary treatment vector indicating whether a unit was randomly assigned to the
#    job program or not
# y: outcome vector representing a unit's earnings in 1978.

# Computing observed values of the test statistics
tauhat = t.test(y[z==1], y[z==0], var.equal=TRUE)$statistic
tauhat
# t = 2.835321

student = t.test(y[z==1], y[z==0], var.equal=FALSE)$statistic
student
# t = 2.674146

W = wilcox.test(y[z==1], y[z==0])$statistic
W
# W = 27402.5

D = ks.test(y[z==1], y[z==0])$statistic
D
# D = 0.1321206

# By randomly permuting the treatment vector, 
# we can obtain the Monte Carlo approximation of the randomization distributions
# of the test statistics
MC = 10^4
Tauhat = rep(0, MC)
Student = rep(0, MC)
Wilcox = rep(0, MC)
Ks = rep(0, MC)

for (mc in 1:MC)
{
  zperm = sample(z)
  Tauhat[mc] = t.test(y[zperm==1], y[zperm==0],
                      val.equal=TRUE)$statistic
  Student[mc] = t.test(y[zperm==1], y[zperm==0],
                       val.equal=FALSE)$statistic
  Wilcox[mc] = wilcox.test(y[zperm==1], y[zperm==0])$statistic
  Ks[mc] = ks.test(y[zperm==1], y[zperm==0])$statistic
}  

# One-sided p-values based on FRT are all smaller than 0.05
exact.pv = c(mean(Tauhat >= tauhat),
             mean(Student >= student),
             mean(Wilcox >= W),
             mean(Ks >= D))
round(exact.pv, 3)
# 0.001 0.002 0.004 0.034


# Without using Monte Carlo, we can compute the asymptotic p-values
# all smaller than 0.05
asym.pv = c(t.test(y[z==1], y[z==0], var.equal = TRUE)$p.value,
            t.test(y[z==1], y[z==0], var.equal = FALSE)$p.value,
            wilcox.test(y[z==1], y[z==0])$p.value,
            ks.test(y[z==1], y[z==0])$p.value)
round(asym.pv, 3)
# [1] 0.005 0.008 0.011 0.046


# “when possible, report a Fisher-exact p-value
# and display its underlying null randomization distribution.”







