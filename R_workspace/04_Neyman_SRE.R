rm(list=ls())

# Neyman point and variance estimators under the SRE

Neyman_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K = length(xlevels)
  PiK = rep(0, K)
  TauK = rep(0, K)
  varK = rep(0, K)
  for (k in 1:K)
  {
    xk = xlevels[k]
    zk = z[x == xk]
    yk = y[x == xk]
    PiK[k] = length(zk) / length(z)
    TauK[k] = mean(yk[zk==1]) - mean(yk[zk==0])
    varK[k] = var(yk[zk==1]) / sum(zk) + var(yk[zk==0]) / sum(1-zk)
  }
  
  return (c(sum(PiK*TauK), sum(PiK^2*varK)))
}

# Simulation K=5, each stratum has 80 units.
# TauHat & VarHat: point & variance estimators over 10^4 simulations
K = 5
n = 80
n1 = 50
n0 = 30
x = rep(1:K, each = n)
y0 = rexp(n*K, rate=x)
y1 = y0 + 1
zb = c(rep(1, n1), rep(0, n0))
MC = 10^4

TauHat = rep(0, MC)
VarHat = rep(0, MC)

for (mc in 1:MC)
{
  z = replicate(K, sample(zb))
  z = as.vector(z)
  y = z*y1 + (1-z)*y0
  est = Neyman_SRE(z, y, x)
  TauHat[mc] = est[1]
  VarHat[mc] = est[2]
}

par(mfrow = c(2, 1), mai = c(1, 0.1, 0.2, 0.2))
hist(TauHat, xlab=expression(hat(tau)[S]),
     ylab="", main="a few large strata", border=FALSE, col='grey',
     breaks = 30, yaxt = 'n', xlim=c(0.8, 1.2))
abline(v=1)

var(TauHat)
# [1] 0.002651365
mean(VarHat)
# [1] 0.002595814

# The average value of variance estimator is almost identical to variance of estimators
# because individual causal effects are constant.

# Simulation K = 50, each stratum has 8 units.
K = 50
n = 8
n1 = 5
n0 = 3
x = rep(1:K, each=n)
y0 = rexp(n*K, rate = log(x + 1))
y1 = y0 + 1
zb = c(rep(1, n1), rep(0, n0))
MC = 10^4

TauHat = rep(0, MC)
VarHat = rep(0, MC)

for (mc in 1:MC)
{
  z = replicate(K, sample(zb))
  z = as.vector(z)
  y = z*y1 + (1-z)*y0
  est = Neyman_SRE(z, y, x)
  TauHat[mc] = est[1]
  VarHat[mc] = est[2]
}

hist(TauHat, xlab = expression(hat(tau)[S]),
     ylab = "", main = "many small strata",
     border=FALSE, col="grey",
     breaks=30, yaxt='n',
     xlim = c(0.8, 1.2))
abline(v=1)

var(TauHat)
mean(VarHat)
# [1] 0.001660618
# [1] 0.001683868

# histogram of the point estimator is symmetric and bell-shaped around true parameter.

# Neyman inference in SRE
# Penn Bonus Experiment

penndata = read.table(("Penn46_ascii.txt"))
z = penndata$treatment
y = log(penndata$duration)
block = penndata$quarter

est = Neyman_SRE(z, y, block)
est[1]
# [1] -0.08990646
sqrt(est[2])
# [1] 0.03079775

# So the job training program significantly shortens 
# the log of the duration time before employment.




















