setwd("C:\\Users\\boris\\바탕 화면\\R_workspace")
getwd()


# Fisher Randomization Test in the Stratified Randomized Experiment

# Penn Bonus experiment to illustrate the FRT in the SRE.
# Data.set is from a job training program stratified on quarter
# outcome being the duration before employment.

penndata = read.table("Penn46_ascii.txt")

z = penndata$treatment
y = log(penndata$duration)
block = penndata$quarter

table(penndata$treatment, penndata$quarter)
#     0   1   2   3   4   5
# 0 234  41 687 794 738 860
# 1  87  48 757 866 811 461

# compute tau_hat_s and W_s
stat_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K = length(xlevels)
  PiK = rep(0, K)
  TauK = rep(0, K)
  WK = rep(0, K)
  for (k in 1:K)
  {
    xk = xlevels[k]
    zk = z[x == xk]
    yk = y[x == xk]
    PiK[k] = length(zk) / length(z)
    TauK[k] = mean(yk[zk==1]) - mean(yk[zk==0])
    WK[k] = wilcox.test(yk[zk==1], yk[zk==0])$statistic
  }
  
  return (c(sum(PiK*TauK), sum(WK/PiK)))
}

# function generates a random treatment assignment in SRE
# based on the observed data
zRandomSRE = function(z, x)
{
  xlevels = unique(x)
  K = length(xlevels)
  zrandom=z
  for (k in 1:K)
  {
    xk = xlevels[k]
    zrandom[x == xk] = sample(z[x == xk])
  }
  return (zrandom)
}

# Simulate the randomization distributions of the test statistics & compute p-values
stat.obs = stat_SRE(z, y, block)
MC = 10^3
statSREMC = matrix(0, MC, 2)
for (mc in 1:MC)
{
  zrandom = zRandomSRE(z, block)
  statSREMC[mc, ] = stat_SRE(zrandom, y, block)
}
mean(statSREMC[, 1] <= stat.obs[1])
# [1] 0.004
mean(statSREMC[, 2] <= stat.obs[2])
# [1] 0.001
# calculate the p-values based on left-tail probabilities 
# because the treatment has a negative effect on the outcome.
