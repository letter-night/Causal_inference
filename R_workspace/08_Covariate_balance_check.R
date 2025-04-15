

## balance check for K=5
library(ggplot2)

Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  q.pscore = quantile(pscore, (1:4)/5)
                  ps.strata = cut(pscore, breaks = c(0, q.pscore, 1),
                                  labels=1:5)
                  Neyman_SRE(z, x[, px], ps.strata)
                })

dat_balance = data.frame(est = Bcheck[1, ],
                         upper = Bcheck[1, ] + 1.96*Bcheck[2, ],
                         lower = Bcheck[1, ] - 1.96*Bcheck[2, ],
                         cov = factor(1:11))

ggplot(dat_balance) +
  geom_errorbar(aes(x = cov,
                    ymin=lower,
                    ymax=upper),
                alpha=0.6) +
  geom_point(aes(x=cov,
                 y=est),
             alpha=0.6) + 
  geom_hline(aes(yintercept=0),
             alpha=0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()) +
  xlab("balance check based on stratification with K=5")
ggsave("balance_PSstratification5.pdf", height=3, width=8)

## balance check based on Hajek
Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  ipw.boot(z, x[, px], x)[2, ]
                })

dat_balance = data.frame(est = Bcheck[1, ],
                         upper = Bcheck[1, ] + 1.96*Bcheck[2, ],
                         lower = Bcheck[1, ] - 1.96*Bcheck[2, ],
                         cov = factor(1:11))
ggplot(dat_balance) + 
  geom_errorbar(aes(x = cov,
                    ymin = lower,
                    ymax = upper),
                alpha = 0.6) + 
  geom_point(aes(x = cov,
                 y = est),
             alpha = 0.6) +
  geom_hline(aes(yintercept = 0),
             alpha = 0.3) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank()) +
  xlab("balance check based on weighting")
ggsave("balance_PShajek.pdf", height = 3, width = 8)