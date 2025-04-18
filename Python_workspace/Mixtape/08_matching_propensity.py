import numpy as np 
import pandas as pd 
import statsmodels.api as sm 
import statsmodels.formula.api as smf 
from itertools import combinations 
import plotnine as p

# read data
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
def read_data(file): 
	return pd.read_stata("https://github.com/scunning1975/mixtape/raw/master/" + file)

nsw_dw = read_data('nsw_mixtape.dta')

# print(nsw_dw.head(9))

mean1 = nsw_dw[nsw_dw.treat==1].re78.mean()
mean0 = nsw_dw[nsw_dw.treat==0].re78.mean()

ate = np.unique(mean1 - mean0)[0]
# print("The experimental ATE estimate is {:.2f}".format(ate))
# The experimental ATE estimate is 1794.34

# Prepare data for logit 
nsw_dw_cpscontrol = read_data('cps_mixtape.dta')

nsw_dw_cpscontrol = pd.concat((nsw_dw_cpscontrol, nsw_dw))
nsw_dw_cpscontrol['u74'], nsw_dw_cpscontrol['u75'] = 0, 0
nsw_dw_cpscontrol.loc[nsw_dw_cpscontrol.re74==0, 'u74'] = 1
nsw_dw_cpscontrol.loc[nsw_dw_cpscontrol.re75==0, 'u75'] = 1

# estimating propensity score
logit_nsw = smf.glm(formula="""treat ~ age + I(age**2) + I(age**3) + educ + I(educ**2) + \
	marr + nodegree + black + hisp + re74 + re75 + u74 + u75 + educ*re74""",
	family=sm.families.Binomial(),
	data=nsw_dw_cpscontrol).fit()

nsw_dw_cpscontrol['pscore'] = logit_nsw.predict(nsw_dw_cpscontrol)

propensity_score_estimate = nsw_dw_cpscontrol.groupby('treat')['pscore'].mean()

p.ggplot(nsw_dw_cpscontrol, p.aes(x='pscore')) + p.geom_histogram(bins=50) + p.facet_wrap('treat', scales='free')


print(propensity_score_estimate)
