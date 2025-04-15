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

## Simple Difference in Outcomes
titanic = read_data("titanic.dta")

titanic['d'] = 0
titanic.loc[titanic['class']=='1st class', 'd'] = 1

titanic['sex_d'] = 0
titanic.loc[titanic['sex'] == 'man', 'sex_d'] = 1

titanic['age_d'] = 0
titanic.loc[titanic['age'] == 'adults', 'age_d'] = 1

titanic['survived_d'] = 0
titanic.loc[titanic['survived'] == 'yes', 'survived_d'] = 1

ey0 = titanic.loc[titanic['d'] == 0, 'survived_d'].mean()
ey1 = titanic.loc[titanic['d'] == 1, 'survived_d'].mean()

sdo = ey1 - ey0
# print("The simple difference in outcomes is {:.2f} %".format(sdo))
# The simple difference in outcomes is 0.35 %

titanic['s'] = 0
titanic.loc[(titanic.sex_d == 0) & (titanic.age_d == 1), 's'] = 1	# old female
titanic.loc[(titanic.sex_d == 0) & (titanic.age_d == 0), 's'] = 2	# young female
titanic.loc[(titanic.sex_d == 1) & (titanic.age_d == 1), 's'] = 3 	# old male
titanic.loc[(titanic.sex_d == 1) & (titanic.age_d == 0), 's'] = 4 	# young male

obs = titanic.loc[titanic.d == 0].shape[0] # ~ first class 

def weighted_avg_effect(df):
	diff = df[df.d==1].survived_d.mean() - df[df.d==0].survived_d.mean() # first class survived - ~first class survived
	weight = df[df.d==0].shape[0] / obs
	return diff * weight

wate = titanic.groupby('s').apply(weighted_avg_effect).sum()

print('The weighted average treatment effect estimate is {:.2f}%'.format(wate))
# The weighted average treatment effect estimate is 0.19%