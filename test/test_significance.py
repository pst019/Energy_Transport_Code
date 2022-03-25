#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 11:32:46 2022

@author: pst019
"""

import numpy as np



def statistic(x, y):
    return np.mean(x) - np.mean(y)

from scipy.stats import permutation_test

res = permutation_test((x, y), statistic, vectorized=True,
                       n_resamples=np.inf, alternative='less')


####own testing

def mean_diff(ds1, ds2, dim=['time', 'Member']):
    return ds1.mean(dim) - ds2.mean(dim)


distr_1= dstot.where(ds["time.year"]>= split_2, drop= True)
distr_2= dstot.where(ds["time.year"]<= split_1, drop= True)

d1= distr_1.stack(x= ("time", "Member"))
d2= distr_2.stack(x= ("time", "Member"))



###diffmean
diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') 
           - dstot.where(ds["time.year"]<= split_1).mean(dim='time') 
           ).mean(dim='Member')


ds1= dstot.where(ds["time.year"]>= split_2, drop= True)
ds2= dstot.where(ds["time.year"]<= split_1, drop= True)

dm= mean_diff(ds1, ds2, dim=['time', 'Member'])

permutation_test((ds1, ds2), mean_diff)


####
def mean_diff(ds1, ds2, dim=('time', 'Member')):
    d1= np.array(ds1.stack(x= dim))
    d2= np.array(ds2.stack(x= dim))

    return d1.mean(axis=0) - d2.mean(axis= -1)

ds1= dstot.where(ds["time.year"]>= split_2, drop= True).stack(x= ("time", "Member") )
ds2= dstot.where(ds["time.year"]<= split_1, drop= True).stack(x= ("time", "Member") )

###
ads1= np.array(dstot.where(ds["time.year"]>= split_2, drop= True).stack(x= ("time", "Member")).isel(lat= 0))
ads2= np.array(dstot.where(ds["time.year"]<= split_1, drop= True).stack(x= ("time", "Member")).isel(lat= 0))


ads1= np.array(dstot.where(ds["time.year"]>= split_2, drop= True).stack(x= ("time", "Member")) )
ads2= np.array(dstot.where(ds["time.year"]<= split_1, drop= True).stack(x= ("time", "Member")) )

res= permutation_test((ads1, ads2), statistic)
res= permutation_test((ads1, ads2), statistic, axis= -1, n_resamples= 1000)

def statistic(x, y):
    return np.mean(x, axis= -1) - np.mean(y, axis= -1)

#test1
#is the ensemble-mean of the standdard deviation the same
# as the standard deviation of the whole
# almost

a= dstot.mean(dim='time').mean(dim='Member')
b= dstot.mean(dim=['time', 'Member'])
a-b

a= dstot.std(dim='time').mean(dim='Member')
b= dstot.std(dim=['time', 'Member'])
a-b



"""significance t-test"""
# from scipy import stats
# distr_1= dstot.where(ds["time.year"]>= split_2, drop= True)
# distr_2= dstot.where(ds["time.year"]<= split_1, drop= True)
# d1= distr_1.stack(x= ("time", "Member"))
# d2= distr_2.stack(x= ("time", "Member"))
# p_value= stats.ttest_ind(d1, d2, axis= 1)[1]
# axs[ip].plot(ds.lat, fac*diffmean.where(p_value >.05, np.nan), label= 'Total', color= 'k')

# axs[ip].plot(ds.lat, fac*diffmean, label= 'Total', color= 'k')
