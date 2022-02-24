#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup/'

else:
    Mediadir= '/run/media/pst019/Backup/'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob as glob
import datetime

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}

year = 1999
month= 1
Member = 1

fignr= 2

var= 'vEmc'
var= 'vEse'
#var= 'vEte'
#var= 'vEsetot'



file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

#ds = xr.open_mfdataset(file_dir+var +'.*.nc')
#ds = xr.merge([xr.open_dataset(f) for f in glob.glob(file_dir+var +'.*.nc')])

ds= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc')

time= pd.date_range(str(year)+'-'+str(month), periods= 1)

#ds['time']= time

ds = ds.expand_dims('time').assign_coords(time=('time', time))



ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').expand_dims('time').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))


ds3= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

syear= 1950
for year in range(syear, 2100):
    for month in range(1,13):
        if year == syear and month== 1:
            ds= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
    
        else:
            ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
            ds2['lat']= ds0.lat #to correct for some file issues

            ds= xr.concat([ds, ds2], dim= 'time')

#ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month+1).zfill(2)+ '.nc')
#
#
#ds0= xr.merge([ds, ds2])

#times
fig= plt.figure(fignr)
fignr+=1
plt.clf()

plt.plot(ds.lat, ds[var].mean(dim='time'))