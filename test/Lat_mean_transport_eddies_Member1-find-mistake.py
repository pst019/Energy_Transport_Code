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


save= True
save= False
imp= False
imp= True

Member = 1
syear= 1980
eyear= 1990
#eyear= 2100
#eyear= 2100

fignr= 3

var= 'vEmc'
var= 'vEse'
var= 'vEte'

varlist= ['vEse']
#varlist= ['vEmc', 'vEse', 'vEte']

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

fig= plt.figure(fignr)
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

axs=fig.subplots(2, 1, sharex= 'col')

if imp:
    for var in varlist:
        print(var)
        for year in range(syear, eyear+1):
            print(year)
            for month in range(1,13):
                print(month)
                if year == syear and month== 1:
                    ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    print(len(ds0.lat))
                else:
                    ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    ds2['lat']= ds0.lat # sometimes the latitudes are differently interpreted by xarray since there is a difference how the files are saves, check: ncdump -h vEse.1990.01.nc and ncdump -h vEse.1989.12.nc. This solves the strange issue.  
#                     print(ds2.isel(lat= -32))
                    print(len(ds2.lat))
                    
                    ds0= xr.concat([ds0, ds2], dim= 'time')
#                    print(ds0.isel(lat= -32))
                    print(len(ds0.lat))
                    
                """the last loop step to merge the different variables in one ds"""
                if year == eyear and month == 12:
    #                print(ds)
                    if var== varlist[0]:
                        ds= ds0
                    else:
                        ds= xr.merge([ds, ds0])


axs[0].plot(ds.lat, np.sum(ds[var].mean(dim='time') for var in varlist), label= 'Total', color= 'k')
#axs[1].plot(ds.lat, np.sum(ds[var].where(ds["time.year"]>= 2050).mean(dim='time') for var in varlist)
#                    - np.sum(ds[var].where(ds["time.year"]<= 2000).mean(dim='time') for var in varlist), label= 'Total', color= 'k')

for var in varlist:    
    axs[0].plot(ds.lat, ds[var].mean(dim='time'), label= var)

#    axs[1].plot(ds.lat, ds[var].where(ds["time.year"]>= 2050).mean(dim='time') - ds[var].where(ds["time.year"]<= 2000).mean(dim='time'), label= var)



#ds.where(ds["time.year"]< 2000, drop= True)
    
axs[0].legend(ncol= 2)

axs[0].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
#axs[1].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-90, 90)
plt.xlabel('Latitude')
axs[0].set_ylabel('Energy transport')
#axs[1].set_ylabel('2050-2100 - 1950-2000')


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Global/'
    
    file_name= 'Eddies_M'+str(Member) +'_' +str(syear)+'-'+str(eyear)+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 