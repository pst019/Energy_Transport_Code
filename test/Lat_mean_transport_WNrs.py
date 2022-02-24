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


save= False
save= True
imp= False
imp= True

Member = 3
syear= 1950
#eyear= 1960
eyear= 2100

fignr= 1

#var= 'vEtot'
var= 'vQtot'


file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

fig= plt.figure(fignr)
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

axs=fig.subplots(2, 1, sharex= 'col')

if imp:
    for year in range(syear, eyear+1):
        if year//10 == 0: print(year)
        for month in range(1,13):

            if year == syear and month== 1:
                ds= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds.time), freq='6H').to_timestamp()
        
            else:
                ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
                ds2['lat']= ds.lat #to correct for some file issues
                ds= xr.concat([ds, ds2], dim= 'time')
    

#
#for WN in range(5):    
#    axs[0].plot(ds.lat, ds[var].sel(WaveNumb= WN).mean(dim='time'), label= WN)

axs[0].plot(ds.lat, ds[var].sum(dim= "WaveNumb").mean(dim='time'), label= 'Total', c= 'k')
#axs[0].plot(ds.lat, ds[var].where(ds.WaveNumb >=1, drop=True).sum(dim= "WaveNumb").mean(dim='time'), label= 'Total se', c= 'k')

axs[0].plot(ds.lat, ds[var].sel(WaveNumb= 0).mean(dim='time'), label= '0 = mc')
axs[0].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb').mean(dim='time'), label= '1-4') 
axs[0].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb').mean(dim='time'), label= '5-9') 
axs[0].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 10, ds.WaveNumb <= 22), drop= True ).sum(dim='WaveNumb').mean(dim='time'), label= '10+') 



axs[1].plot(ds.lat, ds[var].sum(dim="WaveNumb").where(ds["time.year"]> 2050).mean(dim='time')
                    - ds[var].sum(dim= "WaveNumb").where(ds["time.year"]< 2000).mean(dim='time'), label= 'Total', c= 'k')
axs[1].plot(ds.lat, ds[var].sel(WaveNumb= 0).where(ds["time.year"]> 2050).mean(dim='time')
                    - ds[var].sel(WaveNumb= 0).where(ds["time.year"]< 2000).mean(dim='time'), label= '0')
axs[1].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
                    - ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), label= '1-4') 
axs[1].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
                    - ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), label= '1-4') 
axs[1].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 10, ds.WaveNumb <= 22), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
                    - ds[var].where(np.logical_and(ds.WaveNumb >= 10, ds.WaveNumb <= 22), drop= True ).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), label= '1-4') 



axs[0].legend(ncol= 2)
#
axs[0].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[1].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
#
#
plt.xlim(-90, 90)
plt.xlabel('Latitude')
axs[0].set_ylabel('Energy transport')
axs[1].set_ylabel('2050-2100 - 1950-2000')


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Global/'
    
    file_name= 'WaveNr_M'+str(Member) +'_' + str(syear)+'-'+str(eyear)+'_'
    file_name += var
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
    
axs[0].set_title ('Member '+ str(Member))