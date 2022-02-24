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
#save= True
#imp= False
imp= True

Member = 3
syear= 1950
eyear= 1955
#eyear= 2100

fignr= 1

var= 'vEtot'

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
                ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()
#                if Member == Memberlist[0]: dslat= ds0.lat
#                else:
                ds0['lat']= dslat        
            else:
                ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()
                ds1['lat']= dslat #to correct for some file issues

                ds0= xr.concat([ds0, ds1], dim= 'time')
    

ds= ds0

a= 6370
LatExtent= 2*np.pi* a* np.cos(np.deg2rad(ds.lat))

Split1= 8000
#Split2= 4000

WaveSplit1= LatExtent/Split1

#
#for WN in range(5):    
#    axs[0].plot(ds.lat, ds[var].sel(WaveNumb= WN).mean(dim='time'), label= WN)

axs[0].plot(ds.lat, ds[var].sum(dim= "WaveNumb").mean(dim='time'), label= 'Total', c= 'k')
#axs[0].plot(ds.lat, ds[var].where(ds.WaveNumb >=1, drop=True).sum(dim= "WaveNumb").mean(dim='time'), label= 'Total se', c= 'k')

axs[0].plot(ds.lat, ds[var].sel(WaveNumb= 0).mean(dim='time'), label= '0 = mc')

#Planetary
#axs[0].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= WaveSplit1)).sum(dim='WaveNumb').mean(dim='time') , label= '$\lambda$ >'+str(Split1))
axs[0].plot(ds.lat, np.nansum([ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb').mean(dim='time'),
                     + WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).mean(dim='time')  ], axis= 0) , label= '$\lambda$ >'+str(Split1)) 


axs[0].plot(ds.lat, np.nansum([ds[var].where(ds.WaveNumb > WaveSplit1//1).sum(dim='WaveNumb').mean(dim='time'), 
                     (1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).mean(dim='time') ], axis= 0)
, label= str(Split1)+ '> $\lambda$') 

axs[1].plot(ds.lat, WaveSplit1)

#axs[1].plot(ds.lat, ds[var].sum(dim="WaveNumb").where(ds["time.year"]> 2050).mean(dim='time')
#                    - ds[var].sum(dim= "WaveNumb").where(ds["time.year"]< 2000).mean(dim='time'), label= 'Total', c= 'k')
#axs[1].plot(ds.lat, ds[var].sel(WaveNumb= 0).where(ds["time.year"]> 2050).mean(dim='time')
#                    - ds[var].sel(WaveNumb= 0).where(ds["time.year"]< 2000).mean(dim='time'), label= '0')
#axs[1].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= WaveSplit1)).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
#                    - ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= WaveSplit1)).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), label= '$\lambda$ >'+str(Split1)) 
#axs[1].plot(ds.lat, ds[var].where(np.logical_and(ds.WaveNumb > WaveSplit1, ds.WaveNumb <= WaveSplit2)).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
#                    - ds[var].where(np.logical_and(ds.WaveNumb > WaveSplit1, ds.WaveNumb <= WaveSplit2)).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), str(Split1) + '> $\lambda$ >'+str(Split2)) 
#axs[1].plot(ds.lat, ds[var].where(ds.WaveNumb > WaveSplit2).sum(dim='WaveNumb').where(ds["time.year"]> 2050).mean(dim='time')
#                    - ds[var].where(ds.WaveNumb > WaveSplit2).sum(dim='WaveNumb').where(ds["time.year"]< 2000).mean(dim='time'), label= str(Split2)+ '> $\lambda$') 



axs[0].legend(ncol= 2)
#
axs[0].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
#
#
plt.xlim(-90, 90)
plt.xlabel('Latitude')
axs[0].set_ylabel('Energy transport')

#axs[1].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
#axs[1].set_ylabel('2050-2100 - 1950-2000')


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Global/'
    
    file_name= 'WaveLenghts_M'+str(Member) +'_' + str(syear)+'-'+str(eyear)+'_'
    file_name += var
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
    
axs[0].set_title ('Member '+ str(Member))