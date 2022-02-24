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
import cartopy.crs as ccrs
from scipy import fftpack

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save= False
#save= True
imp= False
imp= True

Member = 2
syear= 1950
#eyear= 1970
eyear= 2100


fignr= 4

#var, varfile= 'Z', 'Z850'
#var, varfile= 'STR', 'STR'
#var, varfile= 'T2M', 'T2M'
var, varfile= 'vQtot', 'vQtot'
#var, varfile = 'vEtotLLd', 'vEtotLL'

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

fig= plt.figure(fignr)
fignr+=1
plt.clf()

#axs= plt.subplot(2, 1, 1)

axs=fig.subplots(2, 1)

if imp:
    for year in range(syear, eyear+1):
        if year%10 == 0: print(year)
        for month in range(1,13):

            if year == syear and month== 1:
                if varfile == 'vEtotLL':
                    ds= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                    ds['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds.time), freq='1D').to_timestamp()
                elif varfile == 'vQtot':
                    ds= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds.time), freq='6H').to_timestamp()
                else: ds= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                
                
                if varfile in ['Z850']: ds= ds.isel(plev= 0)
                ds= ds.resample(time='1M').mean()
            else:
                if varfile == 'vEtotLL':
                    ds2= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                    ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='1D').to_timestamp()
                elif varfile == 'vQtot':
                    ds2= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
                else: ds2= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')

                if varfile in ['Z850']:  ds2= ds2.isel(plev= 0)
                ds2= ds2.resample(time='1M').mean()
#                 ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
#                ds2['lat']= ds.lat #to correct for some file issues

                ds= xr.concat([ds, ds2], dim= 'time')
    
#    ds= ds.resample(time='1M').mean()
#    print('Finished import')

#    ds= ds.where(ds["time.year"] < 2097, drop= True)

if var == 'T2M':
    weights= np.cos(np.deg2rad(ds.lat))
    #weights.name = 'weights'
    
    lat= 60
    ds0= ds.where(ds.lat > lat, drop= True)
    time_series= ds0[var].weighted(weights).mean(("lon", "lat"))



if var == 'vQtot':
    ds0= ds
    lat = 70
    time_series= ds[var].sel(lat=lat, method='nearest').where(np.logical_and(ds.WaveNumb>0, ds.WaveNumb < 5), drop= True ).sum(dim='WaveNumb')

axs[0].plot(ds0.time, time_series, label='Monthly mean' )

axs[0].plot(ds0.time, time_series.rolling(time= 12, center= True).mean(), label='1 year running mean' )
axs[0].plot(ds0.time, time_series.rolling(time= 120, center= True).mean(), label='10 year running mean' )







annual =time_series.resample(time='1Y').mean()
z= np.polyfit(annual["time.year"], annual.values, 1)
p = np.poly1d(z)
axs[0].plot(annual.time, p(annual["time.year"]), label= 'Linear trend')

print('Trend [unit per 100 year]:', 100* z[0])
print('Trend [% per 100 years]:', 100* z[0]/annual.mean().values * 100)


axs[0].legend()
axs[0].set_xlabel('Year')
axs[0].set_ylabel(var)
axs[0].set_xlim(ds['time'][0], ds['time'][-1])





fft = fftpack.fft(time_series.values)

#f= 1/12 #frequency
#n = int(len(ds.time))
#Fs= 12 #sampling rate
#fr= Fs/2* np.linspace(0,1, n//2)
#y_m = 2/n * abs(fft[0:np.size(fr)])
#axs[1].stem(fr[1:], y_m[1:])




f_s= 12 #sampling rate (measurements per year)
freqs = fftpack.fftfreq(len(time_series.values)) * f_s

n = len(ds.time)

#axs[1].stem(freqs[1:n//2], 2/n *np.abs(fft[1:n//2]))
#axs[1].set_xlabel('Frequency [1/year]')
#axs[1].set_ylabel('Signal amplityde')
#axs[1].set_xscale('log')


axs[1].stem(1/freqs[1:n//2], 2/n *np.abs(fft[1:n//2]))
axs[1].set_xlabel('Period [year]')
axs[1].set_ylabel('Signal amplitude')

axs[1].set_xscale('log')


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Time_Series/'
    
    file_name= 'TimeSeries_M'+str(Member) +'_' + str(syear)+'-'+str(eyear)+'_'
    file_name += var
    if var == 'T2M': file_name += '_Latbound'+str(lat)
    if var == 'vQtot': file_name += '_Lat'+str(lat)
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
