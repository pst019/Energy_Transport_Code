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

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save= False
save= True
imp= False
#imp= True

Member = 3
syear= 1950
#eyear= 1970
#syear= 2080
eyear= 2100

#split_1= 2000 #bound last years
#split_2= 2047 #bound first years

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 4

#var, varfile= 'Z', 'Z850'
#var, varfile= 'STR', 'STR'
#var, varfile= 'T2M', 'T2M'
#var, varfile= 'vEtotLLd', 'vEtotLL'
var, varfile= 'vQtotLLd', 'vQtotLL'


file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

fig= plt.figure(fignr)
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

#axs=fig.subplots(3, 1, sharex= 'col')

if imp:
    for year in range(syear, eyear+1):
        print(year)
        if year//10 == 0: print(year)
        for month in range(1,13):

            if year == syear and month== 1:
                if varfile in ['vEtotLL', 'vQtotLL']:
                    ds= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                    ds['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds.time), freq='1D').to_timestamp()
                else: ds= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                
                if varfile in ['Z850']: ds= ds.isel(plev= 0)
                ds= ds.resample(time='1M').mean()
            else:
                if varfile in ['vEtotLL', 'vQtotLL']:
                    ds2= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                    ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='1D').to_timestamp()

                else: ds2= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')

                if varfile in ['Z850']:  ds2= ds2.isel(plev= 0)
                ds2= ds2.resample(time='1M').mean()
#                 ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
                ds= xr.concat([ds, ds2], dim= 'time')
    
#    ds= ds.resample(time='1M').mean()
#    print('Finished import')

#    ds= ds.where(ds["time.year"] < 2097, drop= True)



d= ds[var].sum(dim= 'WaveNumb').mean(dim='time')
d_extreme = 0.5* np.max([d.max(), -d.min()])

#from matplotlib.colors import Normalize
#
#class MidpointNormalize(Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        Normalize.__init__(self, vmin, vmax, clip)
#
#    def __call__(self, value, clip=None):
#        # I'm ignoring masked values and all kinds of edge cases to make a
#        # simple example...
#        extreme= np.max([self.vmax, -self.vmin])
#        x, y = [-extreme, self.midpoint, extreme], [0, 0.5, 1]
##        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y))
#
#
#norm = MidpointNormalize(midpoint=0)


if var in ['vEtotLLd', 'vQtotLLd']:
    for WNr in range(6):
#        norm = MidpointNormalize(midpoint=0)

        ax = fig.add_subplot(3, 3, WNr+1, projection=ccrs.NorthPolarStereo())
        ax.set_title('WNr '+str(WNr))
        ax.coastlines('110m', alpha=0.5)    
        ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
        
        
        c = ax.contourf(ds.lon, ds.lat, ds[var].sel(WaveNumb= WNr).mean(dim='time'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
        fig.colorbar(c, ax=ax)
    
    
    ax = fig.add_subplot(3, 3, 7, projection=ccrs.NorthPolarStereo())
    ax.set_title('Total')
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    
    c = ax.contourf(ds.lon, ds.lat, ds[var].sum(dim= 'WaveNumb').mean(dim='time'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -2* d_extreme, vmax= 2*d_extreme)
    fig.colorbar(c, ax=ax)
    
    
else:
       
    ax = fig.add_subplot(2, 2, 1, projection=ccrs.NorthPolarStereo())
    ax.set_title(str(syear)+'-'+str(eyear))
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    
    c = ax.contourf(ds.lon, ds.lat, ds[var].mean(dim='time'), transform=ccrs.PlateCarree() )
    fig.colorbar(c, ax=ax)
    
    ax = fig.add_subplot(2, 2, 2, projection=ccrs.NorthPolarStereo())
    ax.set_title(str(split_2)+'-'+str(eyear))
    
    ax.set_global()
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    c = ax.contourf(ds.lon, ds.lat, ds[var].where(ds["time.year"] > split_2).mean(dim='time'), transform=ccrs.PlateCarree() )
    fig.colorbar(c, ax=ax)
    
    
    ax = fig.add_subplot(2, 2, 4, projection=ccrs.NorthPolarStereo())
    ax.set_title(str(syear)+'-'+str(split_1))
    
    ax.set_global()
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    c = ax.contourf(ds.lon, ds.lat, ds[var].where(ds["time.year"] < split_1).mean(dim='time'), transform=ccrs.PlateCarree() )
    fig.colorbar(c, ax=ax)
    
    
    ax = fig.add_subplot(2, 2, 3, projection=ccrs.NorthPolarStereo())
    ax.set_title('Last - First')
    
    ax.set_global()
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    c = ax.contourf(ds.lon, ds.lat, ds[var].where(ds["time.year"] > split_2).mean(dim='time')- ds[var].where(ds["time.year"] < split_1).mean(dim='time'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r' )
    fig.colorbar(c, ax=ax)
    
    
    #
    ## And black line contours.
    #line_c = ax.contour(x, y, z, levels=filled_c.levels,
    #                    colors=['black'],
    #                    transform=ccrs.PlateCarree())


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
    
    file_name= 'Map_M'+str(Member) +'_' + str(syear)+'-'+str(eyear)+'_'
    file_name += var
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
