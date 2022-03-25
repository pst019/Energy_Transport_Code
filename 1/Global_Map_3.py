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
#save= True
imp= False
imp= True

Member = 3
#syear= 1950
#eyear= 1970
syear= 2080
eyear= 2100

#split_1= 2000 #bound last years
#split_2= 2047 #bound first years

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 4


#var, varfile= 'vEtotLLd', 'vEtotLL'
var, varfile= 'vQtotLLd', 'vQtotLL'


file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'



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
d_extreme =  np.max([d.max(), -d.min()])



syear_list=[syear, syear, split_2]
eyear_list=[eyear, split_1, eyear]


for iy in range(len(syear_list)):
    fig= plt.figure(fignr, (10,4))
    fignr+=1
    plt.clf()
    
    dsn= ds.where(np.logical_and(ds["time.year"] >= syear_list[iy], ds["time.year"] <= eyear_list[iy]  ), drop=True ).mean(dim='time')
    if var in ['vEtotLLd', 'vQtotLLd']:
    
 
        
        ax = fig.add_subplot(1, 3, 1, projection=ccrs.NorthPolarStereo())
        ax.set_title('Plan (WNr 4-)')
        ax.coastlines('110m', alpha=0.5)    
        ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
        
        c = ax.contourf(dsn.lon, dsn.lat, dsn[var].where(dsn['WaveNumb'] < 5, drop= True).sum(dim= 'WaveNumb'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
        

        ax = fig.add_subplot(1, 3, 2, projection=ccrs.NorthPolarStereo())
        ax.set_title('Syno (WNr 5+)')
        ax.coastlines('110m', alpha=0.5)
        ax.gridlines()
        ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
        
        c = ax.contourf(dsn.lon, dsn.lat, dsn[var].sel(WaveNumb= 5), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
        
        
        ax = fig.add_subplot(1, 3, 3, projection=ccrs.NorthPolarStereo())
        ax.set_title('Total')
        ax.coastlines('110m', alpha=0.5)    
        ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
        
        c = ax.contourf(dsn.lon, dsn.lat, dsn[var].sum(dim= 'WaveNumb'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme)

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
        cbar= fig.colorbar(c, cax=cbar_ax)        
        cbar.set_label(varfile)

   
    if save:
        savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
        
        file_name= 'Map_PS_M'+str(Member) +'_' + str(syear_list[iy])+'-'+str(eyear_list[iy])+'_'
        file_name += var
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 


""" difference"""

fig= plt.figure(fignr, (10,4))
fignr+=1
plt.clf()

dsn1= ds.where(np.logical_and(ds["time.year"] >= syear_list[-1], ds["time.year"] <= eyear_list[-1]  ), drop=True ).mean(dim='time')
dsn0= ds.where(np.logical_and(ds["time.year"] >= syear_list[-2], ds["time.year"] <= eyear_list[-2]  ), drop=True ).mean(dim='time')

dsn= dsn1-dsn0


d= dsn[var].sum(dim= 'WaveNumb')
d_extreme = np.max([d.max(), -d.min()])

if var in ['vEtotLLd', 'vQtotLLd']:

    
    ax = fig.add_subplot(1, 3, 1, projection=ccrs.NorthPolarStereo())
    ax.set_title('Plan (WNr 4-)')
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    c = ax.contourf(dsn.lon, dsn.lat, dsn[var].where(dsn['WaveNumb'] < 5, drop= True).sum(dim= 'WaveNumb'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
    

    ax = fig.add_subplot(1, 3, 2, projection=ccrs.NorthPolarStereo())
    ax.set_title('Syno (WNr 5+)')
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    c = ax.contourf(dsn.lon, dsn.lat, dsn[var].sel(WaveNumb= 5), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
    
    
    ax = fig.add_subplot(1, 3, 3, projection=ccrs.NorthPolarStereo())
    ax.set_title('total')
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    
    
    c = ax.contourf(dsn.lon, dsn.lat, dsn[var].sum(dim= 'WaveNumb'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme)
    
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
    cbar= fig.colorbar(c, cax=cbar_ax)
    cbar.set_label(varfile)
    
   
if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
        
    
    file_name= 'Map_PS_M'+str(Member) +'_diff' + str(syear_list[-1])+'-'+str(eyear_list[-1])+ 'mi' + str(syear_list[-2])+'-'+str(eyear_list[-2])+'_'
    file_name += var
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 