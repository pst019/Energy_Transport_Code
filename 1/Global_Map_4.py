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

fignr= 1


var, varfile= 'vEtotLLd', 'vEtotLL'
#var, varfile= 'vQtotLLd', 'vQtotLL'



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





"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)





for i in range(2):
    fig= plt.figure(fignr, (10,4))
    fignr+=1
    plt.clf()
    
#    dsn= ds.where(np.logical_and(ds["time.year"] >= syear_list[iy], ds["time.year"] <= eyear_list[iy]  ), drop=True ).mean(dim='time')
    if i == 0: dsn= ds.mean(dim='time')

#    if i == 0: dsn= ds.resample(time='1Y').mean().std(dim='time') #annual mean standard deviation

    if i == 1: 
        dsn1= ds.where(ds["time.year"] >= split_2, drop=True ).mean(dim='time')
        dsn0= ds.where(ds["time.year"] <= split_1, drop=True ).mean(dim='time')

        dsn= dsn1-dsn0
    
    
    d= dsn[var].sum(dim= 'WaveNumb')
    d_extreme =  np.max([d.max(), -d.min()])
    
    dsn['Plan'] = dsn[var].where(dsn['WaveNumb'] < 5, drop= True).sum(dim= 'WaveNumb')
    dsn['Syno'] = dsn[var].sel(WaveNumb= 5)
    dsn['Total'] = dsn[var].sum(dim= 'WaveNumb')
    catlist= ['Plan', 'Syno', 'Total']
    
    
    for ci, cat in enumerate(catlist):
        if var in ['vEtotLLd', 'vQtotLLd']:
        
            
            ax = fig.add_subplot(1, len(catlist), ci+1, projection=ccrs.NorthPolarStereo())
            
            ax.coastlines('110m', alpha=0.5)    
            ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
            ax.gridlines(ylocs=np.arange(0, 90, 10) )
            ax.set_boundary(circle, transform=ax.transAxes)
            
            c = ax.contourf(dsn.lon, dsn.lat, dsn[cat], transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
            ax.set_title(cat)
    

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
    cbar= fig.colorbar(c, cax=cbar_ax)        
    cbar.set_label(varfile)

   
    if save:
        savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
        
        file_name= 'Map_PS_M'+str(Member) +'_' 
        if i == 0: file_name += str(syear)+'-'+str(eyear)+'_'
        if i == 1: file_name += 'Diff' + str(split_2)+'-'+str(eyear)+ 'mi' + str(syear)+'-'+str(split_1)+'_'

        
        file_name += varfile
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 

