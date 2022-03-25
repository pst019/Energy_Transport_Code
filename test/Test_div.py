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
syear= 1950
eyear= 1960
#syear= 2090
#eyear= 2100

#split_1= 2000 #bound last years
#split_2= 2047 #bound first years

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 1


#var, varfile= 'vEtotLLd', 'vEtotLL'
#var, varfile= 'vQtotLLd', 'vQtotLL'


varlist=[ 'uQtotLLd', 'vQtotLLd'] #important that the u component is first
varfilelist=[ 'uQtotLL', 'vQtotLL']

varlist=[ 'uEtotLLd', 'vEtotLLd'] #important that the u component is first
varfilelist=[ 'uEtotLL', 'vEtotLL']

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'



#ax1= plt.subplot(2, 1, 1)

#axs=fig.subplots(3, 1, sharex= 'col')

if imp:
    for vi, var in enumerate(varlist):
        print(var)
        varfile= varfilelist[vi]
    
        for year in range(syear, eyear+1):
            print(year)
            if year//10 == 0: print(year)
            for month in range(1,13):
    
                if year == syear and month== 1:
                    if varfile in ['vEtotLL', 'vQtotLL','uEtotLL', 'uQtotLL']:
                        ds0= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                        ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='1D').to_timestamp()
                    else: ds0= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                    
                    if varfile in ['Z850']: ds0= ds0.isel(plev= 0)
                    ds0= ds0.resample(time='1M').mean()
                else:
                    if varfile in ['vEtotLL', 'vQtotLL','uEtotLL', 'uQtotLL']:
                        ds1= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                        ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='1D').to_timestamp()
    
                    else: ds1= xr.open_dataset(file_dir + varfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
    
                    if varfile in ['Z850']:  ds1= ds1.isel(plev= 0)
                    
                    ds1= ds1.resample(time='1M').mean()
    #                 ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
                    ds0= xr.concat([ds0, ds1], dim= 'time')


        """the last loop step to merge the different variables in one ds"""
        if year == eyear and month == 12:
            if var== varlist[0]:
                ds2= ds0
            else:
                ds2= xr.merge([ds2, ds0])


ds= ds2





"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)




dvar= 'div' #the display variable


fig= plt.figure(fignr, (10,7))
fignr+=1
plt.clf()


dsn= ds.mean(dim='time').sum(dim= 'WaveNumb')
    
    ###calc the divergence:

dx= (float(dsn.lon[1])- float(dsn.lon[0]) ) *110E3 * np.cos(np.deg2rad(dsn.lat.values))
#dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[np.newaxis, :, np.newaxis] #axis 0 is latitude
dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[:, np.newaxis] #axis 1 is longitude


dy= (float(dsn.lat[1])- float(dsn.lat[0])) *110E3
dvdy= np.gradient(dsn[varlist[1]], axis= 0)/dy #axis 0 is latitude

dsn['div']= (dsn[varlist[0]].dims, dudx+ dvdy)

#    if i == 0: dsn= ds.resample(time='1Y').mean().std(dim='time') #annual mean standard deviation



#d= dsn[display_var].sum(dim= 'WaveNumb')
#d_extreme =  np.max([dsn[dvar].max(), -dsn[dvar].min()])





ax = fig.add_subplot(2, 3, 1, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)

c = ax.contourf(dsn.lon, dsn.lat, dsn[varlist[1]], transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.set_title(varlist[1])
fig.colorbar(c, ax=ax, shrink= 0.5)


ax = fig.add_subplot(2, 3, 2, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)

c = ax.contourf(dsn.lon, dsn.lat, dvdy, transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.set_title('d(vQ) dy')
fig.colorbar(c, ax=ax, shrink= 0.5)
   


ax = fig.add_subplot(2, 3, 4, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)

c = ax.contourf(dsn.lon, dsn.lat, dsn[varlist[0]], transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.set_title(varlist[0])
fig.colorbar(c, ax=ax, shrink= 0.5)


ax = fig.add_subplot(2, 3, 5, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)

c = ax.contourf(dsn.lon, dsn.lat, dudx, transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.set_title('d(uQ) dx')
fig.colorbar(c, ax=ax, shrink= 0.5)
   

ax = fig.add_subplot(2, 3, 6, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)

c = ax.contourf(dsn.lon, dsn.lat, dsn[dvar], transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.set_title('Divergence')
fig.colorbar(c, ax=ax, shrink= 0.5)

#        c = ax.contourf(dsn.lon, dsn.lat, dsn[display_var], transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)



#fig.subplots_adjust(right=0.85)
#cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
#cbar= fig.colorbar(c, cax=cbar_ax)        
#cbar.set_label(varfile)

   
#if save:
#    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
#    
#    file_name= 'Map_PS_M'+str(Member) +'_' 
#    if i == 0: file_name += str(syear)+'-'+str(eyear)+'_'
#    if i == 1: file_name += 'Diff' + str(split_2)+'-'+str(eyear)+ 'mi' + str(syear)+'-'+str(split_1)+'_'
#
#    
#    file_name += varfile
#    savefile= savedir+ file_name
#    print(savefile)
#    plt.savefig(savefile , bbox_inches='tight') 

