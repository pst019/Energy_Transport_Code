#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    # Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/EC_Earth'
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/ERA5'

elif user=='media':
    # Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/EC_Earth'
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/ERA5'

else:
    Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd
import cartopy.crs as ccrs
import somoclu



save= True
# save= False





syear= 2010
eyear= 2010



fignr= 1





var= 'vQtot'
latcircle= 70
catlist = ['Syno', 'Plan', 'total']

clevels= []

#cvar, cvarfile, clabel, cunit= 'T2M', 'T2M', 'Temperature', 'K'
cvar, clabel, cunit, clevels= 'z', 'Geopotential height', 'm', np.arange(-200, 201, 40)


latbound= 50 #for the import

SOM_latlow = 65
SOM_lathigh = 75

# monthlist= [9, 10, 11] #[6, 7, 8] #[3, 4, 5] #[12, 1, 2] #
monthlist= np.arange(1,13)
nrows = 4
ncols = 5


print('Import ...')




cds= xr.open_dataset(Mediadir +f'/Johanne/1979-2018.30-90N.res0.5.daymean.Z850.nc')

cds= cds.isel(plev= 0)
cds= cds.sel(lat = cds.lat[cds.lat > latbound][::2])
cds= cds.sel(lon = cds.lon[::2])

cds['z'] /= 9.81



#reduce the time cds= cds[]
cds= cds.sel( time = cds.time[(cds.time.dt.year >= syear) & (cds.time.dt.year < eyear+1)] )




varlist=[ 'uQtotLL', 'vQtotLL']

# varlist=[ 'uEtotLL', 'vEtotLL']

# file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'
file_dir= Mediadir+ '/EnergySplit/Waves/'



#ax1= plt.subplot(2, 1, 1)

#axs=fig.subplots(3, 1, sharex= 'col')

for vi, var in enumerate(varlist):
    print(var)
    varfile= varlist[vi]

    for year in range(syear, eyear+1):
        print(year)
        for month in  monthlist:
            print(month)
            if year == syear and month== monthlist[0]:
                ds0= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()
                
                ds0= ds0.sel(lat = ds0.lat[ds0.lat > latbound])
                ds0= ds0.resample(time='1D').mean()
            else:
                ds1= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()

                ds1= ds1.sel(lat = ds1.lat[ds1.lat > latbound])                
                ds1= ds1.resample(time='1D').mean()
#                 ds2['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds2.time), freq='6H').to_timestamp()
                ds0= xr.concat([ds0, ds1], dim= 'time')
                # print(ds0)

    """the last loop step to merge the different variables in one ds"""
    if year == eyear and month ==monthlist[-1]:
        if var== varlist[0]:
            ds2= ds0
        else:
            ds2= xr.merge([ds2, ds0])


ds= ds2



dx= (float(ds.lon[1])- float(ds.lon[0]) ) *110E3 * np.cos(np.deg2rad(ds.lat.values))
#dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[np.newaxis, :, np.newaxis] #axis 0 is latitude
dudx= np.gradient(ds[varlist[0]], axis= 1)/dx[:, np.newaxis] #axis 1 is longitude


dy= (float(ds.lat[1])- float(ds.lat[0])) *110E3
dvdy= np.gradient(ds[varlist[1]], axis= 0)/dy #axis 0 is latitude


ds['dudx']= (ds[varlist[0]].dims, dudx)
ds['dvdy']= (ds[varlist[0]].dims, dvdy)


ds['div']= (ds[varlist[0]].dims, dudx+ dvdy)


# var= 'div'



"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)






    # """energy transport"""
    # if catlist == ['Syno', 'Plan', 'total']:
    
    #     ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 3), drop= True ).sum(dim='WaveNumb')
    #     ds['Syno']= ds[var].where(ds.WaveNumb >= 4, drop= True ).sum(dim='WaveNumb')
    
    
    
    
ds['total']= ds[var].sum(dim='WaveNumb')

# ds= ds.drop(var)

    
ds= ds.sel(lat = ds.lat[ds.lat<= 85])    


# for cat in catlist[-1:]:
    
fig= plt.figure(fignr, figsize= (10,7) )
fignr+=1
plt.clf()
            
    #fig, axs = plt.subplots(nrows,ncols, figsize = (10,10))
    #d = np.reshape(d,(dim_time, dim_x, dim_y))


time_step= -1

""""z 850"""
ax = fig.add_subplot(2, 3, 1, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    


c = ax.contourf(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('z 850')
fig.colorbar(c, ax=ax, shrink= 0.5)

 
"""uQtotLL"""
ax = fig.add_subplot(2, 3, 2, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



ds_field= ds.uQtotLL.mean(dim='time').sum(dim='WaveNumb')
v_extr= np.max([ds_field.max(), -1* ds_field.min()])
c = ax.contourf(ds.lon, ds.lat, ds_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('uQtot')
fig.colorbar(c, ax=ax, shrink= 0.5)


"""vQtotLL"""
ax = fig.add_subplot(2, 3, 3, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



ds_field= ds.vQtotLL.mean(dim='time').sum(dim='WaveNumb')
v_extr= np.max([ds_field.max(), -1* ds_field.min()])
c = ax.contourf(ds.lon, ds.lat, ds_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('vQtot')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""d vQ dy"""
ax = fig.add_subplot(2, 3, 6, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



ds_field=-1* ds.dvdy.mean(dim='time').sum(dim='WaveNumb')
v_extr= np.max([ds_field.max(), -1* ds_field.min()])
c = ax.contourf(ds.lon, ds.lat, ds_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('-d(vQ) dy')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""d uQ dx"""
ax = fig.add_subplot(2, 3, 5, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



ds_field= -1* ds.dudx.mean(dim='time').sum(dim='WaveNumb')
v_extr= np.max([ds_field.max(), -1* ds_field.min()])
c = ax.contourf(ds.lon, ds.lat, ds_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('-d(uQ) dx')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""conv Q"""
ax = fig.add_subplot(2, 3, 4, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



ds_field= -1* ds.div.mean(dim='time').sum(dim='WaveNumb')
v_extr= np.max([ds_field.max(), -1* ds_field.min()])
c = ax.contourf(ds.lon, ds.lat, ds_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('conv')
fig.colorbar(c, ax=ax, shrink= 0.5)



print('have to get the units right')