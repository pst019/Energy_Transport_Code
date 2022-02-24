#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import matplotlib.path as mpath
import pandas as pd
import matplotlib as mpl
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
user = os.getcwd().split('/')[2]

if user == 'pst019':
    Mediadir = '/media/'+user+'/Backup/'

else:
    Mediadir = '/run/media/pst019/Backup/'


#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save = False
# save = True
imp = False
imp = True

# syear = 2000
# eyear = 2014

selyear= 2007
running_mean_years= 7
syear = selyear - running_mean_years//2
eyear = selyear + running_mean_years//2

monthlist = np.arange(1, 13)


latextend= 70


fignr = 3

var = 'SIC'
varfile = 'ci'

file_dir = Mediadir + 'data/Energy_Transport/ERA5/' 


mvar= 'LWS'
mvarname= 'var177'



#ax1= plt.subplot(2, 1, 1)

#axs=fig.subplots(3, 1, sharex= 'col')

if imp:

    for year in range(syear, eyear+1):
        print(year)
        for month in monthlist:
            if year == syear and month == monthlist[0]:
                ds0 = xr.open_dataset(
                    file_dir + var + '/'+ 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc')#.resample(time='1M').mean()


                ds0 = ds0.where(ds0.lat >= latextend, drop= True)
                weights = np.cos(np.deg2rad(ds0.lat - 0.5))
                ds0 = ds0.weighted(weights).sum(("lon", "lat")) * 111**2
                
                
            else:
                ds1 = xr.open_dataset(
                    file_dir + var + '/'+ 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc')#.resample(time='1M').mean()
                
                ds1 = ds1.where(ds1.lat >= latextend, drop= True)
                ds1 = ds1.weighted(weights).sum(("lon", "lat")) * 111**2                
                
                
                ds0 = xr.concat([ds0, ds1], dim='time')

    ds0 = ds0.rename_vars({varfile: var})

    
    #to get monthly files
    for year in range(syear, eyear+1):
        if year == syear:
            mds0 = xr.open_dataset(file_dir + mvar + '_m/'+ mvar +'.'+str(year)+'.nc')

            # ds0 = ds0.where(ds0.lat >= latextend, drop= True)
            # weights = np.cos(np.deg2rad(ds0.lat - 0.5))
            # ds0 = ds0.weighted(weights).sum(("lon", "lat")) * 111**2
                    
        else:
            mds1 = xr.open_dataset(file_dir + mvar + '_m/'+ mvar +'.'+str(year)+'.nc')

                        # ds1 = ds1.where(ds1.lat >= latextend, drop= True)
            # ds1 = ds1.weighted(weights).sum(("lon", "lat")) * 111**2                          
            
            mds0 = xr.concat([mds0, mds1], dim='time')

    # ds0 = ds0.rename_vars({varfile: var})

    # to get hourly files
    # for year in range(syear, eyear+1):
    #     for month in monthlist:
    #         if year == syear and month == monthlist[0]:        
    #             dds0 = xr.open_dataset(
    #                 file_dir + mvar + '_d/'+ mvar +'.'+str(year)+'.'+ str(month).zfill(2)+'.nc')
    
    
                # ds0 = ds0.where(ds0.lat >= latextend, drop= True)
                # weights = np.cos(np.deg2rad(ds0.lat - 0.5))
                # ds0 = ds0.weighted(weights).sum(("lon", "lat")) * 111**2
                
                
    #         else:
    #             ds1 = xr.open_dataset(
    #                 file_dir + var + '/'+ 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc')#.resample(time='1M').mean()
                
    #             ds1 = ds1.where(ds1.lat >= latextend, drop= True)
    #             ds1 = ds1.weighted(weights).sum(("lon", "lat")) * 111**2                
            
            
    #             ds0 = xr.concat([ds0, ds1], dim='time')

    # ds0 = ds0.rename_vars({varfile: var})



ds = ds0

ds[var+'_mean'] = ds[var].rolling(time=30, center=True).mean()

ds = ds.where(ds["time.dayofyear"] != 366, drop= True) #remove schaltjahr




ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 366)),names=('year', 'dayofyear'))

#Finally, replace the time index in the Dataset by this new index, and then unstack its levels to create the two required dimensions.

ds = ds.assign(time=ind).unstack('time')

ds[var+'_mean2'] = ds[var+'_mean'].rolling(year=running_mean_years, center=True).mean()



# dsr = dsr.rename({'new_time':'time'})


# ds= ds.where(ds['ci']< 0.3, 1) # no

# assume that the grid cells are located at the upper corner
# weights = np.cos(np.deg2rad(ds.lat - 0.5))


# time_series = ds[var].weighted(weights).sum(
#     ("lon", "lat")) * 111**2  # assume grid cell spacing of 110 km
import matplotlib.dates as mdates


fig= plt.figure(fignr)
fignr += 1
plt.clf()
ax= fig.subplots()

timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# plt.plot( timearray, ds[var].sel(year= year) , label='Values')
    
plt.plot( timearray, ds[var+'_mean'].sel(year= selyear) , label=str(selyear) ) #30day running mean')

plt.plot(timearray, ds[var+'_mean2'].sel(year= selyear) , label='Climatology' )# (9 year running mean)')

plt.ylabel('Sea ice area [km$^2$]')
plt.xlabel('Time')

plt.legend()
plt.gca().set_ylim(bottom=0)

ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))


fig= plt.figure(fignr)
fignr += 1
plt.clf()
ax= fig.subplots()

timearray= (np.asarray(year, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# plt.plot( timearray, ds[var].sel(year= selyear) , label='Values')
    
plt.plot(timearray, ds[var+'_mean'].sel(year= selyear)- ds[var+'_mean2'].sel(year= selyear) , label=var +'anomaly' )# (9 year running mean)')

plt.ylabel('Sea ice area [km$^2$]')
plt.xlabel('Time')

plt.legend()

ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))


# plt.xlim(ds['time'][0], ds['time'][-1])
# make tikes every year

# print('it is somewhat off')


if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Evolution_' + varfile + '_' +str(year)

    savefile= savedir+ file_name 
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight')



# """to make the map a circle"""
# theta = np.linspace(0, 2*np.pi, 100)
# center, radius = [0.5, 0.5], 0.5
# verts = np.vstack([np.sin(theta), np.cos(theta)]).T
# circle = mpath.Path(verts * radius + center)


# fig= plt.figure(fignr)
# fignr += 1
# plt.clf()

# latextend= 70

# ds = ds.where(ds.lat > latextend, drop= True)



# levels= np.arange(0, 1.01, 0.1)
# date = str(year) + "-09-30"

# ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

# ax.coastlines('110m', alpha=0.5)
# ax.set_extent([-180, 180, latextend+1, 90], crs=ccrs.PlateCarree())
# ax.gridlines(ylocs=np.arange(0, 90, 10) )
# ax.set_boundary(circle, transform=ax.transAxes)


# norm = mpl.colors.BoundaryNorm(boundaries=levels, ncolors=256)#, clip=True) #extend= 'both')
# c = ax.contourf(ds.lon, ds.lat, ds[var].sel(time = date) 
#                 , transform=ccrs.PlateCarree(), cmap= 'Blues_r', levels= levels, norm=norm) #norm= norm)


# fig.subplots_adjust(right=0.85)
# cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
# cbar= fig.colorbar(c, cax=cbar_ax)
# # cbar.set_label(var)
# cbar.set_label('Sea-ice cover')


# if save:
#     savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'

#     file_name = 'Map_' + varfile + '_Sep_'+str(year)

#     savefile= savedir+ file_name
#     print(savefile)
#     plt.savefig(savefile , bbox_inches='tight')





# fig= plt.figure(fignr)
# fignr += 1
# plt.clf()


# ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

# ax.coastlines('110m', alpha=0.5)
# ax.set_extent([-180, 180, latextend+1, 90], crs=ccrs.PlateCarree())
# ax.gridlines(ylocs=np.arange(0, 90, 10) )
# ax.set_boundary(circle, transform=ax.transAxes)


# norm = mpl.colors.BoundaryNorm(boundaries=levels, ncolors=256)#, clip=True) #extend= 'both')
# c = ax.contourf(ds.lon, ds.lat, ds[var].rolling(time=9, center=True).mean().sel(time = date) 
#                 , transform=ccrs.PlateCarree(), cmap= 'Blues_r', levels= levels, norm=norm) #norm= norm)


# fig.subplots_adjust(right=0.85)
# cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
# cbar= fig.colorbar(c, cax=cbar_ax)
# # cbar.set_label(var)
# cbar.set_label('Sea-ice cover')

# if save:
#     savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'

#     file_name = 'Map_' + varfile + 'running-mean-state_Sep_'+str(year)

#     savefile= savedir+ file_name
#     print(savefile)
#     plt.savefig(savefile , bbox_inches='tight')



# fig= plt.figure(fignr)
# fignr += 1
# plt.clf()


# ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

# ax.coastlines('110m', alpha=0.5)
# ax.set_extent([-180, 180, latextend+1, 90], crs=ccrs.PlateCarree())
# ax.gridlines(ylocs=np.arange(0, 90, 10) )
# ax.set_boundary(circle, transform=ax.transAxes)

# levels= np.arange(-1, 1.01, 0.1)

# norm = mpl.colors.BoundaryNorm(boundaries=levels, ncolors=256)#, clip=True) #extend= 'both')
# c = ax.contourf(ds.lon, ds.lat, ds[var].sel(time = date) - ds[var].rolling(time=9, center=True).mean().sel(time = date) 
#                 , transform=ccrs.PlateCarree(), cmap= 'RdBu', levels= levels, norm=norm) #norm= norm)


# fig.subplots_adjust(right=0.85)
# cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
# cbar= fig.colorbar(c, cax=cbar_ax)
# # cbar.set_label(var)
# cbar.set_label('Sea-ice cover anomaly')

# if save:
#     savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'

#     file_name = 'Map_' + varfile + 'anomaly_Sep_'+str(year)

#     savefile= savedir+ file_name
#     print(savefile)
#     plt.savefig(savefile , bbox_inches='tight')



