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
save = True
imp = False
imp = True

syear = 1979
eyear = 2020

monthlist = [9]


fignr = 5

var = 'ci'
varfile = 'SIC'

file_dir = Mediadir + 'data/Energy_Transport/ERA5/' + varfile + '/'


#ax1= plt.subplot(2, 1, 1)

#axs=fig.subplots(3, 1, sharex= 'col')

if imp:

    for year in range(syear, eyear+1):
        print(year)
        for month in monthlist:
            if year == syear and month == monthlist[0]:
                ds0 = xr.open_dataset(
                    file_dir + 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc').resample(time='1M').mean()

            else:
                ds1 = xr.open_dataset(
                    file_dir + 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc').resample(time='1M').mean()
                ds0 = xr.concat([ds0, ds1], dim='time')


ds = ds0

# ds= ds.where(ds['ci']< 0.3, 1) # no

ds= ds.assign_coords(lon=ds.lon + 0.5)
# ds= ds.assign_coords(lat=ds.lat + 0.5)

# assume that the grid cells are located at the upper corner
weights = np.cos(np.deg2rad(ds.lat- 0.5))

time_series = ds[var].weighted(weights).sum(
    ("lon", "lat")) * 111**2  # assume grid cell spacing of 110 km


plt.figure(fignr)
fignr += 1
plt.clf()
plt.plot(time_series.time.dt.year, time_series.values/1E6, c= 'k', label='September')
plt.ylabel('Arktisk sjøis [Millioner km$^2$]')
plt.xlabel('År')


plt.plot(time_series.time.dt.year, time_series.rolling(
    time=9, center=True).mean()/1E6, '--', c= 'k', label='September klimatologi')

plt.legend()
plt.gca().set_ylim(bottom=0)
plt.xlim(ds['time'].dt.year[0], ds['time'].dt.year[-1])
# make tikes every year

print('it is somewhat off')


if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Timeseries_' + varfile

    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 300)



"""to make the map a circle"""
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


fig= plt.figure(fignr)
fignr += 1
plt.clf()

latextend= 68

ds = ds.where(ds.lat > latextend, drop= True)



levels= np.arange(0, 1.01, 0.1)
year = 2012
date = str(year) + "-09-30"

ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)
ax.set_extent([-180, 180, latextend+1, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)


norm = mpl.colors.BoundaryNorm(boundaries=levels, ncolors=256)#, clip=True) #extend= 'both')
c = ax.contourf(ds.lon, ds.lat, ds[var].sel(time = date) 
                , transform=ccrs.PlateCarree(), cmap= 'Blues_r', levels= levels, norm=norm) #norm= norm)

ax.contour(ds.lon, ds.lat, ds[var].sel(time = date) 
                , transform=ccrs.PlateCarree(), colors= 'r', levels= [0.15], linestyles= 'dashed')

ax.contour(ds.lon, ds.lat, ds[var].rolling(time=9, center=True).mean().sel(time = date) 
                , transform=ccrs.PlateCarree(), colors= 'r', levels= [0.15])

fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
cbar= fig.colorbar(c, cax=cbar_ax)
# cbar.set_label(var)
cbar.set_label('Sjøis konsentrasjon')


if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'

    file_name = 'Map_' + varfile + '_Sep_'+str(year)

    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 300)





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



