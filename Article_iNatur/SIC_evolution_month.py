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

selyear= 2012
running_mean_years= 7
syear = selyear - running_mean_years//2
eyear = selyear + running_mean_years//2

monthlist = np.arange(1, 13)


latextend= 70


fignr = 1

var = 'SIC'
varname = 'ci'

file_dir = Mediadir + 'data/Energy_Transport/ERA5/' 


mvar_list= ['SWS', 'SWSD', 'LWS', 'LWSD', 'SSH', 'SLH']
mvarname_list= ['var176', 'var169', 'var177', 'var175', 'var146', 'var147']

mvar= 'SWS'
mvarname= 'var176'

mvar= 'SSH'
mvarname= 'var146'

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

    ds0 = ds0.rename_vars({varname: var})

    
    #to get monthly files
    for mi, mvar in enumerate(mvar_list):
        print(mvar)
        mvarname= mvarname_list[mi]
        
        for year in range(syear, eyear+1):
            if year == syear:
                mds0 = xr.open_dataset(file_dir + mvar + '_m/'+ mvar +'.'+str(year)+'.nc')
                mds0 = mds0.isel(time= np.arange(0, len(mds0.time), 2) )# to select every second time step, since the months are double
                
                # if selyear== 2007: mds0 = mds0.where(np.logical_and(mds0.lat >= 75, mds0.lat <= 85) & np.logical_or(mds0.lon >= 120, mds0.lon <= -150) , drop= True)
                # elif selyear== 2012: mds0 = mds0.where(np.logical_and(mds0.lat >= 75, mds0.lat <= 85) & np.logical_or(mds0.lon >= 0, mds0.lon <= -120) , drop= True)
                # else: 
                mds0 = mds0.where(mds0.lat >= latextend, drop= True)

                weights = np.cos(np.deg2rad(mds0.lat - 0.5))
                mds0 = mds0.weighted(weights).mean(("lon", "lat"))/ (24*60**2)
                # ds0 = ds0.where(ds0.lat >= latextend, drop= True)
                # weights = np.cos(np.deg2rad(ds0.lat - 0.5))
                # ds0 = ds0.weighted(weights).sum(("lon", "lat")) * 111**2
                        
            else:
                mds1 = xr.open_dataset(file_dir + mvar + '_m/'+ mvar +'.'+str(year)+'.nc')
                mds1 = mds1.isel(time= np.arange(0, len(mds1.time), 2) )
    
                # if selyear== 2007: mds1= mds1.where(np.logical_and(mds1.lat >= 75, mds1.lat <= 85) & np.logical_or(mds1.lon >= 120, mds1.lon <= -150) , drop= True)
                # elif selyear== 2012: mds1 = mds1.where(np.logical_and(mds1.lat >= 75, mds1.lat <= 85) & np.logical_or(mds1.lon >= 0, mds1.lon <= -120) , drop= True)
                # else: 
                mds1 = mds1.where(mds1.lat >= latextend, drop= True)

                mds1 = mds1.weighted(weights).mean(("lon", "lat"))/ (24*60**2)
                
                # ds1 = ds1.where(ds1.lat >= latextend, drop= True)
                # ds1 = ds1.weighted(weights).sum(("lon", "lat")) * 111**2                          
                
                mds0 = xr.concat([mds0, mds1], dim='time')
    
        mds0 = mds0.rename_vars({mvarname: mvar})

        """the last loop step to merge the different variables in one ds"""
        if year == eyear: # and month == 12:
            if mvar== mvar_list[0]:
                mds2= mds0
            else:
                mds2= xr.merge([mds2, mds0])
    
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

ds[var+'_smth'] = ds[var].rolling(time=30, center=True).mean()
ds = ds.where(ds["time.dayofyear"] != 366, drop= True) #remove schaltjahr

ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 366)),names=('year', 'dayofyear'))
ds = ds.assign(time=ind).unstack('time')

ds[var+'_clim'] = ds[var+'_smth'].rolling(year=running_mean_years, center=True).mean()


mds = mds2

# mds = mds.where(mds.lat >= latextend, drop= True)
# weights = np.cos(np.deg2rad(mds.lat - 0.5))
# mds = mds.weighted(weights).mean(("lon", "lat"))/ (24*60**2)

ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 13)),names=('year', 'month'))
mds = mds.assign(time=ind).unstack('time')

for mvar in mvar_list:
    mds[mvar+'_clim'] = mds[mvar].rolling(year=running_mean_years, center=True).mean()

# dsr = dsr.rename({'new_time':'time'})


# ds= ds.where(ds['ci']< 0.3, 1) # no

# assume that the grid cells are located at the upper corner
# weights = np.cos(np.deg2rad(ds.lat - 0.5))


# time_series = ds[var].weighted(weights).sum(
#     ("lon", "lat")) * 111**2  # assume grid cell spacing of 110 km
import matplotlib.dates as mdates


fig= plt.figure(fignr, (5, 10))
fignr += 1
plt.clf()
axs= fig.subplots(1+ len(mvar_list),1)

timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

    
axs[0].plot( timearray, ds[var+'_smth'].sel(year= selyear) , label=str(selyear) ) #30day running mean')

axs[0].plot(timearray, ds[var+'_clim'].sel(year= selyear) , label='Climatology' )# (9 year running mean)')

axs[0].set_ylabel('Sea ice area [km$^2$]')
# axs[0].set_xlabel('Time')

axs[0].legend()

# axs[0].gca().set_ylim(bottom=0)

# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

for mi, mvar in enumerate(mvar_list):
    axs[mi+1].plot(mds.month, mds[mvar].sel(year= selyear) , label=str(selyear))
    axs[mi+1].plot(mds.month, mds[mvar+'_clim'].sel(year= selyear) , label='Climatology')
    axs[mi+1].set_ylabel(mvar)


fig= plt.figure(fignr, (5, 10))
fignr += 1
plt.clf()
axs= fig.subplots(1+ len(mvar_list),1)

# timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# axs[0].plot( timearray, ds[var].sel(year= selyear) , label='Values')
    
axs[0].plot(timearray, ds[var+'_smth'].sel(year= selyear)- ds[var+'_clim'].sel(year= selyear) , label=var +' anomaly' )# (9 year running mean)')

axs[0].set_ylabel('Sea ice area [km$^2$]')
# axs[0].set_xlabel('Time')

axs[0].legend()

# axs[1].plot(mds.month, mds[mvar].sel(year= selyear) - mds[mvar+'_clim'].sel(year= selyear) , label='anomaly')

# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))


# axs[0].xlim(ds['time'][0], ds['time'][-1])
# make ticks every year

# print('it is somewhat off')


for mi, mvar in enumerate(mvar_list):
    axs[mi+1].plot(mds.month,mds[mvar].sel(year= selyear) - mds[mvar+'_clim'].sel(year= selyear) , label='anomaly' )
    axs[mi+1].set_ylabel(mvar)



fig= plt.figure(fignr, (6,6))
fignr += 1
plt.clf()
axs= fig.subplots(3,1)

# timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# axs[0].plot( timearray, ds[var].sel(year= selyear) , label='Values')
    
axs[0].plot(timearray, ds[var+'_smth'].sel(year= selyear)- ds[var+'_clim'].sel(year= selyear) , label=var +' anomaly' )# (9 year running mean)')

axs[0].set_ylabel('Sea ice area [km$^2$]')
# axs[0].set_xlabel('Time')

axs[0].legend()
# axs[0].plot([1,12], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


axs[1].plot(mds.month, mds['LWS'].sel(year= selyear) - mds['LWS_clim'].sel(year= selyear)
            + mds['SSH'].sel(year= selyear) - mds['SSH_clim'].sel(year= selyear) 
            + mds['SLH'].sel(year= selyear) - mds['SLH_clim'].sel(year= selyear) 
            , label='LWS + SSH + SLH' )
axs[1].plot(mds.month, mds['SWS'].sel(year= selyear) - mds['SWS_clim'].sel(year= selyear), label='SWS' )
axs[1].set_ylabel('Flux anomaly [W m$^{-2}$]')
axs[1].legend(ncol= 2)
axs[1].plot([1,12], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


axs[2].plot(mds.month, mds['LWSD'].sel(year= selyear) - mds['LWSD_clim'].sel(year= selyear), label='LWSD' )
axs[2].plot(mds.month, mds['SWSD'].sel(year= selyear) - mds['SWSD_clim'].sel(year= selyear), label='SWSD' )
axs[2].plot(mds.month, mds['SSH'].sel(year= selyear) - mds['SSH_clim'].sel(year= selyear), label='SSH' )
axs[2].plot(mds.month, mds['SLH'].sel(year= selyear) - mds['SLH_clim'].sel(year= selyear), label='SlH' )
axs[2].set_ylabel('Flux anomaly [W m$^{-2}$]')
axs[2].legend(ncol= 4)
axs[2].plot([1,12], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[2].set_xlabel('Month')    


if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Evolution_' + var + 'heating_' +str(selyear)

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



