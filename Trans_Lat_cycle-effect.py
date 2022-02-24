#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup/EC_Earth'

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



save= True
#save= False
imp= False
imp= True
#
#Member = 1
syear= 2056 #1950
eyear= 2096 # 1960

#syear= 1950
#eyear= 2000 #70

#split_1 = int(2/3*syear + 1/3*eyear)
#split_2 = int(2/3*eyear + 1/3*syear)

fignr= 4


#Memberlist= [3]#,2,3,4]
Member= 3







var= 'vQtot'
latcircle= 70
catlist = ['Syno', 'Plan', 'total']
perc= 90


cvar, cvarfile, clabel, cunit= 'T2M', 'T2M', 'Temperature', 'K'
#cvar, cvarfile, clabel, cunit= 'Z', 'Z850', 'Geopotential height', 'm'

latbound= 50




if imp:
    print('Import ...')
#    for Member in Memberlist:
    print('Member', Member)
    file_dir= Mediadir + '/Member'+str(Member) +'/'
    
#    for var in varlist:
    print(var)
    for year in range(syear, eyear+1):
        if year%10 == 0: print(year)
        for month in range(1,13):
#                print(month)
            if year == syear and month== 1:
                

                ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()
                ds0= ds0.sel(lat= 70, method= 'nearest')
                ds0= ds0.resample(time='1D').mean()
                    
#                if Member == Memberlist[0]: 
                dslat= ds0.lat
#                else: ds0['lat']= dslat
                
                
                cds0= xr.open_dataset(file_dir + cvarfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                cds0 = cds0.where(cds0.lat > latbound, drop= True)
                cds0= cds0.resample(time='1D').mean()
                if cvar in ['Z']:
                    cds0= cds0.isel(plev= 0)
                    cds0[cvar]/= 9.81

            else:
                ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()
                ds1= ds1.sel(lat= latcircle, method= 'nearest')
                ds1= ds1.resample(time='1D').mean()
                                       
                ds1['lat']= dslat
                
                ds0= xr.concat([ds0, ds1], dim= 'time')


                cds1= xr.open_dataset(file_dir + cvarfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                cds1 = cds1.where(cds0.lat > latbound, drop= True)
                cds1= cds1.resample(time='1D').mean()
                if cvar in ['Z']: 
                    cds1= cds1.isel(plev= 0)
                    cds1[cvar]/= 9.81

                
                cds0= xr.concat([cds0, cds1], dim= 'time')



                
            """the last loop step to merge the different variables in one ds"""
            if year == eyear and month == 12:
#                print(ds)
#                    if var== varlist[0]:
                ds2= ds0
                cds2= cds0

#                    else:
#                        ds2= xr.merge([ds2, ds0])
#                        cds2= xr.merge([cds2, cds0])
                            
#        if Member == Memberlist[0]:
#            ds= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
#        else:
#            ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
#            ds= xr.merge([ds, ds2])



"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


"""computation of the categories"""
print('compute')

ds= ds2
cds= cds2


nyears = eyear - syear + 1
ds = ds.sel(time= ds["time.dayofyear"] != 366)
cds = cds.sel(time= cds["time.dayofyear"] != 366)




if catlist == ['Syno', 'Plan', 'total']:

    ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb')
    ###    ds['>5']= ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb')
    ds['Syno']= ds[var].where(ds.WaveNumb >= 5, drop= True ).sum(dim='WaveNumb')
    ds['total']= ds[var].sum(dim='WaveNumb')
    
    ds= ds.drop(var)

#ds= ds.sel(time=~((ds.time.dt.month == 2) & (ds.time.dt.day == 29))) #remove 29 feb such that every year has same amount of days
#a running mean


"""computation of the deseasonalized data"""
for cat in catlist:
    dr =ds[cat].rolling(time=31, center=True).mean().dropna("time")
    dr = dr.groupby("time.dayofyear").mean()
    dr =np.tile(dr, nyears)
    #dr =np.tile(dr, (nyears,1))
    
    
    cata= cat+'_ano'
    ds[cata]= ds[cat] - dr
     

 
cdr =cds[cvar].rolling(time=31, center=True).mean().dropna("time")
cdr = cdr.groupby("time.dayofyear").mean()
cdr =np.tile(cdr, (nyears,1,1) )

cvara= cvar+'_ano'
cds[cvara]= cds[cvar] - cdr




"""time series"""
fignr= 4
fig= plt.figure(fignr, (12, 7))
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

axs=fig.subplots(4, 2, sharex= 'col')


for ic, cat in enumerate(catlist):
    cata= cat+'_ano'
    
    
    axs[ic, 0].plot(ds['time'], ds[cat], 'o', markersize= 1)
    axs[ic, 0].set_ylabel(var + ' ' + cat)
    #axs[0].plot(ds['time'], ds[cat].sel(WaveNumb= 1))
    
    axs[ic,1].plot(ds['time'], ds[cat+'_ano'], 'o', markersize= 1)
    axs[ic,1].set_ylabel(var + ' ' + cat+'_ano')


    hthresh= np.percentile(ds[cata], perc)
    lthresh= np.percentile(ds[cata], (100-perc) )
    
    
    htimes= ds['time'].where(ds[cata] > hthresh, drop= True) #high times
    ltimes= ds['time'].where(ds[cata] < lthresh, drop= True) #high times
    
    axs[ic, 1].plot(htimes, ds[cata].sel(time= htimes), 'o', markersize= 1)
    axs[ic, 1].plot(ltimes, ds[cata].sel(time= ltimes), 'o', markersize= 1)#, zorder= 1)
#    axs[ic,1].plot(ds['time'].isel(time= [0, -1]), 2*[hthresh])
#    axs[ic,1].plot(ds['time'].isel(time= [0, -1]), 2*[lthresh])

    axs[ic, 0].plot(htimes, ds[cat].sel(time= htimes), 'o', markersize= 1)#, zorder= 1)
    axs[ic, 0].plot(ltimes, ds[cat].sel(time= ltimes), 'o', markersize= 1)#, zorder= 1)

axs[-1, 0].plot(cds['time'], cds[cvar].sel(lat= 70, method='nearest').mean(dim='lon'), 'o', markersize= 1 )
axs[-1, 0].set_ylabel(cvar+  ' ['+cunit+']')
#axs[0].plot(ds['time'], ds[cat].sel(WaveNumb= 1))

axs[-1,1].plot(cds['time'], cds[cvara].sel(lat= 70, method='nearest').mean(dim='lon'), 'o', markersize= 1, zorder=1)
axs[-1,1].set_ylabel(cvara+  ' ['+cunit+']')

#if cat =='total': #this is for the last cat
axs[-1, 0].plot(htimes, cds[cvar].sel(lat= 70, method='nearest').mean(dim='lon').sel(time= htimes), 'o', markersize= 1, zorder=2 )
axs[-1, 1].plot(htimes, cds[cvara].sel(lat= 70, method='nearest').mean(dim='lon').sel(time= htimes), 'o', markersize= 1, zorder=2 )

axs[-1, 0].plot(ltimes, cds[cvar].sel(lat= 70, method='nearest').mean(dim='lon').sel(time= ltimes), 'o', markersize= 1, zorder=2 )
axs[-1, 1].plot(ltimes, cds[cvara].sel(lat= 70, method='nearest').mean(dim='lon').sel(time= ltimes), 'o', markersize= 1, zorder=2 )


axs[-1,1].set_xlim(ds['time'].isel(time= [0, -1]))




#axs[1].scatter(ltimes, len(ltimes)*[lthresh], zorder= 1)


for cat in catlist:
    cata= cat+'_ano'

    hthresh= np.percentile(ds[cata], perc)
    lthresh= np.percentile(ds[cata], (100-perc) )
    
    htimes= ds['time'].where(ds[cata] > hthresh, drop= True) #high times
    ltimes= ds['time'].where(ds[cata] < lthresh, drop= True) #high times



    """lat - leadtime plot"""
    lackday= np.arange(-20, 21, 2)
    
    for it, times in enumerate([htimes, ltimes]):
        if it == 0: timeslabel= 'high'
        if it == 1: timeslabel= 'low'
    
        lds =xr.DataArray([cds[cvara].sel(time= times.values+ np.timedelta64(ld, 'D'), method= 'nearest').mean(dim=['time', 'lon']).values for ld in lackday], coords=[lackday, cds.lat], dims=['lackday', 'lat'])
        
        fig= plt.figure(fignr)
        fignr+=1
        plt.clf()
        
        extr= np.max([lds.max(), -lds.min()])
        c = plt.contourf(lds.lat, lds.lackday, lds, cmap= 'RdBu_r', vmin= -extr, vmax= extr)
        cbar= plt.colorbar(c)        
        cbar.set_label(clabel +' anomaly ['+cunit+']')
        
        plt.xlabel('Latitude')
        plt.ylabel('Time lag [days]')
        
        if save:
            savedir= '../Figs/Trans-lat-effect/'
            if not os.path.exists(savedir): os.makedirs(savedir)
            
            file_name= 'LeadLag_M' 
        #        for Member in Memberlist:
            file_name += str(Member)
            file_name += '_' +str(syear)+'-'+str(eyear)+'_'
            file_name += var + str(latcircle) + '-'+ cat + '-'+ timeslabel
            file_name += '-'+ cvar
            savefile= savedir+ file_name
            print(savefile)
            plt.savefig(savefile , bbox_inches='tight') 
    
    
    """map of temperature"""
    lacktimes= [-5, 0, 5, 10, 15, 20]
    
    
    for it, times in enumerate([htimes, ltimes]):
        if it == 0: timeslabel= 'high'
        if it == 1: timeslabel= 'low'
    
        fig= plt.figure(fignr, (10, 7))
        fignr+=1
        plt.clf()
        
        for li, ld in enumerate(lacktimes):
            ax = fig.add_subplot(2, 3, li +1, projection=ccrs.NorthPolarStereo())
            #ax.set_title('WNr '+str(WNr))
            ax.coastlines('110m', alpha=0.5)    
            ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
            ax.gridlines(ylocs=np.arange(0, 90, 10) )
            ax.set_boundary(circle, transform=ax.transAxes)
            
            levels= np.round(np.linspace(-extr, extr, 10), 2)
            c = ax.contourf(cds.lon, cds.lat, cds[cvara].sel(time= times.values+ np.timedelta64(ld, 'D'), method= 'nearest').mean(dim='time')
            , transform=ccrs.PlateCarree(), cmap= 'RdBu_r', levels= levels, extend= 'both' ) #vmin= -1.5* extr, vmax= 1.5* extr) #, vmin= -d_extreme, vmax= d_extreme) #norm= norm)
            ax.set_title(str(ld) +' days')
        
        
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
        cbar= fig.colorbar(c, cax=cbar_ax)        
        cbar.set_label(clabel+' anomaly ['+cunit+']')
    
    
        if save:
            savedir= '../Figs/Trans-lat-effect/'
            if not os.path.exists(savedir): os.makedirs(savedir)
            
            file_name= 'Map_M' 
    #        for Member in Memberlist:
            file_name += str(Member)
            file_name += '_' +str(syear)+'-'+str(eyear)+'_'
            file_name += var + str(latcircle) + '-'+ cat +'-'+ timeslabel
            file_name += '-'+ cvar
            savefile= savedir+ file_name
            print(savefile)
            plt.savefig(savefile , bbox_inches='tight') 
    

        fig.suptitle(var + str(latcircle) +' '+ cat +' '+ timeslabel +' '+ cvar, fontsize= 14)