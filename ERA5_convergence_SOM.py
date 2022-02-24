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
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd
import cartopy.crs as ccrs
import somoclu
from f_carto import *


save= True
# save= False





syear= 1979
eyear= 2018



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
nrows = 2
ncols = 3


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
                
                ds0= ds0.sel(lat = ds0.lat[ds0.lat > latbound][::2])
                ds0= ds0.sel(lon= ds0.lon[::4])
                ds0= ds0.resample(time='1D').mean()
            else:
                ds1= xr.open_dataset(file_dir + varfile+'.'+str(year)+'.'+str(month).zfill(2)+'.WN6.nc')
                ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()

                ds1= ds1.sel(lat = ds1.lat[ds1.lat > latbound][::2])
                ds1= ds1.sel(lon= ds1.lon[::4])
                
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

ds = ds.sel(time= ds["time.dayofyear"] != 366)
cds = cds.sel(time= cds["time.dayofyear"] != 366)

ds= ds.sel(lat = ds.lat[ds.lat<= 85])    


dsn= ds.mean(dim='time').sum(dim= 'WaveNumb')


dx= (float(dsn.lon[1])- float(dsn.lon[0]) ) *110E3 * np.cos(np.deg2rad(dsn.lat.values))
# dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[np.newaxis, :, np.newaxis] #axis 0 is latitude
dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[:, np.newaxis] #axis 1 is longitude


dy= (float(dsn.lat[1])- float(dsn.lat[0])) *110E3
dvdy= np.gradient(dsn[varlist[1]], axis= 0)/dy #axis 0 is latitude


dsn['dudx']= (dsn[varlist[0]].dims, dudx)
dsn['dvdy']= (dsn[varlist[0]].dims, dvdy)


dsn['div']= (dsn[varlist[0]].dims, dudx+ dvdy)


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
    
    
    
    
# ds['total']= ds[var].sum(dim='WaveNumb')

# ds= ds.drop(var)

    


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



dsn_field= dsn.uQtotLL
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('uQtot')
fig.colorbar(c, ax=ax, shrink= 0.5)


"""vQtotLL"""
ax = fig.add_subplot(2, 3, 3, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



dsn_field= dsn.vQtotLL
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('vQtot')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""d vQ dy"""
ax = fig.add_subplot(2, 3, 6, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



dsn_field=-1* dsn.dvdy
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('-d(vQ) dy')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""d uQ dx"""
ax = fig.add_subplot(2, 3, 5, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



dsn_field= -1* dsn.dudx
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('-d(uQ) dx')
fig.colorbar(c, ax=ax, shrink= 0.5)



"""conv Q"""
ax = fig.add_subplot(2, 3, 4, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    



dsn_field= -1* dsn.div
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
ax.set_title('conv')
fig.colorbar(c, ax=ax, shrink= 0.5)



print('have to get the units right')



"""compressed version for total"""


fig= plt.figure(fignr, figsize= (6,5) )
fignr+=1
plt.clf()
            

ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

ax.coastlines('110m', alpha=0.5)    
ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
ax.gridlines(ylocs=np.arange(0, 90, 10) )
ax.set_boundary(circle, transform=ax.transAxes)    

""""z 850"""
ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')

"""div"""
dsn_field= -1* dsn.div
v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
fig.colorbar(c, ax=ax, shrink= 0.5)




def PlotWind(dsn, uvar, vvar, ax, nx= 15, ny=20, alen=15):
    """ winds with barbs
    alen= arrow length: the lenght of the arrows can be specified by alen
    nx, ny, specifies how many arrows there are in the x and y direction
    """

    lon=dsn.lon.values
    lat= dsn.lat.values
    u= dsn[uvar].values
    v= dsn[vvar].values
    alen= 2E8


    if len(np.shape(lon))==1: #for ERA data lon and lat are put into 2D array
        lon, lat= np.meshgrid(lon, lat)
               
    everyx= np.shape(u)[0]//nx #put a barb in every nth box, such that you have 20 barbs in the longer dim
    everyy= np.shape(u)[1]//ny #20 for a quadratical map , 50 for era
    u = u[0::everyx, 0::everyy]
    v = v[0::everyx, 0::everyy]
    lon= lon[0::everyx, 0::everyy]
    lat= lat[0::everyx, 0::everyy]
        
    # u[lat > 70][:, ::2] = 0

    Q= ax.quiver(lon,lat, (u/np.cos(lat /180 * np.pi)),v, pivot= 'mid', scale=20* alen,
                 headwidth= 4, width=0.005, transform=ccrs.PlateCarree())
    
    qk = plt.quiverkey(Q, 0.5, 0.95, alen, str(alen)+'', coordinates='figure', labelpos='W')


PlotWind(dsn.sel(lat = dsn.lat[dsn.lat >70]), 'uQtotLL', 'vQtotLL', ax, nx= 5, ny=10, alen=2E8)
PlotWind(dsn.sel(lat = dsn.lat[dsn.lat <=70]), 'uQtotLL', 'vQtotLL', ax, nx= 7, ny=20, alen=2E8)



"""total, planetary, synoptic"""


fig= plt.figure(fignr, figsize= (10,5) )
fignr+=1
plt.clf()
            
for ic, cat in enumerate(['Total', 'Plan', 'Syno']):

    ax = fig.add_subplot(1, 3, ic+1, projection=ccrs.NorthPolarStereo())
    
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
    ax.gridlines(ylocs=np.arange(0, 90, 10) )
    ax.set_boundary(circle, transform=ax.transAxes)    
    
    """"z 850"""
    ax.contour(cds.lon, cds.lat, cds.z.mean(dim='time'), transform=ccrs.PlateCarree(), colors= 'k')
    
    
    if cat== 'Total':
        dsn= ds.mean(dim='time').sum(dim= 'WaveNumb')
    elif cat== 'Plan':
        dsn= ds.mean(dim='time').sel(WaveNumb= ds.WaveNumb[ds.WaveNumb <= 3]).sum(dim= 'WaveNumb')
    elif cat== 'Syno':
        dsn= ds.mean(dim='time').sel(WaveNumb= ds.WaveNumb[ds.WaveNumb >= 4]).sum(dim= 'WaveNumb')

    dx= (float(dsn.lon[1])- float(dsn.lon[0]) ) *110E3 * np.cos(np.deg2rad(dsn.lat.values))
    dudx= np.gradient(dsn[varlist[0]], axis= 1)/dx[:, np.newaxis] #axis 1 is longitude   
    dy= (float(dsn.lat[1])- float(dsn.lat[0])) *110E3
    dvdy= np.gradient(dsn[varlist[1]], axis= 0)/dy #axis 0 is latitude
       
    dsn['dudx']= (dsn[varlist[0]].dims, dudx)
    dsn['dvdy']= (dsn[varlist[0]].dims, dvdy)    
    dsn['div']= (dsn[varlist[0]].dims, dudx+ dvdy)

       
    """div"""
    dsn_field= -1* dsn.div
    # if ic == 0: v_extr= np.max([dsn_field.max(), -1* dsn_field.min()])
    # c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), vmin= -v_extr, vmax= v_extr, cmap= 'RdBu_r')
    divlevels= [-500, -300, -200, -100, -50, -25, 25, 50, 100, 200, 300, 500 ]

    # cmap, norm = mcolors.from_levels_and_colors(levels, plt.get_cmap('RdBu_r')(np.linspace(0, 1, len(levels) +1)))
    norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=256, extend='both')
    c = ax.contourf(dsn.lon, dsn.lat, dsn_field , transform=ccrs.PlateCarree(), levels=divlevels, cmap= 'RdBu_r', norm= norm, extend='both')


    PlotWind(dsn.sel(lat = dsn.lat[dsn.lat >70]), 'uQtotLL', 'vQtotLL', ax, nx= 5, ny=10, alen=2E8)
    PlotWind(dsn.sel(lat = dsn.lat[dsn.lat <=70]), 'uQtotLL', 'vQtotLL', ax, nx= 7, ny=20, alen=2E8)

    ax.set_title(cat)

fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
fig.colorbar(c, cax=cbar_ax)




"""a SOM"""
print('compute SOMs')
    

    
ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 0, ds.WaveNumb <= 3), drop= True ).sum(dim='WaveNumb')
ds['Syno']= ds[var].where(ds.WaveNumb >= 4, drop= True ).sum(dim='WaveNumb')
ds['total']= ds[var].sum(dim='WaveNumb')

nyears = eyear - syear + 1




m_cds= cds['z'].weighted(np.cos(np.deg2rad(cds.lat))).mean(("lat", "lon"))
cds[cvar+ '_geo_ano'] =cds['z'] - m_cds



d= cds[cvar+ '_geo_ano']
d= d.where((d.lat > SOM_latlow) & (d.lat < SOM_lathigh), drop= True )

dim_time= len(d.time)
dim_x= len(d.lat)
dim_y= len(d.lon)
 
d_flat = np.reshape(d.data, (dim_time, dim_x*dim_y ) )
print(d_flat.shape)

# Creating instance of SOM
som = somoclu.Somoclu(ncols,nrows ) #,compactsupport= False)

# Training SOM on data.
som.train(d_flat, epochs=500)



fig= plt.figure(fignr, figsize= (2.3*ncols,2.7*nrows +1) )
fignr+=1
plt.clf()
        
#fig, axs = plt.subplots(nrows,ncols, figsize = (10,10))
#d = np.reshape(d,(dim_time, dim_x, dim_y))


for row in range(nrows):
    temp_row = 0
    for col in range(ncols):
        bool_array = [np.array_equal(b,[col,row]) for b in som.bmus]
        print(sum(bool_array))
        t_dat = np.mean(cds[cvar+ '_geo_ano'][bool_array,:,:],axis = 0)

        # for cat in catlist:
        #     print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(ds[cat].where(ds['time.year'] <= split_1)[bool_array])/1E6 ), float(np.std(ds[cat].where(ds['time.year'] <= split_1)[bool_array])/1E6 ) )
        #     print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(ds[cat].where(ds['time.year'] >= split_2)[bool_array])/1E6 ), float(np.std(ds[cat].where(ds['time.year'] >= split_2)[bool_array])/1E6 ) )


        if (len(clevels) < 1) & (row + col== 0): 
            extr= np.max([t_dat.max(), -t_dat.min()])
            clevels= np.round(np.linspace(-extr, extr, 10))



        ax = fig.add_subplot(nrows, ncols, col+1 + row* ncols, projection=ccrs.NorthPolarStereo())
        ax.coastlines('110m', alpha=0.5)    
        ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
        ax.gridlines(ylocs=np.arange(0, 90, 10) )
        ax.set_boundary(circle, transform=ax.transAxes)

#        levels= np.round(np.linspace(-extr, extr, 10), 2)
        ax.contour(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(), levels= clevels, colors= 'k', linewidths= 1)
        

        c = ax.contourf(ds.lon, ds.lat, np.mean(ds['total'][bool_array,:,:], axis= 0), transform=ccrs.PlateCarree(),
                        cmap= 'RdBu_r', extend= 'both', levels= divlevels) #vmin= -1* extr, vmax= extr)         

        ax.set_title(str(row) +'/'+ str(col) + ' ('+ str(np.round(100*sum( bool_array)/len(bool_array) ,1 ) ) +'%) '+
                                  'Plan: ' + str(np.round(float(np.mean(ds['Plan'][bool_array])/1E6 ),1)) + 
                                  ', Syno: ' + str(np.round(float(np.mean(ds['Syno'][bool_array])/1E6 ),1)) ,
                    fontsize= 10 )
        


        
# fig.subplots_adjust(bottom=0.11)
# cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.015])
# cbar= fig.colorbar(c, cax=cbar_ax, orientation="horizontal")
# cbar.set_label(clabel +' anomaly ['+cunit+']')



# if save:
#     savedir= '../Figs/ERA5/SOM_trans-lat-effect/'
#     if not os.path.exists(savedir): os.makedirs(savedir)
    
#     file_name= 'SOM'+str(ncols)+'-'+str(nrows) +'_'
#     if len(monthlist) != 12: file_name += f'month{monthlist[0]}-{monthlist[-1]}_'
#     file_name +=  cvar+ '-'+ str(SOM_latlow)+'-'+ str(SOM_lathigh)
#     file_name += '_' +str(syear)+'-'+str(split_1)+'_'+str(split_2)+'-'+str(eyear)+'_'
#     file_name += var + str(latcircle)
#     savefile= savedir+ file_name
#     print(savefile)
#     plt.savefig(savefile , bbox_inches='tight') 

