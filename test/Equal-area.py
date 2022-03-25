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
elif user=='media':
    Mediadir= '/run/media/pst019/Backup/EC_Earth'
else:
    Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI'
    
    
import matplotlib.pyplot as plt #appears not to work with version 3.1.1 but not with 3.3
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd
import cartopy.crs as ccrs
import somoclu



save= True
save= False
imp= False
imp= True
#
compute= False
compute= True

#Member = 1
#syear= 2076 #1950
#eyear= 2096 # 1960

syear= 1950
eyear= 1951 #70

#split_1 = int(2/3*syear + 1/3*eyear)
#split_2 = int(2/3*eyear + 1/3*syear)

fignr= 5


#Memberlist= [3]#,2,3,4]
Member= 3




var= 'vQtot'
latcircle= 70
catlist = ['Syno', 'Plan', 'total']
perc= 90

clevels= []

#cvar, cvarfile, clabel, cunit= 'T2M', 'T2M', 'Temperature', 'K'
cvar, cvarfile, clabel, cunit, clevels= 'Z', 'Z850', 'Geopotential height', 'm', np.arange(-200, 201, 40)


latbound= 50 #for the import

SOM_latlow = 60
SOM_lathigh = 80

nrows = 4
ncols = 5


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




to_proj = ccrs.AlbersEqualArea(central_longitude=0, central_latitude=80)

lon, lat= np.meshgrid(cds2.lon.values, cds2.lat.values)

xp, yp, _= to_proj.transform_points(ccrs.Geodetic(), lon, lat).T


slpgridx, slpgridy, slp = interpolate_to_grid(x_masked, y_masked, pres, interp_type='cressman',
                                              minimum_neighbors=1, search_radius=400000,
                                              hres=100000)

#to_proj.transform_point(cds2.lon.values, cds2.lat.values, ccrs.Geodetic())

#
#"""to make the map a circle"""
#import matplotlib.path as mpath
#theta = np.linspace(0, 2*np.pi, 100)
#center, radius = [0.5, 0.5], 0.5
#verts = np.vstack([np.sin(theta), np.cos(theta)]).T
#circle = mpath.Path(verts * radius + center)
#
#
#if compute:
#    """computation of the categories"""
#    print('compute')
#    
#    ds= ds2
#    cds= cds2
#    
#    
#    nyears = eyear - syear + 1
#    ds = ds.sel(time= ds["time.dayofyear"] != 366)
#    cds = cds.sel(time= cds["time.dayofyear"] != 366)
#    
#    
#    m_cds= cds['Z'].weighted(np.cos(np.deg2rad(cds.lat))).mean(("lat", "lon"))
#    
#    cds[cvar+ '_geo_ano'] =cds['Z'] - m_cds
#    
#
#    
#    d= cds[cvar+ '_geo_ano']
#    d= d.where((d.lat > SOM_latlow) & (d.lat < SOM_lathigh), drop= True )
#    
#    dim_time= len(d.time)
#    dim_x= len(d.lat)
#    dim_y= len(d.lon)
#     
#    d_flat = np.reshape(d.data, (dim_time, dim_x*dim_y ) )
#    print(d_flat.shape)
#    
#    # Creating instance of SOM
#    som = somoclu.Somoclu(ncols,nrows ) #,compactsupport= False)
#    
#    # Training SOM on data.
#    som.train(d_flat, epochs=300)
#
#
#
#
#    """energy transport"""
#    if catlist == ['Syno', 'Plan', 'total']:
#    
#        ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb')
#        ds['Syno']= ds[var].where(ds.WaveNumb >= 5, drop= True ).sum(dim='WaveNumb')
#        ds['total']= ds[var].sum(dim='WaveNumb')
#        
#        ds= ds.drop(var)
#    
#    
#    
#
#
#
#fig= plt.figure(fignr, figsize= (2.3*ncols,2.7*nrows +1) )
#fignr+=1
#plt.clf()
#        
##fig, axs = plt.subplots(nrows,ncols, figsize = (10,10))
##d = np.reshape(d,(dim_time, dim_x, dim_y))
#
#
#
#for row in range(nrows):
#    temp_row = 0
#    for col in range(ncols):
#        bool_array = [np.array_equal(b,[col,row]) for b in som.bmus]
#        print(sum(bool_array))
##        t_dat = np.mean(d[bool_array,:,:],axis = 0)
#        t_dat = np.mean(cds[cvar+ '_geo_ano'][bool_array,:,:],axis = 0)
#
#        for cat in catlist:
#            print(str(row) +'/'+ str(col) +' '+ cat , float(np.mean(ds[cat][bool_array])/1E6 ), float(np.std(ds[cat][bool_array])/1E6 ) )
#
#
#        if (len(clevels) < 1) & (row + col== 0): 
#            extr= np.max([t_dat.max(), -t_dat.min()])
#            clevels= np.round(np.linspace(-extr, extr, 10))
##        ax1= plt.subplot(ncols, nrows, (col*ncols + row+1 ) )
#        
##        cf = ax1.contourf(d.x,d.y,t_dat, cmap = 'seismic',extend = 'both')
##        cf = axs[row,col].contourf(t_dat.lon, t_dat.lat, t_dat, cmap = 'RdBu_r',extend = 'both')
##        cf = axs[row,col].contourf(d.x,d.y,t_dat, levels= levels, cmap = 'RdBu_r',extend = 'both')
#
#
#
#        ax = fig.add_subplot(nrows, ncols, col+1 + row* ncols, projection=ccrs.NorthPolarStereo())
#        ax.coastlines('110m', alpha=0.5)    
#        ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
#        ax.gridlines(ylocs=np.arange(0, 90, 10) )
#        ax.set_boundary(circle, transform=ax.transAxes)
#
##        levels= np.round(np.linspace(-extr, extr, 10), 2)
#        c = ax.contourf(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(),
#                        cmap= 'RdBu_r', extend= 'both', levels= clevels) #vmin= -1* extr, vmax= extr) 
#        ax.contour(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(), levels= clevels, colors= 'k', linewidths= 1)
#
#        ax.set_title(str(row) +'/'+ str(col) + ' ('+ str(np.round(100*sum(bool_array)/len(bool_array) ,1 ) ) +'%)\n'+
#                                 'Plan:' + str(np.round(float(np.mean(ds['Plan'][bool_array])/1E6 ),1)) +r'$\pm$' +
#                                       str(np.round(float(np.std(ds['Plan'][bool_array])/1E6 ),1))  +'\n'+
#                                 'Syno:' + str(np.round(float(np.mean(ds['Syno'][bool_array])/1E6 ),1)) +r'$\pm$' +
#                                       str(np.round(float(np.std(ds['Syno'][bool_array])/1E6 ),1)),
#                    fontsize= 10 )
#        
#
#
#        
#fig.subplots_adjust(bottom=0.11)
#cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.015])
#cbar= fig.colorbar(c, cax=cbar_ax, orientation="horizontal")
#cbar.set_label(clabel +' anomaly ['+cunit+']')
#
#
#
#if save:
#    savedir= '../Figs/Trans-lat-effect/'
#    if not os.path.exists(savedir): os.makedirs(savedir)
#    
#    file_name= 'SOM'+str(ncols)+'-'+str(nrows) +'_'
#    file_name +=  cvar+ '-'+ str(SOM_latlow)+'-'+ str(SOM_lathigh)
##        for Member in Memberlist:
#    file_name += 'M'+str(Member)
#    file_name += '_' +str(syear)+'-'+str(eyear)+'_'
#    file_name += var + str(latcircle)
#    savefile= savedir+ file_name
#    print(savefile)
#    plt.savefig(savefile , bbox_inches='tight') 
#
#
