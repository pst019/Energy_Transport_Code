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

#
compute= False
compute= True



syear= 2010
eyear= 2010

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

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
                print(ds0)

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

ds['div']= (ds[varlist[0]].dims, dudx+ dvdy)


var= 'div'



"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


if compute:
    """computation of the categories"""
    print('compute SOMs')
    
    ds= ds2
    
    #get one season 1 - DJF, 2 - MAM, 3 - JJA, 4 - SON
    # cds.where(cds["time.month"]%12 // 3 + 1 == 1, drop= True)
    
    nyears = eyear - syear + 1
    ds = ds.sel(time= ds["time.dayofyear"] != 366)
    cds = cds.sel(time= cds["time.dayofyear"] != 366)
    
    
    
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
    som.train(d_flat, epochs=30)




    """energy transport"""
    if catlist == ['Syno', 'Plan', 'total']:
    
        ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 3), drop= True ).sum(dim='WaveNumb')
        ds['Syno']= ds[var].where(ds.WaveNumb >= 4, drop= True ).sum(dim='WaveNumb')
        ds['total']= ds[var].sum(dim='WaveNumb')
        
        ds= ds.drop(var)
    
    
    

fignr= 1

for cat in catlist[-1:]:
    
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
    #        t_dat = np.mean(d[bool_array,:,:],axis = 0)
            t_dat = np.mean(cds[cvar+ '_geo_ano'][bool_array,:,:],axis = 0)
            var_dat = np.mean(ds[cat][bool_array,:,:],axis = 0) #check!!!
    
    
            # for cat in catlist:
            #     print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(ds[cat].where(ds['time.year'] <= split_1)[bool_array])/1E6 ), float(np.std(ds[cat].where(ds['time.year'] <= split_1)[bool_array])/1E6 ) )
            #     print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(ds[cat].where(ds['time.year'] >= split_2)[bool_array])/1E6 ), float(np.std(ds[cat].where(ds['time.year'] >= split_2)[bool_array])/1E6 ) )
    
    
            if (len(clevels) < 1) & (row + col== 0): 
                extr= np.max([t_dat.max(), -t_dat.min()])
                clevels= np.round(np.linspace(-extr, extr, 10))
    #        ax1= plt.subplot(ncols, nrows, (col*ncols + row+1 ) )
            
    #        cf = ax1.contourf(d.x,d.y,t_dat, cmap = 'seismic',extend = 'both')
    #        cf = axs[row,col].contourf(t_dat.lon, t_dat.lat, t_dat, cmap = 'RdBu_r',extend = 'both')
    #        cf = axs[row,col].contourf(d.x,d.y,t_dat, levels= levels, cmap = 'RdBu_r',extend = 'both')
    
    
    
            ax = fig.add_subplot(nrows, ncols, col+1 + row* ncols, projection=ccrs.NorthPolarStereo())
            ax.coastlines('110m', alpha=0.5)    
            ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
            ax.gridlines(ylocs=np.arange(0, 90, 10) )
            ax.set_boundary(circle, transform=ax.transAxes)
    
    #        levels= np.round(np.linspace(-extr, extr, 10), 2)
            c = ax.contourf(var_dat.lon, var_dat.lat, var_dat, transform=ccrs.PlateCarree(), cmap= 'RdBu_r') #vmin= -1* extr, vmax= extr) 
            ax.contour(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(), levels= clevels, colors= 'k', linewidths= 1)
    
            ax.set_title(str(row) +'/'+ str(col) ,#+ ' ('+ str(np.round(100*sum((ds['time.year'] <= split_1).values & bool_array)/sum((ds['time.year'] <= split_1).values) ,1 ) ) +'% / '+
                                      #                 str(np.round(100*sum((ds['time.year'] >= split_2).values & bool_array)/sum((ds['time.year'] >= split_2).values)  ,1 ) ) +'%)\n'+
                                      # 'Plan: ' + str(np.round(float(np.mean(ds['Plan'].where(ds['time.year'] <= split_1)[bool_array])/1E6 ),1)) + ' / ' +
                                      #       str(np.round(float(np.mean(ds['Plan'].where(ds['time.year'] >= split_2)[bool_array])/1E6 ),1))    +'\n'+
                                      # 'Syno: ' + str(np.round(float(np.mean(ds['Syno'].where(ds['time.year'] <= split_1)[bool_array])/1E6 ),1)) + ' / ' +
                                      #       str(np.round(float(np.mean(ds['Syno'].where(ds['time.year'] >= split_2)[bool_array])/1E6 ),1)) ,
                        fontsize= 10 )
            
    
    
            
    fig.subplots_adjust(bottom=0.11)
    cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.015])
    cbar= fig.colorbar(c, cax=cbar_ax, orientation="horizontal")
    # cbar.set_label(clabel +' anomaly ['+cunit+']')
    
    
    
    if save:
        savedir= '../Figs/ERA5/SOM_trans-lat-effect/'
        if not os.path.exists(savedir): os.makedirs(savedir)
        
        file_name= 'SOM'+str(ncols)+'-'+str(nrows) +'_'
        if len(monthlist) != 12: file_name += f'month{monthlist[0]}-{monthlist[-1]}_'
        file_name +=  cvar+ '-'+ str(SOM_latlow)+'-'+ str(SOM_lathigh)
        # file_name += '_' +str(syear)+'-'+str(split_1)+'_'+str(split_2)+'-'+str(eyear)+'_'
        # file_name += var #+ str(latcircle)
        file_name += '_'+cat
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 
    
    
