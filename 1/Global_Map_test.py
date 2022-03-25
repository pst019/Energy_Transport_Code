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
import matplotlib as mpl

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save= False
# save= True
imp= False
imp= True

Member = 3
year= 2000



fignr= 5


#var, varfile= 'vEtotLLd', 'vEtotLL'
#var, varfile= 'vQtotLLd', 'vQtotLL'


# varlist=[ 'uQtotLLd', 'vQtotLLd'] #important that the u component is first
# varfilelist=[ 'uQtotLL', 'vQtotLL']
# dvar= 'conv_Q' #the display variable


varlist=['uEtotLLd', 'vEtotLLd'] #important that the u component is first
varfilelist=[ 'uEtotLL', 'vEtotLL']

dvar= 'conv_E'

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'



unit='PerM'
# unit='LatCycl'

if imp:
    for vi, vari in enumerate(varlist):
        print(vari)
        varfile= varfilelist[vi]
    
        print(year)
        if year//10 == 0: print(year)
        for month in range(1,13):

            if year == year and month== 1:
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
        if year == year and month == 12:
            if vari== varlist[0]:
                ds2= ds0
            else:
                ds2= xr.merge([ds2, ds0])


ds= ds2




zvarlist= ['vEtot'] #zonal mean file (zonalds)
# zvarlist= ['vQtot'] #zonal mean file (zonalds)


monthlymean = False



if imp:
    print('Import ...')
    file_dir= Mediadir + 'data/Energy_Transport/'+ model +'/'
    file_dir += 'Member'+str(Member) +'/'            
    
    for vi, zvar in enumerate(zvarlist):
        print(zvar)
        # monthlymean = monthlymeanlist[vi]

        for month in range(1,13):
#                print(month)
            if year == year and month== 1:
                print(file_dir+ zvar+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc')
                if monthlymean:
                    ds0= xr.open_dataset(file_dir + zvar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                else:
                    ds0= xr.open_dataset(file_dir + zvar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds0= ds0.mean(dim = 'time') #calculate the monthly mean
                    ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))


                dslat= ds0.lat
                
                
            else:
                if monthlymean:
                    ds1= xr.open_dataset(file_dir + zvar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                else:
                    ds1= xr.open_dataset(file_dir + zvar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds1= ds1.mean(dim = 'time')
                    ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

                
                ds1['lat']= dslat

                
                ds0= xr.concat([ds0, ds1], dim= 'time')

                
            """the last loop step to merge the different zvariables in one ds"""
            if year == year and 12:
#                print(ds)
                if zvar== zvarlist[0]:
                    ds2= ds0
                else:
                    ds2= xr.merge([ds2, ds0])
                        



zonalds= ds2

















if unit =='PerM':
    fac= 1
    zfac= 1
elif unit == 'LatCycl':
    a= 6371E3
    LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))
    fac= LatCirc

    LatCirc = 2* np.pi * a * np.cos(np.deg2rad(zonalds.lat))
    zfac= LatCirc


"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


ds = ds.where(ds.lat < 85, drop= True)





dx= (float(ds.lon[1])- float(ds.lon[0]) ) *110E3 * np.cos(np.deg2rad(ds.lat.values))
dudx= np.gradient(ds[varlist[0]], axis= -1)/ dx[np.newaxis, :, np.newaxis] #axis -1 is longitude


dy= (float(ds.lat[1])- float(ds.lat[0])) *110E3
dvdy= np.gradient(ds[varlist[1]], axis= -2)/ dy #axis -2 is latitude

ds[dvar]= (ds[varlist[0]].dims, -1* (dudx+ dvdy) )


zonaldy= (float(zonalds.lat[1])- float(zonalds.lat[0])) *110E3

dvdy= np.gradient(zonalds[zvarlist[0]], axis= -1)/ zonaldy #axis -2 is latitude

zonalds[dvar]= (zonalds[zvarlist[0]].dims, -1* dvdy)


"""zonal-mean annual-mean meridional transport"""
fig= plt.figure(fignr)
fignr+=1
plt.clf()

ds_mean= ds.mean(dim= ['time', 'lon', 'WaveNumb'])

plt.plot(ds.lat, fac* ds.mean(dim= ['time', 'lon']).sum(dim= 'WaveNumb')[varlist[1]], label='Total', c= 'k')
plt.plot(zonalds.lat, zfac* zonalds.mean(dim= ['time']).sum(dim= 'WaveNumb')[zvarlist[0]], label='Total', c= 'k', linestyle='dashed')

plt.plot(ds.lat, fac* ds.sel(WaveNumb= [0,1]).mean(dim= ['time', 'lon']).sum(dim= 'WaveNumb')[varlist[1]], label='0-1', c= 'b')
plt.plot(zonalds.lat, zfac* zonalds.sel(WaveNumb= [0,1]).mean(dim= ['time']).sum(dim= 'WaveNumb')[zvarlist[0]], label='0-1', c= 'b', linestyle='dashed')

plt.plot(ds.lat, fac* ds.sel(WaveNumb= [2,3]).mean(dim= ['time', 'lon']).sum(dim= 'WaveNumb')[varlist[1]], label='2-3', c='g')
plt.plot(zonalds.lat, zfac* zonalds.sel(WaveNumb= [2,3]).mean(dim= ['time']).sum(dim= 'WaveNumb')[zvarlist[0]], label='2-3', c= 'g', linestyle='dashed')

plt.plot(ds.lat, fac* ds.sel(WaveNumb= [4,5]).mean(dim= ['time', 'lon']).sum(dim= 'WaveNumb')[varlist[1]], label='4-5')


plt.ylabel(zvarlist[0])

plt.legend()
plt.title('To show that the LLfiles produce similar results to zonal mean')




"""zonal-mean annual-mean energy convergence"""
fig= plt.figure(fignr)
fignr+=1
plt.clf()

ds_mean= ds.mean(dim= ['time', 'lon', 'WaveNumb'])

plt.plot(ds.lat, ds.mean(dim= ['time', 'lon']).sum(dim= 'WaveNumb')[dvar], label='Total' , c= 'k')
plt.plot(zonalds.lat, zonalds.mean(dim= ['time']).sum(dim= 'WaveNumb')[dvar], label='Total' , c= 'k', linestyle= 'dashed')

# plt.plot(ds.lat, ds.sel(WaveNumb= [0,1]).mean(dim= ['time', 'lon', 'WaveNumb'])[dvar], label='0-1')
# plt.plot(ds.lat, ds.sel(WaveNumb= [2,3]).mean(dim= ['time', 'lon', 'WaveNumb'])[dvar], label='2-3')
# plt.plot(ds.lat, ds.sel(WaveNumb= [4,5]).mean(dim= ['time', 'lon', 'WaveNumb'])[dvar], label='4-5')

plt.ylabel(dvar)


plt.legend()


#catlist= ['1', '2', '3', '4', 'Syno', 'Total']
#figsize= (10,7)
#figx, figy= 3, 2


catlist= ['Plan', 'Syno', 'Total']
figsize= (10,4)
figx, figy= 3, 1

fig= plt.figure(fignr, figsize)
fignr+=1
plt.clf()
levels=[]

#    dsn= ds.where(np.logical_and(ds["time.year"] >= syear_list[iy], ds["time.year"] <= eyear_list[iy]  ), drop=True ).mean(dim='time')
dsn= ds.mean(dim='time')
if dvar == 'conv_Q':
    levels= [-400, -200, -100, -50, -20, 20, 50, 100, 200, 400] #the first and last are not plotted by necessary for the norm

levels= np.arange(-230, 231, 20)
#    if i == 0: dsn= ds.resample(time='1Y').mean().std(dim='time') #annual mean standard deviation

 

d= dsn[dvar].sum(dim= 'WaveNumb')
d_extreme =  np.max([d.max(), -d.min()])


"""make the categories to plot"""
if catlist == ['Plan', 'Syno', 'Total']:
    dsn['Plan'] = dsn[dvar].where(dsn['WaveNumb'] < 5, drop= True).sum(dim= 'WaveNumb')
    dsn['Syno'] = dsn[dvar].sel(WaveNumb= 5)
    dsn['Total'] = dsn[dvar].sum(dim= 'WaveNumb')
    catlist= ['Plan', 'Syno', 'Total']   

if catlist == ['1', '2', '3', '4', 'Syno', 'Total']:
    dsn['1'] = dsn[dvar].sel(WaveNumb= 1)
    dsn['2'] = dsn[dvar].sel(WaveNumb= 2)
    dsn['3'] = dsn[dvar].sel(WaveNumb= 3)
    dsn['4'] = dsn[dvar].sel(WaveNumb= 4)
    dsn['Syno'] = dsn[dvar].sel(WaveNumb= 5)
    dsn['Total'] = dsn[dvar].sum(dim= 'WaveNumb')



for ci, cat in enumerate(catlist):
   
        
    ax = fig.add_subplot(figy, figx, ci+1, projection=ccrs.NorthPolarStereo())
    
    ax.coastlines('110m', alpha=0.5)    
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    ax.gridlines(ylocs=np.arange(0, 90, 10) )
    ax.set_boundary(circle, transform=ax.transAxes)
    
    if levels==[]:
        c = ax.contourf(dsn.lon, dsn.lat, dsn[cat], transform=ccrs.PlateCarree(), cmap= 'RdBu', vmin= -d_extreme, vmax= d_extreme) #norm= norm)
    else:
        norm = mpl.colors.BoundaryNorm(boundaries=levels, ncolors=256)#, clip=True) #extend= 'both')
        c = ax.contourf(dsn.lon, dsn.lat, dsn[cat], transform=ccrs.PlateCarree(), cmap= 'RdBu', levels= levels[1:-1], norm=norm, extend= 'both') #norm= norm)

    ax.set_title(cat)
    
#        c = ax.contourf(dsn.lon, dsn.lat, dsn[dvar], transform=ccrs.PlateCarree(), cmap= 'RdBu_r', vmin= -d_extreme, vmax= d_extreme) #norm= norm)



fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.87, 0.25, 0.01, 0.5])
cbar= fig.colorbar(c, cax=cbar_ax)        
cbar.set_label(dvar)

   
if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Maps/'
    
    file_name= 'Map_PS_M'+str(Member) +'_' 
    # if i == 0: file_name += str(syear)+'-'+str(eyear)+'_'
    # if i == 1: file_name += 'Diff' + str(split_2)+'-'+str(eyear)+ 'mi' + str(syear)+'-'+str(split_1)+'_'

    file_name += dvar
    for cat in catlist: file_name += '-'+ cat
    
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 



# fig= plt.figure(fignr)
# fignr+=1
# plt.clf()


# plt.plot(ds.lat, ds[dvar].mean(dim= 'lon').mean(dim='time').mean(dim='WaveNumb') )
