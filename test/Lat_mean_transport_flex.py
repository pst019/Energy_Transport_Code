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
#import glob as glob
#import datetime

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save= True
save= False
imp= False
imp= True

#Member = 1
syear= 1950
eyear= 2100
split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 4


Memberlist= [3]#,2,3]


#typ = 'Eddies'
typ = 'Wavelength_smooth'
typ = 'WaveNR'


if typ == 'Eddies':
    #varlist= ['vEmc', 'vEse', 'vEte']
    varlist= ['vQmc', 'vQse', 'vQte']
    catlist= varlist
    monthlymean = True

elif typ == 'Wavelength_smooth':
    varlist= ['vEtot']
    varlist= ['vQtot']
    catlist = ['Meri', 'Plan', 'Syno']#, 'Meso']
    monthlymean = False


elif typ == 'WaveNR':
    varlist= ['vEtot']
    varlist= ['vQtot']
    catlist = ['0', '1-4', '>5']#, '>10']
    monthlymean = False


fig= plt.figure(fignr)
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

axs=fig.subplots(2, 1, sharex= 'col')

if imp:
    print('Import ...')
    for Member in Memberlist:
        print('Member', Member)
        file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

        for var in varlist:
            print(var)
            for year in range(syear, eyear+1):
                if year%10 == 0: print(year)
                for month in range(1,13):
    #                print(month)
                    if year == syear and month== 1:
                        
                        if monthlymean:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
#                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
#                            ds0= ds0.mean(dim = 'time')
#                            ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()
                        
                        if Member == Memberlist[0]: dslat= ds0.lat
                        else: ds0['lat']= dslat
                        
                        
                    else:
                        if monthlymean:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
#                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
#                            ds1= ds1.mean(dim = 'time')
#                            ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                            
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()

                        
                        ds1['lat']= dslat
    
                        
                        ds0= xr.concat([ds0, ds1], dim= 'time')
    
                        
                    """the last loop step to merge the different variables in one ds"""
                    if year == eyear and month == 12:
        #                print(ds)
                        if var== varlist[0]:
                            ds2= ds0
                        else:
                            ds2= xr.merge([ds2, ds0])
                            
        if Member == Memberlist[0]:
            ds= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
        else:
            ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            ds= xr.merge([ds, ds2])


    if not monthlymean: ds= ds.resample(time='1M').mean()


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(dslat))





"""Total transport"""
if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
if typ in ['Wavelength_smooth', 'WaveNR']: dstot = ds[var].sum(dim= "WaveNumb")



mean= dstot.mean(dim='time').mean(dim='Member')
std= dstot.mean(dim='time').std(dim='Member') 

line, = axs[0].plot(ds.lat, LatCirc*mean, label= 'Total', color= 'k')
axs[0].plot(ds.lat, LatCirc*(mean + std), label='', lw= 1, color= line.get_color() )
axs[0].plot(ds.lat, LatCirc*(mean - std), label='', lw= 1, color= line.get_color() )


diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).mean(dim='Member')
diffstd= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).std(dim='Member')

#diffmean= sum(ds[var].where(ds["time.year"]>= split_2).mean(dim='time').mean(dim='Member') for var in varlist)  - sum(ds[var].where(ds["time.year"]<= split_1).mean(dim='time').mean(dim='Member') for var in varlist)
#diffstd= (sum(ds[var].where(ds["time.year"]>= split_2).mean(dim='time') for var in varlist)  - sum(ds[var].where(ds["time.year"]<= split_1).mean(dim='time') for var in varlist) ).std(dim='Member')


line, = axs[1].plot(ds.lat, LatCirc*diffmean, color= 'k')
axs[1].plot(ds.lat, LatCirc*(diffmean +diffstd) , lw= 1, color= line.get_color())
axs[1].plot(ds.lat, LatCirc*(diffmean -diffstd) , lw= 1, color= line.get_color())

#axs[1].plot(ds.lat, np.sum(ds[var].where(ds["time.year"]>= 2050).mean(dim='time') for var in varlist)
#                    - np.sum(ds[var].where(ds["time.year"]<= 2000).mean(dim='time') for var in varlist), label= 'Total', color= 'k')


"""Transport of each category"""
if typ == 'Wavelength_smooth':
    Split1= 8000E3
#    Split2= 2000E3
    
    WaveSplit1= LatCirc/Split1
#    WaveSplit2= LatCirc/Split2
    ds['Meri']= ds[var].sel(WaveNumb= 0)
    ds['Plan']= ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb')+  WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0)
#    ds['Syno']=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit1)//1, ds.WaveNumb < WaveSplit2//1)).sum(dim='WaveNumb') 
#                     + (1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) 
#                     + WaveSplit2%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) )
    ds['Syno']= ( ds[var].where(ds.WaveNumb > WaveSplit1//1).sum(dim='WaveNumb') 
                     +(1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) )



if typ == 'WaveNR':
    ds['0']= ds[var].sel(WaveNumb= 0)
    ds['1-4']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb')
#    ds['>5']= ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb')
    ds['>5']= ds[var].where(ds.WaveNumb >= 5, drop= True ).sum(dim='WaveNumb')

            
for ci, cat in enumerate(catlist):
#    color = next(plt.gca()._get_lines.prop_cycler)['color']
    mean= ds[cat].mean(dim='time').mean(dim='Member')
    std= ds[cat].mean(dim='time').std(dim='Member')
    
    line, = axs[0].plot(ds.lat, LatCirc* mean, label= var)
    axs[0].plot(ds.lat, LatCirc* (mean +std), label='', lw= 1, color= line.get_color() )
    axs[0].plot(ds.lat, LatCirc* (mean -std), label='', lw= 1, color= line.get_color() )


#    axs[1].plot(ds.lat,  ds[var].mean(dim='time').std(dim='Member'), label= var)

    diffmean= (ds[cat].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat].where(ds["time.year"]<= split_1).mean(dim='time')).mean(dim='Member')
    diffstd= (ds[cat].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat].where(ds["time.year"]<= split_1).mean(dim='time')).std(dim='Member')
    
    line, = axs[1].plot(ds.lat, LatCirc*diffmean, label= var)
    axs[1].plot(ds.lat, LatCirc*(diffmean+ diffstd) , lw= 1, color= line.get_color())
    axs[1].plot(ds.lat, LatCirc*(diffmean- diffstd) , lw= 1, color= line.get_color())



#ds.where(ds["time.year"]< 2000, drop= True)
    
axs[0].legend(ncol= 2)

axs[0].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[1].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-90, 90)
plt.xlabel('Latitude')
axs[0].set_ylabel('Energy transport [W]') #energy transport accross a latitude band
#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel(str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Global/'
    
    file_name= typ+'_M' 
    for Member in Memberlist: file_name += str(Member)
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 