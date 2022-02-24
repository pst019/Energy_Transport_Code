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
#save= False
imp= False
imp= True

#Member = 1
syear= 1950
eyear= 2100

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 3

var= 'vEmc'
var= 'vEse'
var= 'vEte'

Memberlist= [1]#,2,3]
#varlist= ['vEse']
#varlist= ['vEmc', 'vEse', 'vEte']
varlist= ['vQmc', 'vQse', 'vQte']


fig= plt.figure(fignr)
fignr+=1
plt.clf()

#ax1= plt.subplot(2, 1, 1)

axs=fig.subplots(2, 1, sharex= 'col')

if imp:
    for Member in Memberlist:
        print('Member', Member)
        file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

        for var in varlist:
            print(var)
            for year in range(syear, eyear+1):
    #            print(year)
                for month in range(1,13):
    #                print(month)
                    if year == syear and month== 1:
                        ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        if Member == Memberlist[0]: dslat= ds0.lat
                        else: ds0['lat']= dslat
                    else:
                        ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        ds1['lat']= dslat #to correct for some file issues
#                        ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
    
                        
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


fac = 1E-15* 2* np.pi *6371E3 * np.cos(np.deg2rad(ds.lat))


"""Total transport"""
mean= sum(ds[var].mean(dim='time').mean(dim='Member') for var in varlist)
std= sum(ds[var].mean(dim='time')for var in varlist).std(dim='Member') 

line, = axs[0].plot(ds.lat, fac*mean, label= 'Total', color= 'k')
axs[0].plot(ds.lat, fac*(mean + std), label='', lw= 1, color= line.get_color() )
axs[0].plot(ds.lat, fac*(mean - std), label='', lw= 1, color= line.get_color() )


diffmean= sum(ds[var].where(ds["time.year"]>= split_2).mean(dim='time').mean(dim='Member') for var in varlist)  - sum(ds[var].where(ds["time.year"]<= split_1).mean(dim='time').mean(dim='Member') for var in varlist)
diffstd= (sum(ds[var].where(ds["time.year"]>= split_2).mean(dim='time') for var in varlist)  - sum(ds[var].where(ds["time.year"]<= split_1).mean(dim='time') for var in varlist) ).std(dim='Member')


line, = axs[1].plot(ds.lat, fac*diffmean, color= 'k')
axs[1].plot(ds.lat, fac*(diffmean+ diffstd) , lw= 1, color= line.get_color())
axs[1].plot(ds.lat, fac*(diffmean- diffstd) , lw= 1, color= line.get_color())

#axs[1].plot(ds.lat, np.sum(ds[var].where(ds["time.year"]>= 2050).mean(dim='time') for var in varlist)
#                    - np.sum(ds[var].where(ds["time.year"]<= 2000).mean(dim='time') for var in varlist), label= 'Total', color= 'k')


#colorlist= ['k','#17becf','#1f77b4','#2ca02c','#d62728','#ff7f0e','#8c564b']
            
for vi, var in enumerate(varlist):
#    color = next(plt.gca()._get_lines.prop_cycler)['color']
    line, = axs[0].plot(ds.lat, fac*ds[var].mean(dim='time').mean(dim='Member'), label= var)
    axs[0].plot(ds.lat, fac*(ds[var].mean(dim='time').mean(dim='Member')+ ds[var].mean(dim='time').std(dim='Member')), label='', lw= 1, color= line.get_color() )
    axs[0].plot(ds.lat, fac*(ds[var].mean(dim='time').mean(dim='Member')- ds[var].mean(dim='time').std(dim='Member')), label='', lw= 1, color= line.get_color() )


#    axs[1].plot(ds.lat,  ds[var].mean(dim='time').std(dim='Member'), label= var)

    diffmean= ds[var].where(ds["time.year"]>= split_2).mean(dim='time').mean(dim='Member') - ds[var].where(ds["time.year"]<= split_1).mean(dim='time').mean(dim='Member')
    diffstd= (ds[var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[var].where(ds["time.year"]<= split_1).mean(dim='time')).std(dim='Member')
    
    line, = axs[1].plot(ds.lat, fac*diffmean, label= var)
    axs[1].plot(ds.lat, fac*(diffmean+ diffstd) , lw= 1, color= line.get_color())
    axs[1].plot(ds.lat, fac*(diffmean- diffstd) , lw= 1, color= line.get_color())



#ds.where(ds["time.year"]< 2000, drop= True)
    
axs[0].legend(ncol= 2)

axs[0].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[1].plot([-90, 90], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-90, 90)
plt.xlabel('Latitude')
axs[0].set_ylabel('Energy transport [PW]') #energy transport accross a latitude band
#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel(str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))


if save:
    savedir= Mediadir+ '/Dropbox/Energy_Transport/Figs/Global/'
    
    file_name= 'Eddies_M' 
    for Member in Memberlist: file_name += str(Member)
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 