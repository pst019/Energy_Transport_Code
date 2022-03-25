#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import os

params = {'axes.labelsize': 13,
          # 'axes.titlesize': 16,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'legend.fontsize': 12}
plt.rcParams.update(params)

user = os.getcwd().split('/')[2]

# model='EC_Earth'
model='ERA5'

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/'
#    Mediadir= '/media/'+user+'/Backup/'
elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/'    
else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/'
#'/nird/projects/NS9063K/Richard_KNMI'

datadir = Mediadir + 'data/Energy_Transport/'+ model +'/'
if model=='ERA5': datadir += 'EnergySplit/res_0.5x0.5/Waves/'


latcut= 89

  
save= True
# save= False


if model == 'EC_Earth':
    syear= 1950
    eyear= 2100
    Memberlist= [1]#, 2, 3, 4, 5]#,2,3,4]


elif model== 'ERA5':
    syear= 1979
    eyear= 2018
    Memberlist=[1]


timeperiod= 'year'


split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)


# Split1= 8000E3
# Split2= 4000E3
# Split3= 2000E3

Splitlist= [8000E3, 4000E3, 2000E3]
Splitlist= [8000E3, 6000E3, 4000E3, 2000E3]

fignr= 1



varlist, energy_label= ['vEtot', 'vEsetot'], 'energy'
# varlist, energy_label= ['vQtot', 'vQsetot'], 'moisture'

monthlymeanlist=[False, True]


intermediate_file_name= Mediadir +'data/Energy_Transport/'+ model +'/intermediate_data/'+timeperiod+'-mean_'
for var in varlist: intermediate_file_name += var+'_'
intermediate_file_name += str(syear)+'-'+str(eyear)
for Member in Memberlist: intermediate_file_name += '_M'+str(Member)
intermediate_file_name +=  '.nc'


if os.path.exists(intermediate_file_name):
    ds= xr.open_dataset(intermediate_file_name)
    
else:
    print('Import ...')
    for Member in Memberlist:
        file_dir= datadir        
        if model== 'EC_Earth': 
            print('Member', Member) 
            file_dir += 'Member'+str(Member) +'/'  
        
        for vi, var in enumerate(varlist):
            monthlymean = monthlymeanlist[vi]
            print(var)
            for year in range(syear, eyear+1):
                if year%10 == 0: print(year)
                for month in range(1,13):
    #                print(month)
                    if year == syear and month== 1:
                        
                        if monthlymean:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds0= ds0.mean(dim = 'time') #calculate the monthly mean
                            ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))


                        if Member == Memberlist[0]: dslat= ds0.lat
                        else: ds0['lat']= dslat
                        
                        
                    else:
                        if monthlymean:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds1= ds1.mean(dim = 'time')
                            ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

                        
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

    ds.to_netcdf(intermediate_file_name)



if latcut != 90: ds= ds.where(np.abs(ds.lat) < latcut, drop= True)


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))




ds= ds.mean(dim=['time', 'Member'])
ds['WaveNumb']= ds['WaveNumb'] -0.5

ds *= np.sign(ds.lat) /1E7

varlist= varlist + [varlist[0][:2]+'trans']
ds[varlist[-1]]= ds[varlist[0]] - ds[varlist[1]]

from matplotlib.ticker import MultipleLocator

for var in varlist:
    fig= plt.figure(fignr, figsize= (9,6))
    fignr+= 1
    plt.clf()

    ax = fig.add_subplot(111)

    # if var== 'vEtot': cbarlabel= 'Mean transport [10$^7$ W m$^{-1}$]'
    # elif var== 'vEsetot': cbarlabel= 'Stationary transport [10$^7$ W m$^{-1}$]'
    # elif var== 'vEtrans': cbarlabel= 'Transient transport [10$^7$ W m$^{-1}$]'
    # else: cbarlabel= var
    cbarlabel= 'Poleward '+energy_label+' transport [10$^7$ W m$^{-1}$]' 
    
    ds[var].plot(y= 'WaveNumb', x= 'lat', vmax= 4., cbar_kwargs= {'label': cbarlabel})


    for Split in Splitlist:
        ###make the lines with the wavelength
        WaveSplit= LatCirc/Split
    
        plt.plot(WaveSplit.lat, WaveSplit, c= 'k')
        plt.text(0, WaveSplit[len(WaveSplit)//2].values-.3, str(int(Split/1E3)) +' km',
                 fontsize=12, ha='center', va= 'top')


    max_transport= ds[var].sum(dim='WaveNumb').max()
    for lat in ds.lat.values:
        latmax= ds[var].sel(lat= lat).argmax()
        transp_latmax= ds[var].sel(lat= lat).max()
        
        if transp_latmax > 0.05* max_transport:
            plt.plot(lat, latmax, 'o', color= 'grey')


    plt.ylabel('Wavenumber')
    plt.yticks(np.arange(0,21, 2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    plt.xlabel('Latitude')

    if var== 'vEtot': cbarlabel= 'Mean transport \n per lenth unit [W m$^{-1}$]'



    if save:
        savedir= Mediadir +'/home/Energy_Transport/Figs/Scale+PowerSpec_3/'
        if not os.path.exists(savedir): os.makedirs(savedir)
        
        file_name= 'Spectrum'
        if model != 'ERA5':
            file_name += '_M'
            for Member in Memberlist: file_name += str(Member)
        file_name += '_' +str(syear)+'-'+str(eyear)+'_'
        file_name += var
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 

