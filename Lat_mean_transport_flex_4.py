#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]


model='EC_Earth'
# model='ERA5'

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/'+ model +'/'

else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/EnergySplit/res_0.5x0.5/Waves/'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd



save= True
save= False
imp= False
imp= True
#
#Member = 1
syear= 1950
eyear= 1960
split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 6


Memberlist= [3]#,2,3,4]


#typ = 'Eddies'
typ = 'Wavelength_smooth'
typ = 'Wavelength_smooth_2'
# typ = 'Wavelength_smooth_3'

#typ = 'WaveNR'

latcut= 87

unit='PerM'
# unit='LatCycl'

if typ == 'Eddies':
    varlist= ['vEmc', 'vEse', 'vEte']
#    varlist= ['vQmc', 'vQse', 'vQte']
    catlist= varlist
    monthlymean = True

elif typ == 'Wavelength_smooth':
    varlist= ['vEtot']
#    varlist= ['vQtot']
    catlist = ['Meri', 'Plan', 'Syno']#, 'Meso']
    
    monthlymean = False

elif typ == 'Wavelength_smooth_2':
    # varlist= ['vEtot', 'vEsetot']
    varlist= ['vQtot', 'vQsetot']
    
    monthlymeanlist=[False, True]
    
#    varlist= ['vQtot']
    catlist = ['Meri', 'Plan', 'Syno', 'Meso']
    monthlymean = False


elif typ == 'Wavelength_smooth_3':
    varlist= ['vEtot']
#    varlist= ['vQtot']
    catlist = ['Meri', 'Plan', 'Syno_l', 'Syno_s', 'Meso']
    
    monthlymean = False


elif typ == 'WaveNR':
    varlist= ['vEtot']
#    varlist= ['vQtot']
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
        file_dir= Mediadir
        if model== 'EC_Earth': 
            print('Member', Member) 
            file_dir += 'Member'+str(Member) +'/'            
        
        for vi, var in enumerate(varlist):
            print(var)
            monthlymean = monthlymeanlist[vi]

            for year in range(syear, eyear+1):
                if year%10 == 0: print(year)
                for month in range(1,13):
    #                print(month)
                    if year == syear and month== 1:
                        print(file_dir+ var+'.'+str(year)+'.'+str(month).zfill(2)+ '.nc')
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


if latcut != 90: ds= ds.where(np.abs(ds.lat) < latcut, drop= True)


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc



"""Total transport"""


for vi, var in enumerate(varlist):
    if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
    if typ in ['Wavelength_smooth','Wavelength_smooth_2','Wavelength_smooth_3', 'WaveNR']: dstot = ds[var].sum(dim= "WaveNumb")

    mean= dstot.mean(dim='time').mean(dim='Member')
    diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).mean(dim='Member')
    
    if vi == 0:
        line, = axs[0].plot(ds.lat, fac*mean, label= 'Total', color= 'k')
        axs[1].plot(ds.lat, fac*diffmean, color= 'k')
    # else:
    #     line, = axs[0].plot(ds.lat, fac*mean, '--' , color= 'k')
    #     axs[1].plot(ds.lat, fac*diffmean, color= 'k')
    
    if len(Memberlist) > 1:
        std= dstot.mean(dim='time').std(dim='Member') 
        axs[0].plot(ds.lat, fac*(mean + std), label='', lw= 1, color= line.get_color() )
        axs[0].plot(ds.lat, fac*(mean - std), label='', lw= 1, color= line.get_color() )
   
        diffstd= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).std(dim='Member')        
        axs[1].plot(ds.lat, fac*(diffmean +diffstd) , lw= 1, color= line.get_color())
        axs[1].plot(ds.lat, fac*(diffmean -diffstd) , lw= 1, color= line.get_color())
    


"""Transport of each category"""
if typ == 'Wavelength_smooth':
    Split1= 8000E3
    
    WaveSplit1= LatCirc/Split1
    ds['Meri']= ds[var].sel(WaveNumb= 0)
    ds['Plan']= ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb')+  WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0)
    ds['Syno']= ( ds[var].where(ds.WaveNumb > WaveSplit1//1).sum(dim='WaveNumb') 
                     +(1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) )
    
if typ == 'Wavelength_smooth_2':
    Split1= 8000E3
    Split2= 2000E3
    
    WaveSplit1= LatCirc/Split1
    WaveSplit2= LatCirc/Split2
    
    for var in varlist:
        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb')+  WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0)
        ds['Syno_'+var]=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit1)//1, ds.WaveNumb < WaveSplit2//1)).sum(dim='WaveNumb') 
                         + (1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) 
                         + WaveSplit2%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) )
        ds['Meso_'+var]= ( ds[var].where(ds.WaveNumb > WaveSplit2//1).sum(dim='WaveNumb') 
                         +(1- WaveSplit2%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) )



if typ == 'Wavelength_smooth_3':
    Split1= 10000E3
    Split2= 4000E3
    Split3= 2000E3
    
    WaveSplit1= LatCirc/Split1
    WaveSplit2= LatCirc/Split2
    WaveSplit3= LatCirc/Split3
    
    ds['Meri']= ds[var].sel(WaveNumb= 0)
    ds['Plan']= ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb')+  WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0)
    ds['Syno_l']=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit1)//1, ds.WaveNumb < WaveSplit2//1)).sum(dim='WaveNumb') 
                     + (1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) 
                     + WaveSplit2%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) )
    ds['Syno_s']=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit2)//1, ds.WaveNumb < WaveSplit3//1)).sum(dim='WaveNumb') 
                     + (1- WaveSplit2%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) 
                     + WaveSplit3%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit3//1).fillna(0) )
    ds['Meso']= ( ds[var].where(ds.WaveNumb > WaveSplit3//1).sum(dim='WaveNumb') 
                     +(1- WaveSplit3%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit3//1).fillna(0) )

if typ == 'WaveNR':
    ds['0']= ds[var].sel(WaveNumb= 0)
    ds['1-4']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 4), drop= True ).sum(dim='WaveNumb')
#    ds['>5']= ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb')
    ds['>5']= ds[var].where(ds.WaveNumb >= 5, drop= True ).sum(dim='WaveNumb')

            
for ci, cat in enumerate(catlist):
    
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    for vi, var in enumerate(varlist):    
        mean= ds[cat+'_'+var].mean(dim='time').mean(dim='Member')
        
        if vi == 0: axs[0].plot(ds.lat, fac* mean, color= color, label= cat)
        else: axs[0].plot(ds.lat, fac* mean, ls= '--', color= color)

        diffmean= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).mean(dim='Member')
        
        if vi == 0:  axs[1].plot(ds.lat, fac*diffmean, label= var, color= color)
        else:  axs[1].plot(ds.lat, fac*diffmean, ls= '--', color= color)

        if len(Memberlist) > 1:        
            std= ds[cat+'_'+var].mean(dim='time').std(dim='Member')        
            axs[0].plot(ds.lat, fac* (mean +std), label='', lw= 1, color= color )
            axs[0].plot(ds.lat, fac* (mean -std), label='', lw= 1, color= color )
          
            diffstd= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).std(dim='Member')
            axs[1].plot(ds.lat, fac*(diffmean+ diffstd) , lw= 1, color= color)
            axs[1].plot(ds.lat, fac*(diffmean- diffstd) , lw= 1, color= color)



#ds.where(ds["time.year"]< 2000, drop= True)
    
axs[0].legend(ncol= 2)

axs[0].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[1].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
if unit =='PerM': axs[0].set_ylabel('Mean energy transport \n per lenth unit [W/m]') #energy transport accross a latitude band
elif unit == 'LatCycl':  axs[0].set_ylabel('Mean energy transport \n across latitude circle [W]')#energy transport accross a latitude band

#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel('Energy transport change \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))


if save:
    savedir= '../Figs/Global/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
