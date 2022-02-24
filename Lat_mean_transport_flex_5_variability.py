#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd

user = os.getcwd().split('/')[2]



model='EC_Earth'
# model='ERA5'


if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/'+ model +'/'
elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/'+ model +'/'    
else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
if model=='ERA5': Mediadir += 'EnergySplit/res_0.5x0.5/Waves/'




save= True
# save= False
imp= False
imp= True
#
#Member = 1

if model == 'EC_Earth':
    syear= 1979 #50
    eyear= 2018 #100

elif model== 'ERA5':
    syear= 1979
    eyear= 2018

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 6


Memberlist= [3]#,2,3,4]

timeperiod= 'year'
# timeperiod= 'DJF'
# timeperiod= 'JJA'


if timeperiod == 'year': monthlist= np.arange(1,13)
elif timeperiod == 'DJF': monthlist= [1,2,12]
elif timeperiod == 'JJA': monthlist= [6,7,8]


#typ = 'Eddies'
# typ = 'Wavelength_smooth' # corrected
typ = 'Wavelength_smooth_2' #corrected
# typ = 'Wavelength_smooth_3' #not corrected

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
    # varlist= ['vEtot', 'vEsetot']
    varlist= ['vQtot', 'vQsetot']
    monthlymeanlist=[False, True]

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
                for month in monthlist:
    #                print(month)
                    if year == syear and month== monthlist[0]:
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
                    if year == eyear and month == monthlist[-1]:
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


if timeperiod== 'year': ds= ds.resample(time='1Y').mean()
print('only for annual timeperiod')
# .groupby('time.season')

a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc




"""Calculate transport of each category"""
if typ == 'Wavelength_smooth':
    Split1= 8000E3
    WaveSplit1= LatCirc/Split1
    
    for var in varlist:

        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
        ds['Syno_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit1//1).sum(dim='WaveNumb') 
                     +(1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    

if typ == 'Wavelength_smooth_2':
    Split1= 8000E3
    Split2= 2000E3
    
    WaveSplit1= LatCirc/Split1
    WaveSplit2= LatCirc/Split2
    WaveSplit2[WaveSplit2 >= 20]= 19.99
    
    for var in varlist:
        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                            +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
        ds['Syno_'+var]=  ( ds[var].where(np.logical_and(ds.WaveNumb > 1+ WaveSplit1//1, ds.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                          + (1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                          + WaveSplit2%1 * ds[var].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) )
        ds['Meso_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit2//1).sum(dim='WaveNumb') 
                          +(1- WaveSplit2%1) * ds[var].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) )




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




"""Plot the Total transport including changes"""
fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col')




for vi, var in enumerate(varlist):
    if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
    if typ in ['Wavelength_smooth','Wavelength_smooth_2','Wavelength_smooth_3', 'WaveNR']: dstot = ds[var].sum(dim= "WaveNumb")

    mean= dstot.mean(dim='time').mean(dim='Member')
    diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).mean(dim='Member')
    
    if vi == 0:
        line, = axs[0].plot(ds.lat, fac*mean, label= 'Total', color= 'k')
        axs[1].plot(ds.lat, fac*diffmean, color= 'k')
        axs[2].plot(ds.lat, fac * dstot.std(dim='time').mean(dim='Member') , color= 'k')
        axs[3].plot(ds.lat, fac * (dstot.where(ds["time.year"]>= split_2).std(dim='time').mean(dim='Member')
                                   - dstot.where(ds["time.year"]<= split_1).std(dim='time').mean(dim='Member') ) , color= 'k')
        
        
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
    

""" Plot the components"""           
for ci, cat in enumerate(catlist):
    
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    for vi, var in enumerate(varlist):    
        mean= ds[cat+'_'+var].mean(dim='time').mean(dim='Member')
        
        if vi == 0: axs[0].plot(ds.lat, fac* mean, color= color, label= cat)
        elif cat != 'Meri': axs[0].plot(ds.lat, fac* mean, ls= '--', color= color, label= cat +' stat')

        diffmean= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).mean(dim='Member')
        
        if vi == 0:  axs[1].plot(ds.lat, fac*diffmean, label= var, color= color)
        elif cat != 'Meri':  axs[1].plot(ds.lat, fac*diffmean, ls= '--', color= color)

        if len(Memberlist) > 1:        
            std= ds[cat+'_'+var].mean(dim='time').std(dim='Member')        
            axs[0].plot(ds.lat, fac* (mean +std), label='', lw= 1, color= color )
            axs[0].plot(ds.lat, fac* (mean -std), label='', lw= 1, color= color )
          
            diffstd= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).std(dim='Member')
            axs[1].plot(ds.lat, fac*(diffmean+ diffstd) , lw= 1, color= color)
            axs[1].plot(ds.lat, fac*(diffmean- diffstd) , lw= 1, color= color)


        #gets the standard deviation
        if vi == 0: axs[2].plot(ds.lat, fac * ds[cat+'_'+var].std(dim='time').mean(dim='Member') , color= color)
        elif cat != 'Meri':  axs[2].plot(ds.lat, fac* ds[cat+'_'+var].std(dim='time').mean(dim='Member') , ls= '--', color= color)
        
        #change in standard deviation
        if vi == 0: axs[3].plot(ds.lat, fac* (ds[cat+'_'+var].where(ds["time.year"]>= split_2).std(dim='time')
                             - ds[cat+'_'+var].where(ds["time.year"]<= split_1).std(dim='time')).mean(dim='Member') , color= color )
        elif cat != 'Meri': axs[3].plot(ds.lat, fac* (ds[cat+'_'+var].where(ds["time.year"]>= split_2).std(dim='time')
                             - ds[cat+'_'+var].where(ds["time.year"]<= split_1).std(dim='time')).mean(dim='Member') , ls= '--', color= color )

#ds.where(ds["time.year"]< 2000, drop= True)

    
axs[0].legend(ncol= 3)
for axnr in range(4):  axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
if unit =='PerM': axs[0].set_ylabel('Mean energy transport \n per lenth unit [W/m]') #energy transport accross a latitude band
elif unit == 'LatCycl':  axs[0].set_ylabel('Mean energy transport \n across latitude circle [W]')#energy transport accross a latitude band

#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel('Energy transport change \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))
axs[2].set_ylabel('Annual variability \n in the energy transport')
axs[3].set_ylabel('Change in annual variability \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))


if save:
    savedir= '../Figs/Global_2/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 




"""Plot the Total transport and variability"""
fig= plt.figure(fignr, figsize= (8, 5.5))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')




for vi, var in enumerate(varlist):
    if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
    if typ in ['Wavelength_smooth','Wavelength_smooth_2','Wavelength_smooth_3', 'WaveNR']: dstot = ds[var].sum(dim= "WaveNumb")

    mean= dstot.mean(dim='time').mean(dim='Member')
    diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).mean(dim='Member')
    
    if vi == 0:
        line, = axs[0].plot(ds.lat, fac*mean, label= 'Total', color= 'k')
        axs[1].plot(ds.lat, fac * dstot.std(dim='time').mean(dim='Member') , color= 'k')
        
        
    # else:
    #     line, = axs[0].plot(ds.lat, fac*mean, '--' , color= 'k')
    #     axs[1].plot(ds.lat, fac*diffmean, color= 'k')
    
    if len(Memberlist) > 1:
        std= dstot.mean(dim='time').std(dim='Member') 
        axs[0].plot(ds.lat, fac*(mean + std), label='', lw= 1, color= line.get_color() )
        axs[0].plot(ds.lat, fac*(mean - std), label='', lw= 1, color= line.get_color() )
   
   

""" Plot the components"""           
for ci, cat in enumerate(catlist):
    
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    for vi, var in enumerate(varlist):    
        mean= ds[cat+'_'+var].mean(dim='time').mean(dim='Member')
        
        if vi == 0: axs[0].plot(ds.lat, fac* mean, color= color, label= cat)
        elif cat != 'Meri': axs[0].plot(ds.lat, fac* mean, ls= '--', color= color, label= cat +' stat')


        if len(Memberlist) > 1:        
            std= ds[cat+'_'+var].mean(dim='time').std(dim='Member')        
            axs[0].plot(ds.lat, fac* (mean +std), label='', lw= 1, color= color )
            axs[0].plot(ds.lat, fac* (mean -std), label='', lw= 1, color= color )
          

        #gets the standard deviation
        if vi == 0: axs[1].plot(ds.lat, fac * ds[cat+'_'+var].std(dim='time').mean(dim='Member') , color= color)
        elif cat != 'Meri':  axs[1].plot(ds.lat, fac* ds[cat+'_'+var].std(dim='time').mean(dim='Member') , ls= '--', color= color)
        

    
axs[0].legend(ncol= 3)
for axnr in range(2):  axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')


plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
if unit =='PerM': axs[0].set_ylabel('Mean energy transport \n per lenth unit [W/m]') #energy transport accross a latitude band
elif unit == 'LatCycl':  axs[0].set_ylabel('Mean energy transport \n across latitude circle [W]')#energy transport accross a latitude band

#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel('Annual variability \n in the energy transport')


if save:
    savedir= '../Figs/Global_2/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Mean+Var_'+ model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 
