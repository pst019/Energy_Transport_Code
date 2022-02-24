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
model='ERA5'


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
    syear= 1950
    eyear= 2100
    Memberlist= [1]#, 2, 3, 4, 5]#,2,3,4]


elif model== 'ERA5':
    syear= 1979
    eyear= 2018
    Memberlist=[1]

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 6



timeperiod= 'year'
# timeperiod= 'DJF'
# timeperiod= 'JJA'


if timeperiod == 'year': monthlist= np.arange(1,13)
elif timeperiod == 'DJF': monthlist= [1,2,12]
elif timeperiod == 'JJA': monthlist= [6,7,8]


#typ = 'Eddies'
typ = 'Wavelength_smooth' # corrected
# typ = 'Wavelength_smooth_2' #corrected
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
    varlist= ['vEtot', 'vEsetot']
    # varlist= ['vQtot', 'vQsetot']
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


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc





if timeperiod== 'year': ds= ds.resample(time='1Y').mean()
print('only for annual timeperiod')
# .groupby('time.season')



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





"""fraction of wave 0"""
fig= plt.figure(fignr, figsize= (8, 2.5))
fignr+=1
plt.clf()

transp = xr.Dataset()

period_list=['all', 'first', 'last']

for period in period_list:
    if period == 'last': dsn= ds.where(ds["time.year"]>= split_2)
    elif period == 'first': dsn= ds.where(ds["time.year"]<= split_1)
    else: dsn = ds
    
    transp['tot_'+period]= dsn[varlist[0]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
    transp['stat_'+period]= dsn[varlist[1]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
    transp['transi_'+period]= transp['tot_'+period]- transp['stat_'+period]
    transp['sum_'+period]= transp['tot_'+period].sum(dim='WaveNumb')
    
for var  in ['tot', 'stat', 'transi']:
    for period in period_list:
        low_transp_filter= transp['sum'+'_'+period]/transp['sum'+'_'+period].max() > 0.05

        transp[var+'_'+period]= transp[var+'_'+period].where(low_transp_filter, np.nan)


color = next(plt.gca()._get_lines.prop_cycler)['color']
for ip, period in enumerate(period_list):
    if ip == 0:  plt.plot(ds.lat,  transp['tot'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], label='Total', color= color)
    elif ip == 1: plt.plot(ds.lat,  transp['tot'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], color= color, ls= '--', label= period)
    else: plt.plot(ds.lat,  transp['tot'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], color= color, ls= ':', label= period)


# plt.plot(ds.lat,  transp['stat'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], label='Stationary')

color = next(plt.gca()._get_lines.prop_cycler)['color']
for ip, period in enumerate(period_list):
    if ip == 0: plt.plot(ds.lat,  transp['transi'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], label= 'Transient', color= color)
    elif ip == 1: plt.plot(ds.lat,  transp['transi'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], color= color, ls= '--')
    else: plt.plot(ds.lat,  transp['transi'+'_'+period].sel(WaveNumb= 0)/transp['sum'+'_'+period], color= color, ls= ':')

# plt.plot(ds.lat,  transp['tot'].sel(WaveNumb= 0)/transp['sum'], label='Total')
# plt.plot(ds.lat,  transp['stat'].sel(WaveNumb= 0)/transp['sum'], label='Stationary')
# plt.plot(ds.lat,  transp['transi'].sel(WaveNumb= 0)/transp['sum'], label= 'Transient')
plt.legend(ncol= 2)
plt.ylabel('Fraction of wave 0 \n on total transport') 
plt.ylim([-1, 1])
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')

if save:
    savedir= '../Figs/Scale+PowerSpec/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Fraction_wave_0_' +model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 


"""mean scale of max transport"""
fig= plt.figure(fignr, figsize= (8, 5))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')


transp = xr.Dataset()

    
transp['tot']= ds[varlist[0]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
transp['stat']= ds[varlist[1]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
transp['transi']= transp['tot']- transp['stat']
transp['sum']= transp['tot'].sum(dim='WaveNumb')


low_transp_filter= transp['sum']/transp['sum'].max() > 0.05
low_transp_filter_transi= transp['transi'].sum(dim='WaveNumb')/transp['sum'].max() > 0.05


axs[0].plot(ds.lat,  transp['tot'].argmax(dim='WaveNumb').where(low_transp_filter, np.nan), label= 'Total' )
axs[0].plot(ds.lat,  transp['stat'].argmax(dim='WaveNumb').where(low_transp_filter, np.nan), label= 'Stationary' )
axs[0].plot(ds.lat,  transp['transi'].argmax(dim='WaveNumb').where((low_transp_filter) & (low_transp_filter_transi), np.nan), label= 'Transient' )

rolling_interval = 10 #in degree latitude
rolling_values= int(rolling_interval/(ds.lat[0].values- ds.lat[1].values))
axs[1].plot(ds.lat,  (LatCirc/transp['tot'].argmax(dim='WaveNumb').where(low_transp_filter, np.nan)).rolling(lat= rolling_values, center= True, min_periods= .8*rolling_values).mean()/1E3 )
axs[1].plot(ds.lat,  (LatCirc/transp['stat'].argmax(dim='WaveNumb').where(low_transp_filter, np.nan)).rolling(lat= rolling_values, center= True, min_periods= .8*rolling_values).mean()/1E3 )
axs[1].plot(ds.lat,  (LatCirc/transp['transi'].argmax(dim='WaveNumb').where((low_transp_filter) & (low_transp_filter_transi), np.nan)).rolling(lat= rolling_values, center= True, min_periods= .8*rolling_values).mean()/1E3 )

axs[0].legend()
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel('Wavenumber')
axs[1].set_ylabel('Wavelength [km]')
axs[1].plot([-latcut, latcut], [8E3, 8E3], '--', c= 'k', lw= 1)
axs[1].set_ylim(bottom=0)

if save:
    savedir= '../Figs/Scale+PowerSpec/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Scale_max_transport_' +model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 



"""change in the scale of max transport"""

fig= plt.figure(fignr, figsize= (8, 5))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')


for period in period_list:
    if period == 'last': dsn= ds.where(ds["time.year"]>= split_2)
    elif period == 'first': dsn= ds.where(ds["time.year"]<= split_1)
    else: dsn = ds
    
    transp['tot_'+period]= dsn[varlist[0]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
    transp['stat_'+period]= dsn[varlist[1]].mean(dim=['time', 'Member']) *np.sign(ds.lat)
    transp['transi_'+period]= transp['tot_'+period]- transp['stat_'+period]
    transp['sum_'+period]= transp['tot_'+period].sum(dim='WaveNumb')
    
# transp['tot']= ds['vEtot'].mean(dim=['time', 'Member']) *np.sign(ds.lat)
# transp['stat']= ds['vEsetot'].mean(dim=['time', 'Member']) *np.sign(ds.lat)
# transp['transi']= transp['tot']- transp['stat']
# transp['sum']= transp['tot'].sum(dim='WaveNumb')


for var  in ['tot', 'stat', 'transi']:
    for period in period_list:
        low_transp_filter= transp['sum'+'_'+period]/transp['sum'+'_'+period].max() > 0.05
        transp[var+'_'+period]= transp[var+'_'+period].where(low_transp_filter, 0)

for var  in ['tot', 'stat', 'transi']:
    for period in period_list:
        transp[var+'_'+period]= transp[var+'_'+period].where(transp['sum'+'_'+period]/transp['sum'+'_'+period].max() > 0.05, 0)

color = next(plt.gca()._get_lines.prop_cycler)['color']
for ip, period in enumerate(period_list[1:]):
    if ip == 0: axs[0].plot(ds.lat,  transp['tot'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= '-', label= 'Total '+ period)
    else: axs[0].plot(ds.lat,  transp['tot'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= ':', label='Total '+ period)    

    if ip == 0: axs[1].plot(ds.lat,   (LatCirc/transp['tot'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= '-', label= 'Total '+ period)
    else: axs[1].plot(ds.lat,   (LatCirc/transp['tot'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= ':', label='Total '+ period)    

color = next(plt.gca()._get_lines.prop_cycler)['color']
for ip, period in enumerate(period_list[1:]):
    if ip == 0: axs[0].plot(ds.lat,  transp['stat'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= '-', label= 'Stationary' + period)
    else: axs[0].plot(ds.lat,  transp['stat'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= ':', label= 'Stationary' + period)    

    if ip == 0: axs[1].plot(ds.lat,   (LatCirc/transp['stat'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= '-', label= 'Stationary ' + period)
    else: axs[1].plot(ds.lat,   (LatCirc/transp['stat'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= ':', label= 'Stationary ' + period)    
    
    
color = next(plt.gca()._get_lines.prop_cycler)['color']
for ip, period in enumerate(period_list[1:]):
    if ip == 1: axs[0].plot(ds.lat,  transp['transi'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= '-', label= 'Transient' + period)
    else: axs[0].plot(ds.lat,  transp['transi'+'_'+period].argmax(dim='WaveNumb'), color= color, ls= ':', label= 'Transient' + period)        

    if ip == 1: axs[1].plot(ds.lat,   (LatCirc/transp['transi'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= '-', label= 'Transient ' + period)
    else: axs[1].plot(ds.lat,   (LatCirc/transp['transi'+'_'+period].argmax(dim='WaveNumb')).rolling(lat= 10, center= True, min_periods= 8).mean()/1E3, color= color, ls= ':', label= 'Transient ' + period)      
    
    
axs[1].legend(ncol= 3)
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel('Wavenumber')
axs[1].set_ylabel('Wavelength [km]')
axs[1].plot([-latcut, latcut], [8E3, 8E3], '--', c= 'k', lw= 1)
axs[1].set_ylim(bottom=0)

if save:
    savedir= '../Figs/Scale+PowerSpec/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Change_scale_max_transport_' +model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 