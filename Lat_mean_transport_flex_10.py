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
from scipy.stats import permutation_test
from scipy import stats

from functions import *
#import functions

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
    # Mediadir= '/media/'+user+'/Backup1/'
    Mediadir= '/media/'+user+'/Backup/'

elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/'    
else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/'
#'/nird/projects/NS9063K/Richard_KNMI'

Mediadir_0= Mediadir    
Mediadir += 'data/Energy_Transport/'+ model +'/'
if model=='ERA5': Mediadir += 'EnergySplit/res_0.5x0.5/Waves/'


fignr= 8

scale= 'sine'
scale=''

fix_y_axis= False
trans_var= False # plot the transient component - not included in many panels - is specified in the loop for the transient transport, needed here for an if statement of the legend
plot_large= False

timeperiod= 'year'
# timeperiod= 'DJF'
# timeperiod= 'JJA'
# timeperiod= 'MAM'
# timeperiod= 'SON'

rolling= True # a rolling filter on the latitude of 5 datapoints - maybe remove - but might be useful for the convergence

typ = 'Wavelength_smooth' # corrected
# typ = 'Wavelength_smooth_2' #corrected
# typ = 'Wavelength_smooth_3' #not corrected
# typ = 'Wavelength_smooth_6tsd'
# typ = 'Wavelength_smooth_7tsd'
# typ = 'Wavelength_smooth_10tsd'
# typ, rolling = 'Wavelength_strict', False # corrected
# typ, rolling = 'WaveNR', False
#typ = 'Eddies'


energytyp= 'E' #total energy
# energytyp= 'E_nostat' #total energy - no stationary contributions
# energytyp= 'Q'
# energytyp= 'D' # dry static

####'Mean' , 'Var', 'Change-Mean', 'Change-Var', 'Conv', 'Change-Conv', 'Change-Fraction', 'Change-Frac-Var'
# pannellist= ['Mean']
# pannellist= ['Mean' , 'Var'] #, 'Change-Mean', 'Change-Var', 'Conv', 'Change-Conv', 'Change-Fraction', 'Change-Frac-Var'
# pannellist= ['Mean' , 'Change-Mean'] #, 'Change-Var', 'Conv', 'Change-Conv'
# pannellist= ['Mean' , 'Conv'] 
# pannellist, fix_y_axis= ['Change-Mean-ttest'], [[-8E6, 12E6]]

# pannellist= ['Change-Fraction-ttest']

# pannellist= ['Change-Conv-ttest']
# pannellist= ['Change-Var-ftest']

# pannellist= [ 'Change-Frac-Var-ftest']

# pannellist= ['Change-Var-Roll']
# pannellist= ['Change-Var-Roll', 'Change-Frac-Var-Roll']
# pannellist= ['Var', 'Change-Var']


# pannellist, fix_y_axis, plot_large= ['Mean'], [2E8], True
# pannellist, timeperiod, fix_y_axis, plot_large = ['Mean'], 'DJF', [[-1.8E8, 2.5E8]], True
# pannellist, timeperiod, fix_y_axis, plot_large = ['Mean'], 'JJA', [[-2.5E8, 1.8E8]], True
# pannellist, fix_y_axis= ['Mean' , 'Var'], [[-2E8, 2E8], [0, 1.15E7]]  #for EC-Earth

# pannellist, fix_y_axis= ['Change-Mean-ttest', 'Change-Fraction-ttest', 'Change-Conv-ttest'], [[-2E7, 2E7], [-.25, .35], [-15, 16]]
# pannellist, fix_y_axis= ['Change-Mean-ttest'], [[-2E7, 2E7]]


pannellist, timeperiod, fix_y_axis= ['Change-Mean-ttest', 'Change-Fraction-ttest'], 'DJF', [[-2.2E7, 4E7], [-.4, .4]]
# pannellist, fix_y_axis=  ['Change-Var-ftest', 'Change-Frac-Var-ftest'], [[-2.5E6, 3E6], [-.35, .5]]

# pannellist, fix_y_axis= ['Var', 'Var-Fraction'], [[0, 1.2E7], [0, .4]]
# pannellist, fix_y_axis= ['Var'], [[0, 1.2E7]]
# pannellist, scale, fix_y_axis= ['Conv'], 'sine', [[-80, 120]]


save= False
save= True
savedir= Mediadir_0 +'/home/Energy_Transport/Figs/Global_5/'



if model == 'EC_Earth':
    syear= 1950
    eyear= 2100
    # syear= 1979
    # eyear= 2018    
    Memberlist= [1, 2, 3, 4, 5]


elif model== 'ERA5':
    syear= 1979
    eyear= 2018
    Memberlist=[1]

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)



if timeperiod == 'year': monthlist= np.arange(1,13)
elif timeperiod == 'DJF': monthlist= [1,2,12]
elif timeperiod == 'JJA': monthlist= [6,7,8]
elif timeperiod == 'MAM': monthlist= [3,4,5]
elif timeperiod == 'SON': monthlist= [9,10,11]


latcut= 89

unit='PerM'
# unit='LatCycl'

if typ == 'Eddies':
    if energytyp== 'E': varlist= ['vEmc', 'vEse', 'vEte']
    elif energytyp== 'Q': varlist= ['vQmc', 'vQse', 'vQte']
    # elif energytyp== 'D': implist= ['vEtot', 'vEsetot', 'vQtot', 'vQsetot']
    
    catlist= varlist
    monthlymean = True

elif 'Wavelength' in typ or typ == 'WaveNR':
    if energytyp== 'E': varlist, monthlymeanlist= ['vEtot', 'vEsetot'], [False, True]
    elif energytyp== 'E_nostat': varlist, monthlymeanlist= ['vEtot'], [False]
    elif energytyp== 'Q': varlist, monthlymeanlist= ['vQtot', 'vQsetot'], [False, True]
    elif energytyp== 'D': implist, monthlymeanlist= ['vEtot', 'vEsetot', 'vQtot', 'vQsetot'], [False, True, False, True]
 
    if typ == 'Wavelength_smooth':    
        catlist = ['Meri', 'Plan', 'Syno']
        Split1= 8000E3

    if typ == 'Wavelength_smooth_2':    
        catlist = ['Meri', 'Plan', 'Syno', 'Meso']
        Split1= 8000E3
        Split2= 2000E3 

    if typ == 'Wavelength_smooth_6tsd':    
        catlist = ['Meri', 'Plan', 'Syno', 'Meso']
        Split1= 6000E3
        Split2= 2000E3        

    if typ == 'Wavelength_smooth_7tsd':    
        catlist = ['Meri', 'Plan', 'Syno', 'Meso']
        Split1= 7000E3
        Split2= 2000E3  
        
    if typ == 'Wavelength_smooth_10tsd':    
        catlist = ['Meri', 'Plan', 'Syno', 'Meso']
        Split1= 10000E3
        Split2= 2000E3        
        
    if typ == 'Wavelength_smooth_3':    
        catlist = ['Meri', 'Plan', 'Syno_l', 'Syno_s', 'Meso']
        Split1= 10000E3
        Split2= 4000E3
        Split3= 2000E3        
        
    if typ == 'Wavelength_strict':    
        catlist = ['Meri', 'Plan', 'Syno']
        Split1= 8000E3
   
    if typ == 'WaveNR':
        # catlist = ['0', '1-4', '>5']#, '>10']
        catlist = ['Meri', '1-3', '>4']#, '>10']





if energytyp in ['E', 'Q', 'E_nostat']:  implist= varlist
elif energytyp== 'D': varlist= ['vDtot', 'vDsetot']

if pannellist == ['Conv']: varlist= varlist[:1]



intermediate_file_name= Mediadir_0 +'data/Energy_Transport/'+ model +'/intermediate_data/'+timeperiod+'-mean_'
for var in implist: intermediate_file_name += var+'_'
intermediate_file_name += str(syear)+'-'+str(eyear)
for Member in Memberlist: intermediate_file_name += '_M'+str(Member)
intermediate_file_name +=  '.nc'


if os.path.exists(intermediate_file_name):
    ds= xr.open_dataset(intermediate_file_name)
    
else:
    print('Import ...')
    for Member in Memberlist:
        file_dir= Mediadir
        if model== 'EC_Earth': 
            print('Member', Member) 
            file_dir += 'Member'+str(Member) +'/'            
        
        for vi, var in enumerate(implist):
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
                        if var== implist[0]:
                            ds2= ds0
                        else:
                            ds2= xr.merge([ds2, ds0])
                            
        if Member == Memberlist[0]:
            ds= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
        else:
            ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            ds= xr.merge([ds, ds2])


    if energytyp== 'D':
        ds['vDtot']= ds['vEtot'] - ds['vQtot']
        ds['vDsetot']= ds['vEsetot'] - ds['vQsetot']
        ds= ds.drop(['vEtot', 'vQtot', 'vEsetot', 'vQsetot'])

    ds.to_netcdf(intermediate_file_name)




if latcut != 90: ds= ds.where(np.abs(ds.lat) < latcut, drop= True)

a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc

if timeperiod== 'year': ds= ds.resample(time='1Y').mean()
print('only for annual timeperiod')
# .groupby('time.season')


####here the sinescale function was before


"""Calculate transport of each category"""
if typ == 'Wavelength_smooth':
    WaveSplit1= LatCirc/Split1
    
    for var in varlist:
        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
        ds['Syno_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit1//1).sum(dim='WaveNumb') 
                     +(1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    

if typ in ['Wavelength_smooth_2', 'Wavelength_smooth_6tsd', 'Wavelength_smooth_7tsd', 'Wavelength_smooth_10tsd']:  
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
    WaveSplit1= LatCirc/Split1
    WaveSplit2= LatCirc/Split2
    WaveSplit3= LatCirc/Split3
    WaveSplit3[WaveSplit3 >= 20]= 19.99

    for var in varlist:
        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                            +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
        ds['Syno_l_'+var]=  ( ds[var].where(np.logical_and(ds.WaveNumb > 1+ WaveSplit1//1, ds.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                          + (1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                          + WaveSplit2%1 * ds[var].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) )
        ds['Syno_s_'+var]=  ( ds[var].where(np.logical_and(ds.WaveNumb > 1+ WaveSplit2//1, ds.WaveNumb < 1+ WaveSplit3//1)).sum(dim='WaveNumb') 
                          + (1- WaveSplit2%1) * ds[var].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) 
                          + WaveSplit3%1 * ds[var].sel(WaveNumb= 1+ WaveSplit3//1).fillna(0) )
        ds['Meso_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit3//1).sum(dim='WaveNumb') 
                          +(1- WaveSplit3%1) * ds[var].sel(WaveNumb= 1+ WaveSplit3//1).fillna(0) )



    # ds['Meri']= ds[var].sel(WaveNumb= 0)
    # ds['Plan']= ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb < WaveSplit1//1)).sum(dim='WaveNumb')+  WaveSplit1%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0)
    # ds['Syno_l']=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit1)//1, ds.WaveNumb < WaveSplit2//1)).sum(dim='WaveNumb') 
    #                  + (1- WaveSplit1%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit1//1).fillna(0) 
    #                  + WaveSplit2%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) )
    # ds['Syno_s']=  ( ds[var].where(np.logical_and(ds.WaveNumb > (WaveSplit2)//1, ds.WaveNumb < WaveSplit3//1)).sum(dim='WaveNumb') 
    #                  + (1- WaveSplit2%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit2//1).fillna(0) 
    #                  + WaveSplit3%1 * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit3//1).fillna(0) )
    # ds['Meso']= ( ds[var].where(ds.WaveNumb > WaveSplit3//1).sum(dim='WaveNumb') 
    #                  +(1- WaveSplit3%1) * ds[var].where(ds.WaveNumb >= 1).sel(WaveNumb= WaveSplit3//1).fillna(0) )

if typ == 'Wavelength_strict':
    WaveSplit1= LatCirc/Split1
    
    for var in varlist:
        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
        ds['Syno_'+var]= ds[var].where(ds.WaveNumb >= 1+ WaveSplit1//1).sum(dim='WaveNumb') 

if typ == 'WaveNR':
    ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
    ds['1-3_'+var]=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 3), drop= True ).sum(dim='WaveNumb')
#    ds['>5']= ds[var].where(np.logical_and(ds.WaveNumb >= 5, ds.WaveNumb <= 9), drop= True ).sum(dim='WaveNumb')
    ds['>4_'+var]= ds[var].where(ds.WaveNumb >= 4, drop= True ).sum(dim='WaveNumb')



if rolling:
    ds= ds.rolling(lat= 5, center= True).mean()
    ds= ds.dropna(dim= 'lat', how= 'all')



def stack_ds(ds, dim= ("time", "Member")):
    return np.array(ds.stack(x= dim) )


def split_stack_ds(ds, split_2, split1, dim= ("time", "Member")):
    return (np.array(ds.where(ds["time.year"]>= split_2, drop= True).stack(x= ("time", "Member")) ),
            np.array(ds.where(ds["time.year"]<= split_1, drop= True).stack(x= ("time", "Member")) ) )



def mean_diff(ds1, ds2):
    return np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)

def conv_diff(ds1, ds2, dslat= ds.lat):
    from scipy import ndimage

    diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
    diff *= np.cos(np.deg2rad(dslat.values)) #.mean(dim='Member')    
    # conv_diff= -1* diff.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(dslat.values))
    conv_diff= np.gradient(diff, axis= 0)/110E3 * 1/np.cos(np.deg2rad(dslat.values))
    
    # conv_diff= ndimage.uniform_filter1d(conv_diff, 5) # a runnig mean filter

    return conv_diff


def frac_diff(ds1, ds2):
    """difference/ mean"""
    m1= np.mean(ds1, axis= -1)
    m2= np.mean(ds2, axis= -1)
    return (m1 - m2) /(0.5* (m1 + m2))

# def moving_average(x, w):
#     return np.convolve(x, np.ones(w), 'valid') / w



def sign_change(ds1, dslat= ds.lat):
    #ds1 is a stacked data variable
    #dslat is just ds.lat
    from scipy.signal import argrelextrema
    from scipy import ndimage


    ds1_s= np.sign(ds1)
    sign_change= (np.roll(ds1_s, 1, axis= 1)- ds1_s != 0).astype(int)
    sign_change[0]= 0
    # return sign_change

    sign_change_mean= sign_change.mean(axis= 1)
    locmax= argrelextrema(sign_change_mean, np.greater)
    sign_change_mean= ndimage.uniform_filter1d(sign_change_mean, 7)
    sign_change_mean[sign_change_mean<1E-7]=0
    return dslat.values[locmax], sign_change_mean

p_level= 0.05




"""Plot the pannels from the pannellist"""
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]


if plot_large: fig= plt.figure(fignr, figsize= (8, 5*len(pannellist)))
else: fig= plt.figure(fignr, figsize= (8, 3.5*len(pannellist)))
fignr+=1
plt.clf()

axs=fig.subplots(len(pannellist), 1, sharex= 'col')

if len(pannellist)== 1: axs= [axs]

"""calculation relevant for all pannels"""
if typ == 'Eddies':
    dstot= sum(ds[var] for var in varlist)
    
if 'Wavelength' in typ or typ in ['WaveNR']:
    var= varlist[0]
    dstot = ds[var].sum(dim= "WaveNumb")


for ip, ptype in enumerate(pannellist):


#pannellist= ['Mean'] #'Var', 'Change-Mean', 'Change-Var', 'Conv', 'Change-Conv'

    if ptype == 'Mean':
        mean= dstot.mean(dim='time').mean(dim='Member')#.rolling(lat= 5, center= True).mean()
        axs[ip].plot(ds.lat, fac*mean, label= 'Total', color= 'k')

        if len(Memberlist) > 1:
            mean_std= dstot.mean(dim='time').std(dim='Member') 
            axs[ip].fill_between(ds.lat, fac*(mean - mean_std),  fac*(mean + mean_std), color= 'k', alpha= 0.2   )


        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                mean= ds[cat+'_'+var].mean(dim='time').mean(dim='Member')
                
                if vi == 0: axs[ip].plot(ds.lat, fac* mean, color= color, label= cat)               
                elif (typ == 'Wavelength_smooth_2' and cat != 'Meri') or cat == 'Plan':
                    axs[ip].plot(ds.lat, fac* mean, ls= '--', color= color, label= cat +' stat')
                # elif cat == 'Plan': axs[ip].plot(ds.lat, fac* mean, ls= '--', color= color, label= cat +' stat')


                if len(Memberlist) > 1:        
                    mean_std= ds[cat+'_'+var].mean(dim='time').std(dim='Member')        
                    axs[ip].fill_between(ds.lat, fac*(mean - mean_std),  fac*(mean + mean_std), color= color, alpha= 0.2   )

        if unit =='PerM': axs[ip].set_ylabel('Mean transport \n per lenth unit [W m$^{-1}$]') #energy transport accross a latitude band
        elif unit == 'LatCycl':  axs[ip].set_ylabel('Mean transport \n across latitude circle [W]')#energy transport accross a latitude band



    if ptype == 'Var':
        variability= dstot.std(dim='time').mean(dim='Member')
        axs[ip].plot(ds.lat, fac * variability , label= 'Total', color= 'k')
 
        if len(Memberlist) > 1:
            var_std= dstot.std(dim='time').std(dim='Member') 
            axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= 'k', alpha= 0.2   )

        
        if len(pannellist) == 1: trans_var = True #plot the transient variability
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                variability= ds[cat+'_'+var].std(dim='time').mean(dim='Member') #gets the standard deviation
                if vi == 0: axs[ip].plot(ds.lat, fac * variability, color= color, label=cat)
                elif cat not in  ['Meri', 'Syno']:  axs[ip].plot(ds.lat, fac* variability , ls= '--', color= color, label= cat +' stat')

                if len(Memberlist) > 1:        
                    var_std= ds[cat+'_'+var].std(dim='time').std(dim='Member')        
                    axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= color, alpha= 0.2   )


            if trans_var and cat not in  ['Meri', 'Syno']:
                variability= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time').mean(dim='Member')
                axs[ip].plot(ds.lat, fac * variability, ls= ':', color= color, label=cat + ' trans')
                if len(Memberlist) > 1:        
                    var_std= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time').std(dim='Member')        
                    axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= color, alpha= 0.2   )


                
        # axs[ip].set_ylabel('Annual variability \n transport [W m$^{-1}$]')
        axs[ip].set_ylabel('Transport variability \n [W m$^{-1}$]')


    if ptype == 'Change-Mean-ttest':
        """significance t-test"""
        ds1, ds2= split_stack_ds(dstot, split_2, split_1, dim= ("time", "Member"))
        
        diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
        pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]
        
        meanfilter= pval < p_level        
        
        axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff), label= 'Total', color= 'k', lw= 2)
        axs[ip].plot(ds.lat, fac* diff, color= 'k', lw= 0.5)
   
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                ds1, ds2= split_stack_ds(ds[cat+'_'+var], split_2, split_1, dim= ("time", "Member"))

                diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
                pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]
                
                if vi == 0:
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(pval > p_level, diff), label= cat, color= color, lw= 2)
                    axs[ip].plot(ds.lat, fac*diff, color= color, lw= 0.5)

                elif cat not in ['Meri', 'Syno']:
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(pval > p_level, diff), ls= '--', label= cat +' stat', color= color, lw= 2)
                    axs[ip].plot(ds.lat, fac* diff, ls= '--', color= color, lw=.5)

                    # axs[ip].plot(ds.lat, fac*diffmean, ls= '--', label= cat +' stat', color= color)


        axs[ip].set_ylabel('Transport change'#' \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1) 
                           +'\n [W m$^{-1}$]')



    if ptype == 'Change-Var-ftest': #the variance from the rolling mean #in test
        # p_level= .1
        anomalie= dstot- dstot.rolling(time= 30, center= True).mean().dropna(dim= 'time')
        ds1, ds2= split_stack_ds(anomalie, split_2, split_1, dim= ("time", "Member"))

        pval= np.array([stats.bartlett(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])])       
        # pval= np.array([stats.levene(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])])

        diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim=('time', 'Member'))
                        - anomalie.where(ds["time.year"]<= split_1).std(dim=('time', 'Member')) )#.mean(dim='Member')

        axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff_variability), label= 'Total', color= 'k', lw= 2)
        axs[ip].plot(ds.lat, fac* diff_variability, color= 'k', lw= 0.5)

    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                #change in standard deviation
                anomalie= ds[cat+'_'+var]- ds[cat+'_'+var].rolling(time= 30, center= True).mean().dropna(dim= 'time')
                ds1, ds2= split_stack_ds(anomalie, split_2, split_1, dim= ("time", "Member"))

                pval= np.array([stats.bartlett(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])]) 
                # pval= np.array([stats.levene(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])])

                diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim=('time', 'Member'))
                                  - anomalie.where(ds["time.year"]<= split_1).std(dim=('time', 'Member')) )
                
                if vi == 0:
                    axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff_variability), color= color, label= cat, lw= 2)
                    axs[ip].plot(ds.lat, fac*  diff_variability, color= color, lw= .5)

                elif cat not in  ['Meri', 'Syno']: #axs[ip].plot(ds.lat, fac* diff_variability , ls= '--', color= color, label= cat +' stat' )
                    axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff_variability), color= color, label= cat +' stat',ls= '--', lw= 2)
                    axs[ip].plot(ds.lat, fac*  diff_variability, color= color, ls= '--',lw= .5)
                # if len(Memberlist) > 1:                 
                #     diff_variability_std= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                #         - anomalie.where(ds["time.year"]<= split_1).std(dim='time') ).std(dim='Member')
                #     axs[ip].fill_between(ds.lat, fac*(diff_variability - diff_variability_std),  fac*(diff_variability + diff_variability_std), color= color, alpha= 0.2   )


        axs[ip].set_ylabel('Variability change\n' #' \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1)
                           +' [W m$^{-1}$]')


    if ptype == 'Conv':
        ####d(transport)/dy = divergence
        # tot= dstot.mean(dim='time') #.mean(dim='Member')      
        # conv= -1* tot.differentiate('lat')/110E3
        # conv= conv.rolling(lat= 9, center= True).mean()  #smooth
        # conv_mean= conv.mean(dim='Member')

        #### div= d(cos(lat) transport)/dy *1/cos(lat) - to also consider the convergence of latitude
        tot= dstot.mean(dim='time')
        tot*= np.cos(np.deg2rad(tot.lat)) #.mean(dim='Member')    
        conv= -1* tot.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(tot.lat)) 
        conv= conv.rolling(lat= 5, center= True).mean()  #smooth
        conv_mean= conv.mean(dim='Member')
        
        axs[ip].plot(ds.lat, conv_mean, label= 'Total', color= 'k')

        if len(Memberlist) > 1: 
            conv_std= conv.std(dim='Member')
            axs[ip].fill_between(ds.lat, (conv_mean - conv_std),  conv_mean + conv_std, color= 'k', alpha= 0.2   )
            
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist[:1]):    
                tot= ds[cat+'_'+var].mean(dim='time')
                tot*= np.cos(np.deg2rad(tot.lat)) #.mean(dim='Member')                    
                conv= -1* tot.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(tot.lat)) #np.gradient(mean)/dy #divergence
                conv= conv.rolling(lat= 5, center= True).mean()
                conv_mean= conv.mean(dim='Member')

                if vi == 0: axs[ip].plot(ds.lat, conv_mean, label= cat, color= color)
                elif cat not in  ['Meri', 'Syno']: axs[ip].plot(ds.lat, conv_mean, ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1: 
                    conv_std= conv.std(dim='Member')                   
                    axs[ip].fill_between(ds.lat, (conv_mean - conv_std),  conv_mean + conv_std,
                                    color= color, alpha= 0.2   )
                
        axs[ip].set_ylabel('Convergence [W m$^{-2}$]')


    if ptype == 'Change-Conv-ttest':
        ds_var= dstot.copy()
        ds_var *= np.cos(np.deg2rad(ds_var.lat))
        ds_var = -1* ds_var.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(ds_var.lat))
        ds_var= ds_var.rolling(lat= 5, center= True).mean() #smooth

        ds1, ds2= split_stack_ds(ds_var, split_2, split_1, dim= ("time", "Member"))
        
        diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
        pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]
        meanfilter= pval < p_level        
        
        axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff), label= 'Total', color= 'k', lw= 2)
        axs[ip].plot(ds.lat, fac* diff, color= 'k', lw= 0.5)


    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist[:1]):    
                ds_var= ds[cat+'_'+var].copy()
                ds_var *= np.cos(np.deg2rad(ds_var.lat))
                ds_var = -1* ds_var.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(ds_var.lat))
                ds_var= ds_var.rolling(lat= 5, center= True).mean() #smooth
        
                ds1, ds2= split_stack_ds(ds_var, split_2, split_1, dim= ("time", "Member"))
                
                diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
                pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]
                meanfilter= pval < p_level        
                
                axs[ip].plot(ds.lat, fac* np.ma.masked_where(pval > p_level, diff), label= cat, color= color, lw= 2)
                axs[ip].plot(ds.lat, fac* diff, color= color, lw= 0.5)

                # if vi == 0: axs[ip].plot(ds.lat, conv_diff_mean, label= cat, color= color)
                # elif cat not in  ['Meri', 'Syno']: axs[ip].plot(ds.lat, conv_diff_mean, ls= '--', color= color, label= cat +' stat' )

        axs[ip].set_ylabel('Convergence change\n [W m$^{-2}$]')



    if ptype== 'Var-Fraction':
        ds_mean= dstot.mean(dim=('time', 'Member'))
        
        lat_change, cross_filter= sign_change_ds(ds_mean)

        variability= dstot.std(dim=('time', 'Member'))
        variability_frac= variability/ np.abs(ds_mean)
        axs[ip].plot(ds.lat, fac* np.ma.masked_where(cross_filter!= 0, variability_frac), color= 'k', label='Total')

        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                ds_mean= ds[cat+'_'+var].mean(dim=('time', 'Member'))
                
                lat_change, cross_filter= sign_change_ds(ds_mean)

                variability= ds[cat+'_'+var].std(dim=('time', 'Member')) #gets the standard deviation
                variability_frac= variability/ np.abs(ds_mean)

                if vi == 0:  axs[ip].plot(ds.lat, fac* np.ma.masked_where(cross_filter!= 0, variability_frac), label= cat, color= color)
                elif cat not in  ['Meri', 'Syno']: axs[ip].plot(ds.lat, fac* np.ma.masked_where(cross_filter!= 0, variability_frac), label= cat, ls= '--', color= color)
                

        axs[ip].set_ylabel('Variability fraction') 



    if ptype== 'Change-Fraction-ttest': # (transport change) / (mean transport)
        filter_level= 0.1* np.abs(dstot.mean(dim='time').mean(dim='Member')).max()
        
        ds_var= dstot.copy()
        
        ds_mean= ds_var.mean(dim=('time', 'Member')) 
        meanfilter= np.abs(ds_mean) < filter_level        

        ds_var/= ds_mean
        ds1, ds2= split_stack_ds(ds_var, split_2, split_1, dim= ("time", "Member"))
        
        diff= np.mean(ds1, axis= 1) - np.mean(ds2, axis= 1)
        pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]
      
        lat_change, cross_filter= sign_change_ds(ds_mean)
        
        FILTER= np.logical_or(pval > p_level, meanfilter)
        FILTER= np.logical_or(FILTER, cross_filter)       
        
        axs[ip].plot(ds.lat, fac* np.ma.masked_where(FILTER!= 0, diff), label= 'Total', color= 'k', lw= 2)
        axs[ip].plot(ds.lat, fac* np.ma.masked_where(cross_filter!= 0, diff), color= 'k', lw= 0.5)


        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist): 
                ds_var= ds[cat+'_'+var].copy()
                
                ds_mean= ds_var.mean(dim=('time', 'Member'))
                meanfilter= np.abs(ds_mean) < filter_level        

                ds_var/= ds_mean
                ds1, ds2= split_stack_ds(ds_var, split_2, split_1, dim= ("time", "Member"))

                diff= np.mean(ds1, axis= -1) - np.mean(ds2, axis= -1)
                pval= stats.ttest_ind(ds1, ds2, axis= 1)[1]

                lat_change, cross_filter= sign_change_ds(ds_mean)
        
                FILTER= np.logical_or(pval > p_level, meanfilter)
                FILTER= np.logical_or(FILTER, cross_filter)     

                if vi == 0:
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(FILTER != 0, diff), label= cat, color= color, lw= 2)
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(cross_filter != 0, diff), color= color, lw= .5)
                    # axs[ip].plot(ds.lat, fac* diffmean, label= cat, color= color)

                elif cat not in ['Meri', 'Syno']:
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(FILTER != 0, diff), ls= '--', label= cat +' stat', color= color, lw= 2)
                    axs[ip].plot(ds.lat, fac*np.ma.masked_where(cross_filter != 0, diff), ls= '--', color= color, lw= .5)
                    # axs[ip].plot(ds.lat, fac*diffmean, ls= '--', label= cat +' stat', color= color)

        axs[ip].set_ylabel('Fraction change') #' \n'+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))



    if ptype== 'Change-Frac-Var-ftest': #(change transport variance) / (mean variance)
        anomalie= dstot- dstot.rolling(time= 30, center= True).mean().dropna(dim= 'time')
        
        ds1, ds2= split_stack_ds(anomalie, split_2, split_1, dim= ("time", "Member"))        
        pval= np.array([stats.bartlett(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])])       

        variability= anomalie.std(dim=('time','Member') )
        filter_level= 0.2* np.abs(variability).max()
        absfilter= np.abs(variability) < filter_level     
        FILTER= np.logical_or(pval > p_level, absfilter)


        diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim=('time', 'Member'))
                        - anomalie.where(ds["time.year"]<= split_1).std(dim=('time', 'Member') )) #.where(varfilter, np.nan)
        variability_frac= diff_variability/variability


        axs[ip].plot(ds.lat, fac * np.ma.masked_where(FILTER, variability_frac), label= 'Total', color= 'k', lw= 2)
        axs[ip].plot(ds.lat, fac *  variability_frac, color= 'k', lw= .5)


        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                anomalie= ds[cat+'_'+var]- ds[cat+'_'+var].rolling(time= 30, center= True).mean().dropna(dim= 'time')
                ds1, ds2= split_stack_ds(anomalie, split_2, split_1, dim= ("time", "Member"))        
                pval= np.array([stats.bartlett(ds1[n], ds2[n]).pvalue for n in range(ds1.shape[0])])   
                
                variability= anomalie.std(dim=('time', 'Member'))
                absfilter= np.abs(variability) < filter_level     
                FILTER= np.logical_or(pval > p_level, absfilter)

                diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim=('time', 'Member'))
                        - anomalie.where(ds["time.year"]<= split_1).std(dim=('time', 'Member')))
                variability_frac= diff_variability/variability

                if vi == 0:
                    axs[ip].plot(ds.lat, fac* np.ma.masked_where(FILTER, variability_frac), color= color, label= cat , lw= 2)
                    axs[ip].plot(ds.lat, fac* variability_frac, color= color, lw= .5)                    
                elif cat not in ['Meri', 'Syno']:
                    axs[ip].plot(ds.lat, fac* np.ma.masked_where(FILTER, variability_frac), color= color, label= cat  + ' stat', ls= '--', lw= 2)
                    axs[ip].plot(ds.lat, fac* variability_frac, color= color, ls= '--', lw= .5)




        axs[ip].set_ylabel('Fraction in variability change\n' #' \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1)
                            +' [W m$^{-1}$]')









    """Common stuff"""
    axs[ip].ticklabel_format(useMathText=True) #to make 10**n on the axis
    
    if type(fix_y_axis) == list:
        if fix_y_axis[ip]:
            if type(fix_y_axis[ip]) == float:
                axs[ip].set_ylim(-1* fix_y_axis[ip], fix_y_axis[ip])
            else:
                axs[ip].set_ylim(fix_y_axis[ip][0], fix_y_axis[ip][1])


if scale== 'sine':
    axs[ip].set_xscale('sine')
axs[ip].set_xticks(np.arange(-80,80.1, 20))
axs[ip].set_xticks(np.arange(-80,80.1, 10), minor=True)


if pannellist in [['Conv'], ['Change-Conv-ttest']] or energytyp== 'E_nostat':  axs[0].legend(ncol= 2)
elif typ in ['Wavelength_smooth_2', 'Wavelength_smooth_6tsd', 'Wavelength_smooth_7tsd', 'Wavelength_smooth_10tsd'] or trans_var: axs[0].legend(ncol= 4)
else:  axs[0].legend(ncol= 3)


for axnr in range(ip+1): 
    axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    # axs[axnr].set_ylim(-.2, .6)


# plt.xlim(-latcut+2, latcut-2)
plt.xlim(ds.lat.min(), ds.lat.max())
plt.xlabel('Latitude')




if save:
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = model +'_'
    for pannel in pannellist: file_name += pannel + '_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    if scale== 'sine': file_name += '_sinescale'    
    savefile= savedir+ file_name
    print('save: \n', savefile)
    plt.savefig(savefile , bbox_inches='tight') 



