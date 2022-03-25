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

params = {'axes.labelsize': 13,
          # 'axes.titlesize': 16,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'legend.fontsize': 12}
plt.rcParams.update(params)

user = os.getcwd().split('/')[2]

model='EC_Earth'
# model='ERA5'

user = os.getcwd().split('/')[2]
media = os.getcwd().split('/')[3]


if user=='pst019':
    Mediadir= '/media/'+user+'/'+media+'/data/Energy_Transport/'+ model +'/'
else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
if model=='ERA5': Mediadir += 'EnergySplit/res_0.5x0.5/Waves/'

# scale= 'sine'
scale=''

fix_y_axis= False
trans_var= False # plot the transient component - not included in many panels - is specified in the loop for the transient transport, needed here for an if statement of the legend


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
# pannellist= ['Mean' , 'Var'] #, 'Change-Mean', 'Change-Var', 'Conv', 'Change-Conv', 'Change-Fraction', 'Change-Frac-Var'
# pannellist= ['Mean' , 'Change-Mean'] #, 'Change-Var', 'Conv', 'Change-Conv'
# pannellist= ['Mean' , 'Conv'] 
# pannellist= ['Change-Mean', 'Change-Conv']
# pannellist= ['Change-Mean', 'Change-Fraction']
# pannellist= ['Change-Mean', 'Change-Fraction', 'Change-Conv']
# pannellist= ['Change-Mean']
# pannellist = ['Mean', 'Change-Mean', 'Change-Fraction', 'Change-Conv']
# pannellist= ['Change-Var-Roll']
# pannellist= ['Change-Var-Roll', 'Change-Frac-Var-Roll']
# pannellist= ['Var', 'Change-Var']


# pannellist, fix_y_axis= ['Mean'], [2E8]
# pannellist, timeperiod, fix_y_axis = ['Mean'], 'DJF', [[-1.8E8, 2.5E8]]
# pannellist, timeperiod, fix_y_axis = ['Mean'], 'JJA', [[-2.5E8, 1.8E8]]
# pannellist, fix_y_axis= ['Mean' , 'Var'], [[-2E8, 2E8], [0, 1.15E7]]  #for EC-Earth

pannellist, fix_y_axis= ['Change-Mean', 'Change-Fraction', 'Change-Conv'], [[-2E7, 2E7], [-.25, .35], [-15, 16]]
# pannellist, timeperiod, fix_y_axis= ['Change-Mean', 'Change-Fraction'], 'SON', [[-2.2E7, 4E7], [-.4, .4]]
# pannellist, fix_y_axis= ['Var', 'Var-Fraction'], [[0, 1.2E7], [0, .4]]
# pannellist, fix_y_axis=  ['Change-Var-Roll', 'Change-Frac-Var-Roll'], [[-3.5E6, 4E6], [-.4, .6]]

# pannellist, fix_y_axis= ['Var'], [[0, 1.2E7]]
# pannellist, scale, fix_y_axis= ['Conv'], 'sine', [[-80, 120]]


save= False
save= True
savedir= '../Figs/Global_4/'

imp= False
imp= True
#
#Member = 1

if model == 'EC_Earth':
    syear= 1950
    eyear= 1960 #2100
    # syear= 1979
    # eyear= 2018    
    Memberlist= [1, 2, 3, 4, 5]#,2,3,4]


elif model== 'ERA5':
    syear= 1979
    eyear= 2018
    Memberlist=[1]

split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 4





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


# elif typ == 'Wavelength_smooth':
#     if energytyp== 'E': varlist, monthlymeanlist= ['vEtot', 'vEsetot'], [False, True]
#     elif energytyp== 'Q': varlist, monthlymeanlist= ['vQtot', 'vQsetot'], [False, True]
#     elif energytyp== 'D': implist, monthlymeanlist= ['vEtot', 'vEsetot', 'vQtot', 'vQsetot'], [False, True, False, True]
    
#     catlist = ['Meri', 'Plan', 'Syno']#, 'Meso']
#     monthlymean = False

# elif typ == 'Wavelength_smooth_2':
#     if energytyp== 'E': varlist, monthlymeanlist= ['vEtot', 'vEsetot'], [False, True]
#     elif energytyp== 'Q': varlist, monthlymeanlist= ['vQtot', 'vQsetot'], [False, True]
#     elif energytyp== 'D': implist, monthlymeanlist= ['vEtot', 'vEsetot', 'vQtot', 'vQsetot'], [False, True, False, True]

#     catlist = ['Meri', 'Plan', 'Syno', 'Meso']
#     monthlymean = False


# elif typ == 'Wavelength_smooth_3':
#     varlist= ['vEtot']
# #    varlist= ['vQtot']
#     catlist = ['Meri', 'Plan', 'Syno_l', 'Syno_s', 'Meso']
    
#     monthlymean = False






if energytyp in ['E', 'Q', 'E_nostat']:  implist= varlist
elif energytyp== 'D': varlist= ['vDtot', 'vDsetot']

if pannellist == ['Conv']: varlist= varlist[:1]




if imp:
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


if latcut != 90: ds= ds.where(np.abs(ds.lat) < latcut, drop= True)


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc


if energytyp== 'D':
    ds['vDtot']= ds['vEtot'] - ds['vQtot']
    ds['vDsetot']= ds['vEsetot'] - ds['vQsetot']


if timeperiod== 'year': ds= ds.resample(time='1Y').mean()
print('only for annual timeperiod')
# .groupby('time.season')


"""make the sine scale"""
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker

class SineScale(mscale.ScaleBase):
    """
    ScaleBase class for generating square root scale.
    """
 
    name = 'sine'
 
    def __init__(self, axis, **kwargs):
        # note in older versions of matplotlib (<3.1), this worked fine.
        # mscale.ScaleBase.__init__(self)

        # In newer versions (>=3.1), you also need to pass in `axis` as an arg
        mscale.ScaleBase.__init__(self, axis)
 
    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter())
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())
 
    def limit_range_for_scale(self, vmin, vmax, minpos):
        return  max(-90., vmin), min(90, vmax)
 
    class SineTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform_non_affine(self, a): 
            # return np.array(a)**0.5
            return np.sin(np.deg2rad(a))
 
        def inverted(self):
            return SineScale.InvertedSineTransform()
 
    class InvertedSineTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform(self, a):
            return np.rad2deg(np.arcsin(a))
 
        def inverted(self):
            return SineScale.SineTransform()
 
    def get_transform(self):
        return self.SineTransform()
 
mscale.register_scale(SineScale)






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









"""Plot the pannels from the pannellist"""
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]


fig= plt.figure(fignr, figsize= (8, 2+3*len(pannellist)))
fignr+=1
plt.clf()

axs=fig.subplots(len(pannellist), 1, sharex= 'col')

if len(pannellist)== 1: axs= [axs]

"""calculation relevant for all pannels"""
if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
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
                elif cat != 'Meri': axs[ip].plot(ds.lat, fac* mean, ls= '--', color= color, label= cat +' stat')

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
                elif cat != 'Meri':  axs[ip].plot(ds.lat, fac* variability , ls= '--', color= color, label= cat +' stat')

                if len(Memberlist) > 1:        
                    var_std= ds[cat+'_'+var].std(dim='time').std(dim='Member')        
                    axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= color, alpha= 0.2   )


            if trans_var and cat != 'Meri':
                variability= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time').mean(dim='Member')
                axs[ip].plot(ds.lat, fac * variability, ls= ':', color= color, label=cat + ' trans')
                if len(Memberlist) > 1:        
                    var_std= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time').std(dim='Member')        
                    axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= color, alpha= 0.2   )


                
        # axs[ip].set_ylabel('Annual variability \n transport [W m$^{-1}$]')
        axs[ip].set_ylabel('Transport variability \n [W m$^{-1}$]')



    if ptype == 'Change-Mean':
        diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).mean(dim='Member')
        axs[ip].plot(ds.lat, fac*diffmean, label= 'Total', color= 'k')

        """significance t-test"""
        from scipy import stats
        distr_1= dstot.where(ds["time.year"]>= split_2, drop= True)
        distr_2= dstot.where(ds["time.year"]<= split_1, drop= True)
        d1= distr_1.stack(x= ("time", "Member"))
        d2= distr_2.stack(x= ("time", "Member"))
        p_value= stats.ttest_ind(d1, d2, axis= 1)[1]
        axs[ip].plot(ds.lat, fac*diffmean.where(p_value >.05, np.nan), label= 'Total', color= 'k')
        
        # axs[ip].plot(ds.lat, fac*diffmean, label= 'Total', color= 'k')

        if len(Memberlist) > 1:
            diffstd= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') 
                  - dstot.where(ds["time.year"]<= split_1).mean(dim='time') ).std(dim='Member')               
            axs[ip].fill_between(ds.lat, fac*(diffmean - diffstd),  fac*(diffmean + diffstd), color= 'k', alpha= 0.2   )

        
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                diffmean= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).mean(dim='Member')
                if vi == 0:  axs[ip].plot(ds.lat, fac*diffmean, label= cat, color= color)
                elif cat != 'Meri':  axs[ip].plot(ds.lat, fac*diffmean, ls= '--', label= cat +' stat', color= color)

                if len(Memberlist) > 1:        
                    diffstd= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time')
                              - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')).std(dim='Member')
                    axs[ip].fill_between(ds.lat, fac*(diffmean - diffstd),  fac*(diffmean + diffstd), color= color, alpha= 0.2   )


        axs[ip].set_ylabel('Transport change'#' \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1) 
                           +'\n [W m$^{-1}$]')


    if ptype == 'Change-Var':
        diff_variability= (dstot.where(ds["time.year"]>= split_2).std(dim='time').mean(dim='Member')
                                    - dstot.where(ds["time.year"]<= split_1).std(dim='time').mean(dim='Member') )

        axs[ip].plot(ds.lat, fac * diff_variability, label= 'Total', color= 'k')

        if len(Memberlist) > 1:
            diff_variability_std= (dstot.where(ds["time.year"]>= split_2).std(dim='time')
               - dstot.where(ds["time.year"]<= split_1).std(dim='time') ).std(dim='Member')
    
            axs[ip].fill_between(ds.lat, fac*(diff_variability - diff_variability_std),  fac*(diff_variability + diff_variability_std), color= 'k', alpha= 0.2   )

    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                #change in standard deviation
                diff_variability= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).std(dim='time')
                                      - ds[cat+'_'+var].where(ds["time.year"]<= split_1).std(dim='time')).mean(dim='Member')
                if vi == 0: axs[ip].plot(ds.lat, fac* diff_variability , color= color, label= cat )
                elif cat != 'Meri': axs[ip].plot(ds.lat, fac* diff_variability , ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1:                 
                    diff_variability_std= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).std(dim='time')
                        - ds[cat+'_'+var].where(ds["time.year"]<= split_1).std(dim='time') ).std(dim='Member')
                    axs[ip].fill_between(ds.lat, fac*(diff_variability - diff_variability_std),  fac*(diff_variability + diff_variability_std), color= color, alpha= 0.2   )


        axs[ip].set_ylabel('Variability change\n' #' \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1)
                           +' [W m$^{-1}$]')


    if ptype == 'Change-Var-Roll': #the variance from the rolling mean

        anomalie= dstot- dstot.rolling(time= 30, center= True).mean()
        diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                        - anomalie.where(ds["time.year"]<= split_1).std(dim='time')).mean(dim='Member')

        axs[ip].plot(ds.lat, fac * diff_variability, label= 'Total', color= 'k')

        if len(Memberlist) > 1:         
            diff_variability_std= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                - anomalie.where(ds["time.year"]<= split_1).std(dim='time')).std(dim='Member')
            
            axs[ip].fill_between(ds.lat, fac*(diff_variability - diff_variability_std),  fac*(diff_variability + diff_variability_std), color= 'k', alpha= 0.2   )

    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                #change in standard deviation
                anomalie= ds[cat+'_'+var]- ds[cat+'_'+var].rolling(time= 30, center= True).mean()

                diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                                  - anomalie.where(ds["time.year"]<= split_1).std(dim='time')).mean(dim='Member')
                
                if vi == 0: axs[ip].plot(ds.lat, fac* diff_variability, color= color, label= cat)
                elif cat != 'Meri': axs[ip].plot(ds.lat, fac* diff_variability , ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1:                 
                    diff_variability_std= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                        - anomalie.where(ds["time.year"]<= split_1).std(dim='time') ).std(dim='Member')
                    axs[ip].fill_between(ds.lat, fac*(diff_variability - diff_variability_std),  fac*(diff_variability + diff_variability_std), color= color, alpha= 0.2   )


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
                elif cat != 'Meri': axs[ip].plot(ds.lat, conv_mean, ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1: 
                    conv_std= conv.std(dim='Member')                   
                    axs[ip].fill_between(ds.lat, (conv_mean - conv_std),  conv_mean + conv_std,
                                    color= color, alpha= 0.2   )
                
        axs[ip].set_ylabel('Convergence [W m$^{-2}$]')


    if ptype == 'Change-Conv':
        difference= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') 
           - dstot.where(ds["time.year"]<= split_1).mean(dim='time')) #.mean(dim='Member')
        difference *= np.cos(np.deg2rad(difference.lat)) #.mean(dim='Member')    
        conv_diff= -1* difference.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(difference.lat))
        conv_diff= conv_diff.rolling(lat= 5, center= True).mean() #smooth
        conv_diff_mean= conv_diff.mean(dim='Member')
    
        axs[ip].plot(ds.lat, conv_diff_mean, label= 'Total', color= 'k')

        if len(Memberlist) > 1: 
            conv_diff_std= conv_diff.std(dim='Member')
            axs[ip].fill_between(ds.lat, (conv_diff_mean - conv_diff_std),  conv_diff_mean + conv_diff_std,
                            color= 'k', alpha= 0.2)
    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist[:1]):    
                difference= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') 
                     - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time')) #.mean(dim='Member')
                difference *= np.cos(np.deg2rad(difference.lat)) #.mean(dim='Member')    

                conv_diff= -1* difference.differentiate('lat')/110E3 * 1/np.cos(np.deg2rad(difference.lat))
                conv_diff= conv_diff.rolling(lat= 5, center= True).mean()
                conv_diff_mean= conv_diff.mean(dim='Member')


                if vi == 0: axs[ip].plot(ds.lat, conv_diff_mean, label= cat, color= color)
                # elif cat != 'Meri': axs[ip].plot(ds.lat, conv_diff_mean, ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1: 
                    conv_diff_std= conv_diff.std(dim='Member')
                    axs[ip].fill_between(ds.lat, (conv_diff_mean - conv_diff_std),  conv_diff_mean + conv_diff_std, color= color, alpha= 0.2)

        axs[ip].set_ylabel('Convergence change\n [W m$^{-2}$]')


    if ptype== 'Var-Fraction': # (variance) / (absolute in the mean)
        mean= np.abs(dstot.mean(dim='time'))
        # filter_level= 0.2* mean.max()
        # meanfilter= mean.mean(dim='Member') > filter_level

        filter_level= 0.2* np.abs(mean.mean(dim='Member').max())
        meanfilter= np.abs(mean.mean(dim='Member')) > filter_level

        print('make same filter other places!')

        variability= dstot.std(dim='time')

        axs[ip].plot(ds.lat, fac* (variability/mean).mean(dim='Member').where(meanfilter, np.nan), color= 'k', label='Total')

        # if len(Memberlist) > 1:
        #     axs[ip].fill_between(ds.lat, ((variability/mean).mean(dim='Member') - (variability/mean).std(dim='Member')).where(meanfilter, np.nan),
        #                     ((variability/mean).mean(dim='Member') + (variability/mean).std(dim='Member')).where(meanfilter, np.nan), color= 'k', alpha= 0.2   )

        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                mean= np.abs(ds[cat+'_'+var].mean(dim='time'))
                meanfilter= mean.mean(dim='Member') > filter_level
                
                
                variability= ds[cat+'_'+var].std(dim='time') #gets the standard deviation

                if vi == 0:  axs[ip].plot(ds.lat, fac* (variability/mean).mean(dim='Member').where(meanfilter, np.nan), label= cat, color= color)
                elif cat not in  ['Meri', 'Syno']: axs[ip].plot(ds.lat, fac*(variability/mean).mean(dim='Member').where(meanfilter, np.nan), label= cat+' stat', ls= '--', color= color)

                # if len(Memberlist) > 1 and cat not in  ['Meri', 'Syno']:                 
                #     axs[ip].fill_between(ds.lat, ((diffmean/mean).mean(dim='Member') - (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan),
                #             ((diffmean/mean).mean(dim='Member') + (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan), color= color, alpha= 0.2   )

                if trans_var and cat != 'Meri': #trans_var - compute the transient component
                    mean= np.abs((ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).mean(dim='time'))
                    meanfilter= mean.mean(dim='Member') > filter_level

                    variability= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time')
                    axs[ip].plot(ds.lat, fac * (variability/mean).mean(dim='Member').where(meanfilter, np.nan), ls= ':', color= color, label=cat + ' trans')
                    # if len(Memberlist) > 1:        
                    #     var_std= (ds[cat+'_'+varlist[0]]-ds[cat+'_'+varlist[1]]).std(dim='time').std(dim='Member')        
                    #     axs[ip].fill_between(ds.lat, fac*(variability - var_std),  fac*(variability + var_std), color= color, alpha= 0.2   )

        axs[ip].set_ylabel('Variability fraction') #' \n'+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))




    if ptype== 'Change-Fraction': # (transport change) / (mean transport)
        mean= dstot.mean(dim='time')
        filter_level= 0.2* np.abs(mean.mean(dim='Member').max())

        meanfilter= np.abs(mean.mean(dim='Member')) > filter_level
        diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time')
               - dstot.where(ds["time.year"]<= split_1).mean(dim='time') )

        axs[ip].plot(ds.lat, (diffmean/mean).mean(dim='Member').where(meanfilter, np.nan), color= 'k', label='Total')

        if len(Memberlist) > 1:
            axs[ip].fill_between(ds.lat, ((diffmean/mean).mean(dim='Member') - (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan),
                            ((diffmean/mean).mean(dim='Member') + (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan), color= 'k', alpha= 0.2   )


        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                mean= ds[cat+'_'+var].mean(dim='time')
                meanfilter= np.abs(mean.mean(dim='Member')) > filter_level
                
                diffmean= (ds[cat+'_'+var].where(ds["time.year"]>= split_2).mean(dim='time') 
                           - ds[cat+'_'+var].where(ds["time.year"]<= split_1).mean(dim='time'))

                if vi == 0:  axs[ip].plot(ds.lat, (diffmean/mean).mean(dim='Member').where(meanfilter, np.nan), label= cat, color= color)
                # elif cat not in  ['Meri', 'Syno']:
                else: axs[ip].plot(ds.lat, (diffmean/mean).mean(dim='Member').where(meanfilter, np.nan), label= cat+' stat', ls= '--', color= color)

                if len(Memberlist) > 1: # and cat not in  ['Meri', 'Syno']:                 
                    axs[ip].fill_between(ds.lat, ((diffmean/mean).mean(dim='Member') - (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan),
                            ((diffmean/mean).mean(dim='Member') + (diffmean/mean).std(dim='Member')).where(meanfilter, np.nan), color= color, alpha= 0.2   )

        axs[ip].set_ylabel('Fraction change') #' \n'+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))


    if ptype== 'Change-Frac-Var': #(change transport variance) / (mean variance)
        print('see _6 around line 900')

    



    if ptype== 'Change-Frac-Var-Roll': #(change transport variance) / (mean variance)
        anomalie= dstot- dstot.rolling(time= 30, center= True).mean()
        variability= anomalie.std(dim='time')
        filter_level= 0.01* np.abs(variability.mean(dim='Member').max())
        varfilter= np.abs(variability.mean(dim='Member')) > filter_level
        
        diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                        - anomalie.where(ds["time.year"]<= split_1).std(dim='time')).where(varfilter, np.nan)
        variability_frac= (diff_variability/variability).mean(dim='Member')



        axs[ip].plot(ds.lat, fac * variability_frac, label= 'Total', color= 'k')

        if len(Memberlist) > 1:                   
            variability_frac_std= (diff_variability/variability).std(dim='Member')
            
            axs[ip].fill_between(ds.lat, fac*(variability_frac - variability_frac_std),
                                 fac*(variability_frac + variability_frac_std), color= 'k', alpha= 0.2   )

    
        for ci, cat in enumerate(catlist):
            color = colors[ci]
            for vi, var in enumerate(varlist):    
                anomalie= ds[cat+'_'+var]- ds[cat+'_'+var].rolling(time= 30, center= True).mean()
                variability= anomalie.std(dim='time')
                varfilter= np.abs(variability.mean(dim='Member')) > filter_level

                diff_variability= (anomalie.where(ds["time.year"]>= split_2).std(dim='time')
                        - anomalie.where(ds["time.year"]<= split_1).std(dim='time')).where(varfilter, np.nan)
                variability_frac= (diff_variability/variability).mean(dim='Member')

                if vi == 0: axs[ip].plot(ds.lat, fac* variability_frac, color= color, label= cat )
                elif cat != 'Meri': axs[ip].plot(ds.lat, fac* variability_frac , ls= '--', color= color, label= cat +' stat' )

                if len(Memberlist) > 1:                 
                    variability_frac_std= (diff_variability/variability).std(dim='Member')

                    axs[ip].fill_between(ds.lat, fac*(variability_frac - variability_frac_std),
                                 fac*(variability_frac + variability_frac_std), color= color, alpha= 0.2   )


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


if pannellist == ['Conv'] or energytyp== 'E_nostat':  axs[0].legend(ncol= 2)
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



