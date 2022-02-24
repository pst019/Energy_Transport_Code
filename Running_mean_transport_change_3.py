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
import scipy.stats as stats

user = os.getcwd().split('/')[2]



model='EC_Earth'
model='ERA5'


if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/'+ model +'/'
elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/'+ model +'/'    
else:
    if model=='EC_Earth': Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI/'
    if model=='ERA5': Mediadir= '/nird/projects/nird/NS9063K/from_stallo/era5/EnergySplit/res_0.5x0.5/Waves/'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    




save= True
# save= False
imp= False
imp= True
#
#Member = 1
syear= 1950
eyear= 2100
split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 6


Memberlist= [1,2,3,4,5]

timeperiod= 'year'
# timeperiod= 'DJF'
# timeperiod= 'JJA'

compare_to = 2000

if timeperiod == 'year': monthlist= np.arange(1,13)
elif timeperiod == 'DJF': monthlist= [1,2,12]
elif timeperiod == 'JJA': monthlist= [6,7,8]


#typ = 'Eddies'
typ = 'Wavelength_smooth' # corrected
# typ = 'Wavelength_smooth_2' #corrected
# typ = 'Wavelength_smooth_3' #not corrected

#typ = 'WaveNR'

across_lat= 70

unit='PerM'
# unit='LatCycl'

if typ == 'Eddies':
    varlist= ['vEmc', 'vEse', 'vEte']
#    varlist= ['vQmc', 'vQse', 'vQte']
    catlist= varlist
    monthlymean = True

elif typ == 'Wavelength_smooth':
    varlist= ['vEtot']#, 'vEsetot']
    # varlist= ['vQtot']#, 'vQsetot']
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
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').sel(lat= across_lat, method='nearest').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').sel(lat= across_lat, method='nearest')
                            ds0= ds0.mean(dim = 'time') #calculate the monthly mean
                            ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))


                        if Member == Memberlist[0]: dslat= ds0.lat
                        else: ds0['lat']= dslat
                        
                        
                    else:
                        if monthlymean:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').sel(lat= across_lat, method='nearest').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').sel(lat= across_lat, method='nearest')
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




a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))


if unit =='PerM': fac= 1
elif unit == 'LatCycl': fac= LatCirc



"""remove years with an outlier months"""
# zscore_thresh= 2.5
value_thresh= 1E6

fignr= 6
fig_raw= plt.figure(fignr, figsize= (10, 5))
fignr+=1
plt.clf()

ax_raw=fig_raw.subplots()


for vi, var in enumerate(varlist):
    if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
    if typ in ['Wavelength_smooth','Wavelength_smooth_2','Wavelength_smooth_3', 'WaveNR']:
        ds['tot'] = ds[var].sum(dim= "WaveNumb")

    for mi, member in enumerate(Memberlist):
        exclude_year_list= []

        if mi == 0: line, = ax_raw.plot(ds.time, fac * ds['tot'].sel(Member= member), label= 'Total') #, ls= '--', lw= 1)
        else: line, = ax_raw.plot(ds.time, fac * ds['tot'].sel(Member= member)) #, ls= '--', lw= 1)
        
        for month in range(1,13):
            ds_month= ds['tot'].sel(Member= member).sel(time = (ds.time.dt.month== month) )
            
            # z_scores = np.abs(stats.zscore(ds_month ) )    # calculate z-scores
            # excludes =(z_scores > zscore_thresh) #.all(axis=1)
            
            # z_scores = stats.zscore(ds_month )    # calculate z-scores
            # excludes =(z_scores < -1* zscore_thresh) #.all(axis=1)         
            
            excludes= ds_month.values < value_thresh
            
            if len(np.where(excludes)[0]) > 0: print(var, 'Member', member, list(pd.to_datetime(ds_month.time[excludes].values).strftime('%Y.%m') ) )   

            ax_raw.plot(ds_month.time[excludes], fac * ds_month.values[excludes], 'o', color= line.get_color() , label=member)
    
            if len(np.where(excludes)[0]) > 0:
                exclude_year_list += list(ds_month.time[excludes].dt.year.values)
                
        for exclude_year in exclude_year_list:
            # print('exclude year', exclude_year)
            ds[var].loc[dict(Member= member, time =(ds.time.dt.year == exclude_year)  )] = np.nan
            ds['tot'].loc[dict(Member= member, time =(ds.time.dt.year == exclude_year)  )] = np.nan


plt.legend()

fig_excl= plt.figure(fignr, figsize= (10, 5))
fignr+=1
plt.clf()

ax_excl=fig_excl.subplots()

exclude_year_list= []

for vi, var in enumerate(varlist):
    # if typ == 'Eddies': dstot= sum(ds[var] for var in varlist)
    # if typ in ['Wavelength_smooth','Wavelength_smooth_2','Wavelength_smooth_3', 'WaveNR']:
    #     ds['tot'] = ds[var].sum(dim= "WaveNumb")

    for mi, member in enumerate(Memberlist):
        if mi == 0: line, = ax_excl.plot(ds.time, fac * ds['tot'].sel(Member= member), label= 'Total') #, ls= '--', lw= 1)
        else: line, = ax_excl.plot(ds.time, fac * ds['tot'].sel(Member= member)) #, ls= '--', lw= 1)
        

"""make annual series"""


# if timeperiod== 'year': ds= ds.resample(time='1Y').mean()
# print('only for annual timeperiod')
# .groupby('time.season')

if timeperiod== 'year': ds= ds.resample(time='1Y').mean()

# .sel(WaveNumb= 1, Member= 1)['vQtot']

"""Calculate transport of each category"""
if typ == 'Wavelength_smooth':
    Split1= 8000E3
    WaveSplit1= LatCirc/Split1
    
    for var in varlist:

        ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
        ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1) )
        ds['Syno_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit1//1).sum(dim='WaveNumb') 
                      +(1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1) )
    

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



"""Plot the Total transport"""
fig0= plt.figure(fignr, figsize= (10, 5))
fignr+=1
plt.clf()

ax0=fig0.subplots()


fig= plt.figure(fignr, figsize= (10, 8))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')


fig2= plt.figure(fignr, figsize= (10, 8))
fignr+=1
plt.clf()

axs2=fig2.subplots(2, 1, sharex= 'col')


for vi, var in enumerate(varlist):
    for mi, member in enumerate(Memberlist):
        if mi == 0: line, = ax0.plot(ds.time, fac * ds['tot'].sel(Member= member), label= 'Total', color= 'k') #, ls= '--', lw= 1)
        else: line, = ax0.plot(ds.time, fac * ds['tot'].sel(Member= member), color= 'k') #, ls= '--', lw= 1)
      
    mean= ds['tot'].rolling(time= 30, min_periods= 30-3, center= True).mean().mean(dim='Member')
    std= ds['tot'].rolling(time= 30, min_periods= 30-3, center= True).std().mean(dim='Member')


    # if compare_to:
    compare_mean = mean.sel(time= f'{compare_to}-12-31')
    compare_std = std.sel(time= f'{compare_to}-12-31')

        # line, = ax0.plot(ds.time, fac * dstot, label= 'Total', color= 'k') #, ls= '--', lw= 1)
    axs[0].plot(ds.time, fac * (mean - compare_mean), label= 'Total', color= 'k')
    axs[1].plot(ds.time, fac* (std - compare_std) , color= 'k')

    axs2[0].plot(ds.time, mean /compare_mean, label= 'Total', color= 'k')
    axs2[1].plot(ds.time, std / compare_std , color= 'k')


    if len(Memberlist) > 1:
        mean_var= ds['tot'].rolling(time= 30, min_periods= 30-7, center= True).mean().std(dim='Member')
        std_var= ds['tot'].rolling(time= 30, min_periods= 30-7, center= True).std().std(dim='Member')

        axs[0].fill_between(ds.time, fac*(mean - compare_mean - mean_var), fac*(mean - compare_mean + mean_var), color= line.get_color(), alpha= 0.3)
        axs[1].fill_between(ds.time, fac*(std - compare_std- std_var), fac*(std - compare_std + std_var), color= line.get_color(), alpha= 0.3)

        # mean_indiv = ds['tot'].rolling(time= 30, min_periods= 30-3, center= True).mean()
        # mean_frac_var= (mean_indiv/mean_indiv.sel(time= f'{compare_to}-12-31')).std(dim='Member')

        # std_indiv = ds['tot'].rolling(time= 30, min_periods= 30-3, center= True).std()
        # std_frac_var= (std_indiv/std_indiv.sel(time= f'{compare_to}-12-31')).std(dim='Member')

        # axs2[0].fill_between(ds.time, mean / compare_mean - mean_frac_var, mean / compare_mean + mean_frac_var, color= line.get_color(), alpha= 0.3)
        # axs2[1].fill_between(ds.time, std / compare_std - std_frac_var, std / compare_std + std_frac_var, color= line.get_color(), alpha= 0.3)

          
     
""" Plot the components"""           
for ci, cat in enumerate(catlist):
    
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    for vi, var in enumerate(varlist):    
        
        ds['var']= ds[cat+'_'+var]
        
        for mi, member in enumerate(Memberlist):
            if mi == 0: ax0.plot(ds.time, fac * ds['var'].sel(Member= member), label= cat, color= color) #, ls= '--', lw= 1)
            else: ax0.plot(ds.time, fac * ds['var'].sel(Member= member), color= color) #, ls= '--', lw= 1)
    
        #     z_scores = np.abs(stats.zscore(ds['var'].sel(Member= member) ) )    # calculate z-scores
        #     excludes =(z_scores > zscore_thresh)
        #     if len(np.where(excludes)[0]) > 0: print(cat, var, 'Member', member, ds.time[excludes].values )

        #     ax0.plot(ds.time[excludes], fac * ds['var'].sel(Member=member)[excludes], 'o', color= color)
        
        # excludes =(np.abs(stats.zscore(ds['var'] )) > zscore_thresh) #.all(axis=1)
        # if len(np.where(excludes)[0]) > 0: ds['var']= ds['var'].where(~excludes)
        
        mean= ds['var'].rolling(time= 30, min_periods= 30-3, center= True).mean().mean(dim='Member')
        std= ds['var'].rolling(time= 30, min_periods= 30-3, center= True).std().mean(dim='Member')

        # if compare_to:
        compare_mean = mean.sel(time= f'{compare_to}-12-31')
        compare_std = std.sel(time= f'{compare_to}-12-31')
        
        axs[0].plot(ds.time, fac* (mean - compare_mean), color= color, label= cat)
        axs[1].plot(ds.time, fac* (std - compare_std), label= var, color= color)

        axs2[0].plot(ds.time, mean / compare_mean, color= color, label= cat)
        axs2[1].plot(ds.time, std / compare_std, label= var, color= color)
        
        if len(Memberlist) > 1:
            mean_var= ds['var'].rolling(time= 30, min_periods= 30-7, center= True).mean().std(dim='Member')
            std_var= ds['var'].rolling(time= 30, min_periods= 30-7, center= True).std().std(dim='Member')
    
            axs[0].fill_between(ds.time, fac*(mean - compare_mean - mean_var), fac*(mean - compare_mean + mean_var), color= color, alpha= 0.3)
            axs[1].fill_between(ds.time, fac*(std - compare_std - std_var), fac*(std  - compare_std + std_var), color= color, alpha= 0.3)

            # mean_indiv = ds['var'].rolling(time= 30, min_periods= 30-3, center= True).mean()
            # mean_frac_var= (mean_indiv/mean_indiv.sel(time= f'{compare_to}-12-31')).std(dim='Member')
    
            # std_indiv = ds['var'].rolling(time= 30, min_periods= 30-3, center= True).std()
            # std_frac_var= (std_indiv/std_indiv.sel(time= f'{compare_to}-12-31')).std(dim='Member')
    
            # axs2[0].fill_between(ds.time, mean / compare_mean - mean_frac_var, mean / compare_mean + mean_frac_var, color= color, alpha= 0.3)
            # axs2[1].fill_between(ds.time, std / compare_std - std_frac_var, std / compare_std + std_frac_var, color= color, alpha= 0.3)
            


ax0.legend(ncol= 2) 
axs[0].legend(ncol= 2)
axs2[0].legend(ncol= 2)

for axnr in range(2):
    axs[axnr].plot([ds.time[0].values, ds.time[-1].values ], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs2[axnr].plot([ds.time[0].values, ds.time[-1].values ], [1,1], c= 'k', linewidth= 0.5, linestyle= '--')


axs[1].set_xlim(ds.time[0], ds.time[-1])
axs[1].set_xlabel('Time')

axs2[1].set_xlim(ds.time[0], ds.time[-1])
axs2[1].set_xlabel('Time')

if unit =='PerM': axs[0].set_ylabel('Mean energy transport \n per lenth unit [W/m]') #energy transport accross a latitude band
elif unit == 'LatCycl':  axs[0].set_ylabel('Mean energy transport \n across latitude circle [W]')#energy transport accross a latitude band
axs[1].set_ylabel('Annual variability \n in the energy transport')

axs2[0].set_ylabel('Fractional change in mean \n energy transport')
axs2[1].set_ylabel('Fractional change in variability \n in the energy transport')

ax0.set_ylabel('Annual energy transport [W/m]')


if save:
    savedir= '../Figs/Running_mean_change_2/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = model +'_'
    file_name += typ+'_'+unit 
    if model == 'EC_Earth': 
        for Member in Memberlist: file_name += '_M'+str(Member)
        
    file_name += '_' +str(syear)+'-'+str(eyear)+'_'+ timeperiod+'_'
    for var in varlist: file_name += var+'+'
    savefile= savedir+ file_name
    print('save: \n', savefile)
    fig.savefig(savefile , bbox_inches='tight') 

    savefile= savedir+ file_name + '_rel_change'
    print('save: \n', savefile)
    fig2.savefig(savefile , bbox_inches='tight') 

    savefile= savedir+ file_name +'_outliers'
    print('save: \n', savefile)
    fig0.savefig(savefile , bbox_inches='tight') 