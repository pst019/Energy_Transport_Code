#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/'

else:
    Mediadir= '/run/media/pst019/Backup1/'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib as mpl
from functions import *
import pandas as pd
import scipy.signal


attrs = {"units": "days since 2000-01-01"}


save= False
save= True

Member = 3
syear_1= 1950
eyear_1= 1979

syear_2= 2070
eyear_2= 2100


fignr= 10



varlist=[ 'T_850', 'vEtot', 'vQtot'] 
fvar= 'T' #the display variable
label='Temperature'
unit='K'
evarlist=['vEtot', 'vQtot', 'vEtot', 'vEtot', 'vEtot']
transp_type_list= ['Total', 'Total', 'Plan', 'Syno', 'Meri']
namelist= [r'E$_{tot}$', r'Q$_{tot}$', r'E$_{plan}$', r'E$_{syno}$', r'E$_{meri}$']
evar_name='Energy'

# varlist=[ 'Q_850', 'vQtot'] 
# fvar= 'Q' #the display variable
# label='Specific humidity'
# unit='g/kg'
# evarlist=['vQtot', 'vQtot', 'vQtot', 'vQtot']
# transp_type_list= ['Total', 'Plan', 'Syno', 'Meri']
# namelist= [r'Q$_{tot}$', r'Q$_{plan}$', r'Q$_{syno}$', r'Q$_{meri}$']
# evar_name='Moisture'


latlim= 80


filter_value= 0.15 #filters the fraction tramsport per temperature gradient where the temp grad is less than filter_value of its meridional maximum

diff_fac= 100E3 #to compute from W/m / K/100km to W/K


file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'



""" import data"""
for vi, var in enumerate(varlist):
    print(var)
    # varfile= varfilelist[vi]
    intermediate_file_name= Mediadir + 'data/Energy_Transport/EC_Earth/intermediate_data/'+var+'_'
    intermediate_file_name += str(syear_1)+'-'+str(eyear_1) + '_' + str(syear_2)+'-'+str(eyear_2)
    intermediate_file_name += '_M'+str(Member)
    intermediate_file_name +=  '.nc'

    if os.path.exists(intermediate_file_name):
        ds0= xr.open_dataset(intermediate_file_name)
        
    else:
        for year in list(np.arange(syear_1, eyear_1+1)) + list(np.arange(syear_2, eyear_2+1)):
            # if year%10 == 0:
            print(year)
            for month in range(1,13):
                print(month)
        
                if year == syear_1 and month== 1:
                    if var in ['vEtot', 'vQtot']:
                        ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')
                    else:
                        ds0= xr.open_dataset(file_dir + var +'/' + var+'_'+str(year)+'_'+ str(month).zfill(2)+ '.nc').mean(dim=['lon', 'time'])
                        
                    ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    
        
                else:
                    if var in ['vEtot', 'vQtot']:
                        ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')
                    else:
                        ds1= xr.open_dataset(file_dir + var +'/' + var+'_'+str(year)+'_'+ str(month).zfill(2)+ '.nc').mean(dim=['lon', 'time'])
        
                    ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    
                    ds0= xr.concat([ds0, ds1], dim= 'time')
        

        if var in ['Z_850', 'Q_850', 'T_850']: ds0= ds0.isel(plev= 0)
        ds0.to_netcdf(intermediate_file_name)


    
    """merge the different variables in one ds"""

    if var== varlist[0]:
        ds2= ds0
    else:
        if np.max(np.abs(ds2.lat.values -ds0.lat.values) ) < .1:
            ds0['lat']= ds2.lat
        else:
            print('lat does not fit')
        ds2= xr.merge([ds2, ds0])
    
    
ds= ds2

latcut= 89
ds= ds.where(np.abs(ds.lat) < latcut, drop= True) #Q has a problem at lat= 89.14


if fvar == 'Q': ds[fvar] *= 1E3


  

ds= ds.reindex(lat= list(reversed(ds.lat)))

a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))
Split1= 8000E3
WaveSplit1= LatCirc/Split1






"""transport difference, gradient +fraction"""


fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col', gridspec_kw={'height_ratios': [1, 1, 1.5, 1.5]})



dlat= ds.lat[1].values- ds.lat[0].values

mean_last= ds[fvar].where(ds["time.year"]>= syear_2).mean(dim='time') 
# mean_last_smth= mean_last.rolling(lat= 5, center= True).mean()
mean_last_grad_sav= scipy.signal.savgol_filter(mean_last, window_length=11, polyorder=1, deriv=1)/dlat *1E2/110 #K per 100km
filter_level= filter_value* np.abs(mean_last_grad_sav).max()
absfilter_sav_last= np.abs(mean_last_grad_sav) < filter_level  

mean_first= ds[fvar].where(ds["time.year"]<= eyear_1).mean(dim='time') 
# mean_first_smth= mean_first.rolling(lat= 5, center= True).mean()
mean_first_grad_sav= scipy.signal.savgol_filter(mean_first, window_length=11, polyorder=1, deriv=1)/dlat *1E2/110 #K per 100km
filter_level= filter_value* np.abs(mean_first_grad_sav).max()
absfilter_sav_first= np.abs(mean_first_grad_sav) < filter_level        

# axs[1].plot(ds.lat, mean )
# axs[1].plot(ds.lat, mean_smth, label=f'{syear}-{eyear}', color= line.get_color() )
line, =axs[0].plot(ds.lat, mean_first, 'k', label=f'{syear_1}-{eyear_1}')
axs[0].plot(ds.lat, mean_last, color= 'k', ls= '--', label=f'{syear_2}-{eyear_2}')

axs[1].plot(ds.lat, mean_first_grad_sav, 'k')
axs[1].plot(ds.lat, mean_last_grad_sav , color= 'k' , ls='--')


for evar, transport_type, varname in zip(evarlist, transp_type_list, namelist):
    if transport_type== 'Total':
        transp = ds[evar].sum(dim= "WaveNumb")
    elif transport_type== 'Meri':
            transp= ds[evar].sel(WaveNumb= 0)
    elif transport_type== 'Plan':        
            transp= (ds[evar].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
                    +  WaveSplit1%1 * ds[evar].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    elif transport_type== 'Syno':        
            transp= ( ds[evar].where(ds.WaveNumb > 1+ WaveSplit1//1).sum(dim='WaveNumb') 
                         +(1- WaveSplit1%1) * ds[evar].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
        
    
    trans_mean_last= transp.where(ds["time.year"]>= syear_2).mean(dim='time')
    trans_mean_first= transp.where(ds["time.year"]<= eyear_1).mean(dim='time')
    
    line, =axs[2].plot(ds.lat, trans_mean_first, label=f'{varname}') # {syear_1}+' )
    axs[2].plot(ds.lat, trans_mean_last, color= line.get_color() , ls='--' ) #, label=f'{varname} {syear_2}+'
    
    # axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_sav_first!= 0, -1* trans_mean_first/mean_first_grad_sav), color= line.get_color())
    # axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_sav_last!= 0, -1* trans_mean_last/mean_last_grad_sav), color= line.get_color(), ls= '--')
    
    if fvar== 'Q': latfilter= 10
    else: latfilter= 23
    
    tropical_filter= np.abs(ds.lat) < latfilter
    axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(tropical_filter!= 0, -1* trans_mean_first/mean_first_grad_sav), color= line.get_color())
    axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(tropical_filter!= 0, -1* trans_mean_last/mean_last_grad_sav), color= line.get_color(), ls= '--')



for axnr in [1, 2,3]:
    axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs[axnr].ticklabel_format(useMathText=True) #to make 10**n on the axis

# axs[1].plot([-latcut, latcut], [filter_level,filter_level], c= 'grey', linewidth= 0.5, linestyle= '--')
# axs[1].plot([-latcut, latcut], [-filter_level, -filter_level], c= 'grey', linewidth= 0.5, linestyle= '--')

axs[3].set_xticks(np.arange(-80,80.1, 10), minor=True)
axs[3].set_xlim(-latlim, latlim)
axs[3].set_xlabel('Latitude')
axs[2].set_ylabel(f'{evar_name} transport \n per lenth unit [W/m]') #energy transport accross a latitude band
axs[0].set_ylabel(f'{label} \n [{unit}]')
axs[1].set_ylabel(f'{label} gradient \n [{unit} (100 km)'+r'$^{-1}$]')
axs[3].set_ylabel(f'Transport per \n{label} gradient  \n[W/{unit}]') #energy transport accross a latitude band
axs[3].set_ylabel(f'Diffusion constant  \n[W/{unit}]') #energy transport accross a latitude band


axs[0].legend(ncol= 2)
if len(namelist) == 4:axs[2].legend(ncol= 4)
else: axs[2].legend(ncol= 3)


if save:
    savedir= Mediadir +'/home/Energy_Transport/Figs/Transport-vs-gradient_2/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    # file_name= fvar +'_' + evar 
    file_name= savedir + varlist[0] + '_'
    for evar, transp_type in zip(evarlist, transp_type_list): file_name += f'{evar}-{transp_type}_'
    file_name += f'M{Member}'
    file_name += '_' +str(syear_1)+'-'+str(eyear_1)+'_'+str(syear_2)+'-'+str(eyear_2)
    print(file_name)
    fig.savefig(file_name , bbox_inches='tight') 

