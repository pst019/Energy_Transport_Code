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

attrs = {"units": "days since 2000-01-01"}


save= False
save= True

Member = 3
syear_1= 1950
eyear_1= 1979

syear_2= 2070
eyear_2= 2100


fignr= 7


varlist=[ 'T2M', 'vEtot'] #important that the u component is first
fvar= 'T2M' #the display variable
label='temperature'
unit='K'
evar='vEtot'
evar_name='energy'


# varlist=[ 'T_850', 'vEtot'] #important that the u component is first
# fvar= 'T_850' #the display variable
# label='temperature 850hPa'
# unit='K'
# evar='vEtot'
# evar_name='energy'

# varlist=['TCWV','vQtot'] #important that the u component is first
# fvar= 'TCWV' #the display variable
# label='water vapor'
# unit=r'kg/m$^2$' #not sure
# evar='vQtot'
# evar_name='latent energy'

# varlist=['Q850','vQtot'] #important that the u component is first
# fvar= 'Q' #the display variable
# label='Specific humidity'
# unit='g/kg'
# evar='vQtot'
# evar_name='latent energy'

latcut= 89

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

# latgrid= np.arange(-85, 85, 1)


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
                        ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')#.interp(lat= latgrid)
                    else:
                        ds0= xr.open_dataset(file_dir + var +'/' + var+'_'+str(year)+'_'+ str(month).zfill(2)+ '.nc').mean(dim=['lon', 'time'])#.interp(lat= latgrid)
        
                    ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    
        
                else:
                    if var in ['vEtot', 'vQtot']:
                        ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')#.interp(lat= latgrid)
                    else:
                        ds1= xr.open_dataset(file_dir + var +'/' + var+'_'+str(year)+'_'+ str(month).zfill(2)+ '.nc').mean(dim=['lon', 'time'])#.interp(lat= latgrid)
        
                    ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                    
                    ds0= xr.concat([ds0, ds1], dim= 'time')
        
        ds0.to_netcdf(intermediate_file_name)


    
    """merge the different variables in one ds"""
    if var in ['Z850', 'Q850']: ds0= ds0.isel(plev= 0)

    if var== varlist[0]:
        ds2= ds0
    else:
        if np.max(np.abs(ds2.lat.values -ds0.lat.values) ) < .1:
            ds0['lat']= ds2.lat
        else:
            print('lat does not fit')
        ds2= xr.merge([ds2, ds0])
    
    
ds= ds2

if fvar == 'Q': ds[fvar] *= 1E3



"""test plot """
# fig= plt.figure(fignr)
# fignr+=1
# plt.clf()

# plt.plot(ds.time, ds[fvar].sel(lat=70, method='nearest'))


import scipy.signal

    

ds= ds.reindex(lat= list(reversed(ds.lat)))

fignr= 8

filter_value= 0.05

"""method transport, gradient +fraction"""
fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col')


dstot = ds[evar].sum(dim= "WaveNumb")
trans_mean= dstot.mean(dim='time')


line, = axs[0].plot(ds.lat, trans_mean, label= 'Total')


mean= ds[fvar].mean(dim='time')
mean_smth= mean.rolling(lat= 5, center= True).mean()
mean_smth= mean_smth.interpolate_na(dim='lat', method= 'linear', fill_value= 'extrapolate')

mean_grad= mean_smth.differentiate('lat')*1E2/110 #K per 100km #np.gradient(mean)/dy #divergence filter_level= filter_value* np.abs(mean_grad).max()
filter_level= filter_value* np.abs(mean_grad).max()

absfilter= np.abs(mean_grad) < filter_level 


mean_grad_smth= mean_grad.rolling(lat= 5, center= True).mean()
mean_grad_smth= mean_grad_smth.interpolate_na(dim='lat', method= 'nearest', fill_value= 'extrapolate')
filter_level= filter_value* np.abs(mean_grad_smth).max()
absfilter_smth= np.abs(mean_grad_smth) < filter_level    
# lat_change, cross_filter= sign_change_ds(mean_grad_smth)


dlat= ds.lat[1].values- ds.lat[0].values

# mean_grad_sav= scipy.signal.savgol_filter(mean_smth, window_length=rolling_intervall, polyorder=2, deriv=1)
mean_grad_sav= scipy.signal.savgol_filter(mean, window_length=11, polyorder=1, deriv=1)/dlat *1E2/110 #K per 100km
# mean_grad_sav = xr.DataArray(mean_grad_sav, dims= ['lat'], coords=dict(lat= ds.lat), name= 'sav' )
# lat_change, cross_filter_sav= sign_change_ds(mean_grad_sav)

filter_level= filter_value* np.abs(mean_grad_sav).max()
absfilter_sav= np.abs(mean_grad_sav) < filter_level        


axs[1].plot(ds.lat, mean, '--', label='mean')
axs[1].plot(ds.lat, mean_smth, label='mean smooth', color= line.get_color() )

axs[2].plot(ds.lat, mean_grad, ':', color= line.get_color(), label= 'grad smooth' )
axs[2].plot(ds.lat, mean_grad_smth, '--', color= line.get_color(), label= 'smooth grad smooth' )
axs[2].plot(ds.lat, mean_grad_sav, '-', color= line.get_color(), label= 'grad savitzky' )

diff_fac= 100E3 #to compute from W/m / K/100km to W/K

axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter!= 0, trans_mean/mean_grad), ':', label= 'grad mean smooth')

# axs[3].plot(ds.lat, np.ma.masked_where(cross_filter!= 0, trans_mean/mean_grad_smth), '--', label= 'smooth grad mean smooth')
axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_smth!= 0, trans_mean/mean_grad_smth), '--', label= 'smooth grad mean smooth', color= line.get_color())

axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_sav!= 0,trans_mean/mean_grad_sav), '-', label= 'savitzky', color= line.get_color())


for axnr in [0,1, 2,3]:
    if axnr != 1: axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs[axnr].legend()

plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'Mean {evar_name} transport \n per lenth unit [W/m]') #energy transport accross a latitude band
axs[1].set_ylabel(f'{label} \n [{unit}]')
axs[2].set_ylabel(f'{label} gradient \n [{unit} (100 km)'+r'$^{-1}$]')
axs[3].set_ylabel(f'Transport per \n{label} gradient  \n[W/m /{unit}/m]') #energy transport accross a latitude band





"""transport, gradient +fraction"""
latcut= 85


fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col')


dstot = ds[evar].sum(dim= "WaveNumb")
trans_mean_last= dstot.where(ds["time.year"]>= syear_2).mean(dim='time')
trans_mean_first= dstot.where(ds["time.year"]<= eyear_1).mean(dim='time')

line2, =axs[0].plot(ds.lat, trans_mean_first, label=f'{syear_1}-{eyear_1}' )
line1, =axs[0].plot(ds.lat, trans_mean_last, label=f'{syear_2}-{eyear_2}' )

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
axs[1].plot(ds.lat, mean_first, label=f'{syear_1}-{eyear_1}' )
axs[1].plot(ds.lat, mean_last, label=f'{syear_2}-{eyear_2}'  )

axs[2].plot(ds.lat, mean_first_grad_sav, '-', color= line2.get_color() )
axs[2].plot(ds.lat, mean_last_grad_sav , '-', color= line1.get_color() )

diff_fac= 100E3 #to compute from W/m / K/100km to W/K

axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_sav_first!= 0, trans_mean_first/mean_first_grad_sav), '-', color= line2.get_color())
axs[3].plot(ds.lat, diff_fac* np.ma.masked_where(absfilter_sav_last!= 0, trans_mean_last/mean_last_grad_sav), '-', color= line1.get_color())


for axnr in [0, 2,3]:
    axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs[axnr].ticklabel_format(useMathText=True) #to make 10**n on the axis

axs[3].set_xticks(np.arange(-80,80.1, 10), minor=True)
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'Mean {evar_name} transport \n per lenth unit [W/m]') #energy transport accross a latitude band
axs[1].set_ylabel(f'{label} \n [{unit}]')
axs[2].set_ylabel(f'{label} gradient \n [{unit} (100 km)'+r'$^{-1}$]')
axs[3].set_ylabel(f'Transport per \n{label} gradient  \n[W/{unit}]') #energy transport accross a latitude band

axs[0].legend()


if save:
    savedir= Mediadir +'/home/Energy_Transport/Figs/Transport-vs-gradient_2/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    file_name= fvar +'_' + evar
    file_name += f'_M{Member}'
    file_name += '_' +str(syear_1)+'-'+str(eyear_1)+'_'+str(syear_2)+'-'+str(eyear_2)
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 



