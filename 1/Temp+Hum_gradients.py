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

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save= False
save= True

Member = 3
syear= 1950
# eyear= 1952
#syear= 2090
eyear= 2096



split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 3



varlist=[ 'T2M', 'vEtot'] #important that the u component is first
fvar= 'T2M' #the display variable
label='temperature'
unit='K'
evar='vEtot'
evar_name='energy'


# varlist=['TCWV','vQtot'] #important that the u component is first
# fvar= 'TCWV' #the display variable
# label='water vapor'
# unit=''
# evar='vQtot'
# evar_name='latent energy'

# varlist=['Q850','vQtot'] #important that the u component is first
# fvar= 'Q' #the display variable
# label='Specific humidity'
# unit='g/kg'
# evar='vQtot'
# evar_name='latent energy'

latcut= 85

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'

latgrid= np.arange(30, 85, 1)


""" import data"""
for vi, var in enumerate(varlist):
    print(var)
    # varfile= varfilelist[vi]

    for year in range(syear, eyear+1):
        if year%10 == 0: print(year)
        for month in range(1,13):
    
            if year == syear and month== 1:
                if var in ['vEtot', 'vQtot']:
                    ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time').interp(lat= latgrid)
                else:
                    ds0= xr.open_dataset(file_dir + var+'_'+str(year)+str(month).zfill(2)+'_daily.nc').mean(dim=['lon', 'time']).interp(lat= latgrid)
    
                ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                
    
            else:
                if var in ['vEtot', 'vQtot']:
                    ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time').interp(lat= latgrid)
                else:
                    ds1= xr.open_dataset(file_dir + var+'_'+str(year)+str(month).zfill(2)+'_daily.nc').mean(dim=['lon', 'time']).interp(lat= latgrid)
    
                ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                
                ds0= xr.concat([ds0, ds1], dim= 'time')
    
    
    """the last loop step to merge the different variables in one ds"""
    if year == eyear and month == 12:
        if var in ['Z850', 'Q850']: ds0= ds0.isel(plev= 0)

        if var== varlist[0]:
            ds2= ds0
        else:
            ds2= xr.merge([ds2, ds0])
    
    
ds= ds2

if fvar == 'Q': ds[fvar] *= 1E3

# print(evar)
#     # varfile= varfilelist[vi]

# for year in range(syear, eyear+1):
#     if year%10 == 0: print(year)
#     for month in range(1,13):

#         if year == syear and month== 1:
#             if evar in ['vEtot', 'vQtot']:
#                 eds0= xr.open_dataset(file_dir + evar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')

#             eds0= eds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
            
#         else:
#             if evar in ['vEtot', 'vQtot']:
#                 eds1= xr.open_dataset(file_dir + evar+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').mean(dim = 'time')

#             eds1= eds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
            
#             eds0= xr.concat([eds0, eds1], dim= 'time')

# eds= eds0



"""plot var"""
fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col')


mean= ds[fvar].mean(dim='time')
mean_grad= mean.differentiate('lat')/110E3 *1E6 #np.gradient(mean)/dy #divergence

mean_last= ds[fvar].where(ds["time.year"]>= split_2).mean(dim='time') 
mean_last_grad= mean_last.differentiate('lat')/110E3 *1E6 #np.gradient(mean)/dy #divergence

mean_first= ds[fvar].where(ds["time.year"]<= split_1).mean(dim='time') 
mean_first_grad= mean_first.differentiate('lat')/110E3 *1E6 #np.gradient(mean)/dy #divergence

diff_mean= (ds[fvar].where(ds["time.year"]>= split_2).mean(dim='time') - ds[fvar].where(ds["time.year"]<= split_1).mean(dim='time') ) #.mean(dim='Member')
diff_mean_grad= diff_mean.differentiate('lat')/110E3 *1E6 #np.gradient(mean)/dy #divergence


line, =axs[0].plot(ds.lat, mean, label=f'{syear}-{eyear}' )
line2, =axs[0].plot(ds.lat, mean_last, label=f'{split_2}-{eyear}' )
line3, =axs[0].plot(ds.lat, mean_first, label=f'{syear}-{split_1}' )

axs[1].plot(ds.lat, mean_grad, '--', color= line.get_color() )
axs[1].plot(ds.lat, mean_grad.rolling(lat= 10, center= True).mean(), color= line.get_color() )
axs[1].plot(ds.lat, mean_last_grad.rolling(lat= 10, center= True).mean(), color= line2.get_color() )
axs[1].plot(ds.lat, mean_first_grad.rolling(lat= 10, center= True).mean(), color= line3.get_color() )

axs[2].plot(ds.lat, diff_mean)
axs[3].plot(ds.lat, diff_mean_grad, '--')
axs[3].plot(ds.lat, diff_mean_grad.rolling(lat= 10, center= True).mean(), color= line.get_color() )


for axnr in [1,3]:  
    axs[axnr].plot([ds.lat[0], ds.lat[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'{label} [{unit}]')
axs[1].set_ylabel(f'{label} gradient \n [{unit} (1000 km)$^{-1}$]')

axs[2].set_ylabel(f'{label} change [{unit}]')
axs[3].set_ylabel(f'Change in \n {label} gradient \n [{unit} (1000 km)$^{-1}$]')

axs[0].legend()#ncol= 3)




"""plot evar"""

fig= plt.figure(fignr, figsize= (8, 5))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')


dstot = ds[evar].sum(dim= "WaveNumb")

trans_mean= dstot.mean(dim='time')
trans_mean_last= dstot.where(ds["time.year"]>= split_2).mean(dim='time')
trans_mean_first= dstot.where(ds["time.year"]<= split_1).mean(dim='time')


trans_diffmean= (dstot.where(ds["time.year"]>= split_2).mean(dim='time') - dstot.where(ds["time.year"]<= split_1).mean(dim='time') )


line, = axs[0].plot(ds.lat, trans_mean, label= 'Total')
line2, =axs[0].plot(ds.lat, trans_mean_last, label=f'{split_2}-{eyear}' )
line3, =axs[0].plot(ds.lat, trans_mean_first, label=f'{syear}-{split_1}' )

axs[1].plot(ds.lat, trans_diffmean, color= line.get_color())


for axnr in range(2):  axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'Mean {evar_name} transport \n per lenth unit [W/m]') #energy transport accross a latitude band

#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel(f'{evar_name} transport change \n '+str(split_2)+'-'+str(eyear)+' - '+str(syear)+'-'+str(split_1))

axs[0].legend()


"""fraction"""

fig= plt.figure(fignr, figsize= (8, 5))
fignr+=1
plt.clf()

axs=fig.subplots(2, 1, sharex= 'col')


line, = axs[0].plot(ds.lat, trans_mean/mean_grad.rolling(lat= 10, center= True).mean(), label= 'Total')
axs[0].plot(ds.lat, trans_mean_last/mean_last_grad.rolling(lat= 10, center= True).mean(), label= f'{split_2}-{eyear}' )
axs[0].plot(ds.lat, trans_mean_first/mean_first_grad.rolling(lat= 10, center= True).mean(), label= f'{syear}-{split_1}')

axs[1].plot(ds.lat, trans_diffmean/diff_mean_grad.rolling(lat= 10, center= True).mean(), color= 'k')

for axnr in range(2):  axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'Transport \n per {label} gradient  \n[W/m /{unit}/1000km]') #energy transport accross a latitude band

#axs[1].set_ylabel('2050-2100 - 1950-2000')
axs[1].set_ylabel('Transport change \n per gradient change')



"""transport, gradient +fraction"""
fig= plt.figure(fignr, figsize= (8, 10))
fignr+=1
plt.clf()

axs=fig.subplots(4, 1, sharex= 'col')


line, = axs[0].plot(ds.lat, trans_mean, label= 'Total')
line2, =axs[0].plot(ds.lat, trans_mean_last, label=f'{split_2}-{eyear}' )
line3, =axs[0].plot(ds.lat, trans_mean_first, label=f'{syear}-{split_1}' )

axs[1].plot(ds.lat, mean, label=f'{syear}-{eyear}' )
axs[1].plot(ds.lat, mean_last, label=f'{split_2}-{eyear}' )
axs[1].plot(ds.lat, mean_first, label=f'{syear}-{split_1}' )

axs[2].plot(ds.lat, mean_grad, '--', color= line.get_color() )
axs[2].plot(ds.lat, mean_grad.rolling(lat= 10, center= True).mean(), color= line.get_color() )
axs[2].plot(ds.lat, mean_last_grad.rolling(lat= 10, center= True).mean(), color= line2.get_color() )
axs[2].plot(ds.lat, mean_first_grad.rolling(lat= 10, center= True).mean(), color= line3.get_color() )

axs[3].plot(ds.lat, trans_mean/mean_grad.rolling(lat= 10, center= True).mean(), label= 'Total')
axs[3].plot(ds.lat, trans_mean_last/mean_last_grad.rolling(lat= 10, center= True).mean(), label= f'{split_2}-{eyear}' )
axs[3].plot(ds.lat, trans_mean_first/mean_first_grad.rolling(lat= 10, center= True).mean(), label= f'{syear}-{split_1}')


for axnr in [0, 2,3]:  axs[axnr].plot([-latcut, latcut], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
plt.xlim(-latcut, latcut)
plt.xlabel('Latitude')
axs[0].set_ylabel(f'Mean {evar_name} transport \n per lenth unit [W/m]') #energy transport accross a latitude band
axs[1].set_ylabel(f'{label} \n [{unit}]')
axs[2].set_ylabel(f'{label} gradient \n [{unit} (1000 km)$^{-1}$]')
axs[3].set_ylabel(f'Transport per \n{label} gradient  \n[W/m /{unit}/1000km]') #energy transport accross a latitude band

axs[0].legend()



if save:
    savedir= '../Figs/Transport-vs-gradient/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    file_name= fvar +'_' + evar
    file_name += f'_M{Member}'
    file_name += '_' +str(syear)+'-'+str(split_1)+'_'+str(split_2)+'-'+str(eyear)
    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 



print('In the future there is somewhat more transport in total energy and moisture per gradient of temperature and humidity, respectively. Hence, the transport becomes somewhat more effectiv. However they scale quite well and the change is rather small')
