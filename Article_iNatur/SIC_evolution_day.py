#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import matplotlib.path as mpath
import pandas as pd
import matplotlib as mpl
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
user = os.getcwd().split('/')[2]

if user == 'pst019':
    Mediadir = '/media/'+user+'/Backup/'

else:
    Mediadir = '/run/media/pst019/Backup/'


#time = pd.date_range("2000-01-01", freq="D", periods=365)

attrs = {"units": "days since 2000-01-01"}


save = False
# save = True
imp = False
imp = True

# syear = 2000
# eyear = 2014

selyear= 2012
running_mean_years= 9
meandays= 30


syear = selyear - running_mean_years//2
eyear = selyear + running_mean_years//2

monthlist = np.arange(1, 13)


latextend= 70


fignr = 1

var = 'SIC'
varname = 'ci'

file_dir = Mediadir + 'data/Energy_Transport/ERA5/' 


# mvar_list= ['SWS', 'SWSD', 'LWS', 'LWSD', 'SSH', 'SLH']
# mvarname_list= ['var176']#, 'var169', 'var177', 'var175', 'var146', 'var147']

mvar_list= ['SWS',  'LWS', 'SSH', 'SLH']
mvarname_list= ['var176', 'var177', 'var146', 'var147']

# mvar_list= ['SWS']
# mvarname_list= ['var176']

# mvar_list= ['LWS']
# mvarname_list= ['var177']

# mvar_list= ['SWS', 'LWS']
# mvarname_list= ['var176', 'var177']


# mvar= 'SSH'
# mvarname= 'var146'



evar= 'vQsum'
elat= 70

evar_list= ['vQsum', 'vEsum']
elat= 70


if imp:

    for year in range(syear, eyear+1):
        print(year)
        for month in monthlist:
            if year == syear and month == monthlist[0]:
                ds0 = xr.open_dataset(
                    file_dir + var + '/'+ 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc')#.resample(time='1M').mean()


                ds0 = ds0.where(ds0.lat >= latextend, drop= True)
                weights = np.cos(np.deg2rad(ds0.lat - 0.5))
                ds0 = ds0.weighted(weights).sum(("lon", "lat")) * 111**2
                
                
            else:
                ds1 = xr.open_dataset(
                    file_dir + var + '/'+ 'SIC.'+str(year)+'.'+str(month).zfill(2)+'.nc')#.resample(time='1M').mean()
                
                ds1 = ds1.where(ds1.lat >= latextend, drop= True)
                ds1 = ds1.weighted(weights).sum(("lon", "lat")) * 111**2                
                
                
                ds0 = xr.concat([ds0, ds1], dim='time')

    ds0 = ds0.rename_vars({varname: var})

    
    #to get monthly files
    for mi, mvar in enumerate(mvar_list):
        print(mvar)
        mvarname= mvarname_list[mi]
    

        for year in range(syear, eyear+1):
        # for year in range(2016, eyear+1):
            print(year)
            for month in monthlist:
                print(month)
                
                if year == syear and month == monthlist[0]:
                    mds0 = xr.open_dataset(file_dir + mvar + '_d/'+ mvar +'.'+str(year)+'.'+str(month).zfill(2)+'.nc')
    
                    mds0 = mds0.where(mds0.lat >= latextend, drop= True)
            
                    weights = np.cos(np.deg2rad(mds0.lat - 0.5))
                    mds0 = mds0.weighted(weights).mean(("lon", "lat"))/ (60**2)
    
                                    
                else:
                    mds1 = xr.open_dataset(file_dir + mvar + '_d/'+ mvar +'.'+str(year)+'.'+str(month).zfill(2)+'.nc')
    
                    mds1 = mds1.where(mds1.lat >= latextend, drop= True)        
                    mds1 = mds1.weighted(weights).mean(("lon", "lat"))/ (60**2)
                            
                            
                    mds0 = xr.concat([mds0, mds1], dim='time')
            
        mds0 = mds0.rename_vars({mvarname: mvar})
        
        """the last loop step to merge the different variables in one ds"""
        if year == eyear: # and month == 12:
            if mvar== mvar_list[0]:
                mds2= mds0
            else:
                mds2= xr.merge([mds2, mds0])


for ei, evar in enumerate(evar_list):

    print(evar)
    for year in range(syear, eyear+1):
        print(year)
        for month in monthlist:
    
            if year == syear and month == monthlist[0]:
                eds0 = xr.open_dataset(file_dir+ '/EnergySplit/res_0.5x0.5/Total/'+ evar +'.'+str(year)+'.'+str(month).zfill(2)+'.nc' )
                eds0 = eds0.sel(lat= 70)
                eds0 = eds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= len(eds0.time), freq='6H') ) )
    
            else:
                eds1 = xr.open_dataset(file_dir+ '/EnergySplit/res_0.5x0.5/Total/'+ evar +'.'+str(year)+'.'+str(month).zfill(2)+'.nc' )
                eds1 = eds1.sel(lat= elat)
                eds1= eds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= len(eds1.time), freq='6H') ) )
    
                eds0 = xr.concat([eds0, eds1], dim='time')

    """the last loop step to merge the different variables in one ds"""
    if year == eyear: # and month == 12:
        if evar== evar_list[0]:
            eds2= eds0
        else:
            eds2= xr.merge([eds2, eds0])

    
file_name= 'SIC_evolution-save-' 
for mvar in mvar_list: file_name += mvar

mds2.to_netcdf(file_dir+ file_name +'.nc')

ds = ds0

ds[var+'_smth'] = ds[var].rolling(time=meandays, center=True).mean()
ds = ds.where(ds["time.dayofyear"] != 366, drop= True) #remove schaltjahr

ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 366)),names=('year', 'dayofyear'))
ds = ds.assign(time=ind).unstack('time')

ds[var+'_clim'] = ds[var+'_smth'].rolling(year=running_mean_years, center=True).mean()


# mds = mds2

# mds = mds.where(mds.lat >= latextend, drop= True)
# weights = np.cos(np.deg2rad(mds.lat - 0.5))
# mds = mds.weighted(weights).mean(("lon", "lat"))/ (24*60**2)


mds= mds2.resample(time='1D').mean()
mds = mds.where(mds["time.dayofyear"] != 366, drop= True) #remove schaltjahr

mds = mds.where(mds["time.year"]<= eyear, drop= True) #to remove the last 6 hours
ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 366)),names=('year', 'dayofyear'))
mds = mds.assign(time=ind).unstack('time')


for mvar in mvar_list:
    mds[mvar+'_clim'] = mds[mvar].rolling(year=running_mean_years, center=True).mean()


eds= eds2.resample(time='1D').mean()
eds['vEsum'] -= eds['vQsum'] #it is the energy transport excluding the moisture transport

eds = eds.where(eds["time.dayofyear"] != 366, drop= True) #remove schaltjahr

ind = pd.MultiIndex.from_product((range(syear, eyear +1) , range(1, 366)),names=('year', 'dayofyear'))
eds = eds.assign(time=ind).unstack('time')

for evar in evar_list:
    eds[evar+'_clim'] = eds[evar].rolling(year=running_mean_years, center=True).mean()

a= 6371E3
Surface_Kugelschnitt= 2*np.pi* a**2*(1- np.cos(np.deg2rad(90- elat)) ) #+ 0.5* (np.sin(np.deg2rad(elat)))**2 )
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(elat))

eds *= LatCirc/Surface_Kugelschnitt

# time_series = ds[var].weighted(weights).sum(
#     ("lon", "lat")) * 111**2  # assume grid cell spacing of 110 km
import matplotlib.dates as mdates


fignr= 1
"""absolute plots"""
fig= plt.figure(fignr)
fignr += 1
plt.clf()

# plt.plot(mds[mvarname])


axs= fig.subplots(2 + len(mvar_list), 1) #+ len(mvar_list),1)

timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(mds['dayofyear'], dtype='timedelta64[D]')-1)

# axs[0].plot( timearray, ds[var].sel(year= selyear) , label=str(selyear)+'unsmth' ) #30day running mean')
    
axs[0].plot( timearray, ds[var+'_smth'].sel(year= selyear) , label=str(selyear) ) #30day running mean')

axs[0].plot(timearray, ds[var+'_clim'].sel(year= selyear) , label='Climatology' )# (9 year running mean)')

axs[0].set_ylabel('Sea ice area [km$^2$]')

axs[0].legend()

axs[0].set_ylim(bottom=0)

import locale
# locale.setlocale(locale.LC_TIME, 'de_DE.UTF-8')
locale.setlocale(locale.LC_TIME, 'norwegian')

axs[0].set_xlim([timearray[0], timearray[-1]])
axs[0].xaxis.set_major_locator(mdates.MonthLocator())
axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))


for mi, mvar in enumerate(mvar_list):
    axs[mi+1].plot(timearray, mds[mvar].sel(year= selyear) , label=str(selyear))
    axs[mi+1].plot(timearray, mds[mvar+'_clim'].sel(year= selyear) , label='Climatology')
    axs[mi+1].set_ylabel(mvar)

    axs[mi+1].set_xlim([timearray[0], timearray[-1]])
    axs[mi+1].xaxis.set_major_locator(mdates.MonthLocator())
    axs[mi+1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

axs[mi+2].plot(timearray, eds[evar].sel(year= selyear) , label=str(selyear))
axs[mi+2].plot(timearray, eds[evar+'_clim'].sel(year= selyear) , label='Climatology')
axs[mi+2].set_ylabel(evar)

axs[mi+2].set_xlim([timearray[0], timearray[-1]])
axs[mi+2].xaxis.set_major_locator(mdates.MonthLocator())
axs[mi+2].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

"""individual anomaly plots"""
fig= plt.figure(fignr) #, (5, 10))
fignr += 1
plt.clf()
axs= fig.subplots(1+ len(mvar_list)+ len(evar_list),1)

# timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# axs[0].plot( timearray, ds[var].sel(year= selyear) , label='Values')
    
axs[0].plot(timearray, ds[var+'_smth'].sel(year= selyear)- ds[var+'_clim'].sel(year= selyear) , label='sjø is anomali' )# (9 year running mean)')
axs[0].plot([timearray[0], timearray[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')

axs[0].set_ylabel('SIC') #('Sea ice area [km$^2$]')

# axs[0].legend()

axs[0].set_xlim([timearray[0], timearray[-1]])
axs[0].xaxis.set_major_locator(mdates.MonthLocator())
axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

for mi, mvar in enumerate(mvar_list):
    # axs[mi+1].plot(timearray ,mds[mvar].sel(year= selyear) - mds[mvar+'_clim'].sel(year= selyear) , label='anomaly' )
    axs[mi+1].plot(timearray , (mds[mvar].sel(year= selyear) - mds[mvar+'_clim'].sel(year= selyear)).rolling(dayofyear=meandays, center=True).mean() , label='anomaly' )
    axs[mi+1].plot([timearray[0], timearray[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs[mi+1].set_ylabel(mvar)


    axs[mi+1].set_xlim([timearray[0], timearray[-1]])
    axs[mi+1].xaxis.set_major_locator(mdates.MonthLocator())
    axs[mi+1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))


for ei, evar in enumerate(evar_list):
    axs[mi+ei+2].plot(timearray , (eds[evar].sel(year= selyear) - eds[evar+'_clim'].sel(year= selyear)).rolling(dayofyear=meandays, center=True).mean() , label='anomaly' )
    axs[mi+ei+2].plot([timearray[0], timearray[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    axs[mi+ei+2].set_ylabel(evar)
    
    
    axs[mi+ei+2].set_xlim([timearray[0], timearray[-1]])
    axs[mi+ei+2].xaxis.set_major_locator(mdates.MonthLocator())
    axs[mi+ei+2].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = str(selyear) + 'anomaly+SIC'
    for nvar in mvar_list + evar_list: file_name += '+'+nvar

    savefile= savedir+ file_name 
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 300)

"""combined anomaly plot"""
fig= plt.figure(fignr, (6,4))
fignr += 1
plt.clf()
ax= fig.subplots() #2,1)

axs = [ax, ax.twinx()]


# timearray= (np.asarray(selyear, dtype='datetime64[Y]')-1970)+(np.asarray(ds['dayofyear'], dtype='timedelta64[D]')-1)

# axs[0].plot( timearray, ds[var].sel(year= selyear) , label='Values')
    
axs[0].plot(timearray, (ds[var+'_smth'].sel(year= selyear)- ds[var+'_clim'].sel(year= selyear))/1E6 , c= 'k', lw= 2, label='Sjøis'  )# (9 year running mean)')

axs[0].plot([timearray[0], timearray[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')

axs[0].set_ylabel('Arktisk sjøis anomali [millioner km$^2$]')

axs[0].legend(loc="lower left")

# axs[0].set_xlim([timearray[0], timearray[-1]])
# axs[0].xaxis.set_major_locator(mdates.MonthLocator())
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))


axs[1].plot(timearray, (mds['LWS'].sel(year= selyear) - mds['LWS_clim'].sel(year= selyear)
            + mds['SSH'].sel(year= selyear) - mds['SSH_clim'].sel(year= selyear) 
            + mds['SLH'].sel(year= selyear) - mds['SLH_clim'].sel(year= selyear) ).rolling(dayofyear=meandays, center=True).mean()
            , label='Varmestråling' )
axs[1].plot(timearray, (mds['SWS'].sel(year= selyear) - mds['SWS_clim'].sel(year= selyear)).rolling(dayofyear=meandays, center=True).mean(), label='Solstråling' )

axs[1].plot(timearray , (eds['vQsum'].sel(year= selyear) - eds['vQsum_clim'].sel(year= selyear)).rolling(dayofyear=meandays, center=True).mean() , label='Fuktighetstransport' )


axs[1].set_ylabel('Stråling og transport anomali [W m$^{-2}$]')
axs[1].legend(ncol= 1)

axs[1].plot([timearray[0], timearray[-1]], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
axs[1].set_xlim([timearray[0], timearray[-1]])
axs[1].xaxis.set_major_locator(mdates.MonthLocator())
axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))


# axs[2].plot(mds.month, mds['LWSD'].sel(year= selyear) - mds['LWSD_clim'].sel(year= selyear), label='LWSD' )
# axs[2].plot(mds.month, mds['SWSD'].sel(year= selyear) - mds['SWSD_clim'].sel(year= selyear), label='SWSD' )
# axs[2].plot(mds.month, mds['SSH'].sel(year= selyear) - mds['SSH_clim'].sel(year= selyear), label='SSH' )
# axs[2].plot(mds.month, mds['SLH'].sel(year= selyear) - mds['SLH_clim'].sel(year= selyear), label='SlH' )
# axs[2].set_ylabel('Flux anomaly [W m$^{-2}$]')
# axs[2].legend(ncol= 4)
# axs[2].plot([1,12], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
# axs[2].set_xlabel('Month')    


#put 0 in the middle
axs[0].set_ylim([-1.2, 1.2])
axs[1].set_ylim([-9, 9])


if save:
    savedir = Mediadir+ '/Dropbox/Energy_Transport/Article_iNatur/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name = 'Evolution_' + var + 'heating_' +str(selyear)

    savefile= savedir+ file_name 
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 300)






