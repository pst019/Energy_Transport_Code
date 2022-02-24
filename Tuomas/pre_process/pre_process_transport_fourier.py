import numpy as np
import netCDF4 as nc4
import glob
import xarray as xr
import bottleneck as bn
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

def moving_average(a, n=3):
    """ Moving average with symmetric extension at left boundary"""
    ret = np.cumsum(a, dtype=float,axis = 0)
    ret = np.concatenate([np.flip(ret[:n-1],axis = 0),ret])
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

## GLOBALS

# Start and end years
startyr = 1979
endyr = 2012
numyrs = endyr - startyr + 1
ndays_yr = 365
years = np.arange(startyr,endyr +1,1)
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
# Earth radius
R = 6371e3
datadir = 'Fourier/Waves/'
fname = 'vQtot_D.1979-2012.nc'
# fname = 'vQtot'
# varname =  'v'
varname = 'vQtot'

d = xr.open_dataset(datadir + fname)
d = d.sel(time=~((d.time.dt.month == 2) & (d.time.dt.day == 29)))
wn = d.wn.data
lat = d.lat.data
#
d2 = np.reshape(d[varname].data,(numyrs,ndays_yr,len(wn),len(lat)))
dm = np.mean(d2,axis=0)
# Number of days for MA
ndays = 30
ma = bn.move_mean(dm,window=ndays,min_count = 1,axis = 0)


temp = d2 - dm
de_seasonalized = np.reshape(temp, (ndays_yr*numyrs,len(wn),len(lat)))
# plt.plot(d.time,np.mean(np.mean(de_seasonalized,axis = 1),axis = 1))
# plt.plot(d.time,moving_average(np.mean(np.mean(de_seasonalized,axis = 1),axis = 1),n=10))
# 2nd number of days for MA
ndays2 = 7
de_seasonalized = bn.move_mean(de_seasonalized, window = ndays2, min_count =1,axis = 0)
de_seasonalized /= 1e15
deses2 = np.zeros((de_seasonalized.shape[0],3,de_seasonalized.shape[2]))
deses2[:,0,:] = de_seasonalized[:,0,:]
deses2[:,1,:] = np.sum(de_seasonalized[:,1:4,:],axis =1)
deses2[:,2,:] = np.sum(de_seasonalized[:,4:,:],axis =1)
scales = np.array([0,1,2])
data = xr.DataArray(deses2, coords = [d.time.data,scales,lat], dims = ['time','scale','lat'])
data.name = varname
data.attrs['units'] = 'PW [10^15 W]'
data.to_netcdf(datadir + 'vQtot_D_deseas.1979-2012.nc')
# plt.show()
#Deseasonalizing data



#
# D = 2*np.pi*R*np.cos(np.deg2rad(lat))
# print(D.shape)
# print(d[varname].data.shape)
# d[varname].data *= D
#
# d2 = d.resample(time='D').mean()
# print(d2)
# d2.to_netcdf(datadir + 'vQtot_D.1979-2012.nc')

#
# first = True
# time = pd.date_range(start = '01/01/1979',end = '01/01/2013',freq='6H')[:-1]
# for year in years:
#     for month in months:
#         nc = nc4.Dataset(f'{datadir}{fname}.{year}.{month}.WN20.nc')
#         if first:
#             t_array = np.array(nc.variables[varname]).transpose(1,0,2)
#             print(t_array.shape)
#             first = False
#         else:
#             t_array = np.concatenate([t_array,np.array(nc.variables[varname]).transpose(1,0,2)])
# lats = np.array(nc.variables['lat'])
# wn = np.array(nc.variables['WaveNumb'])
# # lons = np.array(nc.variables['g0_lon_2'])
# # print(lats.shape)
# # print(lons.shape)
# # print(t_array.shape)
# data = xr.DataArray(t_array, coords = [time,wn,lats], dims =['time','wn','lat'])
# data.name = 'vQtot'
# data.to_netcdf(datadir + 'vQtot_6h.1979-2012.nc')
