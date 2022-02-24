import numpy as np
import netCDF4 as nc4
import glob
import xarray as xr
import pandas as pd
import bottleneck as bn
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

## GLOBALS

# Start and end years
startyr = 1979
endyr = 2017
numyrs = endyr - startyr + 1
ndays_yr = 365
years = np.arange(startyr,endyr +1,1)
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
# Paths for data directories
datadir = '/media/stallo/SAT/'
varname = '2T_GDS0_SFC'
fname = 'T2m.1979-2017.nc'
varname =  'T2m'
# Defining the southernmost latitude of the dataset
s_lat = 30

d = xr.open_dataset(datadir + fname)
# Dropping leap year days
d = d.sel(time=~((d.time.dt.month == 2) & (d.time.dt.day == 29)))
d = d.sel(lat=(d.lat >= s_lat))
# Conversion to celcius
d['T2m'].data -= 273.15
lats = d['lat'].data
lons = d['lon'].data
# d4 = d.groupby('time.dayofyear').mean('time')
# Reshaping for computing daily means
d2 = np.reshape(d['T2m'].data,(numyrs,ndays_yr,len(lats),len(lons)))
lats = lats[lats>=s_lat]
dm = np.mean(d2,axis = 0)
# Number of days for MA
ndays = 30
ma = bn.move_mean(dm,window=ndays,min_count =1,axis=0)
# Rolling back to match days
ma = np.roll(ma,-int(ndays/2),axis =0)
temp = d2 - ma
#plt.plot(np.mean(np.mean(dm,axis=1),axis=1))
#plt.plot(np.mean(np.mean(ma,axis=1),axis=1))
#plt.figure()
#print(temp.shape)
de_seasonalized = np.reshape(temp, (ndays_yr*numyrs,len(lats),len(lons)))
#plt.plot(d.time,np.mean(np.mean(d['T2m'].data,axis =1),axis =1))
#plt.plot(d.time,np.mean(np.mean(de_seasonalized,axis = 1),axis = 1))
# plt.plot(np.mean(np.mean(ma,axis = 1),axis = 1))
# print(d2.transpose(1,0,2,3).shape)
# print(ma.shape)

#Smoothing the deseasonalized data
ndays2 = 10
de_seasonalized =bn.move_mean(de_seasonalized, window=ndays2,min_count=1,axis = 0)

#plt.plot(d.time,np.mean(np.mean(de_seasonalized,axis=1),axis=1))
#plt.show()
data = xr.DataArray(de_seasonalized, coords = [d.time.data,lats,lons], dims = ['time','lat','lon'])
data.name = varname
data.attrs['units'] = 'C'
data.to_netcdf(datadir + 'T2m_D_deseas.1979-2017.nc')
