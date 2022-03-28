#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:52:35 2022

@author: pst019
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from scipy import stats
import xrft
import numpy.fft as npft
import scipy.signal as signal
import dask.array as dsar


params = {'axes.labelsize': 13,
          # 'axes.titlesize': 16,
          'xtick.labelsize': 13,
          'ytick.labelsize': 13,
          'legend.fontsize': 12}
plt.rcParams.update(params)

user = os.getcwd().split('/')[2]

model='ERA5'

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/'
 
datadir = Mediadir + 'data/ERA5_Clim/ERA5_data/'

fignr= 1


year= 2013
month= 1

filedir= datadir + f'levels_500_era5_{year}_'+str(month).zfill(2)+'.nc'

ds= xr.open_dataset(filedir)#.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

ds['plev']= ds['plev']/100





var= 'v'
plev= 500
lat= 50
time= 0


ds1= ds[var].sel(plev= plev, lat= lat).isel(time= 0)

fig= plt.figure(fignr)
fignr+= 1
plt.clf()

ds1.plot()

# fds1 = xrft.fft(ds1, true_phase=True, true_amplitude=True)
fds1= xrft.power_spectrum(ds1)

fig= plt.figure(fignr)
fignr+= 1
plt.clf()

fds1.plot()