#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 09:19:25 2022

@author: pst019
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd

from functions import *
#import functions

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
#    Mediadir= '/media/'+user+'/Backup/'


energytyp= 'E' #total energy
var= 'vEtot'

#energytyp= 'Q'



file_dir= Mediadir + 'data/Energy_Transport/'+ model +'/EnergySplit/res_0.5x0.5/Waves/'

year= 2000
month = 1

ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()


