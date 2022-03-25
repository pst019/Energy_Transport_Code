#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup/'

else:
    Mediadir= '/run/media/pst019/Backup/'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob as glob
import datetime

import pandas as pd
#time = pd.date_range("2000-01-01", freq="D", periods=365)


year = 2000
month= 2
Member = 3




fignr= 1

var= 'vEtot'

file_dir= Mediadir + 'data/Energy_Transport/EC_Earth/Member'+str(Member) +'/'



#attrs = {"units": "6 hours since 2000-01-01"}

ds= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
ds['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds.time), freq='6H').to_timestamp()







