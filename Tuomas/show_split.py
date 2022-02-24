#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 13:41:33 2021

@author: pst019
"""


a= 6371E3
LatCirc = 2* np.pi * a * np.cos(np.deg2rad(ds.lat))

# varlist= ['vEtot', 'vEsetot']
varlist= ['vQtot', 'vQsetot']


Split1= 8000E3
WaveSplit1= LatCirc/Split1

for var in varlist:

    ds['Meri_'+var]= ds[var].sel(WaveNumb= 0)
    ds['Plan_'+var]= (ds[var].where(np.logical_and(ds.WaveNumb > 0, ds.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')
            +  WaveSplit1%1 * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    ds['Syno_'+var]= ( ds[var].where(ds.WaveNumb > 1+ WaveSplit1//1).sum(dim='WaveNumb') 
                 +(1- WaveSplit1%1) * ds[var].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
