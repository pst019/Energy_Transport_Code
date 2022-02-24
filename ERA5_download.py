#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 11:50:04 2021

@author: pst019
"""

import cdsapi
c = cdsapi.Client()








# c.retrieve("reanalysis-era5-pressure-levels",
# {
# "variable": "temperature",
# "pressure_level": "1000",
# "product_type": "reanalysis",
# "year": "2008",
# "month": "01",
# "day": "01",
# "time": "12:00",
# "format": "grib"
# }, "download.grib")









# c.retrieve(
#     'reanalysis-era5-pressure-levels',
#     {
#         'product_type': 'reanalysis',
#         'variable': 'temperature',
#         'pressure_level': '1000',
#         'year': '2008',
#         'month': '01',
#         'day': '01',
#         'time': '12:00',
#         'format': 'netcdf',                 # Supported format: grib and netcdf. Default: grib
#         'area'          : [90, -180, 70, 180], # North, West, South, East.          Default: global
#         'grid'          : [1.0, 1.0],       # Latitude/longitude grid.           Default: 0.25 x 0.25
#     },
#     'era5_temperature_sub_area.nc')         # Output file. Adapt as you wish.


c.retrieve('reanalysis-era5-complete', {
    'class': 'ea',
    'date': '1979-05-01/to/1979-05-31',
    'expver': '1',
    'levtype': 'sfc',
    'param': '38.235',
    'step': '6/12',
    'stream': 'oper',
    'time': '06:00:00/18:00:00',
    'type': 'fc',
    'format': 'netcdf',                 # Supported format: grib and netcdf. Default: grib
    'area'          : [90, -180, 70, 180], # North, West, South, East.          Default: global
    'grid'          : [1.0, 2.0],       # Latitude/longitude grid.           Default: 0.25 x 0.25    
}, 'output.nc')