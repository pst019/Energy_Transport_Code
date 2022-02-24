#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:59:33 2021

@author: pst019
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def PlotWind(lon, lat, u, v, ax, nx= 15, ny=20, arrowtype= 'quiver', alen=None, scale=False, color= 'k'):
    """ winds with barbs
    alen= arrow length: the lenght of the arrows can be specified by alen
    arrowtype= ('quiver', 'barb') - chose one type of arrows, quiver are normal arrows, barbs have a tail
    nx, ny, specifies how many arrows there are in the x and y direction
    """
#    print('plot winds')
    if len(np.shape(lon))==1: #for ERA data lon and lat are put into 2D array
        lon, lat= np.meshgrid(lon, lat)
        
    everyx= np.shape(u)[0]//nx #put a barb in every nth box, such that you have 20 barbs in the longer dim
    everyy= np.shape(u)[1]//ny #20 for a quadratical map , 50 for era
    u = u[0::everyx, 0::everyy]
    v = v[0::everyx, 0::everyy]
    lon= lon[0::everyx, 0::everyy]
    lat= lat[0::everyx, 0::everyy]
    
    
#    if arrowtype== 'barb':
#        map.barbs(Lon,Lat,u,v, np.sqrt(u**2+v**2), cmap= plt.cm.Reds)
        
    if arrowtype== 'quiver':
        if alen != None:
#            Q = ax.quiver(lon,lat,u,v, pivot= 'mid', scale=20* alen, headwidth= 4, width=0.005, transform=ccrs.PlateCarree(), color= color)           
            #this should fix for wrong cartopy arrow direction:
            Q = ax.quiver(lon,lat, u/np.cos(lat /180 * np.pi),v, pivot= 'mid', scale=20* alen,
                          headwidth= 4, width=0.005, transform=ccrs.PlateCarree(), color= color , angles = "xy")
            qk = plt.quiverkey(Q, 0.18, 0.98, alen, str(alen)+' m/s', coordinates='figure', labelpos='W')
        