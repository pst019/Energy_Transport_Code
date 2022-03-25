#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

# model='EC_Earth'
model='ERA5'

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/'+ model +'/' 
elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/'+ model +'/'   
else:
    Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI'
#'/nird/projects/NS9063K/Richard_KNMI'

if model=='ERA5': Mediadir += 'EnergySplit/res_0.5x0.5/Waves/'

    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd



save= True
save= False
imp= False
imp= True
#

if model == 'EC_Earth':
    syear= 1950
    eyear= 2100
    Memberlist= [1]#, 2, 3, 4, 5]#,2,3,4]


elif model== 'ERA5':
    syear= 1979
    eyear= 2018
    Memberlist=[1]


split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)


Split1= 8000E3
Split2= False #4000E3
Split_up= 2000E3



fignr= 7



varlist= ['vEtot', 'vEsetot']
# varlist= ['vQtot', 'vQsetot']

monthlymeanlist=[False, True]



if imp:
    print('Import ...')
    for Member in Memberlist:
        file_dir= Mediadir        
        if model== 'EC_Earth': 
            print('Member', Member) 
            file_dir += 'Member'+str(Member) +'/'  
        
        for vi, var in enumerate(varlist):
            monthlymean = monthlymeanlist[vi]
            print(var)
            for year in range(syear, eyear+1):
                if year%10 == 0: print(year)
                for month in range(1,13):
    #                print(month)
                    if year == syear and month== 1:
                        
                        if monthlymean:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds0= ds0.mean(dim = 'time') #calculate the monthly mean
                            ds0= ds0.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))


                        if Member == Memberlist[0]: dslat= ds0.lat
                        else: ds0['lat']= dslat
                        
                        
                    else:
                        if monthlymean:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc').assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))
                        else:
                            ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                            ds1= ds1.mean(dim = 'time')
                            ds1= ds1.assign_coords(time=('time', pd.date_range(str(year)+'-'+str(month), periods= 1)))

                        
                        ds1['lat']= dslat
    
                        
                        ds0= xr.concat([ds0, ds1], dim= 'time')
    
                        
                    """the last loop step to merge the different variables in one ds"""
                    if year == eyear and month == 12:
        #                print(ds)
                        if var== varlist[0]:
                            ds2= ds0
                        else:
                            ds2= xr.merge([ds2, ds0])
                            
        if Member == Memberlist[0]:
            ds= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
        else:
            ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            ds= xr.merge([ds, ds2])


from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)

a= 6371E3
# LatCirc = 2* np.pi * a * np.cos(np.deg2rad(dslat))



var0= varlist[0]
var1= varlist[1]

for lat in [40]: #np.arange(-80, 85, 10):
    fig= plt.figure(fignr)
    plt.clf()
    
    ax = fig.add_subplot(111)
    
    ax2 = ax.twiny()

    LatCirc = 2* np.pi * a * np.cos(np.deg2rad(lat))
    
    WaveSplit1= LatCirc/Split1
    if Split2: WaveSplit2= LatCirc/Split2
    WaveSplit_up= LatCirc/Split_up
    if WaveSplit_up >= 20: WaveSplit_up= 19.99
    
    dsn= ds.sel(lat= lat, method= 'nearest')
    dsn= dsn.mean(dim=['time', 'Member'])
    
    
    # plt.plot(dsn.WaveNumb[1:], dsn[var0][1:], 'x')
    # # ax.plot(LatCirc/1E6 * 1/dsn.WaveNumb[1:], dsn[var0][1:], 'x')
    
    # plt.plot([1, 20], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    # plt.xlim([1,20])
    # plt.xlabel('Wavenumber')

    # ax.scatter(LatCirc/np.array([Split1, Split2, Split_up]), 3*[0])
    ax.plot([0, 0], [min([0] + list(dsn[var0].sel(WaveNumb= [0,1]).values)), max([0] + list(dsn[var0].sel(WaveNumb= [0,1]).values))], '--', c= 'k', lw= 1)
    ax.plot(LatCirc/np.array([Split1, Split1]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).values], '--', c= 'k', lw= 1)
    if Split2: ax.plot(LatCirc/np.array([Split2, Split2]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).values], '--', c= 'k', lw= 1)
    ax.plot(LatCirc/np.array([Split_up, Split_up]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit_up//1).values], '--', c= 'k', lw= 1)


    # ax.plot([-1, 0, 0] + list(dsn.WaveNumb.values[1:]), 2*[dsn[var0].values[0]] + 2*[dsn[var0].values[1]] + list(dsn[var0].values[2:]), color='k')
    ax.step([-1] +list(dsn.WaveNumb.values), [dsn[var0].values[0]] + list(dsn[var0].values), label= 'Total' )#, color='k')
    # ax.fill_between([-1, 0, 0] + list(dsn.WaveNumb.values[1:]), 0, 2*[dsn[var0].values[0]] + 2*[dsn[var0].values[1]] + list(dsn[var0].values[2:]), label= 'Total', alpha=0.5)

    # ax.plot(dsn.WaveNumb, dsn[var0], label='Total')
    # ax.plot([-1, 0, 0] + list(dsn.WaveNumb.values[1:]), 2*[dsn[var1].values[0]] + 2*[dsn[var1].values[1]] + list(dsn[var1].values[2:]), color= 'k')
    ax.step([-1] +list(dsn.WaveNumb.values), [dsn[var1].values[0]] + list(dsn[var1].values), label= 'Stationary' ) #, color='k')
    # ax.fill_between([-1, 0, 0] + list(dsn.WaveNumb.values[1:]), 0,  2*[dsn[var1].values[0]] + 2*[dsn[var1].values[1]] + list(dsn[var1].values[2:]), label= 'Stationary', alpha=0.5)

    # ax.plot(dsn.WaveNumb, dsn[var1], label='Stationary')
    
    # # ax.plot(LatCirc/1E6 * 1/dsn.WaveNumb[1:], dsn[var0][1:], 'x')
    
    # ax.plot([-.5, 20], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    ax.set_xlim([-0.99,20])
    ax.set_xlabel('Wavenumber')
    ax.plot([-1, 20], [0,0], '--', c= 'k', lw= 1)
    ax.set_xticks(np.arange(0, 20.1))
    ax.legend()
    ax.set_ylabel('Meridional energy transport per length unit [W/m]')
    

    
    vEtotal= float(dsn[var0].sum(dim='WaveNumb'))
    vEmer= float(dsn[var0].sel(WaveNumb= 0))
    ax.text(-0.5, float(dsn[var0].sel(WaveNumb= 0))/2 , str(np.round(100* vEmer/vEtotal, 1))+'%'
            , horizontalalignment='center', verticalalignment='center', rotation='vertical' , weight='bold')    

    # var= 'vEsetot'
    vEseplan= float(dsn[var1].where(np.logical_and(dsn.WaveNumb > 0, dsn.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')+  
                    WaveSplit1%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    
    ax.text(WaveSplit1/2, 0 ,  str(np.round(100* vEseplan/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEtot'
    vEteplan= -1*vEseplan + float(dsn[var0].where(np.logical_and(dsn.WaveNumb > 0, dsn.WaveNumb < 1+ WaveSplit1//1)).sum(dim='WaveNumb')+ 
            WaveSplit1%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) )
    
    # ax.text(WaveSplit1/2, float(dsn[var0].sel(WaveNumb= WaveSplit1//2)),  str(np.round(100* vEteplan/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )
    ax.text(WaveSplit1/2, dsn[var0].max(),  str(np.round(100* vEteplan/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )


    # var= 'vEsetot'
    if Split2: vEsesyno= float( dsn[var1].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var1].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit2%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0)  )
    else: vEsesyno= float( dsn[var1].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit_up//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var1].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit_up%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit_up//1).fillna(0)  )

    if Split2: ax.text( (WaveSplit1 + WaveSplit2)/2, 0,  str(np.round(100* vEsesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )
    else: ax.text( (WaveSplit1 + WaveSplit_up)/2, 0,  str(np.round(100* vEsesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEtot'
    if Split2: vEtesyno= -1*vEsesyno + float( dsn[var0].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit2%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) )
    else: vEtesyno= -1*vEsesyno + float( dsn[var0].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit_up//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit_up%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit_up//1).fillna(0) )

    if Split2: ax.text( (WaveSplit1 + WaveSplit2)/2, float(dsn[var0].sel(WaveNumb= WaveSplit1//1+2)),  str(np.round(100* vEtesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )
    else: ax.text( (WaveSplit1 + WaveSplit_up)/2, dsn[var0].max(),  str(np.round(100* vEtesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )


    # var= 'vEsetot'
    if Split2: vEsesyno_s= float( dsn[var1].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit2//1, dsn.WaveNumb < 1+ WaveSplit_up//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit2%1) * dsn[var1].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) 
                       + WaveSplit_up%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit_up//1).fillna(0)  )

    if Split2: ax.text( (WaveSplit2 + WaveSplit_up)/2 , 0,  str(np.round(100* vEsesyno_s/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEtot'
    if Split2: vEtesyno_s= -1*vEsesyno_s + float( dsn[var0].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit2//1, dsn.WaveNumb < 1+ WaveSplit_up//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit2%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) 
                       + WaveSplit_up%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit_up//1).fillna(0) )

    if Split2: ax.text( (WaveSplit2 + WaveSplit_up)/2, float(dsn[var0].sel(WaveNumb= (WaveSplit2)//1+3)),  str(np.round(100* vEtesyno_s/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )



    # var= 'vEtot'
    vEmeso= float(dsn[var0].where(dsn.WaveNumb > 1+ WaveSplit_up//1).sum(dim='WaveNumb')
                    +  (1-WaveSplit_up%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit_up//1).fillna(0) )
    ax.text( min([WaveSplit_up +2, 19]), 0 ,  str(np.round(100* vEmeso/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )


    
    new_tick_locations = np.array(np.arange(2, 20.1, 3))
    
    def tick_function(x):
        V= LatCirc/1E3 * 1/x
        # V = 1/X
        return ["%.0f" % z for z in V]
    
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel("Wavelength [km]")
    
    
    
    if save:
        savedir= '../Figs/Scale+PowerSpec/'
        if not os.path.exists(savedir): os.makedirs(savedir)
        
        file_name= 'PowerSpec'
        if model != 'ERA5':
            file_name += '_M'
            for Member in Memberlist: file_name += str(Member)
        file_name += '_lat'+str(lat)
        file_name += '_' +str(syear)+'-'+str(eyear)+'_'
        for var in varlist: file_name += var+'+'
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 
