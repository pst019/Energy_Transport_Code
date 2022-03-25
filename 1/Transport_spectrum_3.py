#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup/EC_Earth'
if user=='media':
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/'+ model +'/'   
else:
    Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd



save= True
# save= False
imp= False
imp= True
#
#Member = 1
syear= 1950
eyear= 2100
split_1 = int(2/3*syear + 1/3*eyear)
split_2 = int(2/3*eyear + 1/3*syear)

fignr= 7


Memberlist= [3]#,2,3,4]


##typ = 'Eddies'
#typ = 'Wavelength_smooth'
#typ = 'Wavelength_smooth_2'
#typ = 'Wavelength_smooth_3'

#typ = 'WaveNR'


#unit='PerM'
#unit='LatCycl'
#
#if typ == 'Eddies':
#    varlist= ['vEmc', 'vEse', 'vEte']
##    varlist= ['vQmc', 'vQse', 'vQte']
#    catlist= varlist
#    monthlymean = True
#
#elif typ == 'Wavelength_smooth':
# varlist= ['vEtot', 'vEsetot']
varlist= ['vQtot', 'vQsetot']

monthlymeanlist=[False, True]
##    varlist= ['vQtot']
#    catlist = ['Meri', 'Plan', 'Syno']#, 'Meso']
#    
#    monthlymean = False
#
#
#elif typ == 'Wavelength_smooth_3':
#    varlist= ['vEtot']
##    varlist= ['vQtot']
#    catlist = ['Meri', 'Plan', 'Syno_l', 'Syno_s', 'Meso']
#    
# monthlymean = False
#
#
#elif typ == 'WaveNR':
#    varlist= ['vEtot']
##    varlist= ['vQtot']
#    catlist = ['0', '1-4', '>5']#, '>10']
#    monthlymean = False



# fignr+=1
# plt.clf()
#ax = fig.subplots() #constrained_layout=True)
#ax1= plt.subplot(2, 1, 1)

# fig, ax = plt.subplots(constrained_layout=True)
#fignr+=1
#plt.clf()
#axs=fig.subplots(2, 1, sharex= 'col')

if imp:
    print('Import ...')
    for Member in Memberlist:
        print('Member', Member)
        file_dir= Mediadir + '/Member'+str(Member) +'/'
        
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

Split1= 8000E3
Split2= 4000E3
Split3= 2000E3

var0= varlist[0]
var1= varlist[1]

for lat in np.arange(-80, 85, 10):
    fig= plt.figure(fignr)
    plt.clf()
    
    ax = fig.add_subplot(111)
    
    ax2 = ax.twiny()

    LatCirc = 2* np.pi * a * np.cos(np.deg2rad(lat))
    
    WaveSplit1= LatCirc/Split1
    WaveSplit2= LatCirc/Split2
    WaveSplit3= LatCirc/Split3
    if WaveSplit3 >= 20: WaveSplit3= 19.99
    
    dsn= ds.sel(lat= lat, method= 'nearest')
    dsn= dsn.mean(dim=['time', 'Member'])
    
    
    # plt.plot(dsn.WaveNumb[1:], dsn[var0][1:], 'x')
    # # ax.plot(LatCirc/1E6 * 1/dsn.WaveNumb[1:], dsn[var0][1:], 'x')
    
    # plt.plot([1, 20], [0,0], c= 'k', linewidth= 0.5, linestyle= '--')
    # plt.xlim([1,20])
    # plt.xlabel('Wavenumber')

    # ax.scatter(LatCirc/np.array([Split1, Split2, Split3]), 3*[0])
    ax.plot([0, 0], [min([0] + list(dsn[var0].sel(WaveNumb= [0,1]).values)), max([0] + list(dsn[var0].sel(WaveNumb= [0,1]).values))], '--', c= 'k', lw= 1)
    ax.plot(LatCirc/np.array([Split1, Split1]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).values], '--', c= 'k', lw= 1)
    ax.plot(LatCirc/np.array([Split2, Split2]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).values], '--', c= 'k', lw= 1)
    ax.plot(LatCirc/np.array([Split3, Split3]), [0, dsn[var0].sel(WaveNumb= 1+ WaveSplit3//1).values], '--', c= 'k', lw= 1)


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
    
    ax.text(WaveSplit1/2, float(dsn[var0].sel(WaveNumb= WaveSplit1//2)),  str(np.round(100* vEteplan/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEsetot'
    vEsesyno= float( dsn[var1].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var1].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit2%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0)  )

    ax.text( (WaveSplit1 + WaveSplit2)/2, 0,  str(np.round(100* vEsesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEtot'
    vEtesyno= -1*vEsesyno + float( dsn[var0].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit1//1, dsn.WaveNumb < 1+ WaveSplit2//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit1%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit1//1).fillna(0) 
                       + WaveSplit2%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) )

    ax.text( (WaveSplit1 + WaveSplit2)/2, float(dsn[var0].sel(WaveNumb= WaveSplit1//1+2)),  str(np.round(100* vEtesyno/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )


    # var= 'vEsetot'
    vEsesyno_s= float( dsn[var1].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit2//1, dsn.WaveNumb < 1+ WaveSplit3//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit2%1) * dsn[var1].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) 
                       + WaveSplit3%1 * dsn[var1].sel(WaveNumb= 1+ WaveSplit3//1).fillna(0)  )

    ax.text( (WaveSplit2 + WaveSplit3)/2 , 0,  str(np.round(100* vEsesyno_s/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )

    # var= 'vEtot'
    vEtesyno_s= -1*vEsesyno_s + float( dsn[var0].where(np.logical_and(dsn.WaveNumb > 1+ WaveSplit2//1, dsn.WaveNumb < 1+ WaveSplit3//1)).sum(dim='WaveNumb') 
                      + (1- WaveSplit2%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit2//1).fillna(0) 
                       + WaveSplit3%1 * dsn[var0].sel(WaveNumb= 1+ WaveSplit3//1).fillna(0) )

    ax.text( (WaveSplit2 + WaveSplit3)/2, float(dsn[var0].sel(WaveNumb= (WaveSplit2)//1+3)),  str(np.round(100* vEtesyno_s/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )



    # var= 'vEtot'
    vEmeso= float(dsn[var0].where(dsn.WaveNumb > 1+ WaveSplit3//1).sum(dim='WaveNumb')
                    +  (1-WaveSplit3%1) * dsn[var0].sel(WaveNumb= 1+ WaveSplit3//1).fillna(0) )
    ax.text( min([WaveSplit3 +2, 19]), 0 ,  str(np.round(100* vEmeso/vEtotal,1) )+'%', horizontalalignment='center', verticalalignment='center', weight='bold'  )


    
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
        savedir= '../Figs/PowerSpec/correct/'
        if not os.path.exists(savedir): os.makedirs(savedir)
        
        file_name= 'PowerSpec_M'
        # file_name= typ+'_'+unit+'_M' 
        for Member in Memberlist: file_name += str(Member)
        file_name += '_lat'+str(lat)
        file_name += '_' +str(syear)+'-'+str(eyear)+'_'
        for var in varlist: file_name += var+'+'
        savefile= savedir+ file_name
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 
