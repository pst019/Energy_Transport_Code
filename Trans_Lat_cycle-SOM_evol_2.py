#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:55:54 2020

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup1/data/Energy_Transport/EC_Earth'
elif user=='media':
    Mediadir= '/run/media/pst019/Backup1/data/Energy_Transport/EC_Earth'
else:
    Mediadir= '/nird/projects/nird/NS9063K/Richard_KNMI'
#'/nird/projects/NS9063K/Richard_KNMI'
    
    
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import glob as glob
#import datetime

import pandas as pd
import cartopy.crs as ccrs
import somoclu



save= True
# save= False
imp= False
imp= True
#
compute= False
compute= True


syear_0= 1950
eyear_0= 1980

syear_1= 2066 #1950
eyear_1= 2096 # 1960

# syear_0= 1950
# eyear_0= 1965

# syear_1= 2081 #1950
# eyear_1= 2096 # 1960


#split_1 = int(2/3*syear + 1/3*eyear)
#split_2 = int(2/3*eyear + 1/3*syear)

fignr= 5

SOM_iterations= 1000

# Memberlist= [1,2,3,4,5]
Memberlist= [2,3, 4, 5]

# Member= 3




var= 'vQtot'
latcircle= 70
catlist = ['Syno', 'Plan', 'total']

clevels= []

#cvar, cvarfile, clabel, cunit= 'T2M', 'T2M', 'Temperature', 'K'
cvar, cvarfile, clabel, cunit, clevels= 'Z', 'Z850', 'Geopotential height', 'm', np.arange(-200, 201, 40)


latbound= 50 #for the import

SOM_latlow = 60
SOM_lathigh = 80

# monthlist= [9, 10, 11] #[3, 4, 5] 
# monthlist= [12, 1, 2] #
monthlist= [6, 7, 8]
# monthlist= np.arange(1,13)
nrows = 3
ncols = 4


if imp:
    print('Import ...')
    for Member in Memberlist:
        print('Member', Member)
        file_dir= Mediadir + '/Member'+str(Member) +'/'
        
    #    for var in varlist:
        print(var)
        for year in np.concatenate((np.arange(syear_0, eyear_0+1), np.arange(syear_1, eyear_1+1) )):
            if year%10 == 0: print(year)
            for month in monthlist:
    #                print(month)
                if year == syear_0 and month== monthlist[0]:
                    
    
                    ds0= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds0['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds0.time), freq='6H').to_timestamp()
                    ds0= ds0.sel(lat= 70, method= 'nearest')
                    ds0= ds0.resample(time='1D').mean()
                        
    #                if Member == Memberlist[0]: 
                    dslat= ds0.lat
    #                else: ds0['lat']= dslat
                    
                    
                    cds0= xr.open_dataset(file_dir + cvarfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                    cds0 = cds0.where(cds0.lat > latbound, drop= True)
                    cds0= cds0.resample(time='1D').mean()
                    if cvar in ['Z']:
                        if Member == 1: cds0= cds0.isel(lev= 0)
                        else: cds0= cds0.isel(plev= 0)
                        cds0[cvar]/= 9.81
    
                else:
                    ds1= xr.open_dataset(file_dir + var+'.'+str(year)+'.'+str(month).zfill(2)+ '.WN20.nc')
                    ds1['time']= pd.period_range(start=str(year)+'-'+str(month).zfill(2)+'-01', periods= len(ds1.time), freq='6H').to_timestamp()
                    ds1= ds1.sel(lat= latcircle, method= 'nearest')
                    ds1= ds1.resample(time='1D').mean()
                                           
                    ds1['lat']= dslat
                    
                    ds0= xr.concat([ds0, ds1], dim= 'time')
    
    
                    cds1= xr.open_dataset(file_dir + cvarfile+'_'+str(year)+str(month).zfill(2)+'_daily.nc')
                    cds1 = cds1.where(cds0.lat > latbound, drop= True)
                    cds1= cds1.resample(time='1D').mean()
                    if cvar in ['Z']: 
                        if Member == 1 and year < 2066: cds1= cds1.isel(lev= 0)
                        else: cds1= cds1.isel(plev= 0)
                        cds1[cvar]/= 9.81
    
                    
                    cds0= xr.concat([cds0, cds1], dim= 'time')
    
    
    
                    
            """the last loop step to merge the different variables in one ds"""
            if year == eyear_1 and month == monthlist[-1]:
#                print(ds)
#                    if var== varlist[0]:
                ds2= ds0
                cds2= cds0

#                    else:
#                        ds2= xr.merge([ds2, ds0])
#                        cds2= xr.merge([cds2, cds0])
                                
        if Member == Memberlist[0]:
            ds= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            cds= cds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)

        else:
            ds2= ds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            ds= xr.merge([ds, ds2])
            cds2= cds2.assign_coords(Member= Member).expand_dims('Member', axis= -1)
            cds= xr.merge([cds, cds2])    


cds = cds.drop('plev')
# cds= cds.sel(lon= cds.lon[::2]) #reduce the amount of data

"""to make the map a circle"""
import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


if compute:
    """computation of the categories"""
    
    #get one season 1 - DJF, 2 - MAM, 3 - JJA, 4 - SON
    # cds.where(cds["time.month"]%12 // 3 + 1 == 1, drop= True)
    
    nyears = eyear_1 - syear_1 + eyear_0 - syear_0 + 2
    ds = ds.sel(time= ds["time.dayofyear"] != 366)
    cds = cds.sel(time= cds["time.dayofyear"] != 366)
    
    
    m_cds= cds['Z'].weighted(np.cos(np.deg2rad(cds.lat))).mean(("lat", "lon"))
    
    cds[cvar+ '_geo_ano'] =cds['Z'] - m_cds
    

    
    d= cds[cvar+ '_geo_ano']
    d= d.where((d.lat > SOM_latlow) & (d.lat < SOM_lathigh), drop= True )
    
    # dim_time= len(d.time)
    # dim_x= len(d.lat)
    # dim_y= len(d.lon)
    # n_members= len(d.Member)
     
    d_flat = d.stack(step=('time', 'Member'), cell= ('lat', 'lon')) #np.reshape(d.data, (dim_time*n_members, dim_x*dim_y ) )
    print(d_flat.shape)
    
    print('compute SOMs ', SOM_iterations)

    # Creating instance of SOM
    som = somoclu.Somoclu(ncols,nrows ) #,compactsupport= False)
    
    # Training SOM on data.
    som.train(d_flat, epochs=SOM_iterations)


    # d_flat= d_flat.to_dataset()
    cds= cds.stack(step=('time', 'Member'))
    cds['node_col']=(('step'), som.bmus[:,0])
    cds['node_row']=(('step'), som.bmus[:,1])

    """energy transport"""
    if catlist == ['Syno', 'Plan', 'total']:
    
        ds['Plan']=  ds[var].where(np.logical_and(ds.WaveNumb >= 1, ds.WaveNumb <= 3), drop= True ).sum(dim='WaveNumb')
        ds['Syno']= ds[var].where(ds.WaveNumb >= 4, drop= True ).sum(dim='WaveNumb')
        ds['total']= ds[var].sum(dim='WaveNumb')
        
        ds= ds.drop(var).drop('WaveNumb')
    
    
    
ads= cds.unstack() #all datasets
ads = xr.merge([ads, ds])

fig= plt.figure(fignr, figsize= (2.3*ncols,2.7*nrows +1) )
fignr+=1
plt.clf()
        
#fig, axs = plt.subplots(nrows,ncols, figsize = (10,10))
# d = np.reshape(d,(dim_time, dim_x, dim_y))

def to_perc(value):
    return np.round(value* 100, 1)

def rd6(value): #round to 10-6 
    return np.round(value/ 1E6, 1)

from scipy.stats import ttest_rel
# from matplotlib import rc
# rc('text', usetex=True)


for row in range(nrows):
    temp_row = 0
    for col in range(ncols):
        print('\n ', row, col)
        
        title_text= f'{row}/{col}'
        # bool_array = [np.array_equal(b,[col,row]) for b in som.bmus]

#        t_dat = np.mean(d[bool_array,:,:],axis = 0)
        # t_dat = np.mean(cds[cvar+ '_geo_ano'][bool_array,:,:],axis = 0)
        nodeds= ads.where((ads.node_row== row) & (ads.node_col== col)) #dataset for the SOM node

        nvalues= int((~np.isnan(nodeds.Plan)).sum()) #number of values that are not nan
        print('total values in node:', nvalues)
        
        nvalues_members= (~np.isnan(nodeds.Plan)).sum(dim= 'time').values
        # print('values for the members:', nvalues_members)
        # print('fraction for the members:', to_perc(nvalues_members/len(ads.time) ) )
        
        n_time_1= len(nodeds['time'].where(nodeds['time.year'] <= eyear_0, drop= True))
        nodeds_1 = nodeds.where(nodeds['time.year'] <= eyear_0)
        nvalues_members_1= (~np.isnan(nodeds_1.Plan)).sum(dim= 'time').values

        n_time_2= len(nodeds['time'].where(nodeds['time.year'] >= syear_1, drop= True))
        nodeds_2 = nodeds.where(nodeds['time.year'] >= syear_1)
        nvalues_members_2= (~np.isnan(nodeds_2.Plan)).sum(dim= 'time').values
        
        # print('values for the members and periods (1/2):', nvalues_members_1, nvalues_members_2)        
        # print('fraction for the members and periods (1/2) :', to_perc(nvalues_members_1/n_time_1), to_perc(nvalues_members_2/n_time_2 ) )

        mean_frac_1, mean_frac_2= to_perc(np.mean(nvalues_members_1/n_time_1)), to_perc(np.mean(nvalues_members_2/n_time_2) )
        # print('mean fraction for periods:', mean_frac_1, mean_frac_2)
        # print('mean standard deviation for periods:', to_perc(np.std(nvalues_members_1/n_time_1)), to_perc(np.std(nvalues_members_2/n_time_2) ) )

        tval, pval= ttest_rel( nvalues_members_1/n_time_1, nvalues_members_2/n_time_2)
        if pval < 0.1:
            print('significant')
            title_frac= f' ({mean_frac_1}\% / {mean_frac_2}\%)'
            
            title_frac = r" $\bf{"+title_frac+"}$"
        else:
            print('non significant')
            title_frac= f' ({mean_frac_1}% / {mean_frac_2}%)'

        print(pval)

        title_text += title_frac +'\n'        
        print(title_text)

        Planvalue_members_1= nodeds_1['Plan'].mean(dim= 'time').values
        Planvalue_members_2= nodeds_2['Plan'].mean(dim= 'time').values

        title_values = f'Plan: {rd6(np.mean(Planvalue_members_1))} / {rd6(np.mean(Planvalue_members_2))}'
        tval, pval= ttest_rel( Planvalue_members_1, Planvalue_members_2)
        if pval < 0.1:
            print('Plan: significant')
            title_values = r"$\bf{"+title_values+"}$"
        else:
            print('Plan: non significant')
        print(pval)

        title_text += title_values +'\n'     


        Synovalue_members_1= nodeds_1['Syno'].mean(dim= 'time').values
        Synovalue_members_2= nodeds_2['Syno'].mean(dim= 'time').values

        title_values = f'Syno: {rd6(np.mean(Synovalue_members_1))} / {rd6(np.mean(Synovalue_members_2))}'
        tval, pval= ttest_rel( Synovalue_members_1, Synovalue_members_2)
        if pval < 0.1:
            print('Plan: significant')
            title_values = r"$\bf{"+title_values+"}$"
        else:
            print('Syno: non significant')
        print(pval)


        title_text += title_values        



        # for cat in catlist:
            # print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(nodeds[cat].where(ds['time.year'] <= eyear_0))/1E6 ) )#, float(np.std(ds[cat].where(ds['time.year'] <= eyear_0)[bool_array])/1E6 ) )
            # print(str(row) +'/'+ str(col) +' '+ cat ,float(np.mean(ds[cat].where(ds['time.year'] >= syear_1)[bool_array])/1E6 ), float(np.std(ds[cat].where(ds['time.year'] <= syear_1)[bool_array])/1E6 ) )

        t_dat= nodeds[cvar+ '_geo_ano'].mean(dim=['time', 'Member'])
        
        if (len(clevels) < 1) & (row + col== 0): 
            extr= np.max([t_dat.max(), -t_dat.min()])
            clevels= np.round(np.linspace(-extr, extr, 10))


        ax = fig.add_subplot(nrows, ncols, col+1 + row* ncols, projection=ccrs.NorthPolarStereo())
        ax.coastlines('110m', alpha=0.5)    
        ax.set_extent([-180, 180, latbound, 90], crs=ccrs.PlateCarree())
        ax.gridlines(ylocs=np.arange(0, 90, 10) )
        ax.set_boundary(circle, transform=ax.transAxes)

# #        levels= np.round(np.linspace(-extr, extr, 10), 2)
        c = ax.contourf(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(),
                        cmap= 'RdBu_r', extend= 'both', levels= clevels) #vmin= -1* extr, vmax= extr) 
        ax.contour(t_dat.lon, t_dat.lat, t_dat, transform=ccrs.PlateCarree(), levels= clevels, colors= 'k', linewidths= 1)

        ax.set_title( title_text)
        # ax.set_title(str(row) +'/'+ str(col) + ' ('+ str(np.round(100*sum((ds['time.year'] <= eyear_0).values & bool_array)/sum((ds['time.year'] < eyear_0).values) ,1 ) ) +'% / '+
#                                                   str(np.round(100*sum((ds['time.year'] >= syear_1).values & bool_array)/sum((ds['time.year'] > syear_1).values)  ,1 ) ) +'%)\n'+
#                                   'Plan: ' + str(np.round(float(np.mean(ds['Plan'].where(ds['time.year'] <= eyear_0)[bool_array])/1E6 ),1)) + ' / ' +
#                                         str(np.round(float(np.mean(ds['Plan'].where(ds['time.year'] >= syear_1)[bool_array])/1E6 ),1))    +'\n'+
#                                   'Syno: ' + str(np.round(float(np.mean(ds['Syno'].where(ds['time.year'] <= eyear_0)[bool_array])/1E6 ),1)) + ' / ' +
#                                         str(np.round(float(np.mean(ds['Syno'].where(ds['time.year'] >= syear_1)[bool_array])/1E6 ),1)) ,
#                     fontsize= 10 )
        


        
fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.015])
cbar= fig.colorbar(c, cax=cbar_ax, orientation="horizontal")
cbar.set_label(clabel +' anomaly ['+cunit+']')



if save:
    savedir= '../Figs/SOM_trans-lat-effect_ensembles/'
    if not os.path.exists(savedir): os.makedirs(savedir)
    
    file_name= 'SOM'+str(ncols)+'-'+str(nrows) +'_'
    if len(monthlist) != 12: file_name += f'month{monthlist[0]}-{monthlist[-1]}_'
    file_name +=  cvar+ '-'+ str(SOM_latlow)+'-'+ str(SOM_lathigh) +'_M'
    for Member in Memberlist:
        file_name += str(Member)
    file_name += '_' +str(syear_0)+'-'+str(eyear_0)+'_'+str(syear_1)+'-'+str(eyear_1)
    file_name += var + str(latcircle)
    file_name += f'_iterate-{SOM_iterations}'

    savefile= savedir+ file_name
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 


