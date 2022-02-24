import numpy as np
import netCDF4 as nc4
import glob
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
from reg_funcs import linreg
register_matplotlib_converters()

# Datapaths
### Temperature
t_dir = '/media/stallo/SAT/'
t_fname = 'T2m_D_deseas.1979-2017.nc'
t_vname = 'T2m'
### Latent heat
## Wavelets
# Fourier = False
# vq_dir = '/home/the038/workspace/greenland/regressions/Waves/'
# vq_fname = 'vQtot_D_deseas.1979-2017.nc'
# vq_varname = 'vQtot'
## Fourier
Fourier = True
vq_dir = '/home/the038/workspace/greenland/regressions/Fourier/Waves/'
vq_fname = 'vQtot_D_deseas.1979-2012.nc'
vq_varname = 'vQtot'
### Loading files
T = xr.open_dataset(t_dir + t_fname)
if Fourier:
    T = T.sel(time = (T.time.dt.year < 2013))
# print(T)
vQ = xr.open_dataset(vq_dir + vq_fname)
### Latitude where vQ is computed
s_lat = 70
vQ = vQ.sel(lat = (vQ.lat == s_lat))

# print(T.lat)
# print(vQ.lat[::5])
### Timelags in days
tlag = 20
tl = tlag*2 + 1
tlags = np.linspace(-tlag,tlag,tl)

### Planetary or Synoptic transport
rvar = 'plan'
# rvar = 'syn'
if rvar == 'syn':
    regvar = vQ[vq_varname].data[:,2,0]
elif rvar == 'plan':
    regvar = vQ[vq_varname].data[:,1,0]

### Zonal mean of T
T = T.mean(dim = 'lon')
mT = T[t_vname].data
lats = T.lat.data

b = np.zeros((len(lats),tl))
sign = np.zeros((len(lats),tl))
for i in range(len(lats)):
    for k,t in enumerate(tlags):
        X = regvar[int(tlag+t):-int(tlag-t) -1 ]
        y = mT[int(tlag):-int(tlag)-1,i]
        b[i,k] = linreg(X,y)
print(b.shape)
ll,tt = np.meshgrid(np.flip(lats),tlags)
# levels = np.linspace(-np.max(np.abs(b)),np.max(np.abs(b)),21)
levels = np.linspace(-20,20,41)
c = plt.contourf(ll,tt,(np.flipud(np.fliplr(b.T))),cmap='RdBu_r',levels = levels,extend = 'both')
plt.colorbar(c)
plt.xlim(np.max(lats),np.min(lats))
plt.show()
