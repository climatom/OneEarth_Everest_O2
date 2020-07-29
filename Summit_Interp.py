#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script interpolates quantities to the summit of Everest, and it also 
evaluates the temperature lapse rate. Note that function is slow (but steady!)
"""
import xarray as xa
import numpy as np
import pandas as pd
from numba import jit

# Functions
@jit
def interp_vert(yin,zs,zi):
    nt=yin.shape[0]
    yout=np.zeros(nt)
    for i in range(nt):
        yout[i]=np.interp(zi,zs[i,::-1],yin[i,::-1])
    return yout


# Assign filenames and read in
di="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/"
fin=di+"ERA5i.nc"
fout_p=di+"Scol_press_Era5.csv"
fout_p_bal=di+"Bal_press_Era5.csv"
fout_t=di+"Scol_temp_Era5.csv"
fout_t2=di+"Summit_temp_Era5.csv"
fout_l=di+"Lapse_Era5.csv"
fout_w=di+"Summit_wind_Era5.csv"
data=xa.open_dataset(fin)
scol=7945.
bal=8430.
summit=8848.
g=9.80665

# Extract/manipulate
p=data.level
z=np.squeeze(data.z)/g
nt=z.shape[0]
t=data.t.squeeze()
w=np.squeeze(np.sqrt(data.u**2+data.v**2))

# Interp col pressure and write out
prep=np.tile(p,(nt,1))
pcol=interp_vert(prep,z,scol)
pcol=pd.Series(pcol,index=data.time.values[:])
pcol.to_csv(fout_p)

# Interp temp and write out
tcol=interp_vert(t,z,scol)
tcol=pd.Series(tcol,index=data.time.values[:])
tcol.to_csv(fout_t)

# Repeat temp interpolation for the summit
tsum=interp_vert(t,z,summit)
tsum=pd.Series(tcol,index=data.time.values[:])
tsum.to_csv(fout_t2)

# Interp wind (to the summit) and write out
wsum=interp_vert(w,z,summit)
wsum=pd.Series(wsum,index=data.time.values[:])
wsum.to_csv(fout_w)

# Interp pressure (to the balcony) and write out
pbal=interp_vert(prep,z,bal)
pbal=pd.Series(pbal,index=data.time.values[:])
pbal.to_csv(fout_p_bal)

# Get the lapse rate as dT/dz -- between 300 and 350 hPa surfaces
lapse=(np.squeeze(data["t"])[:,p==400].values-\
       np.squeeze(data["t"])[:,p==300]).values/\
(z[:,p==400.].values-z[:,p==300].values)
lapse=pd.Series(np.squeeze(lapse),index=data.time.values[:])
lapse.to_csv(fout_l)





