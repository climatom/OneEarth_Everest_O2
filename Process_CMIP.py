#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script evaluates the sensitivity of air pressure on the summit of Everest
to changes in global mean air temperature 
"""
import os, xarray as xa, pandas as pd, numpy as np
from scipy.stats import f_oneway

# Functions
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

# Input files etc
dig="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/Gtemps/"
dip="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/VerInterp/"
dio="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/"
figs=[dig + ii for ii in os.listdir(dig) if "field" in ii]
fips=[dip + ii for ii in os.listdir(dip)]
nf=len(fips)

# Extract model names in fips
inst=[]
mod=[]
for m in fips:
    scratch=m.split("/")
    inst.append(scratch[-1].split("_")[0])
    mod.append(scratch[-1].split("_")[1])
out=pd.DataFrame(data={"Inst":inst,"Model":mod}).\
to_csv("/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/mods_used.csv")    

# Other paramaters
w=30 # Years to use for rolling averages

# Preallocate output
outmean=np.zeros((nf,12))
outmin=np.zeros(outmean.shape)
outmax=np.zeros(outmean.shape)
maxdt=np.zeros(21)
count=0
for fp in fips:
    dp=xa.open_dataset(fp)
    dt=xa.open_dataset(fp.replace("VerInterp","Gtemps").replace("rcp85","field"))
    annt=np.squeeze(dt.groupby('time.year').mean("time")) # Annual mean temp
    for m in range(1,13):
        # Identify time
        idxp=dp.time.dt.month==m
        selp=np.squeeze(dp.Pressure[idxp])
        # Compute stats (and convert to Pa)
        annpmean=selp.groupby('time.year').mean("time")/100.
        annpmin=selp.groupby('time.year').min("time")/100.
        annpmax=selp.groupby('time.year').max("time")/100. 
        # Match with annt
        idxp=annpmean.year.isin(annt.year)
        idxt=annt.year.isin(annpmean.year)
        prollmean=moving_average(annpmean.values[idxp],w)
        prollmin=moving_average(annpmin.values[idxp],w)
        prollmax=moving_average(annpmax.values[idxp],w)
        troll=moving_average(annt.tas[idxt],w)
        # Compute the slopes
        outmean[count,m-1]=np.polyfit(troll,prollmean,1)[0]; 
        outmin[count,m-1]=np.polyfit(troll,prollmin,1)[0]; 
        outmax[count,m-1]=np.polyfit(troll,prollmax,1)[0]; 
        
    # Store the temp range in CMIP5
    maxdt[count]=np.max(troll)-np.min(troll)
        
    # Increment counter
    count+=1
        
# Summarise and plot
medmean=np.median(outmean,axis=0); stdmean=np.std(outmean,axis=0)
medmin=np.median(outmin,axis=0);  stdmin=np.std(outmin,axis=0)
medmax=np.median(outmax,axis=0);  stdmax=np.std(outmax,axis=0)  

# Perform ANOVA
fmean,pmean=f_oneway(outmean[:,0],outmean[:,1],outmean[:,2],outmean[:,3],\
                     outmean[:,4],outmean[:,5],outmean[:,6],outmean[:,7],
                     outmean[:,8],outmean[:,9],outmean[:,10],outmean[:,11])
fmin,pmin=f_oneway(outmin[:,0],outmin[:,1],outmin[:,2],outmin[:,3],\
                     outmin[:,4],outmin[:,5],outmin[:,6],outmin[:,7],
                     outmin[:,8],outmin[:,9],outmin[:,10],outmin[:,11])
fmax,pmax=f_oneway(outmax[:,0],outmax[:,1],outmax[:,2],outmax[:,3],\
                     outmax[:,4],outmax[:,5],outmax[:,6],outmax[:,7],
                     outmax[:,8],outmax[:,9],outmax[:,10],outmax[:,11])

# Write out
header=",".join(["%.0f" % i for i in range(1,13)])
np.savetxt(dio+"means.txt",outmean,delimiter=",",header=header)
np.savetxt(dio+"mins.txt",outmin,delimiter=",",header=header)
np.savetxt(dio+"maxs.txt",outmax,delimiter=",",header=header)


