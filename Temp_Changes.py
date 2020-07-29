#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to compute the global mean air temperature using HadCRU data
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats

# Read in properly
fin="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/HadCRU.txt"
nt=2020-1850+1
count=1
row=0
t=np.zeros((nt,14))
with open(fin) as f:
    for l in f:
        if count % 2. == 1:
            t[row]=l.split()
            row+=1
        count+=1
t=t[:-1,:]
yr=t[:,0]; tann=t[:,-1]
fig,ax=plt.subplots(1,1)
ax.plot(yr,tann); ax.grid()

# Cut out the 1981-2010 mean and see how much warmer than PI
ctrl=tann[np.logical_and(yr>1980,yr<2011)]
pi=tann[np.logical_and(yr>1849,yr<1880)]
dpi=np.mean(ctrl)-np.mean(pi); print("Warming so far = %.2fC" % dpi)


# Trend over 1979-2010
years=np.arange(1979,2020)
trend,intercept,lower,upper=\
stats.theilslopes(tann[np.logical_and(yr>1978,yr<2020)], years, 0.95)
print("Trend = %.3fC/decade (%.3f,%.3f)" % (trend*10,lower*10,upper*10))




