#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Extracts ERA5 circulation during the most extreme low pressure events in the 
Everest reconstruction
"""
import cdsapi, os, pandas as pd, numpy as np
fin="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv"
pest_sum=pd.read_csv(fin,parse_dates=True,index_col=0,names=["p",])

# Lowest 1% of all events
n=np.float(len(pest_sum))
ntop=np.int(np.round(n*(1-0.999)))
pest_sum_sort=pest_sum.sort_values(["p"])
top_samp=pest_sum_sort.iloc[:ntop]
top_samp=top_samp.sort_index()
idx=np.ones(len(pest_sum)).astype(np.bool)

# Lowest n events -- avoiding clustering (must be separated by at least k days)
n=20
k=2
dates=[]
j=0
for i in range(n):
    pest_sum_sort=pest_sum_sort.loc[idx].sort_values(["p"])
    dates.append(pest_sum_sort.index[0])
    # This screens out everything <=2 days of the minima
    idx=np.abs((pest_sum_sort.index-pest_sum_sort.index[0]).days)>k
    

c = cdsapi.Client()

direct="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Comps/"
area=['50','30','0','150']
count=0
for d in dates:
    year="%.0f"%np.float(d.year)
    month="%.0f"%np.float(d.month)
    day="%02d"%np.int(d.day)
    time="%02d:00"%np.int(d.hour)
    out=direct+"big_comp_%.0f.nc"%count
    count+=1
    if os.path.isfile(out): print("Already have file: %s"%out); continue

    c.retrieve(

		    'reanalysis-era5-pressure-levels',

		    {
			'product_type': 'reanalysis',
			'format': 'netcdf',
			'variable': [
			    'geopotential', 
			    'temperature', 'u_component_of_wind', 'v_component_of_wind'
			],
			'pressure_level': [
                            '250','300','500'
            ],
			'year': [
			    year,
			],
			'month': [
			    month,
			],
			'day': [
			    day,
			],
			'time': [
			    time,
			],
			'area': area, 

		    },

		    out)

os.system("cdo -b 64 ensmean %sbig_comp*.nc %sbig_ens_comp.nc"%(direct,direct))