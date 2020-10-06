#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script assesses the most extreme low pressure events. 

Note that the getERA_composite must use the same params as in the second half
of this scipt. 
"""

import os, xarray as xa, numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns

from mpl_toolkits.basemap import Basemap
R=6371*1e3
scalar=2.29*10**-11
def wave_speed(lon_ridge,lon_trough,lat_trough,jet):
    
    lon1=np.radians(lon_ridge)
    lon2=np.radians(lon_trough)
    wlen=np.abs((lon2-lon1))*2
    lat=np.radians(lat_trough)
    wlen=wlen*np.cos(lat)*R
    dist=wlen/2.
    beta=scalar*np.cos(lat)
    cd=-beta*(wlen/(2*np.pi))**2
    u=jet+cd
    delta_time=dist/u # time in seconds
    delta_time_days=delta_time/(60**2*24)
    
    return wlen,dist,cd,u,delta_time,delta_time_days

di="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Comps/"
fdir="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/"
fwind="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Summit_wind_Era5.csv"
ftemp="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Summit_temp_Era5.csv"
fp="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv"

# Circulation
# Find all the circulation files that were extracted
fs=[ii for ii in os.listdir(di) if ".nc" in ii and "big" in ii and "ens" not in \
    ii]
# Extract relevant meta information and set params
data_main=xa.open_dataset(di+"big_ens_comp.nc");
lon=data_main.longitude.values[:]
lat=data_main.latitude.values[:]
lon2,lat2=np.meshgrid(lon,lat)
#everlat=27.9881° N, 86.9250° E
elat=27.9881
elon=86.9250
thumbs=True
rows=np.arange(lat2.shape[0])
cols=np.arange(lat2.shape[1])
erow=rows[np.abs(lat-elat)==np.min(np.abs(lat-elat))]
const_idx=np.logical_and(lon<elon,lon>50)
const_idx2=np.logical_and(lon>elon,lon<130)
wind_idx=np.logical_and(lat>=20,lat<=40)
out=np.zeros((len(fs),9))*np.nan

# Plot thumbnails?
if thumbs:
    fig,ax=plt.subplots(5,4)
    clevs=np.linspace(33,80,50)
    clevs2=np.linspace(8500,10000,4)
    clevs2=[8750,9000,9250,9500]
    i=0
    for f in fs:
        data=xa.open_dataset(di+f)
        w=np.sqrt(data.u**2+data.v**2)
        z=np.squeeze(w[0,0,:,:])
        z.values[z.values[:,:]<33]=np.nan
        z2=np.squeeze(data.z[0,1,:,:])/9.80665
        
        # Get delta lons; mean winds, and wave estimates
        lon_ridge=lon[np.squeeze(z2[erow[0],:]==np.max(z2[erow[0],const_idx]))][0]
        lon_ridge2=lon[np.squeeze(z2[erow[0],:]==np.max(z2[erow[0],const_idx]))][0]
        dlon=elon-lon_ridge
        out[i,0]=lon_ridge
        out[i,1]=dlon
        muwind=np.mean(w[0,2,wind_idx,:])
        out[i,2]=muwind
        out[i,3],out[i,4],out[i,5],out[i,6],out[i,7],out[i,8]=\
           wave_speed(lon_ridge,elon,elat,muwind)

                
        m = Basemap(llcrnrlon=lon2.min(),llcrnrlat=lat2.min(),\
                    urcrnrlon=lon2.max(),urcrnrlat=lat2.max(),\
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='l',projection='merc',\
                    lat_0=40.,lon_0=87.,lat_ts=27.,ax=ax.flat[i],fix_aspect=False)  

        # PLot the ridge
        x,y=m(lon_ridge*np.ones(len(lat)),lat)
        m.plot(x,y,color='red',linestyle='--')
        x,y=m((lon_ridge+dlon*2)*np.ones(len(lat)),lat)
        m.plot(x,y,color='red',linestyle='--')

        x,y=m(lon2,lat2)
        c=m.contourf(x,y,z,levels=clevs,cmap="cividis")
#        c2=m.contour(x,y,z2,color=["grey",],levels=clevs2)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries()
        x,y=m(elon,elat)
        m.scatter(x,y,marker="*",s=70,color="red")
        i+=1
plt.subplots_adjust(bottom=0.15)
cax=fig.add_axes([0.15,0.1,0.72,0.02])     
cbar=fig.colorbar(c, cax=cax, orientation='horizontal')   
cbar.set_ticks([35,40,45,50,55,60,65,70,75,80])
cbar.set_label("Wind Speed (m/s)")
fig.set_size_inches(6,7)
fig.savefig(fdir+"Thumbnails.png",dpi=300)

# Calculate some summaries for the individual waves
mu_travel=np.median(out[:,-1])
pclower_travel=np.percentile(out[:,-1],25)
pcupper_travel=np.percentile(out[:,-1],75)
mu_phase=np.median(out[:,-3])
pclower_phase=np.percentile(out[:,-3],25)
pcupper_phase=np.percentile(out[:,-3],75)
lower_ridge=np.percentile(out[:,0],10)
upper_ridge=np.percentile(out[:,0],90)
lower_ridge_2=lower_ridge+2*np.mean(out[:,1])
upper_ridge_2=upper_ridge+2*np.mean(out[:,1])

# Composite analysis/plot
w=np.sqrt(data_main.u**2+data_main.v**2)
z=np.squeeze(w[0,0,:,:])
z.values[z.values[:,:]<33]=np.nan
z2=np.squeeze(data_main.z[0,1,:,:])/9.80665
fig=plt.figure()
fig.set_size_inches(4,8)
ax1=fig.add_axes([0.13,0.45,0.82,0.22])
clevs=np.linspace(33,75,50)
clevs2=np.arange(6000,10100,200)
#clevs2=[8750,9000,9250,9500]
m = Basemap(llcrnrlon=lon2.min(),llcrnrlat=lat2.min(),\
                    urcrnrlon=lon2.max(),urcrnrlat=lat2.max(),\
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='l',projection='merc',\
                    lat_0=40.,lon_0=87.,lat_ts=27.,ax=ax1,fix_aspect=False)  
m.drawcoastlines(linewidth=0.5)
m.drawcountries()
x,y=m(lon2,lat2)
c=m.contourf(x,y,z,levels=clevs,cmap="cividis")

c2=m.contour(x,y,z2,colors=["k",],linestyle=["--",],levels=clevs2)
ax1.clabel(c2, c2.levels, inline=True, fmt="%.0f", fontsize=10)
m.drawparallels(np.arange(10,50.,15.),labels=[True,False,False,False],fontsize=8)
m.drawmeridians(np.arange(30.,150.,15.),labels=[True,True,True,True],fontsize=8)
x,y=m(elon,elat)
m.scatter(x,y,marker="*",s=70,color="red")
dummy=np.linspace(10,50,100)


# Event analysis -- note that this should produce the same 'n' events. 
n=20
k=2
data=pd.read_csv\
("/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv",\
 index_col=0,parse_dates=True,names=["p"])

# Convert all pressures into anomalies
ind=data.index.dayofyear+data.index.hour/24.
uind=np.unique(ind)
mu=np.zeros(len(uind))
anom=np.zeros(len(data))
clim=np.zeros(len(data))
for i in range(len(uind)):
    mu[i]=np.mean(data.loc[ind==uind[i]])

for i in range(len(uind)):
    anom[ind==uind[i]]=data["p"].loc[ind==uind[i]]-mu[uind==uind[i]]
    clim[ind==uind[i]]=mu[uind==uind[i]]
anom=pd.Series(anom,index=data.index)

# Lowest n events -- avoiding clustering (must be separated by at least k days)
dates=[]
pest_sum_sort=data.sort_values(["p"])
idx=np.ones(len(pest_sum_sort)).astype(np.bool)
for i in range(n):
    pest_sum_sort=pest_sum_sort.loc[idx].sort_values(["p"])
    dates.append(pest_sum_sort.index[0])
    idx=np.abs((pest_sum_sort.index-pest_sum_sort.index[0]).days)>k
   
# Composite the anomaly
nhours=240
companom=np.zeros((nhours*2+1,n))    
for d in range(len(dates)):
    dt=np.abs((data.index-dates[d]).total_seconds()/3600.)
    companom[:,d]=anom[dt<=nhours]
x=np.arange(-240,241)/24.
y=np.mean(companom,axis=1)
ystd=np.std(companom,axis=1)
ax2=fig.add_axes([0.13,0.77,0.82,0.22])
ax2.fill_between(x,y-ystd,y+ystd,color='k',alpha=0.2)
ax2.plot(x,y,color='k')
ax2.grid()

#plt.subplots_adjust(bottom=0.2,left=0.05,right=0.9)
cax=fig.add_axes([0.13,0.38,0.82,0.03])
cbar=fig.colorbar(c, cax=cax, orientation='horizontal')
cbar.set_ticks([35,40,45,50,55,60,65,70,75])
cbar.set_label("Wind Speed (m/s)")
ax2.set_xlim(-10,10)
ax2.set_xlabel("Time since minima (days)")
ax2.set_ylabel("Anomaly (hPa)")
ax2.plot([mu_travel,mu_travel],[-13,5],color='red')
ax2.plot([-mu_travel,-mu_travel],[-13,5],color='red')
ax2.fill_between([pclower_travel,pcupper_travel],[-13,-13],[5,5],\
                 color='red',alpha=0.2)
ax2.fill_between([-pclower_travel,-pcupper_travel],[-13,-13],[5,5],\
                 color='red',alpha=0.2)
ax2.set_ylim(-13,5)

# Read in wind, pressure, and temp to query relationship between pressure
# and other met vars
wind=pd.read_csv(fwind,index_col=0,parse_dates=True,names=["u",]); 
temp=pd.read_csv(ftemp,index_col=0,parse_dates=True,names=["t",]); 
wind["p"]=data["p"]
temp["p"]=data["p"]
wint_idx=np.logical_or(wind.index.month==12,wind.index.month<3)
spr_idx=wind.index.month==5
samp=np.zeros((n,3))
ndays=1
for d in range(len(dates)):
    samp[d,0]=data.loc[data.index==dates[d]].values[0]
    dt=(wind.index-dates[d]).total_seconds()/(3600*24)
    idx_wind=np.logical_and(dt>=-ndays,dt<=ndays) 
    samp[d,1]=wind.loc[idx_wind]["u"].mean()
    samp[d,2]=wind.loc[idx_wind]["u"].max()
    
# Scatter plots
ax3=fig.add_axes([0.13,0.09,0.82,0.22])
a1=ax3.scatter(wind.loc[wint_idx]["p"],wind.loc[wint_idx]["u"],\
            color="k",s=0.2,alpha=0.1,label='')
sns.kdeplot(wind.loc[wint_idx]["p"],wind[wint_idx]["u"],ax=ax3,color="white")
a2=ax3.scatter(samp[:,0],samp[:,1],color="red",label="mean")
a3=ax3.scatter(samp[:,0],samp[:,2],color="Purple",label="max")
ax3.set_xlabel("P$_{s}$ (hPa)")
ax3.set_ylabel("Wind speed (m s$^{-1}$)")
ax3.legend(loc=1)
ax3.axhline(wind.loc[wint_idx]["u"].mean(),linestyle='-',color='k')
ax3.axvline(wind.loc[wint_idx]["p"].mean(),linestyle='-',color='k')
rwint=wind.loc[wint_idx].corr()
rtwint=temp.loc[wint_idx].corr()
rspr=wind.loc[spr_idx].corr()
ax3.text(332,10,"r = %.2f" % rwint["p"]["u"])
ax3.set_yticks([15,30,45,60,75])
ax3.grid()
fig.savefig(fdir+"Fig3.pdf",dpi=300)

print("Mean travel time = %.2f (%.2f - %.2f) days" % (mu_travel,\
      pclower_travel,pcupper_travel))