#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script: 
    - 1. Figures out dp/dz based reconstruction
    - 2. Bias corrects pressure, as necessary, and adjusts to the summit. 
    - 3. Compute annual and day of year stats (then write out)
    - 4. TBC -- briefly -- illustrating the "variable height" of Everest
    - 5. Convert to Inspired Oxygen and Vo2 max (then write out)
"""

import pandas as pd, matplotlib.pyplot as plt, numpy as np, GeneralFunctions \
as GF, seaborn as sns
from scipy.ndimage import gaussian_filter1d
from statsmodels.distributions.empirical_distribution import ECDF

# Constants
a=29.3 # Rd/g

# Functions
def adjust_p(p1,t1,z1,z2,lapse):
    
    """
    Adjust pressure using the hypsometric equation. See Stull (2015, p. 17)
    
    In: 
        - p1: pressure at z1
        - z1: altitude at p1
        - z2: altitude we want to interpolate to
        - t1: temperature at z1
        - lapse: lapse rate
        
    Out:
        
        - p2: pressure at z2
        - t2: temperature at z2
        - mut: mean temperature in layer between z1 and z2
        
    """
    t2=t1+lapse*(z2-z1)
    mut=0.5*(t1+t2)
    p2=p1*np.exp((z1-z2)/(a*mut))
    
    
    return p2,[t2,mut]

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

# Params
z1=7945. # South Col
z2=8430. # Balcony
z3=8848. # Summit
lowpc=1 # Percentile for annual stats
frac=0.2095 # Volume fraction of O2 in the atmopshere 
scalar=1.33322 # mmHg --> hPa

# Filenames
eraf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Scol_press_Era5.csv"
eraf_bal="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Bal_press_Era5.csv"
eralf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Lapse_Era5.csv"
scolf="/home/lunet/gytm3/Everest2019/AWS/Logging/south_col.csv"
eratf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Scol_temp_Era5.csv"
balf="/home/lunet/gytm3/Everest2019/AWS/Logging/balcony.csv"
outsummit="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv"
outannf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/ann_summit.txt"
outdoyf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/doy_summit.txt"
outsumf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Summit_press_Era5.csv"
outgrad="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/gradient.csv"
outaer="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Aerobic.csv"
outfigobs="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/Fig1.pdf"

# Read in
era_p=pd.read_csv(eraf,parse_dates=True,index_col=0,names=["PRESS",])
era_p_bal=pd.read_csv(eraf_bal,parse_dates=True,index_col=0,names=["PRESS",])
era_l=pd.read_csv(eralf,parse_dates=True,index_col=0,names=["LAPSE",])
era_t=pd.read_csv(eratf,parse_dates=True,index_col=0,names=["T",],dtype=np.float)
scol=pd.read_csv(scolf,parse_dates=True,index_col=0);\
scol_p=scol["PRESS"].astype(np.float)
scol_t=scol["T_HMP"].astype(np.float)+273.15
bal=pd.read_csv(balf,parse_dates=True,index_col=0);\
bal_p=bal["PRESS"].astype(np.float)

# [1] dPdz-based reconstruction
# (a) BIAS CORRECT P at the South Col
pctls=np.arange(1,100,1)
era_idx=era_p.index.isin(scol.index)
scol_idx=scol_p.index.isin(era_p.index)
sub_era_p=np.squeeze(era_p.loc[era_idx])
sub_scol_p=scol_p.loc[scol_idx]
# Plot the obs and ECDFs
scol_ec=ECDF(np.squeeze(sub_scol_p.values[:]))
era_ec=ECDF(np.squeeze(sub_era_p.values[:]))
# Correlations
ridx=np.logical_and(~np.isnan(sub_scol_p.values[:]),\
                    ~np.isnan(sub_era_p.values[:]))
r=np.corrcoef(sub_scol_p[ridx],sub_era_p[ridx])
fig = plt.figure()
fig.set_size_inches(7,6)
ax1 = fig.add_axes([0.08,0.55,0.9,0.4])
ax2 = fig.add_axes([0.08,0.08,0.4,0.4])
scol_p.plot(ax=ax1,color="black")
#bal_p.plot(ax=ax1,color="black",linestyle="--")
ax2.plot(scol_ec.x,scol_ec.y,color="black",label="AWS",linewidth=3)
ax2.plot(era_ec.x,era_ec.y,color="grey",label="ERA5",linewidth=1)
ax2.set_ylim(0,1); ax2.set_xlim(355,387)
ax1.set_ylabel("South Col [hPa]")
ax1.grid()
ax2.grid()
ax2.set_xlabel("South Col [hPa]")
ax2.set_ylabel("Non-exceedance probability")
# Note -- place "sub_" in second-last place to evaluate how well corr works
scol_corr_p=GF.QQmatch3(sub_scol_p,sub_era_p,pctls,np.squeeze(era_p),extrap=False)
scol_corr_p_plot=GF.QQmatch3(sub_scol_p,sub_era_p,pctls,np.squeeze(sub_era_p),\
                             extrap=False)
era_ec_plot=ECDF(np.squeeze(scol_corr_p_plot))
ax2.plot(era_ec_plot.x,era_ec_plot.y,color="red",linestyle="--",alpha=1,\
         label="ERA5-corrected")
ax2.legend(loc=2)
scol_corr_p=pd.Series(scol_corr_p,index=era_p.index)
scol_idx=scol_corr_p.index.isin(bal_p.index)
scol_sub=scol_corr_p.loc[scol_idx]
bal_idx=bal_p.index.isin(scol_sub.index)
bal_sub=bal_p.loc[bal_idx]
era_idx=era_p.index.isin(scol_sub.index)
# Get dpdT
for m in range(1,13):
    scol_day=scol_p.resample("D").mean()
    dp=np.abs(scol_day.values[1:]-scol_day[:-1]); dp=pd.Series(dp,index=scol_day.index[1:])
    pc99=np.nanpercentile(dp.loc[dp.index.month==m],99)
    maxv=np.max(dp.loc[dp.index.month==m])
    print "In month %.0f pc99 was: %.2f" % (m,pc99)

# Get range of pressures in DJF and JAS (2019)
wint_months=[1,2,12]
mons_months=[7,8,9]
djf=scol_p.loc[scol_p.index.month.isin(wint_months)]
mons=scol_p.loc[np.logical_and(scol_p.index.month.isin(mons_months),\
                               scol_p.index.year==2019)]
djf_range=djf.max()-djf.min(); print("DJF range = %.2f" % djf_range)
mons_range=mons.max()-mons.min(); print("Mons range = %.2f" % mons_range)
rat=djf_range/mons_range; print("Wint/Mons = %.2f" % rat)

# (b) estimate pressure at the Balcony
dz=8430.-7945.
#muT_bal=0.5*(era_t.values[era_idx]+era_t.values[era_idx]+era_l.values[era_idx]*dz)
# Note that we must use the mean temp Scol--> Summit for the regression to work
muT_bal=0.5*(era_t.values[era_idx]+era_t.values[era_idx]+era_l.values[era_idx]*(8848-7945))
dpdz=(np.log(bal_p[bal_idx])-np.log(scol_corr_p[scol_idx]))/dz
denom=1./(a*muT_bal)
x=np.squeeze(denom[~np.isnan(dpdz)]); y=np.squeeze(dpdz[~np.isnan(dpdz)].values[:])
r=np.corrcoef(x,y)
ps=np.polyfit(x,y,1) # Enables (log(p1)-log(p2))/dz to be related to a*muTv
f=np.polyval(ps,denom) # gradient (in log space)
p2=scol_sub.values[:]*np.exp(dz*np.squeeze(f))
# Now add plot for the agreement between modelled and observed Balcony
ax3 = fig.add_axes([0.58,0.08,0.4,0.4])
ax3.scatter(p2,bal_sub,s=5,color='k')
xi=np.linspace(bal_sub.min(),bal_sub.max(),100)
ax3.plot(xi,xi,color='red'); 
ax3.grid()
ax3.set_xlim(xi.min(),xi.max())
ax3.set_ylim(xi.min(),xi.max())
ax3.set_xlabel("Estimated Balcony [hPa]")
ax3.set_ylabel("Observed Balcony [hPa]")
# Get MAE
mae=np.nanmean(np.abs(p2-bal_sub)); print ("MAE [method BC] = %.3f" % mae)
mae2=np.nanmean(np.abs(era_p_bal[era_idx].values[:]-bal_sub.values[:])); 
print ("MAE [method direct] = %.3f" % mae2)
ax3.text(335,357,"MAE = %.2f hPa" %mae)
fig.savefig(outfigobs,dpi=300)


# [2] Reconstruct at summit
dz=8848-7945.
muT_sum=(2*era_t.values[:]+era_l.values[:]*dz)*0.5
denom=1/(a*muT_sum)
f=np.polyval(ps,denom) # This is the gradient (dlog(p)/dz)
pest_sum=scol_corr_p.values[:]*np.exp(dz*np.squeeze(f)) 
pest_sum=pd.Series(pest_sum,index=scol_corr_p.index)
pest_sum.to_csv(outsummit)
summit_day=pest_sum.resample("D").mean()
for m in range(1,13):
    dp=np.abs(summit_day.values[1:]-summit_day[:-1]); \
    dp=pd.Series(dp,index=summit_day.index[1:])
    pc99=np.nanpercentile(dp.loc[dp.index.month==m],99)
    maxv=np.max(dp.loc[dp.index.month==m])
    print "In month %.0f pc99 was: %.2f" % (m,pc99)
 
# Write out the gradient (dlog(p)/dz)-- also needed later on
f=pd.Series(np.squeeze(f),index=era_t.index,name="grad")
f.to_csv(outgrad)

#plt.psd(pest_sum.loc[pest_sum.index.month==1].values[:],24,detrend="constant")

# [3] Compute annual and DoY stats
year=pest_sum.index.year
doy=pest_sum.index.dayofyear
uyear=np.unique(year)
udoy=np.unique(doy)
doymin=gaussian_filter1d(np.squeeze(pest_sum.groupby(pest_sum.index.dayofyear).min()),\
                         sigma=7,mode="wrap")
doymax=gaussian_filter1d(np.squeeze(pest_sum.groupby(pest_sum.index.dayofyear).max()),\
                         sigma=7,mode="wrap")
doymean=gaussian_filter1d(np.squeeze(pest_sum.groupby(pest_sum.index.dayofyear).mean()),\
                         sigma=7,mode="wrap")
fig,ax=plt.subplots(1,1)
ax.plot(doymin)
ax.plot(doymean)
ax.plot(doymax)

# For ann percentile and min, loop over the years 
annmin=np.zeros((len(uyear),1))
annpc=np.zeros((len(uyear),1))
for i in range(len(uyear)):
    annmin[i]=np.min(pest_sum.loc[year==uyear[i]])
    annpc[i]=np.percentile(pest_sum.loc[year==uyear[i]],lowpc)
    
outann=np.column_stack((uyear,annmin,annpc))
outdoy=np.column_stack((udoy,doymin,doymean,doymax))
np.savetxt(outannf,outann,header="year,min,%.0fpc"%lowpc,fmt=\
           ["%.0f","%.2f","%.2f"])
np.savetxt(outdoyf,outdoy,header="doy,min,mean,max",fmt=\
           ["%.0f","%.2f","%.2f","%.2f"])

# [4] -- get the difference in climbing height from the South Col upward -- on
# min and max summit pressure days
mu_may_col=scol_corr_p.loc[scol_corr_p.index.month==5].mean()
mu_may_f=f.loc[f.index.month==5].mean()
max_sum=pest_sum.max(); min_sum=pest_sum.min()
dz_easy=(np.log(max_sum)-np.log(mu_may_col))/mu_may_f
dz_tough=(np.log(min_sum)-np.log(mu_may_col))/mu_may_f
delta_diff=(dz_tough)/dz_easy

# [5] Convert pressure to inspired oxygen and then write out
satvp=GF.satVpBolton(37.0)/100.
pio=frac*(pest_sum-satvp) # partial pressure of oxygen in inspired air
vo2max=(np.log(pio*(1/scalar))-3.25)/0.0308
aerobic=pd.DataFrame(data={"ps":pest_sum,"pio":pio,"vo2":vo2max},index=\
                     pest_sum.index)
# Reminder of the threshold for climbing?
thresh_met=3.5*3.5 # ml/kg/min
thresh_pio=scalar*np.exp((0.0308*thresh_met+3.2500)) # hPa
thresh_pest=satvp+thresh_pio/frac
amax=aerobic["vo2"].loc[aerobic["vo2"]==aerobic["vo2"].max()]
amin=aerobic["vo2"].loc[aerobic["vo2"]==aerobic["vo2"].min()]
delta_diff2=amax.values[0]/amin.values[0]-1





