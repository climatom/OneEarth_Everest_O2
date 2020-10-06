#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Compute the observed trends in summit pressure as a function of statistic
and month; also assess CMIP5 changes. 
"""

import pymannkendall as mk, pandas as pd, matplotlib.pyplot as plt, numpy as np
import datetime, seaborn as sb, GeneralFunctions as GF
from scipy import stats
from matplotlib.ticker import FormatStrFormatter

def p2o(p):
    
    """
    Function to convert between air pressure and VO2 max -- using
    the regression equation of Bailey (2001)
    """
    frac=0.2095 # Volume fraction of O2 in the atmopshere 
    scalar=1.33322 # mmHg --> hPa
    satvp=GF.satVpBolton(37.0)/100.
    pio=frac*(p-satvp) # partial pressure of oxygen in inspired air
    vo2max=(np.log(pio*(1/scalar))-3.25)/0.0308
    
    return vo2max

# Import
data=pd.read_csv\
("/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv",\
 index_col=0,parse_dates=True,names=["p"])
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

# Preallocate 
out_p=np.zeros((12,3))*np.nan
out_slope=np.zeros((12,3))*np.nan

# End date
end_date=datetime.datetime(year=2020,month=1,day=1)
data=data.loc[data.index<end_date]
for i in range(1,13):
    
    # Mins
    v=data.loc[data.index.month==i].resample("Y").min()
    trend, h, out_p[i-1,0], z, Tau, s, var_s, out_slope[i-1,0],\
    intercept = mk.original_test(v)
    
    # Means
    v=data.loc[data.index.month==i].resample("Y").mean()
    trend, h, out_p[i-1,1], z, Tau, s, var_s, out_slope[i-1,1],\
    intercept = mk.original_test(v)  
    
    # Maxs
    v=data.loc[data.index.month==i].resample("Y").max()
    trend, h, out_p[i-1,2], z, Tau, s, var_s, out_slope[i-1,2],\
    intercept = mk.original_test(v)  
    
theta=np.array([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.])+14  
fig1=plt.figure()
ax=fig1.add_subplot(111, projection='polar')
theta=2*np.pi/366.*theta
for i in range(len(theta)):
    ax.plot([theta[i],theta[i]],[0,0.1],color='k')
idx_sig_min=out_p[:,0]<0.05
idx_sig_mu=out_p[:,1]<0.05
idx_sig_max=out_p[:,2]<0.05
ax.scatter(theta[idx_sig_min], out_slope[idx_sig_min,0], color='b',s=80)
ax.scatter(theta[~idx_sig_min], out_slope[~idx_sig_min,0], color='b',alpha=0.3,s=80)
ax.scatter(theta[idx_sig_mu], out_slope[idx_sig_mu,1], color="k",s=80)
ax.scatter(theta[~idx_sig_mu], out_slope[~idx_sig_mu,1], color="k",alpha=0.3,s=80)
ax.scatter(theta[idx_sig_max], out_slope[idx_sig_max,2], color="red",s=80)
ax.scatter(theta[~idx_sig_max], out_slope[~idx_sig_max,2], color="red",alpha=0.3,s=80)
ax.set_ylim(0,0.1)
ax.set_xticks(theta)
months=["Jan","Feb","Mar","Apr","May","June","July","Aug","Sep",\
                     "Oct","Nov","Dec"]
ax.set_xticklabels(months)

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

fig,ax=plt.subplots(12,1)
sublen=10
for i in range(1,13):
    scratch=data.loc[data.index.month==i]
    early=scratch.loc[scratch.index.year<(1979+sublen)]
    late=scratch.loc[scratch.index.year>(2019-sublen)]
    sb.distplot(early,color='blue',ax=ax.flat[i-1],label="1979-1988")
    sb.distplot(late,color='red',ax=ax.flat[i-1],label="2010-2019")
    ax.flat[i-1].set_xlim(310,342)
    ax.flat[i-1].set_yticklabels([])
    ax.flat[i-1].set_ylabel(months[i-1])
    if i<12: ax.flat[i-1].set_xticklabels([])
    if i==7: ax.flat[i-1].legend(loc=2)
ax.flat[i-1].set_xlabel("P$_{s}$ (hPa)")
fig.set_size_inches(9,10)

# Compute the annual mins, means, maxs
mins=data.resample("Y").min()
mus=data.resample("Y").mean()
maxs=data.resample("Y").max()
years=np.arange(1979,2020)

# Trends as f(t)
# Min
beta_min,int_min,beta_min_lower,beta_min_upper=\
stats.theilslopes(mins, years, 0.95)
dummy_min=years*beta_min+int_min
dummy_min_lower=(years-np.mean(years))*beta_min_lower+np.mean(dummy_min)
dummy_min_upper=(years-np.mean(years))*beta_min_upper+np.mean(dummy_min)

# Mean
beta_mu,int_mu,beta_mu_lower,beta_mu_upper=\
stats.theilslopes(mus, years, 0.95)
dummy_mu=years*beta_mu+int_mu
dummy_mu_lower=(years-np.mean(years))*beta_mu_lower+np.mean(dummy_mu)
dummy_mu_upper=(years-np.mean(years))*beta_mu_upper+np.mean(dummy_mu)

# Max
beta_max,int_max,beta_max_lower,beta_max_upper=\
stats.theilslopes(maxs, years, 0.95)
dummy_max=years*beta_max+int_max
dummy_max_lower=(years-np.mean(years))*beta_max_lower+np.mean(dummy_max)
dummy_max_upper=(years-np.mean(years))*beta_max_upper+np.mean(dummy_max)


# # # From Review
max_resid=np.array([ii[0] for ii in maxs.values])-dummy_max
mu_resid=np.array([ii[0] for ii in mus.values])-dummy_mu
rmumax=np.corrcoef(max_resid,mu_resid)

# Convert to f(Temp) (hard coded from tempChanges.py)
beta_max_temp=beta_max*10*1/0.169
beta_max_lower_temp=beta_max_lower*10*1/0.169
beta_max_upper_temp=beta_max_upper*10*1/0.169

beta_mu_temp=beta_mu*10*1/0.169
beta_mu_lower_temp=beta_mu_lower*10*1/0.169
beta_mu_upper_temp=beta_mu_upper*10*1/0.169

beta_min_temp=beta_min*10*1/0.169
beta_min_lower_temp=beta_min_lower*10*1/0.169
beta_min_upper_temp=beta_min_upper*10*1/0.169

print("Max: dmax/dT = %.2f/decade (%.3f-%.3f)" % \
      (beta_max_temp,beta_max_lower_temp,beta_max_upper_temp))
print("Mean: dmean/dT = %.2f/decade (%.3f-%.3f)" % \
      (beta_mu_temp,beta_mu_lower_temp,beta_mu_upper_temp))
print("Max: dmin/dT = %.2f/decade (%.3f-%.3f)" % \
      (beta_min_temp,beta_min_lower_temp,beta_min_upper_temp))

# Adjust and plot
plt.subplots_adjust(right=0.45)
ax=fig.add_axes([0.47,0.76,0.3,0.2])
ax.plot(years,mins,color='k')
ax.plot(years,dummy_min,color="r",linestyle="-")
ax.plot(years,dummy_min_lower,color="r",linestyle="--")
ax.plot(years,dummy_min_upper,color="r",linestyle="--")
ax.yaxis.tick_right()
ax.set_xticklabels([])
ax.text(1979.2,317,"Ann. Min")
ax.text(1979.2,308.2, "Trend = %.2f (%.2f-%.2f)"%\
        (beta_min*10,beta_min_lower*10,beta_min_upper*10))
ax.set_ylim(308,318)
ax.set_ylabel("P$_{s}$ (hPa)")
ax.grid()
ax.yaxis.set_label_position("right")
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax=fig.add_axes([0.47,0.54,0.3,0.2])
ax.plot(years,mus,color='k')
ax.plot(years,dummy_mu,color="r",linestyle="-")
ax.plot(years,dummy_mu_lower,color="r",linestyle="--")
ax.plot(years,dummy_mu_upper,color="r",linestyle="--")
ax.yaxis.tick_right()
ax.set_xticklabels([])
ax.text(1979.2,332.5,"Ann. Mean")
ax.text(1979.2,329.2, "Trend = %.2f (%.2f-%.2f)"%\
        (beta_mu*10,beta_mu_lower*10,beta_mu_upper*10))
ax.set_ylim(329,333)
ax.set_ylabel("P$_{s}$ (hPa)")
ax.yaxis.set_label_position("right")
ax.grid()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax=fig.add_axes([0.47,0.32,0.3,0.2])
ax.plot(years,maxs,color='k')
ax.plot(years,dummy_max,color="r",linestyle="-")
ax.plot(years,dummy_max_lower,color="r",linestyle="--")
ax.plot(years,dummy_max_upper,color="r",linestyle="--")
ax.yaxis.tick_right()
ax.set_xlabel("Year")
ax.text(1979.2,342.5,"Ann. Max")
ax.text(1979.1,339.2, "Trend = %.2f (%.2f-%.2f)"%\
        (beta_max*10,beta_max_lower*10,beta_max_upper*10))
ax.set_ylim(339,343)
ax.set_ylabel("P$_{s}$ (hPa)")
ax.yaxis.set_label_position("right")
ax.grid()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
fig.savefig\
("/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/obs_dists.png",dpi=300)


#### CMIP5 analysis 

# Opts and params
pclower=5
pcupper=95
pcmed=50
frac=0.2095 # Volume fraction of O2 in the atmopshere 
scalar=1.33322 # mmHg --> hPa

# Filenames
di="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/"
sumf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv"
maxf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/maxs.txt"
minf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/mins.txt"
meanf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/CMIP5/means.txt"
climbsf="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/o2_climbs.csv"
figdir="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/"

# Read them in
ps=pd.read_csv(sumf,index_col=0,parse_dates=True,names=["ps",])
maxs=np.loadtxt(maxf,delimiter=",")
mins=np.loadtxt(minf,delimiter=",")
means=np.loadtxt(meanf,delimiter=",")

# Get observed stats over 1981-2010
ctrl=ps.loc[np.logical_and(ps.index.year>=1981,ps.index.year<=2010)]
cmean=np.array([np.mean(ps.loc[ps.index.month==i].resample("Y").mean()) \
                for i in range(1,13)])
cmin=np.array([np.mean(ps.loc[ps.index.month==i].resample("Y").min()) \
                for i in range(1,13)])
cmax=np.array([np.mean(ps.loc[ps.index.month==i].resample("Y").max()) \
                for i in range(1,13)])

# Now iterate over temps and constuct scenarios: -5 + 5 
tscen=np.linspace(-6,6,100)
muscen=np.zeros((len(tscen),3))
minscen=np.zeros((len(tscen),3))
maxscen=np.zeros((len(tscen),3))

for t in range(len(tscen)):
    # Mins
    minscen[t,0]=np.min(\
                [cmin[i]+np.percentile(mins[:,i],pclower)*tscen[t] \
                 for i in range(12)])
    minscen[t,1]=np.min(\
                [cmin[i]+np.percentile(mins[:,i],pcmed)*tscen[t] \
                 for i in range(12)])
    minscen[t,2]=np.min(\
                [cmin[i]+np.percentile(mins[:,i],pcupper)*tscen[t] \
                 for i in range(12)])
    
    # Means
    muscen[t,0]=np.mean(\
                [cmean[i]+np.percentile(means[:,i],pclower)*tscen[t] \
                 for i in range(12)])
    muscen[t,1]=np.mean(\
                [cmean[i]+np.percentile(means[:,i],pcmed)*tscen[t] \
                 for i in range(12)])
    muscen[t,2]=np.mean(\
                [cmean[i]+np.percentile(means[:,i],pcupper)*tscen[t] \
                 for i in range(12)])    

    # Maxs
    maxscen[t,0]=np.max(\
                [cmax[i]+np.percentile(maxs[:,i],pclower)*tscen[t] \
                 for i in range(12)])
    maxscen[t,1]=np.max(\
                [cmax[i]+np.percentile(maxs[:,i],pcmed)*tscen[t] \
                 for i in range(12)])
    maxscen[t,2]=np.max(\
                [cmax[i]+np.percentile(maxs[:,i],pcupper)*tscen[t] \
                 for i in range(12)])    
 
# Summarise what the sensitivities are
beta_min,int_min=np.polyfit(tscen,minscen[:,1],1)
beta_min_lower,int_min_lower=np.polyfit(tscen,minscen[:,0],1)
beta_min_upper,int_min_upper=np.polyfit(tscen,minscen[:,2],1)
str_min="%.2f (%.2f-%.2f)"%(beta_min,beta_min_lower,beta_min_upper)
print("Mins: %s" % str_min)

beta_mu,int_mu=np.polyfit(tscen,muscen[:,1],1)
beta_mu_lower,int_mu_lower=np.polyfit(tscen,muscen[:,0],1)
beta_mu_upper,int_mu_upper=np.polyfit(tscen,muscen[:,2],1)
str_mu="%.2f (%.2f-%.2f)"%(beta_mu,beta_mu_lower,beta_mu_upper)
print("Means: %s" % str_mu)

beta_max,int_max=np.polyfit(tscen,maxscen[:,1],1)
beta_max_lower,int_max_lower=np.polyfit(tscen,maxscen[:,0],1)
beta_max_upper,int_max_upper=np.polyfit(tscen,maxscen[:,2],1)
str_max="%.2f (%.2f-%.2f)"%(beta_max,beta_max_lower,beta_max_upper)
print("Maxs: %s" % str_max)
    
# Now plot
# Plot combined with CMIP5 sensitivities
medmean=np.median(means,axis=0)
medmin=np.median(mins,axis=0)
medmax=np.median(maxs,axis=0)
 
plt.subplots_adjust(bottom=0.38,top=0.96)
ax=fig.add_axes([0.07,0.05,0.4,0.25], polar=True)
theta=np.array([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,1])+14  
theta=2*np.pi/366.*theta
months=["Jan","Feb","Mar","Apr","May","June","July","Aug","Sep",\
                     "Oct","Nov","Dec"]

ax.set_ylim(1.8,2.7)
ax.set_xticks(theta)
months=["Jan","Feb","Mar","Apr","May","June","July","Aug","Sep",\
                     "Oct","Nov","Dec"]
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_xticklabels(months)
ax.plot(theta,np.append(medmean,medmean[0]),color="k",marker=".",markersize=10,label="Mean")
ax.plot(theta,np.append(medmin,medmin[0]),color="blue",marker=".",markersize=10,label="Min")
ax.plot(theta,np.append(medmax,medmax[0]),color="red",marker=".",markersize=10,label="Max")
ax.legend(loc=8)
plt.subplots_adjust(right=0.45)
ax=fig.add_axes([0.47,0.07,0.3,0.2])
ax.yaxis.tick_right()
ax.grid()
ax.plot(tscen,minscen[:,1],color="blue")
ax.fill_between(tscen,np.min(minscen,axis=1),np.max(minscen,axis=1),color='blue',\
                alpha=0.2)
ax.text(-4,300,str_min,color='blue')
ax.plot(tscen,muscen[:,1],color='k')
ax.fill_between(tscen,np.min(muscen,axis=1),np.max(muscen,axis=1),color='k',\
                alpha=0.2)
ax.plot(tscen,maxscen[:,1],color='red')
ax.fill_between(tscen,np.min(maxscen,axis=1),np.max(maxscen,axis=1),color='red',\
                alpha=0.2)
ax.text(-4,345,str_max,color='red')
ax.text(-4,319,str_mu,color='k')
ax.set_xlabel("$\Delta$T$_{g}$ ($^{\circ}$C)")
ax.set_ylabel("P$_{s}$ (hPa)")
ax.yaxis.set_label_position("right")
ax.set_xlim([-6,6])

# Read in the temp sensitivity and plot 
climbs=pd.read_csv(climbsf,index_col=0,parse_dates=True,names=["p"])
mu=np.mean(climbs); mx=np.max(climbs); mn=np.min(climbs)
ax.scatter(0,mu,color='green')
ax.plot([0,0],[mn,mx],color='green')
fig.savefig(figdir+"Fig5.pdf",dpi=300)

# What is the % change in VO2 max? (For t change of 1C)
# Use the regression to summarise this. 
p1_min=int_min
p2_min=int_min+beta_min
dvo2_min=(p2o(p2_min)/p2o(p1_min)-1)*100.; print("dVo2_min (+1C) = %.2f" % dvo2_min)

p1_mu=int_mu
p2_mu=int_mu+beta_mu
dvo2_mu=(p2o(p2_mu)/p2o(p1_mu)-1)*100.; print("dVo2_mu (+1C) = %.2f" % dvo2_mu)

# Repeat for cooler climate
p1_min=int_min
p2_min=int_min-beta_min
dvo2_min=(p2o(p2_min)/p2o(p1_min)-1)*100.; print("dVo2_min (-1C) = %.2f" % dvo2_min)

p1_mu=int_mu
p2_mu=int_mu-beta_mu
dvo2_mu=(p2o(p2_mu)/p2o(p1_mu)-1)*100.; print("dVo2_mu (-1C) = %.2f" % dvo2_mu)

# How much warming until the min pressures breach the 1979-2019 mean?
dwarm=(data.mean()-p1_min)/beta_min
dcool=(302-int_min)/beta_min
print("It would take %.1fC warming until min catches up to mean" % dwarm)
print("It would take %.1fC cooling until min drops to unclimbable" % dcool)