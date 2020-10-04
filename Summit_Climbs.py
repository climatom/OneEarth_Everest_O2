#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script extracts successufl ascents of Everest completed without 
supplemental oxygen, and: 
    1. Computes the mean (& stdev) for the time of day of summits
    2. Extracts the (reconstructed) oxygen availability at the same time
    3. Plots the seasonal cycle 
"""
import datetime, pandas as pd, numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy import stats
import GeneralFunctions as GF
from matplotlib.ticker import FormatStrFormatter

font = {'size'   : 8}

plt.rc('font', **font)

# Function definitions 
def himDate(datecol,hourcol,verb=False):
    
    """
    Simply iterates over datestrings (day/month/year)
    and returns a vector of datetimes
    
    In:
        - datecol     : of form "day/month/year"
        - hourcol     : decimal form
        - verb        : if verb, print outcome from except
        
    Out:
    
        - uc (m/s)    : datetime.datetime
        
    """
    date=[]; years=[]; months=[]; days=[]; hours=[]
    for i in range(len(datecol)):
        year=np.nan; month=np.nan; day=np.nan; day=np.nan
#       print datecol.values[i]
        day=np.int(datecol.values[i][:2])
        month=np.int(datecol.values[i][3:5])
        year=np.int(datecol.values[i][6:10])
#       print hourcol.values[i]
        hour=np.int(np.round(hourcol.values[i])); 
        if hour >23: hour=0; day+=1
        if hour == 0: hour=np.nan; d=np.nan 
        else: d=datetime.datetime(year=year,month=month,day=day,hour=hour)
        date.append(d)
        years.append(year); months.append(month); days.append(day)
        hours.append(hour)
        
        
    return date,[years,months,days,hours]

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


def dp2dz(p1,p2,grad,z1):
    
    """
    Given pressure at lower elevation z1, and at higher elevation z(?), 
    return the required vertical elevation change (dz) to account for the
    difference in pressure -- assuming a gradient of (d log(p) / dz). 
    Also return the absolute height z2
    """
    
    dz=(np.log(p2)-np.log(p1))/grad
    z2=z1+dz
    
    return dz, z2


def dz2p(grad,p1,dz):
    
    """
    Given the pressure at p1, the vertical gradient (grad), and the dz, return
    p2
    """
    
    p2=p1*np.exp(grad*dz)
    
    return p2

def icao(zm):
    
    """
    Uses the ICAO standard atmosphere as defined In Ward et al. (2000) p. 26
    to model the air pressure at height zm (m)
    
    Input: height in metres
    Output: pressure in hPa
    """
    
    zf=zm*3.28084
    p0=760. # mmHg
    p=p0*(288./(288.-1.98*zf/1000.))**-5.256
    p=p*1.33322
    return p

def west93(zm):
 
    """
    Uses the equation provided by Ward et al. (2000) p. 29 to compute air 
    pressure as f(zm). Note: based on results. 
    
    Input: height in metres
    Output: pressure in hPa
    """
    zm=zm/1000.
    p=np.exp(6.63268-0.1112*zm-0.00149*zm**2)*1.33322
    return p
    
    
    
    
# Filenames
finc="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/members_and_hired_everest.csv"
fino="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/summit_recon.csv"
out="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/o2_climbs.csv"
si="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/SI_brief.csv"
si_full="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/SI_full.csv"
fstat="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/doy_summit.txt"
figo="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/summit_o2.pdf"
figo2="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/summit_o2.png"
fwind="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/Summit_wind_Era5.csv"
fgrad="/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/gradient.csv"
figo3="/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/climbs_details.png"

# Read in 
climbs=pd.read_csv(finc)
o2=pd.read_csv(fino,parse_dates=True,index_col=0,names=["p",])
statdata=np.loadtxt(fstat); doy=statdata[:,0]; mins=statdata[:,1]; mus=statdata[:,2]; \
maxs=statdata[:,3]
wind=pd.read_csv(fwind,index_col=0,parse_dates=True,names=["u",]); 
grad=pd.read_csv(fgrad,index_col=0,parse_dates=True,names=["f"])

# Correlation between winds and pressure
wind["p"]=o2["p"]
rall=wind.corr()
idx_wint=np.logical_or(wind.index.month<3,wind.index.month==12)
rwint=wind.loc[idx_wint].corr()

# Pull out no-O2
no=climbs.loc[np.logical_and(climbs["o2none"]==1,climbs["mperhighpt"]==8850)]; 
no2=no.loc[np.logical_and(np.logical_and(no["o2none"]==1,no["smtcnt"]==2),\
                          no["mperhighpt"]==8850)]
no3=no.loc[np.logical_and(np.logical_and(no["o2none"]==1,no["smtcnt"]==3),\
                          no["mperhighpt"]==8850)]
# With o2
wo=climbs.loc[np.logical_and(climbs["o2none"]==0,climbs["mperhighpt"]==8850)] 
wo2=wo.loc[np.logical_and(np.logical_and(wo["o2none"]==0,wo["smtcnt"]==2),\
                          wo["mperhighpt"]==8850)]
wo3=wo.loc[np.logical_and(np.logical_and(wo["o2none"]==0,wo["smtcnt"]==3),\
                          wo["mperhighpt"]==8850)]

    
# Get datetimes
dates1,meta1=himDate(no["msmtdate1"],no["dsmttime1"],verb=False)
dates2,meta2=himDate(no2["msmtdate2"],no2["dsmttime2"],verb=False)
dates3,meta3=himDate(no3["msmtdate3"],no3["dsmttime3"],verb=False)
dates=np.concatenate((dates1,dates2,dates3))
no_comb=pd.concat([no,no2,no3])  
years=np.concatenate((meta1[0],meta2[0]))
months=np.concatenate((meta1[1],meta2[1],meta3[1]))
days=np.concatenate((meta1[2],meta2[2],meta3[2]))
hours=np.concatenate((meta1[3],meta2[3],meta3[3]))

# Ditto for o2
dates1_o,meta1_o=himDate(wo["msmtdate1"],wo["dsmttime1"],verb=False)
dates2_o,meta2_o=himDate(wo2["msmtdate2"],wo2["dsmttime2"],verb=False)
dates3_o,meta3_o=himDate(wo3["msmtdate2"],wo3["dsmttime3"],verb=False)
dates_o=np.concatenate((dates1_o,dates2_o,dates3_o))
years_o=np.concatenate((meta1_o[0],meta2_o[0],meta3_o[0]))
months_o=np.concatenate((meta1_o[1],meta2_o[1],meta3_o[1]))
days_o=np.concatenate((meta1_o[2],meta2_o[2],meta3_o[2]))
hours_o=np.concatenate((meta1_o[3],meta2_o[3],meta3_o[3]))

# Interrogate summit hour (mean and stdev)
hours=np.array([ii.hour for ii in dates if type(ii) == datetime.datetime])
muhour=np.mean(hours)
stdhour=np.std(hours)
#minhour=np.min(hours)
#maxhour=np.max(hours)

# Extract O2 during climbs
o2_climb=np.zeros(len(dates))*np.nan

# Note that here we convert O2 UTC times to NPT
ref_times=o2.index+datetime.timedelta(hours=6)

# Create a secondary time vector -- to store the resultant date
res_date=[]
estimated=0
no_meta=pd.DataFrame()
for i in range(len(dates)):
    idx = ref_times==dates[i]
    if idx.any():
        # For those cases *with* an hour -- straightforward extract
        o2_climb[i]=o2.values[idx]
        res_date.append(dates[i])
        no_meta=no_meta.append(no_comb.iloc[i])
    # For those *without* an hour -- we figure out the *row* in the o2 dataset, 
    # and use this to cut out the *day*, before applying a Gaussian filter 
    # (stdev = 6 hours) and taking the 12th value
    else:
        idx = np.logical_and(np.logical_and(ref_times.year==years[i],\
                             ref_times.month==months[i]),\
                             ref_times.day==days[i])
        date_est=datetime.datetime(year=years[i],month=months[i],day=days[i],\
                                   hour=12)
        if idx.any():
            scratch=np.squeeze(o2.values[idx])
            scratch_hours=np.squeeze(ref_times.hour[idx])
            filt=gaussian_filter1d(scratch,sigma=3.5,mode="constant",cval=\
                                   np.mean(scratch))
            o2_climb[i]=filt[scratch_hours==12] 
            res_date.append(date_est)
            no_meta=no_meta.append(no_comb.iloc[i])
            estimated+=1
            
        else: print dates[i]; continue
no_meta["datetime"]=res_date
no_meta.set_index("datetime",inplace=True)
print("Estimated %.0f" % estimated) 


# Convert to series and write out
o2_climb=pd.Series(o2_climb[~np.isnan(o2_climb)],index=res_date)      
o2_climb.to_csv(out)   
date_mask=o2_climb.isin(dates)
no_meta["Press"]=o2_climb.values[:]
sorted_o2=no_meta.sort_values(by=["Press"])

# Repeat for w_o2
wo2_climb=np.zeros(len(dates_o))*np.nan
res_date_o=[]
estimated_o=0
for i in range(len(dates_o)):
    idx = ref_times==dates_o[i]
    if idx.any():
        # For those cases *with* an hour -- straightforward extract
        wo2_climb[i]=o2.values[idx]
        res_date_o.append(dates_o[i])
    # For those *without* an hour -- we figure out the *row* in the o2 dataset, 
    # and use this to cut out the *day*, before applying a Gaussian filter 
    # (stdev = 6 hours) and taking the 12th value
    else:
        idx = np.logical_and(np.logical_and(ref_times.year==years_o[i],\
                             ref_times.month==months_o[i]),\
                             ref_times.day==days_o[i])
        date_est=datetime.datetime(year=years_o[i],month=months_o[i],day=days_o[i],\
                                   hour=12)
        if idx.any():
            scratch=np.squeeze(o2.values[idx])
            scratch_hours=np.squeeze(ref_times.hour[idx])
            filt=gaussian_filter1d(scratch,sigma=3.5,mode="constant",cval=\
                                   np.mean(scratch))
            wo2_climb[i]=filt[scratch_hours==12] 
            res_date_o.append(date_est)
            
            estimated+=1
            
        else: continue
      
# Convert to series (o2 and no o2)
wo2_climb=pd.Series(wo2_climb[~np.isnan(wo2_climb)],index=res_date_o)      
wo2_climb.to_csv(out.replace(".csv","_with_o2.csv"))

# Print out the total number of climbs in this period (o2 + wo2)
nwith=len(wo2_climb)
nwithout=len(o2_climb)
ntot=nwith+nwithout
frac_no_o2=np.float(len(o2_climb))/ntot*100.
print("Over period of reconstruction (1979-2019), there were %.0f climbs" % \
      ntot + "(%.0f with o2, %.0f without)" % (nwith,nwithout))
print("... %.2f%% were without o2" % frac_no_o2)


# plot the summit O2 -- seasonal cycle in min/max/mean... and add the
# summit climbs 
theta=2*np.pi/366.*doy
fig1=plt.figure()
ax=fig1.add_subplot(111, projection='polar')
# Repeat the 5th points to avoid artefact
ax.plot(np.append(theta,theta[0]), np.append(mus,mus[0]), color="blue",\
        label="This study")
ax.fill_between(np.append(theta,theta[0]),\
                np.append(mins,mins[0]), np.append(maxs,maxs[0]),\
                color="blue",alpha=0.2,linewidth=0)
theta_points=np.array([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,1])+14  
theta_points=2*np.pi/366.*theta_points
months=["Jan","Feb","Mar","Apr","May","June","July","Aug","Sep",\
                     "Oct","Nov","Dec"]
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_ylim(302,np.max(maxs)+5)
ax.set_yticks([310,320,330,340,])
ax.set_xticks(theta_points)
ax.set_xticklabels(months)
doys2=o2_climb.index.dayofyear
theta2=np.pi*2/366.*doys2
ax.scatter(theta2,o2_climb.values[:],s=10,color="red",alpha=0.25)
mus2=np.ones(len(doy))*o2_climb.mean()
mins2=np.ones(len(doy))*o2_climb.min()
maxs2=np.ones(len(doy))*o2_climb.max()
# Get LTM from the entire dataset -- inc. full years only
muall=np.mean(o2.loc[o2.index.year<2020])
ax.plot(theta,np.ones(len(theta))*muall[0],color='k',linestyle="-",label="LTM")
ax.set_title("Summit Pressure (hPa)",pad=15)

##** For review -- add the ICAO mean and the West mean
ic=icao(8850.)
west=west93(8850.)
ax.plot(theta,np.ones(len(theta))*ic,color='k',linestyle="dashed",label="ICAO")
ax.plot(theta,np.ones(len(theta))*west,color='k',linestyle="dotted",label=\
        "West '96")
# Read in the digitized West data
west_stdv=2.024 # hPa (see Excel workbook on OneDrive)
west_dig=pd.read_csv(\
    "/home/lunet/gytm3/Everest2019/Research/OneEarth/Data/west1993.csv",\
    index_col=0)
ax.plot(theta_points,west_dig["hPa"].values[:],linestyle="dashdot",color='k',\
        label="West '83")
#ax.plot(theta_points[:-1],west_dig["hPa"].values[:])
#ax.errorbar(theta_points,west_dig["hPa"].values[:],west_stdv,color='green',marker='s')
ax.legend(bbox_to_anchor=(0.555, 0.82, 0.2, 0.2),ncol=2)

plt.tight_layout()
fig1.savefig(figo,dpi=300)
fig1.savefig(figo2,dpi=300)


# ---------------------------------------------
# Plot the main climbing/height/o2 comparison
# ---------------------------------------------

fig, ax = plt.subplots()
fig.set_size_inches(6,3.5)

# Twin the x-axis twice to make independent y-axes.
axes = [ax, ax.twinx(), ax.twinx()]
axes[-1].spines['right'].set_position(('axes', 1.19))

# Get monthly stats -- mean, max, min
# (make sure --> 2019 only)
p=o2.loc[o2.index.year<2020]
mon_means=p.groupby(p.index.month).mean()

# Get may mean
may_mean_p=mon_means.loc[mon_means.index==5]

# Get the mean may gradient in log P
may_mean_grad=grad.loc[grad.index.month==5].mean()

# Get the mean may o2 availability
may_mean_vo2=p2o(p.loc[p.index.month==5]).mean()

# Convert all the pressures to dz/zs under May-like conditions
dzs,zs=dp2dz(may_mean_p.values[0],p,may_mean_grad.values[0],8850.)

### From review
# Change in elevation from 2C warming
delta=2.57*2.0
p_min=p.min()
p_new=p.min()+delta
zcurrent=dp2dz(may_mean_p.values[0],p_min,may_mean_grad.values[0],8850.)
znew=dp2dz(may_mean_p.values[0],p_new,may_mean_grad.values[0],8850.)

# Again under climbing-only conditions
dzs_climb,zs_climb=dp2dz(may_mean_p.values[0],o2_climb.values[:],\
                         may_mean_grad.values[0],8850.)
zs_climb=pd.Series(zs_climb,index=o2_climb.index)
assert (zs_climb.index == no_meta.index).all()
no_meta["Perceived Elevation (m)"]=zs_climb.values[:]
assert (o2_climb.index == no_meta.index).all()
no_meta["Air Pressure (hPa)"]=o2_climb.values[:]

# Biggest diff?
max_diff_z=dzs.max()-dzs.min()
print("Max diff in apparent z is: %.2f" % max_diff_z.values[0])
# Biggest diff in o2?
vo2_all=p2o(p)
vo2_climb=p2o(o2_climb)
max_diff_o=100-vo2_all.min()/vo2_all.max()*100.
print("Max reduction in  O2 is [100-olow/ohigh*100]is: %.2f%%" % max_diff_o.values[0])
print("Reduction in O2 non-climb/climb: %.2f%%"%\
      (100-vo2_all.min()/vo2_climb.min()*100.))
print("Max O2 climbs [(ohigh-maymean)/maymean] is: %.2f%%" % \
      ((vo2_climb.max()-may_mean_vo2)/may_mean_vo2*100))
print("Min O2 climbs [(olow-maymean)/maymean] is: %.2f%%" % \
      ((vo2_climb.min()-may_mean_vo2)/may_mean_vo2*100))
print("Max O2 all [(ohigh-maymean)/maymean] is: %.2f%%" % \
      ((vo2_all.max()-may_mean_vo2)/may_mean_vo2*100))
print("Min O2 all [(olow-maymean)/maymean] is: %.2f%%" % \
      ((vo2_all.min()-may_mean_vo2)/may_mean_vo2*100))
# Remind when we had min/max pressures (irrespective of climbs)
print("Rememeber -- min press was: ", p["p"].loc[p["p"]==p["p"].min()])
print("Rememeber -- max press was: " ,p["p"].loc[p["p"]==p["p"].max()])


# Take care of ordering and plotting
zs=pd.Series(np.squeeze(zs),index=p.index)
ref=np.arange(1,13)
mudz=zs.groupby(zs.index.month).mean()
order=np.argsort(mudz).values[:][::-1]
x=ref[order]
pseq=[dzs.loc[dzs.index.month == i].values[:] for i in x]
medianprops = dict(linestyle='-', linewidth=0)
labs=np.array(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct",\
               "Nov","Dec"])
axes[0].boxplot(pseq,whis=10,meanline=True,showmeans=True,medianprops=medianprops,\
           labels=labs[order])


# Sort meta on perceived elevation and write out
no_meta=no_meta.sort_values(by="Perceived Elevation (m)",ascending=False)
# Create rank and copy key info across 
rank=[]; yr_out=[]; mon_out=[]; day_out=[]; hr_out=[]
rank_i=1
for i in range(len(no_meta["Perceived Elevation (m)"])):
    yr_out.append(no_meta.index[i].year)
    mon_out.append(no_meta.index[i].month)
    day_out.append(no_meta.index[i].day)
    hr_out.append(no_meta.index[i].hour)
    if i >0:
        if no_meta["Perceived Elevation (m)"].values[i]<\
        no_meta["Perceived Elevation (m)"].values[i-1]: rank_i+=1
    rank.append(rank_i)
no_meta["Rank"]=rank; no_meta["Year"]=yr_out; no_meta["Month"]=mon_out
no_meta["Day"]=day_out; no_meta["Hour"]=hr_out

# Output for SI (no names)
no_meta_out=no_meta[["Rank","name","Year","Month","Day","Hour","Air Pressure (hPa)",\
                     "Perceived Elevation (m)"]]
no_meta_out.to_csv(si)
# Output for SI (names)
no_meta_out_full=no_meta_out[["Rank","name","Year","Month","Day","Hour","Air Pressure (hPa)",\
                     "Perceived Elevation (m)"]]
no_meta_out_full["dz"]=no_meta_out["Perceived Elevation (m)"]-8850.
no_meta_out.to_csv(si_full)

# Now plot the dzs for the climbers
for i in range(len(order)):
    scratch=o2_climb.loc[o2_climb.index.month==(order[i]+1)]
    nc=len(scratch)
    print("Index = %.0f, Month = %.0f, nc = %.0f" % (order[i],order[i]+1,nc))
    if nc==0: 
        continue
    elif nc==1:
        
        ymu,dump=dp2dz(may_mean_p.values[0],scratch,\
                       may_mean_grad.values[0],8850.)
        axes[0].scatter(i+1,ymu,color="red",marker="o",s=25)

        
    elif nc==2:
        
         ymin,dump=dp2dz(may_mean_p.values[0],scratch.max(),\
                       may_mean_grad.values[0],8850.)    
    
         ymax,dump=dp2dz(may_mean_p.values[0],scratch.min(),\
                       may_mean_grad.values[0],8850.)  
         axes[0].plot([i+1.01,i+1.01],[ymax,ymin],color="red")

    else:
        
         ymin,dump=dp2dz(may_mean_p.values[0],scratch.max(),\
                       may_mean_grad.values[0],8850.)    
    
         ymu,dump=dp2dz(may_mean_p.values[0],scratch.mean(),\
                       may_mean_grad.values[0],8850.)    
    
         ymax,dump=dp2dz(may_mean_p.values[0],scratch.min(),\
                       may_mean_grad.values[0],8850.)        
        
         axes[0].scatter(i+1,ymu,color="red",marker="o",s=25)
         axes[0].plot([i+1.01,i+1.01],[ymax,ymin],color="red")

# # # Added from review -- compare the different methods for estimating pressure
zrefs=np.linspace(8500,10000,100)
dzrefs=zrefs-8850
me_p=dz2p(may_mean_grad.values[0],may_mean_p.values[0],dzrefs)
west_p=west93(zrefs)
fig,ax=plt.subplots(1,1)
ax.plot(zrefs,me_p,color="black")
ax.plot(zrefs,west_p,color="grey")
ax.grid()

# Now figure out the remaining axes
# Set the second y-axis to be abs height
yticks=[-200,-100,0,100,200,300,400,500,600]
axes[0].axhline(0,color='k')
axes[1].set_ylim(8850+yticks[0],8850+yticks[-1])
axes[1].set_yticks(np.array(yticks)+8850)

# Set the third y-axis to be %o2 -- use the labels of the main y-axis to 
# determine o2
yticks=[-200,-100,0,100,200,300,400,500,600]
axes[0].set_ylim(yticks[0],yticks[-1])
axes[0].set_yticks(yticks)
# This converts the ticks to absolute pressures
ypress=dz2p(may_mean_grad.values[:],may_mean_p.values[:],yticks)
# This converts those absolute pressures to VO2
yvo2=p2o(np.squeeze(ypress))
# And then to a standardized deviation (%) from the May mean
yvo2=(yvo2-may_mean_vo2.values[0])/may_mean_vo2.values[0]*100
axes[2].set_ylim(yvo2[0],yvo2[-1])
axes[2].set_yticks(yvo2)

# Formatting
axes[2].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axes[0].grid()
axes[0].set_ylabel("$\Delta$z (m)")
axes[1].set_ylabel("z (m)")
axes[2].set_ylabel("$\Delta$VO$_{2}$ max (%)")
axes[0].legend(loc=1)
plt.subplots_adjust(right=0.75)
fig.savefig(figo3,dpi=300)
# ---------------------------------------------
# End main plot
# ---------------------------------------------


# Some more summaries
print("Mean pressure during climbs = %.2f hPa" % o2_climb.mean())
deltas=np.zeros(len(o2_climb))
for i in range(len(o2_climb)):
    deltas[i]=o2_climb.values[i]-mus[doy==o2_climb.index.dayofyear[i]]
print("Mean anomaly during climbs = %.2f hPa" % np.mean(deltas))

# Max-min O2 difference? 
dp=o2_climb.max()-o2_climb.min()
#rato=p2o(o2_climb.max())/p2o(o2_climb.min())
rato=p2o(o2_climb.min())/p2o(o2_climb.max())
do=100-p2o(o2_climb.min())/p2o(o2_climb.max())*100. 


print("Delta pressure = %.2f hPa" % dp)
print("... reduction in o2 of %.2f%%" % do)
print("... ratio in pressure of %.3f%%" % rato)

# Dates of min/max?
print("Max O2 on:" , o2_climb.loc[o2_climb==o2_climb.max()].index)
print("Min O2 on:" , o2_climb.loc[o2_climb==o2_climb.min()].index)

# What % of Jan Days exceed the mean may pressure
frac_jan_may=np.sum(o2.loc[o2.index.month==1]>o2.loc[o2.index.month==5].mean())/\
np.float(np.sum(o2.index.month==1))*100.
jan_pctl=100-stats.percentileofscore(o2.loc[o2.index.month==5],\
                                 o2.loc[o2.index.month==1].max()[0])

# Pressure at time of Ang Rita?
ang=o2_climb.loc[o2_climb.index.month==12]
may=o2.loc[o2.index.month==5]
# higher than what fraction of May climbs?
fra_ang_may=100-np.sum(may>ang[0])/np.float(len(may))*100
print("%.2f%% of hours in May < Ang Rita" % fra_ang_may)

# Maximum excursion below mean?
max_ex=mus-mins; print("Max excursion below mean = %.2f on day %.02f" % \
                       (np.max(max_ex),doy[max_ex==np.max(max_ex)]))

# Create Table 1 
table=np.zeros((12,9))
mu_all=o2.loc[np.logical_and(o2.index.year>1978,o2.index.year<2020)].mean()["p"]
mu_climb=o2_climb.loc[np.logical_and(o2_climb.index.year>1978,\
                                     o2_climb.index.year<2020)].mean()
for m in range(1,13):
    table[m-1,0]=m
    idx_climb=np.logical_and(np.logical_and(\
    o2_climb.index.month==m,o2_climb.index.year>1978),o2_climb.index.year<2020)
    
    idx_all=np.logical_and(np.logical_and(\
    o2.index.month==m,o2.index.year>1978),o2.index.year<2020)
    
    n=np.float(np.sum(idx_all))
    
    # Quantile
    mu_mon_climb=o2_climb.loc[idx_climb].mean()
    mon_quant=stats.percentileofscore(o2[idx_all],mu_mon_climb)
    
    table[m-1,1]=np.sum(idx_climb)
    table[m-1,2]=table[m-1,1]/np.float(len(o2_climb))*100
    table[m-1,3]=o2.loc[idx_all].mean()
    table[m-1,4]=o2.loc[idx_all].max()-o2.loc[idx_all].min()
    table[m-1,5]=np.sum(o2.loc[idx_all]>mu_all)/n*100
    table[m-1,6]=np.sum(o2.loc[idx_all]>mu_climb)/n*100
    table[m-1,7]=mu_mon_climb
    table[m-1,8]=mon_quant
    

# For review -- illustrate the curvilinear relationship between Pio and VO2 max
ref=np.linspace(200,1000,100)
vo2=p2o(ref)
frac=0.2095 # Volume fraction of O2 in the atmopshere 
satvp=GF.satVpBolton(37.0)/100.
pio=frac*(ref-satvp) # partial pressure of oxygen in inspired air
fig,ax=plt.subplots(1,1)
ax.plot(pio,vo2,color='k')
ax.grid()
ax.set_xlabel("Partial pressure of inspired oxygen (hPa)")
ax.set_ylabel("VO$_{2}$max (ml kg$^{-1}$ min$^{-1}$)")
fig.savefig("/home/lunet/gytm3/Everest2019/Research/OneEarth/Figures/FigR1.tif",\
            dpi=300)