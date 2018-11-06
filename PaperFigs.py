### Import useful routines
import numpy as np
import string
import glob
import os
from collections import Counter
import scipy.stats as stats 
### Import CDAT routines ###
import MV2 as MV
import cdms2 as cdms
import genutil
import cdutil
import cdtime
from eofs.cdms import Eof
from seasonal_cycle_utils import mask_data
import Plotting 


### Import scipy routines for smoothing, interpolation
from scipy.interpolate import interp1d
from scipy.optimize import brentq,fminbound
import scipy.ndimage as ndimag

#From python-utils repo
import CMIP5_tools as cmip5
import DA_tools as da
import Plotting


### Import plotting routines
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap  
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib import mpl

### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

import pdsi_functions as b
import soilmoisture_functions as soil

## Data class: organizes individual drought atlases ###
class DATA():
    def __init__(self):
        self.NHDA = b.DroughtAtlas("ALL_OBS_plus_mex")
        self.NHDA.name="NHDA"
        self.NoMEX = b.DroughtAtlas("ALL_OBS")
        self.NADA = b.DroughtAtlas("NADA2.5",cutoff='1100-8-1')
        self.NADA.name="NADA"
        self.OWDA = b.DroughtAtlas("OWDA2.5",cutoff='1100-8-1')
        self.OWDA.name="OWDA"
        self.MADA = b.DroughtAtlas("MADA2.5")
        self.MADA.name="MADA"
        self.MXDA = b.DroughtAtlas("MXDA2.5")
        self.MXDA.name="MXDA"
        self.NoNADA=b.DroughtAtlas("ALL_NONADA")
        self.ALL = b.DroughtAtlas("ALL_ANZDA")
        self.ALL.name="GDA"
        self.ANZDA = b.DroughtAtlas("ANZDA2.5")
        self.ANZDA.name="ANZDA"
        

def colorregions(region):
    """Standardize colors for each region"""
    d={}
    d["ALL"]="k"
    d["NHDA"]=cm.gray(.5)
    d["NADA"]=cm.Purples(.5)
    d["OWDA"]=cm.Blues(.5)
    d["MXDA"]=cm.PiYG(.1)
    d["ANZDA"]=cm.PiYG(.8)
    d["MADA"]=cm.Oranges(.5)
    d["GDA"]="k"
    return d[region]

def get_dataset_color(dataset,depth=None):
    """Standardize colors for each dataset"""
    dataset=string.lower(dataset)
    d={}
    d["dai"]=cm.Blues(.5)
    d["tree"]=cm.summer(.3)
    d["cru"]=cm.Blues(.9)

    #models
    d["picontrol"]=cm.Purples(.8)
    d["h85"]=cm.Oranges(.8)
    d["tree_noise"]=cm.Greens(.8)
    
    #Soil moisture
    
    d["merra2"]={}
    d["merra2"]["30cm"]=cm.copper(.3)
    d["merra2"]["2m"]=cm.copper(.3)
    d["gleam"]={}
    d["gleam"]["30cm"]=cm.Reds(.3)
    d["gleam"]["2m"]=cm.Reds(.7)
    if depth is None:
        return d[dataset]
    else:
        return d[dataset][depth]
    
    
def sm_label(depth):
    if depth == "30cm":
        return "Surface"
    elif depth == "2m":
        return "Root Zone"
    else:
        raise TypeError("must be 30cm or 2m")

    
def SM_PDSI_histogram(D,start_time,stop_time,SM=None):
    """Make histograms of piControl and observed noise and forced trends"""
    if SM == None:
        SM=soil.SoilMoisture(D.ALL.obs[0].mask)

    plt.subplot(222)
    if stop_time.year<=2005:
        include_trees=True
    else:
        include_trees=False
    D.ALL.obs_SN(start_time,stop_time=stop_time,include_piControl=True,include_dai=True,include_cru=True,include_trees=include_trees)
    plt.legend(fontsize=8,ncol=3)
    plt.title("(b): D&A Results for PDSI")
    plt.subplot(221)
    pdsi_time_series(D,start_time,stop_time)
    plt.legend(fontsize=8)
    plt.title("(a): PDSI projections onto fingerprint")    
    plt.subplot(223)
    if start_time.year>=1981:
        SM.obs_SN(start_time,stop_time,"30cm",overlapping=True)
    else:
        SM.histograms(start_time,stop_time,"30cm",overlapping=True)
    plt.title("(c): 30cm soil moisture")
    plt.legend(fontsize=8)
    plt.subplot(224)
    if start_time.year>=1981:
        SM.obs_SN(start_time,stop_time,"2m",overlapping=True)
    else:
        SM.histograms(start_time,stop_time,"2m",overlapping=True)
    plt.title("(d): 2m soil moisture")
    plt.legend(fontsize=8)

    return SM


    

def pdsi_time_series(D,start_time,stop_time,SM=None,best_fit=True,aerosols=False,use_tree=True):
    if aerosols:
    
        aerosol_start = cdtime.comptime(1950,1,1)
        aerosol_stop = cdtime.comptime(1975,12,31)
        cru_p=D.ALL.project_cru_on_solver(start=start_time,solver=truncated_solver(D,aerosol_start,aerosol_stop,include_cru=True))
        dai_p=D.ALL.project_dai_on_solver(start=start_time,solver=truncated_solver(D,aerosol_start,aerosol_stop,include_dai=True))
        cdutil.setTimeBoundsYearly(dai_p)
        tree_p=D.ALL.get_tree_ring_projection(solver=truncated_solver(D,aerosol_start,aerosol_stop))
    else:    
        cru_p=D.ALL.project_cru_on_solver(start=start_time)
        dai_p=D.ALL.project_dai_on_solver(start=start_time)
        cdutil.setTimeBoundsYearly(dai_p)
        tree_p=D.ALL.projection
    x=np.arange(start_time.year,stop_time.year+1)
    if use_tree:
        cdutil.setTimeBoundsYearly(tree_p)
        ytree=tree_p(time=(start_time,stop_time,"oob")).anom()
        Plotting.time_plot(ytree,color=get_dataset_color("tree"),label="Tree ring reconstructions",lw=3)
        xtree=x[:len(ytree)]
        if best_fit:
            p=np.polyfit(xtree,ytree.asma(),1)
            plt.plot(xtree,np.polyval(p,xtree),"--",color=get_dataset_color("tree"),lw=1)

    cdutil.setTimeBoundsYearly(cru_p)
    if start_time.year>=1900:
        ycru=cru_p(time=(start_time,stop_time,"oob")).anom()
        Plotting.time_plot(ycru,color=get_dataset_color("cru"),label ="CRU",lw=3)
        xcru=x[:len(ycru)]
        if best_fit:
            p=np.polyfit(xcru,ycru.asma(),1)
            print "CRU slope is"+str(p[0]*10.)
            plt.plot(xcru,np.polyval(p,xcru),"--",color=get_dataset_color("cru"),lw=1)

   
   
        ydai=dai_p(time=(start_time,stop_time,"oob")).anom()
        Plotting.time_plot(ydai,color=get_dataset_color("dai"),label="Dai",lw=3)
        xdai=x[:len(ydai)]
        if best_fit:
            p=np.polyfit(xdai,ydai.asma(),1)
            print "DAI slope is"+str(p[0]*10.)
            plt.plot(xdai,np.polyval(p,xdai),"--",color=get_dataset_color("dai"),lw=1)
    
        #Plotting.time_plot(dai_p(time=(start_time,stop_time)).anom())
    plt.xlabel("Time")
    plt.ylabel("Projection")

    if SM is not None:
    
        SM.project_soilmoisture("MERRA2")
        SM.project_soilmoisture("GLEAM")
        for depth in ["30cm","2m"]:
            for dataset in ["GLEAM","MERRA2"]:
                ymoist=SM.OBS_PROJECTIONS[dataset][depth](time=(start_time,stop_time))
                xmoist=x[:len(ymoist)]
                Plotting.time_plot(ymoist,color=get_dataset_color(dataset,depth),lw=3,label=dataset+": "+sm_label(depth))
                if best_fit:
                    p=np.polyfit(xmoist,ymoist.asma(),1)
                    plt.plot(xmoist,np.polyval(p,xmoist),"--",color=get_dataset_color(dataset,depth),lw=1)

    if aerosols:
        return ytree,ycru,ydai
    
def plot_sn_data(data,i):
    noise=data["noise"].asma()
    models=data["modslopes"].asma()
    y=np.zeros_like(noise)+i
    plt.plot(noise,y,".",color=cm.Greens(.8))
    plt.plot(models,np.zeros_like(models)+i+.2,color=get_dataset_color("H85"))
    if "tree_rings" in data.keys():
        plt.plot([data["tree_rings"]],[i+.2],"o",color=get_dataset_color("tree"))
    if "cru" in data.keys():
        plt.plot([data["cru"]],[i+.2],"o",color=get_dataset_color("cru"))
    if "dai" in data.keys():
        plt.plot([data["dai"]],[i+.2],"o",color=get_dataset_color("dai"))
    if "picontrol" in data.keys():
        mnoise=data["picontrol"].asma()
        plt.plot(mnoise,np.zeros_like(mnoise)+i+0.05,".",color=get_dataset_color("picontrol"))
import seaborn as sns
import pandas as pd
def violinplot(D):
    gettime=lambda year: cdtime.comptime(year,1,1)
     
    getdata=lambda year: D.ALL.for_figure_4(gettime(year),gettime(year).add(30,cdtime.Years),include_piControl=True,include_dai=True,include_cru=True)
    myyears=[1945,1975]
    
    data1915=getdata(1915)
    ns=data1915["noise"].asma()
    ms=data1915["modslopes"].asma()
    
    labels=np.concatenate((np.repeat(["30 year pre-industrial trends"],len(ns)),np.repeat(["1915-1944"],len(ms))))
    slopes = np.concatenate((ns,ms))
    for year in myyears:
        data=getdata(year)
        ms=data["modslopes"].asma()
        lab=str(year)+"-"+str(year+29)
        labels=np.append(labels,np.repeat(lab,len(ms)))
        slopes = np.append(slopes,ms)
    
    d={"labels":labels,"slopes":slopes}
    df2=pd.DataFrame(d)
    sns.violinplot(y=df2["slopes"],x=df2["labels"])
    
    
#depths=["30cm","2m"]    
def soil_SN_figure(SM,depth,start_time=None):
    noise=[]
    signal_gleam=[]
    signal_merra2=[]
    
    if start_time is None:
        start_time=cdtime.comptime(1981,1,1)
    SM.project_soilmoisture("MERRA2")
    SM.project_soilmoisture("GLEAM")
    SM.project_piControl_on_solver(depth)
    SM.model_projections(depth)
    stop_time=cdtime.comptime(2017,12,31)
    nt=stop_time.year-start_time.year
    nmodel=SM.P[depth].shape[0]
    H85=np.ma.zeros((nmodel,nt))
    t=start_time.add(1,cdtime.Years)
    i=0
    while t.cmp(stop_time)<0:
        
        L=t.year-start_time.year+1
        modslopes,noiseterm = SM.sn_at_time(start_time,L,depth)
        H85[:,i] = modslopes
        noise += [np.std(noiseterm)]
        signal_gleam += [float(cmip5.get_linear_trends(SM.OBS_PROJECTIONS["GLEAM"][depth](time=(start_time,t))))]
        signal_merra2 += [float(cmip5.get_linear_trends(SM.OBS_PROJECTIONS["MERRA2"][depth](time=(start_time,t))))]
        t=t.add(1,cdtime.Years)
        i+=1
    timex=np.arange(start_time.year+1,start_time.year+1+nt)

    for i in range(nmodel):
        plt.plot(timex,H85[i]/np.array(noise),c="k",lw=1,alpha=.2)
    plt.plot(timex,np.array(signal_gleam)/np.array(noise),label="GLEAM",color=get_dataset_color("GLEAM",depth),lw=3)
    plt.plot(timex,np.array(signal_merra2)/np.array(noise),label="MERRA2",color=get_dataset_color("MERRA2",depth),lw=3)
    likely=stats.norm.interval(.66)[1]
    vlikely=stats.norm.interval(.90)[1]
    certain=stats.norm.interval(.99)[1]
    plt.axhline(likely,lw=3,color=cm.Reds(.33),label="Likely",alpha=.5)
    plt.axhline(vlikely,lw=3,color=cm.Reds(.66),label="Very Likely",alpha=.5)
    plt.axhline(certain,lw=3,color=cm.Reds(.99),label="Virtually Certain",alpha=.5)
    plt.xlabel("Last year of trend beginning in 1981")
    plt.ylabel("Signal to Noise Ratio")
   

    
    
    return noise,signal_gleam,signal_merra2,H85
        
        
                                   
def pdsi_SN_figure(D,start_time=None,stop_time=None,use_dai=True):
    noise_cru=[]
    noise_dai=[]
    noise_tree=[]
    noise=[]
    signal_cru=[]
    signal_dai=[]
    signal_tree=[]
    
    if start_time is None:
        start_time=cdtime.comptime(1981,1,1)
    if stop_time is None:
        stop_time=cdtime.comptime(2017,12,31)
    
    stop_cru=cdtime.comptime(2017,12,31)
    stop_dai=cdtime.comptime(2014,12,31)
    start_cru = cdtime.comptime(1901,1,1)
    start_dai=cdtime.comptime(1901,1,1)
    pcru=D.ALL.project_cru_on_solver(start=start_cru)
    pdai=D.ALL.project_dai_on_solver(start=start_dai)
    start_tree=cdtime.comptime(1400,1,1)
    stop_tree=cdtime.comptime(1975,12,31)
    nt=stop_cru.year-start_time.year
    nmodel=D.ALL.P.shape[0]
    H85=np.ma.zeros((nmodel,nt))
    t=start_time.add(1,cdtime.Years)
    i=0
    cru_time=[]
    tree_time=[]
    dai_time=[]
    while t.cmp(stop_time)<0:
        
        L=t.year-start_time.year+1
        modslopes,noiseterm = D.ALL.sn_at_time(start_time,L)
        H85[:,i] = modslopes
        noise+=[np.std(noiseterm)]
        if (t.cmp(stop_cru)<=0) and (t.cmp(start_cru)>0):
            signal_cru += [float(cmip5.get_linear_trends(pcru(time=(start_time,t))))]
            cru_time+=[t.year]
            noise_cru += [np.std(noiseterm)]
        if (t.cmp(stop_dai)<=0) and (t.cmp(start_dai)>0):
            signal_dai += [float(cmip5.get_linear_trends(pdai(time=(start_time,t))))]
            dai_time+=[t.year]
            noise_dai += [np.std(noiseterm)]
        if t.cmp(stop_tree)<=0:
            signal_tree += [float(cmip5.get_linear_trends(D.ALL.projection(time=(start_time,t))))]
            tree_time +=[t.year]
            noise_tree += [np.std(noiseterm)]
        t=t.add(1,cdtime.Years)
        i+=1
    timex=np.arange(start_time.year+1,start_time.year+1+nt)
    

    #for i in range(nmodel):
     #   plt.plot(timex,H85[i]/np.array(noise),c="k",lw=1,alpha=.2)
    plt.plot(cru_time,np.array(signal_cru)/np.array(noise_cru),label="CRU",color=get_dataset_color("cru"),lw=3)
    if use_dai:
        plt.plot(dai_time,np.array(signal_dai)/np.array(noise_dai),label="Dai",color=get_dataset_color("dai"),lw=3)
    plt.plot(tree_time,np.array(signal_tree)/np.array(noise_tree),label="Tree Rings",color=get_dataset_color("tree"),lw=3)
   # likely=stats.norm.interval(.66)[1]
   # vlikely=stats.norm.interval(.90)[1]
    #certain=stats.norm.interval(.99)[1]
    #plt.axhline(likely,lw=3,color=cm.Reds(.33),label="Likely",alpha=.5)
    #plt.axhline(vlikely,lw=3,color=cm.Reds(.66),label="Very Likely",alpha=.5)
    #plt.axhline(certain,lw=3,color=cm.Reds(.99),label="Virtually Certain",alpha=.5)

    #if start_time.year <=1960:
     #   plt.axhline(-1*likely,lw=3,color=cm.Reds(.33),alpha=.5)
      #  plt.axhline(-1*vlikely,lw=3,color=cm.Reds(.66),alpha=.5)
       # plt.axhline(-1*certain,lw=3,color=cm.Reds(.99),alpha=.5)
    #plt.axhline(0,c="k",ls=":",lw=1)

        
    #plt.xlabel("Last year of trend beginning in "+str(start_time.year))
    #plt.ylabel("Signal to Noise Ratio")
   

    
    
    #return noise,signal_gleam,signal_merra2,H85        
        
        
        
        
def projection_figure(D,fortalk=False):
    start = cdtime.comptime(1975,8,1)
    trends = cmip5.get_linear_trends(D.ALL.obs(time=(start,'2005-12-31')))
    if fortalk:
        plt.figure()
    else:
        plt.subplot(211)
    m=b.landplot(trends,vmin=-2,vmax=2)
    m.drawcoastlines(color="gray")
    plt.colorbar(orientation='horizontal',label="Trend (PDSI/decade)")
    plt.title("(a): 1975-2005 GDA trends")
    if fortalk:
        plt.figure()
    else:
        plt.subplot(212)
    Plotting.Plotting.time_plot(D.ALL.projection(time=('1100-1-1','1850-12-31')),color=cm.Greens(.9),lw=3,label="pre-industrial noise")
    Plotting.Plotting.time_plot(D.ALL.projection(time=('1851-1-1','1974-12-31')),c=cm.Greys(.5),lw=3)
    Plotting.Plotting.time_plot(D.ALL.projection(time=('1975-1-1','2005-12-31')),c=cm.Blues(.8),label="Target period",lw=3)
    plt.xlabel("Year")
    plt.ylabel("Projection")
    plt.title("(b): GDA projection on fingerprint")
    plt.xlim(1400,2010)
    plt.legend()
     
     
def average_vs_projection(data):
    """Show the difference between PDSI regional average and projections"""
    #start = cdtime.comptime(1975,8,1)
    #stop = cdtime.comptime(2005,12,31)
    start = cdtime.comptime(1900,1,1)
    stop = cdtime.comptime(1949,12,31)
    ax1=plt.subplot(211)
    ax2=plt.subplot(212)
    for thing in sorted(["MXDA","OWDA","NADA","MADA","ANZDA"])+["ALL"]:
        X=getattr(data,thing)
        pdsi_av = cdutil.averager(X.obs,axis='xy')(time=(start,stop))
        c=colorregions(X.name)
        Plotting.time_plot(pdsi_av,lw=3,color=c,label=X.name,ax=ax1)
        Plotting.time_plot(X.projection(time=(start,stop)),lw=3,color=c,label=X.name,ax=ax2)
        if thing =="ALL":
           
            x=np.arange(1975,2006)
            y=X.projection(time=(start,stop))
            p=np.polyfit(x,y,1)
            ax2.plot(x,np.polyval(p,x),"k--",lw=1)
            y=pdsi_av
            p=np.polyfit(x,y,1)
            ax1.plot(x,np.polyval(p,x),"k--",lw=1)
            
            
    ax1.set_title("(a): Regional mean PDSI")
    ax1.set_ylabel("PDSI")
    ax2.set_title("(b): Projection on fingerprint")
   # ax2.set_ylim(-22,22)
    ax2.set_ylabel("Projection")
    ax1.set_ylim(-1.6,1.6)
    plt.legend(ncol=4,fontsize=8)

def time_of_emergence_figure(data,noisestart="1400-8-1"):
    start = cdtime.comptime(1981,7,1)
    for thing in sorted(["MXDA","OWDA","NADA","MADA","ANZDA"])+["ALL"]:
        X=getattr(data,thing)
        c=colorregions(X.name)
        X.time_of_emergence(start,plot=True,color=c,noisestart=noisestart,shade=True)
def noisefigure(data,L=50,starttime="1400-8-1",showhist=True):
    for thing in sorted(["MXDA","OWDA","NADA","MADA","ANZDA"])+["ALL"]:
        
        X=getattr(data,thing)
        noiseterm = b.bootstrap_slopes(X.noise(time=(starttime,"2018-12-31")),L)
        if showhist:
            plt.hist(noiseterm,color=colorregions(X.name),normed=True,alpha=.4,histtype="bar",edgecolor="w")
        da.fit_normals_to_data(noiseterm,color=colorregions(X.name),label=X.name,lw=3)
        plt.xlabel(str(L)+"-year slope (PDSI/decade)")
        plt.ylabel("Normalized frequency")
        

def SignalToNoise(D,fortalk=False):
    if fortalk:
        plt.figure()
    else:
        plt.subplot(211)
    time_of_emergence_figure(D,noisestart=cmip5.start_time(D.ALL.obs))
    
    plt.title("(a): Time of emergence for PDSI signal") 
    plt.xlim(1985,2050)
    plt.legend(ncol=2)
    if fortalk:

        plt.figure()
        ax=plt.subplot(111)
    else:
        ax=plt.subplot(212)
    start_time = cdtime.comptime(1981,1,1)
    ALL_SM = soil.SoilMoisture(D.ALL.obs.mask[0])
    ALL_SM.time_of_emergence(start_time,"30cm",ax=ax,color=cm.Set1(1/2.))
    ALL_SM.time_of_emergence(start_time,"2m",ax=ax,color=cm.Set1(2/2.))
    plt.xlim(1985,2050)
    plt.title("(b): Times of emergence for soil moisture metrics")
    #noisefigure(D)
    plt.legend()
    plt.title("(b): Preindustrial \"noise\" terms")
def compare_pre_post_1100_noise(X,L=31,latbounds=None):
    time1=('1100-1-1','1399-12-31')
    c1=cm.Purples(.8)
    time2=('1400-1-1','2005-12-31')
    if latbounds is not None:
        obs=X.obs(latitude=latbounds)
        mma = MV.average(X.model(latitude=latbounds),axis=0)
        mma = mask_data(mma,obs[0].mask)
        solver = Eof(mma)
        obs = mask_data(obs,solver.eofs()[0].mask)
        truncnoise=solver.projectField(obs)[:,0]*da.get_orientation(solver)
        noisy1=truncnoise(time=time1)
        noisy2=truncnoise(time=time2)
    else:
        noisy1=X.noise(time=time1)
        noisy2=X.noise(time=time2)
    c2=cm.viridis(.1)
    plt.subplot(121)
    Plotting.Plotting.time_plot(noisy1,c=c1)
    Plotting.Plotting.time_plot(noisy2,c=c2)
    plt.ylabel("Projection")
    plt.title("(a): Noise time series")
    plt.subplot(122)
   
    plt.hist(b.bootstrap_slopes(noisy1,L),color=c1,normed=True,alpha=.5)
    da.fit_normals_to_data(b.bootstrap_slopes(noisy1,L),c=c1,label="1100-1400")
    plt.hist(b.bootstrap_slopes(noisy2,L),color=c2,normed=True,alpha=.5)
    da.fit_normals_to_data(b.bootstrap_slopes(noisy2,L),c=c2,label="1400-2005")
    plt.legend()
    plt.title("(b): 31-year trend distributions")
    return np.std(b.bootstrap_slopes(noisy1,L)),np.std(b.bootstrap_slopes(noisy2,L))
    
    
           


class ALL_SM():
    def __init__(self,D,regions):
        depths = ["30cm","2m"]
        start_time = cdtime.comptime(1975,8,1)
        ax1 = plt.subplot(211)
        ax2=plt.subplot(212)
        axes={}
        axes["30cm"]=ax1
        axes["2m"]=ax2
        for region in regions:
            ALL = getattr(D,region)
            mask = ALL.obs[0].mask
            SM = b.SoilMoisture(mask)
            
            for depth in depths:
                ax=axes[depth]
                SM.time_of_emergence(start_time,depth,ax=ax)
            setattr(self,region,SM)
            plt.savefig("JuneFigs/TOEsoilmoisture_"+region+".png")
            
    
        
def fingerprint_agreement_percentages(D,SM=None):
    pdsi_eof=da.get_orientation(D.ALL.solver)*D.ALL.solver.eofs()[0]
    mask = D.ALL.obs[0].mask
    if SM is None:
        SM = b.SoilMoisture(mask)
    SM30_eof=da.get_orientation(SM.solvers["30cm"] )*SM.solvers["30cm"].eofs()[0]
    SM2m_eof=da.get_orientation(SM.solvers["2m"] )*SM.solvers["2m"].eofs()[0]
    samesign=cmip5.cdms_clone(np.sign(pdsi_eof)*np.sign(SM30_eof),pdsi_eof)
    aw=cdutil.area_weights(pdsi_eof)
    test_area=np.ma.sum(MV.absolute(samesign)*aw)
    samesign_area=np.ma.sum(MV.masked_where(samesign<1,samesign)*aw)
    print "PDSI and 30cm have same sign in "+str(samesign_area/test_area*100)+"% of area"

    samesign=cmip5.cdms_clone(np.sign(pdsi_eof)*np.sign(SM2m_eof),pdsi_eof)
    samesign_area=np.ma.sum(MV.masked_where(samesign<1,samesign)*aw)
    print "PDSI and 2m have same sign in "+str(samesign_area/test_area*100)+"% of area"

    samesign=cmip5.cdms_clone(np.sign(SM30_eof)*np.sign(SM2m_eof),pdsi_eof)
    samesign_area=np.ma.sum(MV.masked_where(samesign<1,samesign)*aw)
    print "30cm and 2m have same sign in "+str(samesign_area/test_area*100)+"% of area"
    

    
    
#MAINTEXT FIGURES
def Figure1(D,fortalk=False):
    mask = D.ALL.obs[0].mask
    soil.soilmoisture_fingerprints(mask,fortalk=fortalk)
    ax=plt.gcf().axes[-1]
    ax.set_xlim(1900,2050)
    ax.set_ylim(-0.4,0.4)

def Figure2(D,fortalk=False):
    if fortalk:
        c="w"
        
    else:
        c="k"
    plt.subplot(121)
    Plotting.time_plot(D.ALL.get_noise(),color=cm.Greens(.8),label="Pre-industrial Noise",lw=1)
    Plotting.time_plot(D.ALL.projection(time=('1850-1-1','1975-12-31')),c=c,lw=1)
    plt.ylabel("Projection (temporal amplitude)")
    plt.xlabel("Year")
    plt.legend(fontsize=8)
    plt.title("(a): Global Drought Atlas Projection Onto Fingerprint")
    plt.subplot(122)
    t1900 = cdtime.comptime(1900,7,1)
    times = np.arange(1,201)
    

    likely=stats.norm.interval(.66)[1]
    vlikely=stats.norm.interval(.90)[1]
    certain=stats.norm.interval(.99)[1]
    plt.axhline(likely,lw=1,color=cm.Reds(.33),label="Likely",alpha=.5)
    plt.axhline(vlikely,lw=1,color=cm.Reds(.66),label="Very Likely",alpha=.5)
    plt.axhline(certain,lw=1,color=cm.Reds(.99),label="Virtually Certain",alpha=.5)
    D.ALL.time_of_emergence(t1900,times=times,color="k",lw=3)
   
    plt.title("(b): Model-predicted time of emergence")

   # plt.subplot(133)
    pdsi_SN_figure(D,cdtime.comptime(1900,1,1),use_dai=True)
    plt.legend(loc=0,ncol=2,fontsize=8)


def old_Figure3(D,SM=None):
    start_time = cdtime.comptime(1901,1,1)
    stop_time = cdtime.comptime(1949,12,31)
    SM_PDSI_histogram(D,start_time,stop_time,SM=SM)
def Figure3(D,SM):
    recent_start = cdtime.comptime(1981,1,1)
    recent_stop=cdtime.comptime(2017,12,31)
    #historical_stop = cdtime.comptime(2005,12,31)

    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)

    hist_start = cdtime.comptime(1900,1,1)
    hist_stop = cdtime.comptime(1949,12,31)

    plt.subplot(421)
    pdsi_time_series(D,hist_start,hist_stop,best_fit=True)
    plt.title("(a): 1900-1949 projections")
    L=plt.legend(loc=0,ncol=3,fontsize=8)
    
    L.set_frame_on(False)

    plt.subplot(422)
    D.ALL.obs_SN(hist_start,stop_time=hist_stop,include_piControl=True,include_dai=True,include_cru=True,include_trees=True)
    plt.title("(b): D&A Results: 1900-1949")
    L=plt.legend(loc=0,ncol=1,fontsize=8)
    L.get_texts()[-1].set_text("50-year trends in PiControl simulations")
    L.set_frame_on(False)

    plt.subplot(423)
    pdsi_time_series(D,aerosol_start,aerosol_stop,best_fit=True)
    plt.title("(c): 1950-1975 projections")
    L=plt.legend(loc=0,ncol=3,fontsize=8)
    L.set_frame_on(False)

    plt.subplot(424)
    D.ALL.obs_SN(aerosol_start,stop_time=aerosol_stop,include_piControl=True,include_dai=True,include_cru=True,include_trees=True)
    plt.title("(d): D&A Results: 1950-1975")
    L=plt.legend(loc=0,ncol=1,fontsize=8)
    L.get_texts()[-1].set_text("26-year trends in PiControl simulations")
    L.set_frame_on(False)
    #plt.legend().set_visible(False)  

    plt.subplot(425)
    pdsi_time_series(D,recent_start,recent_stop,best_fit=False,SM=SM,use_tree=False)
    plt.title("(e): 1981-2005 projections")
    L=plt.legend(loc=0,ncol=3,fontsize=8)
    L.set_frame_on(False)
    
    plt.subplot(426)
    D.ALL.obs_SN(recent_start,stop_time=recent_stop,include_piControl=True,include_dai=False,include_cru=True,include_trees=False)
    plt.title("(f): D&A Results: 1981-2017 ")
    L=plt.legend(loc=0,ncol=1,fontsize=8)
    L.get_texts()[-1].set_text("37-year trends in PiControl simulations")
    L.set_frame_on(False)
    #plt.legend().set_visible(False)

    ax=plt.subplot(427)
    recent_start = cdtime.comptime(1981,1,1)
   
    recent_stop = cdtime.comptime(2017,12,31)
    #plt.subplot(121)
    SM.obs_SN(recent_start,recent_stop,"30cm")
    L=ax.legend(fontsize=8,loc=2)
    L.get_texts()[0].set_text("37 year trends in PiControl simulations")
    L.get_texts()[1].set_text("1981-2017 model projections")
    L.set_frame_on(False)
    plt.title("(g): D&A Results: 1981-2017 (surface soil moisture)")
    plt.subplot(428)
    SM.obs_SN(recent_start,recent_stop,"2m")
    L=ax.legend(fontsize=8,loc=2)
    L.get_texts()[0].set_text("37 year trends in PiControl simulations")
    L.get_texts()[1].set_text("1981-2017 model projections")
    plt.title("(h): D&A Results: 1981-2017 (root zone soil moisture)")
    L.set_frame_on(False)

    
    

    
def old_Figure4(D,SM):
    axes=[]
    axes+=[plt.subplot(221)]
    t1981=cdtime.comptime(1981,1,1)
    t2017=cdtime.comptime(2017,12,31)
    pdsi_time_series(D,t1981,t2017,SM=SM,best_fit=False)
    
    axes+=[plt.subplot(222)]
    pdsi_SN_figure(D)
    plt.title("(b): PDSI")
    plt.legend(fontsize=8)
    axes+=[plt.subplot(223)]
    noise,signal_gleam,signal_merra2,H85=soil_SN_figure(SM,"30cm")
    plt.legend(fontsize=8)
    plt.title("(c): Surface soil moisture")
    axes+=[plt.subplot(224)]
    noise2,signal_gleam2,signal_merra22,H852=soil_SN_figure(SM,"2m")
    plt.title("(d): Root zone soil moisture")
    plt.legend(fontsize=8)

    h,l=axes[1].get_legend_handles_labels()
    axes[1].legend(h[:3],l[:3],fontsize=8)

    h,l=axes[2].get_legend_handles_labels()
    axes[2].legend(h[:2],l[:2],fontsize=8)

    h,l=axes[3].get_legend_handles_labels()
    axes[3].legend(h[:2],l[:2],fontsize=8)
    axes[0].legend(ncol=2,fontsize=8)
    

    
   
#SUPPORTING FIGURES
# def S_SM_DA(SM):
#     recent_start = cdtime.comptime(1981,1,1)
#     # recent_stop = cdtime.comptime(2005,12,31)
#     # plt.subplot(221)
#     # SM.obs_SN(recent_start,recent_stop,"30cm")
#     # plt.title("(a): Surface")
#     # plt.subplot(222)
#     # SM.obs_SN(recent_start,recent_stop,"2m")
#     # plt.title("(b): Root Zone")

#     recent_stop = cdtime.comptime(2017,12,31)
#     plt.subplot(121)
#     SM.obs_SN(recent_start,recent_stop,"30cm")
#     plt.title("(a): Surface")
#     plt.subplot(122)
#     SM.obs_SN(recent_start,recent_stop,"2m")
#     plt.title("(b): Root Zone")
    
    
def Smodel_trends(D):
    m=b.landplot(cmip5.get_linear_trends(D.ALL.mma))
def Saerosol(D):
    pdsi_SN_figure(D,start_time=cdtime.comptime(1950,1,1))
   
    
def Spdsi_sm(D,SM):
    pdsi_time_series(D,t1981,t2005,SM=SM)
    plt.legend(ncol=3,fontsize=10)

def STOE(D):
    time_of_emergence_figure(D)
    plt.legend(ncol=3,fontsize=12)
    vlikely=stats.norm.interval(.9)[1]
    plt.axhline(vlikely,lw=3,color=cm.Reds(.66),label="Very Likely",alpha=.5)
    
    
def Snoise(D):
    noisefigure(D)
    plt.legend(ncol=3,fontsize=10)

def Savproj(D):
    average_vs_projection(D)

def Saerosol_fingerprint(D):
    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)
    aerosolsolver=Eof(D.ALL.mma(time=(aerosol_start,aerosol_stop)),weights='area')
    fac=da.get_orientation(aerosolsolver)
    plt.subplot(221)
    m=b.landplot(fac*aerosolsolver.eofs()[0],vmin=-.1,vmax=.1)
    varex= str(int(100*np.round(aerosolsolver.varianceFraction()[0],2)))
    plt.title("(a): 1950-1975 historical fingerprint ("+varex+"% of variance explained)")
    m.drawcoastlines(color='gray')
    plt.colorbar(orientation='horizontal',label='EOF loading')
    plt.subplot(222)
    Plotting.time_plot(fac*aerosolsolver.pcs()[:,0],color="k",lw=3)
    plt.title("(b): Associated PC")
    plt.ylabel("Temporal amplitude")

    plt.subplot(223)

    target_obs,cru_proj,dai_proj=pdsi_time_series(D,aerosol_start,aerosol_stop,aerosols=True)
    plt.title("(c): Projections on fingerprint")
    plt.subplot(224)

   # target_obs = D.ALL.get_tree_ring_projection(solver = aerosolsolver)(time=(aerosol_start,aerosol_stop))
    L=len(target_obs)
    modslopes,noiseterm = D.ALL.sn_at_time(aerosol_start,L,overlapping=True,solver=aerosolsolver)
    ns=np.std(noiseterm)
    signal = float(cmip5.get_linear_trends(target_obs))
    plt.hist(modslopes/ns,20,normed=True,color=get_dataset_color("h85"),alpha=.5)
    lab = str(aerosol_start.year)+"-"+str(aerosol_stop.year)
    da.fit_normals_to_data(modslopes/ns,color=get_dataset_color("h85"),lw=3,label=lab+" Model projections")

    plt.hist(noiseterm/ns,20,normed=True,color=get_dataset_color("tree_noise"),alpha=.5)
    da.fit_normals_to_data(noiseterm/ns,color=get_dataset_color("tree_noise"),lw=3,label="Pre-1850 tree-ring reconstructions")
    percentiles=[]
    plt.axvline(signal/ns,color=get_dataset_color("tree"),lw=3,label=lab+" Tree-ring reconstructions")
    
    noise_percentile=stats.percentileofscore(noiseterm.tolist(),signal)
    h85_percentile=stats.percentileofscore(modslopes.tolist(),signal)
    percentiles += [noise_percentile,h85_percentile]


    daitrend = cmip5.get_linear_trends(dai_proj)
    print "DAI slope is "+str(daitrend)
    daisignal = daitrend/ns
    
    plt.axvline(daisignal,color=get_dataset_color("dai"),lw=3,label=lab+" Dai dataset")
    print "DAI signal/noise is "+str(daisignal)

    
    
    crutrend = cmip5.get_linear_trends(cru_proj)
    print "CRU slope is "+str(crutrend)
    crusignal = crutrend/ns
    
    plt.axvline(crusignal,color=get_dataset_color("cru"),lw=3,label=lab+" Cru dataset")
    print "CRU signal/noise is "+str(crusignal)

   
            
       
    plt.legend(loc=0)
    plt.xlabel("S/N")
    plt.ylabel("Normalized Frequency")
    plt.title("(d): Detection and Attribution Results")
    
       
        
        
def for_figure_4(self,start_time,stop_time=None,overlapping=True,include_trees=True,include_dai=False,include_cru=False,include_piControl=False,noisestart=None,solver=None):
        data={}
        if stop_time is None:
            stop_time=cmip5.stop_time(self.get_tree_ring_projection())
        target_obs = self.get_tree_ring_projection()(time=(start_time,stop_time))
        L=len(target_obs)
        modslopes,noiseterm = self.sn_at_time(start_time,L,overlapping=True,noisestart=noisestart,solver=solver)
        ns=np.std(noiseterm)
        signal = float(cmip5.get_linear_trends(target_obs))
        data["noise"]=noiseterm
        data["modslopes"]=modslopes
        data["tree_rings"]=signal
        
        if include_dai:
            dai_proj = self.project_dai_on_solver(start=start_time,solver=solver)
            daitrend = cmip5.get_linear_trends(dai_proj(time=(start_time,stop_time)))
            data["dai"]=daitrend
           
        if include_cru:
            cru_proj = self.project_cru_on_solver(start=start_time,solver=solver)
            crutrend = cmip5.get_linear_trends(cru_proj(time=(start_time,stop_time)))
            data["cru"]=crutrend
            
        if include_piControl:
            p=self.project_piControl_on_solver(solver=solver)
            noiseterm_mod=bootstrap_slopes(p,L)
            data["picontrol"]=noiseterm_mod
    
    
#supporting table
    
def sensitivity_tests(D):
    t1949=cdtime.comptime(1949,12,31)
    t1900=cdtime.comptime(1900,1,1)
    for thing in sorted(["MXDA","OWDA","NADA","MADA","ANZDA"]):
        X = getattr(D,thing)
        data = X.for_figure_4(t1900,t1949)
        sn = str(np.round(data["tree_rings"]/np.std(data["noise"]),2))
        print X.name +" & "+sn
    X = getattr(D,"NoMEX")
    data = X.for_figure_4(t1900,t1949)
    sn = str(np.round(data["tree_rings"]/np.std(data["noise"]),2))
    print "NADA+OWDA+MADA (1100 start)"+" & "+sn
    X=getattr(D,"NoNADA")
    data = X.for_figure_4(t1900,t1949)
    sn = str(np.round(data["tree_rings"]/np.std(data["noise"]),2))
    print "No NADA" +" & "+sn
    X=getattr(D,"ALL")
    data = X.for_figure_4(t1900,t1949)
    sn = str(np.round(data["tree_rings"]/np.std(data["noise"]),2))
    print "GDA" +" & "+sn



def aerosol_solver(D,aerosol_start=None,aerosol_stop=None,include_cru=False):
    if aerosol_start is None:
        aerosol_start = cdtime.comptime(1950,1,1)
    if aerosol_stop is None:
        aerosol_stop=cdtime.comptime(1975,1,1)
    aerosoldata = D.ALL.model(time=(aerosol_start,aerosol_stop))
    aerosol_mma=MV.average(cmip5.ensemble2multimodel(aerosoldata),axis=0)
    aerosolsolver = Eof(aerosol_mma,weights="area")
    return aerosolsolver

def early_solver(D,early_start=None,early_stop=None):
    if early_start is None:
        early_start = cdtime.comptime(1900,1,1)
    if early_stop is None:
        early_stop=cdtime.comptime(1949,12,31)
    earlydata = D.ALL.model(time=(early_start,early_stop))
    early_mma=MV.average(cmip5.ensemble2multimodel(earlydata),axis=0)
    earlysolver = Eof(early_mma,weights="area")
    return earlysolver
def rcp_solver(D,early_start=None,early_stop=None):
    if early_start is None:
        early_start = cdtime.comptime(2006,1,1)
    if early_stop is None:
        early_stop=cdtime.comptime(2099,12,31)
    earlydata = D.ALL.model(time=(early_start,early_stop))
    early_mma=MV.average(cmip5.ensemble2multimodel(earlydata),axis=0)
    earlysolver = Eof(early_mma,weights="area")
    return earlysolver

def truncated_solver(D,the_start=None,the_stop=None,include_cru=False,include_dai=False):
   
    thedata = D.ALL.model(time=(the_start,the_stop))
    the_mma=MV.average(cmip5.ensemble2multimodel(thedata),axis=0)
    if include_cru:
        f = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/CRU_selfcalibrated.nc")
        cru_jja=f("pdsi")
        f.close()
        cru_jja_mask = mask_data(cru_jja,D.ALL.obs[0].mask)(time=(the_start,'2018-12-31'))
        newmask = np.prod(~cru_jja_mask.mask,axis=0)
        cru_jja_mask = mask_data(cru_jja_mask,newmask==0)
        
        thesolver = Eof(mask_data(D.ALL.mma(time=(the_start,the_stop)),newmask==0),weights='area')
    elif include_dai:
        f = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/DAI_selfcalibrated.nc")
        dai_jja=f("pdsi")
        f.close()
        dai_jja_mask = mask_data(dai_jja,D.ALL.obs[0].mask)(time=(the_start,'2018-12-31'))
        newmask = np.prod(~dai_jja_mask.mask,axis=0)
        dai_jja_mask = mask_data(dai_jja_mask,newmask==0)
        
        thesolver = Eof(mask_data(D.ALL.mma(time=(the_start,the_stop)),newmask==0),weights='area')
    else:
        thesolver = Eof(the_mma,weights="area")
    return thesolver

import csv
def forTable(D):
    
    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)

    hist_start = cdtime.comptime(1900,1,1)
    hist_stop = cdtime.comptime(1949,12,31)

    with open("../DroughtAtlasPaper/FIGS/ForPaper/SI/Table1.csv","w") as fw:
        csvwriter=csv.writer(fw,delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        d={}
        d
        for X in ["ANZDA","MADA","MXDA","NADA","OWDA","NHDA","NoNADA","NoMEX"]:
            towrite=getattr(D,X).obs_SN(hist_start,stop_time=hist_stop,plot=False)+getattr(D,X).obs_SN(aerosol_start,stop_time=aerosol_stop,plot=False)
            towrite2 = [X]+(map(str,towrite))
            csvwriter.writerow(towrite2)
   
        
        
