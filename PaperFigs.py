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


    

class DATA():
    """ Data and methods for regional and global drought atlases"""
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
    """Set the colors to ensure a uniform scheme for each region """
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
    """Set the colors to ensure a uniform scheme for each dataset """
    dataset=string.lower(dataset)
    d={}
    d["dai"]=cm.Blues(.5)
    d["tree"]=cm.summer(.3)
    d["cru"]=cm.Blues(.9)

    #models
    d["picontrol"]=cm.Purples(.8)
    d["h85"]="k"
    d["tree_noise"]=cm.PiYG(.2)
    
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


   
def pdsi_time_series(D,start_time,stop_time,SM=None,best_fit=True,aerosols=False,use_tree=True):
    """ plot time series of PDSI data between start_time and stop_time"""
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
        Plotting.time_plot(ytree,color=get_dataset_color("tree"),label="Tree ring reconstructions",lw=1)
        xtree=x[:len(ytree)]
        if best_fit:
            p=np.polyfit(xtree,ytree.asma(),1)
            plt.plot(xtree,np.polyval(p,xtree),"--",color=get_dataset_color("tree"),lw=1)

    cdutil.setTimeBoundsYearly(cru_p)
    if start_time.year>=1900:
        ycru=cru_p(time=(start_time,stop_time,"oob")).anom()
        Plotting.time_plot(ycru,color=get_dataset_color("cru"),label ="CRU",lw=1)
        xcru=x[:len(ycru)]
        if best_fit:
            p=np.polyfit(xcru,ycru.asma(),1)
            print "CRU slope is"+str(p[0]*10.)
            plt.plot(xcru,np.polyval(p,xcru),"--",color=get_dataset_color("cru"),lw=1)

   
   
        ydai=dai_p(time=(start_time,stop_time,"oob")).anom()
        Plotting.time_plot(ydai,color=get_dataset_color("dai"),label="Dai",lw=1)
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
                Plotting.time_plot(ymoist,color=get_dataset_color(dataset,depth),lw=1,label=dataset+": "+sm_label(depth))
                if best_fit:
                    p=np.polyfit(xmoist,ymoist.asma(),1)
                    plt.plot(xmoist,np.polyval(p,xmoist),"--",color=get_dataset_color(dataset,depth),lw=1)

    if aerosols:
        return ytree,ycru,ydai
    
   

                                   
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
def NatureRevisionsFigure2(D,fortalk=False):
    mask = D.ALL.obs[0].mask
    soil.soilmoisture_fingerprints(mask,fortalk=fortalk)
    ax=plt.gcf().axes[-1]
    ax.set_xlim(1900,2050)
    ax.set_ylim(-0.4,0.4)
    if not fortalk:
    #cosmetic changes for Nature
        fig=plt.gcf()
        fig.axes[0].set_title("(a)",fontsize=6)
        fig.axes[2].set_title("(b)",fontsize=6)
        fig.axes[4].set_title("(c)",fontsize=6)
        fig.axes[6].set_title("(d)",fontsize=6)
        pcax=fig.axes[-1]
        pcax.yaxis.set_label_position('right')
        pcax.yaxis.set_ticks_position('right')
        plt.setp(pcax.yaxis.get_label(),fontsize=6)
        plt.setp(pcax.xaxis.get_label(),fontsize=6)
        plt.setp(pcax.get_yticklabels(),fontsize=6)
        plt.setp(pcax.get_xticklabels(),fontsize=6)
        cbaxes=fig.axes[1::2]
        for cax in cbaxes:
            ticklabels=["-0.1","","-0.05","","0","","0.05","","0.1"]
            cax.set_yticklabels(ticklabels)
            plt.setp(cax.yaxis.get_ticklabels(),fontsize=6)
            plt.setp(cax.yaxis.get_label(),fontsize=6)
   
  
def NatureRevisions_Figure4(D,SM):

    recent_start = cdtime.comptime(1981,1,1)
    recent_stop=cdtime.comptime(2017,12,31)
    #historical_stop = cdtime.comptime(2005,12,31)

    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)

    hist_start = cdtime.comptime(1900,1,1)
    hist_stop = cdtime.comptime(1949,12,31)

    plt.subplot(421)
    pdsi_time_series(D,hist_start,hist_stop,best_fit=True)
    plt.title("(a)",fontsize=6)
    L=plt.legend(loc=0,ncol=3,fontsize=6)
    
    L.set_frame_on(False)

    plt.subplot(422)
    D.ALL.obs_SN(hist_start,stop_time=hist_stop,include_piControl=True,include_dai=True,include_cru=True,include_trees=True)
    plt.title("(b)",fontsize=6)
    L=plt.legend(loc=0,ncol=1,fontsize=6)
    #L.get_texts()[-1].set_text("50-year trends in PiControl simulations")
    L.set_frame_on(False)

    plt.subplot(423)
    pdsi_time_series(D,aerosol_start,aerosol_stop,best_fit=True)
    plt.title("(c)",fontsize=6)
    L=plt.legend(loc=0,ncol=3,fontsize=6)
    L.set_frame_on(False)

    plt.subplot(424)
    D.ALL.obs_SN(aerosol_start,stop_time=aerosol_stop,include_piControl=True,include_dai=True,include_cru=True,include_trees=True)
    plt.title("(d)",fontsize=6)
    L=plt.legend(loc=0,ncol=1,fontsize=6)
    #L.get_texts()[-1].set_text("26-year trends in PiControl simulations")
    L.set_frame_on(False)
    #plt.legend().set_visible(False)  

    plt.subplot(425)
    pdsi_time_series(D,recent_start,recent_stop,best_fit=False,SM=SM,use_tree=False)
    plt.title("(e)",fontsize=6)
    L=plt.legend(loc=0,ncol=3,fontsize=6)
    L.set_frame_on(False)
    
    plt.subplot(426)
    D.ALL.obs_SN(recent_start,stop_time=recent_stop,include_piControl=True,include_dai=False,include_cru=True,include_trees=False)
    plt.title("(f)",fontsize=6)
    L=plt.legend(loc=0,ncol=1,fontsize=6)
    #L.get_texts()[-1].set_text("37-year trends in PiControl simulations")
    L.set_frame_on(False)
    #plt.legend().set_visible(False)

    ax=plt.subplot(427)
    recent_start = cdtime.comptime(1981,1,1)
   
    recent_stop = cdtime.comptime(2017,12,31)
    #plt.subplot(121)
    SM.obs_SN(recent_start,recent_stop,"30cm")
    L=ax.legend(fontsize=6,loc=2)
    L.get_texts()[0].set_text("37 year trends in PiControl simulations")
    L.get_texts()[1].set_text("1981-2017 model projections")
    L.set_frame_on(False)
    plt.title("(g)",fontsize=6)
    plt.subplot(428)
    SM.obs_SN(recent_start,recent_stop,"2m")
    L=ax.legend(fontsize=8,loc=2)
    L.get_texts()[0].set_text("37 year trends in PiControl simulations")
    L.get_texts()[1].set_text("1981-2017 model projections")
    plt.title("(h)",fontsize=6)
    L.set_frame_on(False)
    for ax in plt.gcf().axes:
        plt.setp(ax.xaxis.get_label(),fontsize=6)
        plt.setp(ax.yaxis.get_label(),fontsize=6)
        plt.setp(ax.get_xticklabels(),fontsize=6)
        plt.setp(ax.get_yticklabels(),fontsize=6)
    

    
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
    m.fillcontinents(color="gray",zorder=0)
    plt.colorbar(orientation="horizontal",label="1900-2099 trend (PDSI/decade)")
    
    
    plt.ylim(-60,90)
   
    
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

def NatureRevisions_Figure5(D):
    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)
    aerosolsolver=Eof(D.ALL.mma(time=(aerosol_start,aerosol_stop)),weights='area')
    fac=da.get_orientation(aerosolsolver)
    plt.subplot(221)
    m=b.landplot(fac*aerosolsolver.eofs()[0],vmin=-.1,vmax=.1)
    m.fillcontinents(color="gray",zorder=0)
    
    varex= str(int(100*np.round(aerosolsolver.varianceFraction()[0],2)))
    plt.title("(a)")#: 1950-1975 historical fingerprint ("+varex+"% of variance explained)",fontsize=8)
    m.drawcoastlines(color='gray')
    plt.ylim(-60,90)
    plt.colorbar(orientation='horizontal',label='EOF loading')
    plt.subplot(222)
    Plotting.time_plot(fac*aerosolsolver.pcs()[:,0],color=cm.Greys(.8),lw=1)
    plt.title("(b)")#: Associated PC",fontsize=8)
    plt.ylabel("Temporal amplitude")

    plt.subplot(223)

    target_obs,cru_proj,dai_proj=pdsi_time_series(D,aerosol_start,aerosol_stop,aerosols=True)
    plt.legend(fontsize=6)
    plt.title("(c)")#: Projections on fingerprint",fontsize=8)
    plt.subplot(224)

   # target_obs = D.ALL.get_tree_ring_projection(solver = aerosolsolver)(time=(aerosol_start,aerosol_stop))
    L=len(target_obs)
    modslopes,noiseterm = D.ALL.sn_at_time(aerosol_start,L,overlapping=True,solver=aerosolsolver)
    ns=np.std(noiseterm)
    signal = float(cmip5.get_linear_trends(target_obs))
    plt.hist(modslopes/ns,20,normed=True,color=get_dataset_color("h85"),alpha=.5)
    lab = str(aerosol_start.year)+"-"+str(aerosol_stop.year)
    da.fit_normals_to_data(modslopes/ns,color=get_dataset_color("h85"),lw=1,label="H85")

    plt.hist(noiseterm/ns,20,normed=True,color=get_dataset_color("tree_noise"),alpha=.5)
    da.fit_normals_to_data(noiseterm/ns,color=get_dataset_color("tree_noise"),lw=1,label="Pre-1850 tree rings")
    percentiles=[]
    plt.axvline(signal/ns,color=get_dataset_color("tree"),lw=1,label=lab+" GDA trend")
    
    noise_percentile=stats.percentileofscore(noiseterm.tolist(),signal)
    h85_percentile=stats.percentileofscore(modslopes.tolist(),signal)
    percentiles += [noise_percentile,h85_percentile]


    daitrend = cmip5.get_linear_trends(dai_proj)
    print "DAI slope is "+str(daitrend)
    daisignal = daitrend/ns
    
    plt.axvline(daisignal,color=get_dataset_color("dai"),lw=1,label="Dai")
    print "DAI signal/noise is "+str(daisignal)

    
    
    crutrend = cmip5.get_linear_trends(cru_proj)
    print "CRU slope is "+str(crutrend)
    crusignal = crutrend/ns
    
    plt.axvline(crusignal,color=get_dataset_color("cru"),lw=1,label="CRU")
    print "CRU signal/noise is "+str(crusignal)

   
            
       
    plt.legend(loc=0,fontsize=8)
    plt.xlabel("S/N")
    plt.ylabel("Normalized Frequency")
    plt.title("(d)")#: Detection and Attribution Results",fontsize=8)
    fig=plt.gcf()
    for ax in fig.axes:
        plt.setp(ax.xaxis.get_label(),fontsize=6)
        plt.setp(ax.yaxis.get_label(),fontsize=6)
        plt.setp(ax.get_xticklabels(),fontsize=6)
        plt.setp(ax.get_yticklabels(),fontsize=6)
    ax=fig.axes[0]
    ax.set_title("(a)",fontsize=6)
    ax=fig.axes[2]
    ax.set_title("(b)",fontsize=6)
    ax=fig.axes[3]
    ax.set_title("(c)",fontsize=6)
    ax=fig.axes[4]
    ax.set_title("(d)",fontsize=6)
    leg=ax.legend(fontsize=6,ncol=1,loc=2)
    leg.set_frame_on(False)
    cax=fig.axes[1]
    ticklabels=["-0.1","","-0.05","","0","","0.05","","0.1"]
    cax.set_xticklabels(ticklabels)
    plt.setp(cax.xaxis.get_ticklabels(),fontsize=6)
    plt.setp(cax.xaxis.get_label(),fontsize=6)
       
        
        
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
   

#nature figure sizes (cm)
onecol=8.9
oneandahalfcol=12.0
twocol=18.3
height=24.7
def cm2inch(value):
    return value/2.54
def naturefig(colsize=1,heightfrac=1):
    if colsize==1:
        col=8.9
    elif colsize==1.5:
        col=12.0
    elif colsize==2:
        col=18.3
    else:
        print "col must be one of 1,1.5,2"
        
        
    fig=plt.figure(figsize=(cm2inch(col),cm2inch(24.7*heightfrac)))
    return fig

            
def NatureRevisions_Figure1(D):
    #fig=naturefig(1.5,heightfrac=1) 
    i=1
    
    
    letters =["a","b","c","d","e"]
    letters2=["f","g","h","i","j"]
    for X in ["ANZDA","MADA","MXDA","NADA","OWDA"]:
        solver = getattr(D,X).solver
        fac=da.get_orientation(solver)
        calib_period=('1920-1-1','2005-12-31')
        eof1 = fac*solver.eofs()[0]
        
        
       
          
       # plt.subplot(2,5,i)
        plt.subplot(5,2,2*i-1)
        
        m=b.plot_regional(eof1,X,vmin=-.15,vmax=.15)
        
        m.drawcoastlines(color='gray')
            
      
        plt.title("("+letters[i-1]+")",fontsize=6)#: "+X+" fingerprint",fontsize=8)
        if letters[i-1]=="a":
            print "AAAA"
        plt.subplot(5,2,2*i)
        #if i==1:
         #   fig = plt.gcf()
          #  cb_ax = fig.add_axes([0.1, 0.55, 0.83, 0.02])
           # plt.colorbar(label="EOF Loading",cax=cb_ax)
        
        
       # i+=1
        #plt.subplot(2,5,i+5)
        
        
        cru=getattr(D,X).project_cru_on_solver(start='1901-1-1')
        Plotting.time_plot(cru-np.ma.average(cru(time=calib_period)),lw=1,color=get_dataset_color("CRU"),label="CRU")
        dai=getattr(D,X).project_dai_on_solver(start='1901-1-1')
        Plotting.time_plot(dai-np.ma.average(dai(time=calib_period)),lw=1,color=get_dataset_color("DAI"),label="DAI")
        trees = getattr(D,X).projection(time=('1900-1-1','1975-1-1'))
        Plotting.time_plot(trees-np.ma.average(trees(time=calib_period)),lw=1,color=get_dataset_color("tree"),label=X)
        pc1= fac*solver.pcs()[:,0](time=('1900-1-1','2050-12-31'))
        Plotting.time_plot(pc1-np.ma.average(pc1(time=calib_period)),lw=1,color='k',label="PC1")
        plt.ylim(-.3,.3)
        plt.setp(plt.gca().get_yticklabels(),fontsize=6)
        
            
        plt.ylabel("Temporal Amplitude",fontsize=6)
        if i!=5:
            plt.xticks([])
            
        else:
            plt.legend(loc=0,ncol=2,fontsize=6)
            plt.setp(plt.gca().get_xticklabels(),fontsize=6)
            plt.xlabel("Year",fontsize=6)
        
        
        plt.title("("+letters2[i-1]+")",fontsize=6)#: "+X+" PC1 and projections",fontsize=8)
        i+=1
    #colorbar kludge
    if 1:
        fig = plt.gcf()
        #fig.subplots_adjust(left=-.15)
        #fig.subplots_adjust(bottom=.075)
        axes = plt.gcf().axes
        ax = axes[0]
        #left=fig.subplotpars.left+fig.subplotpars.hspace/2.
        left=0.03
        height = 0.01
        #bottom=fig.subplotpars.bottom-height*1.5
        bottom=.05
        width=.22
        #width=0.5-fig.subplotpars.hspace*1.5
    
        
        
        cb_ax = fig.add_axes([left,bottom,width,height])    
        ax.clear()
        eof1= D.ANZDA.solver.eofs()[0]*da.get_orientation(D.ANZDA.solver)
        m=b.plot_regional(eof1,"ANZDA",vmin=-.15,vmax=.15,ax=ax,cax=cb_ax,orientation='horizontal')
        m.drawcoastlines(color='gray',ax=ax)
        #cb_ax.yaxis.set_ticks_position('left')
        #cb_ax.yaxis.set_label_position("left")
        cb_ax.set_xticklabels(["-0.15","","","0","","","0.15"])
        plt.setp(cb_ax.get_xticklabels(),fontsize=6)
        plt.setp(cb_ax.xaxis.get_label(),fontsize=6)
        fig.axes[0].set_title("(a)",fontsize=6) # for some reason this disappears

        return cb_ax

        

def NatureRevisions_Figure3(D,start_year=2019):
    plt.subplot(311)
    Plotting.time_plot(D.ALL.get_noise(),color=get_dataset_color("tree"),label="Pre-industrial Noise",lw=1)
    Plotting.time_plot(D.ALL.projection(time=('1850-1-1','1975-12-31')),c="k",lw=1)
    plt.ylabel("Projection (temporal amplitude)",fontsize=8)
    plt.xlabel("Year")
    plt.legend()
    plt.title("(a)")#: Global Drought Atlas Projection Onto Fingerprint",fontsize=8)
    plt.subplot(312)
    t1900 = cdtime.comptime(1900,7,1)
    times = np.arange(1,201)
    

    likely=stats.norm.interval(.66)[1]
    vlikely=stats.norm.interval(.90)[1]
    certain=stats.norm.interval(.99)[1]
    plt.axhline(likely,lw=1,color=cm.Reds(.33),label="Likely",alpha=.5)
    plt.axhline(vlikely,lw=1,color=cm.Reds(.66),label="Very Likely",alpha=.5)
    plt.axhline(certain,lw=1,color=cm.Reds(.99),label="Virtually Certain",alpha=.5)
    D.ALL.time_of_emergence(t1900,times=times,color="k",lw=3)
    pdsi_SN_figure(D,cdtime.comptime(1900,1,1),use_dai=True)
   
    plt.title("(b)")#: Model-predicted time of emergence",fontsize=8)
    
    plt.subplot(313)
    hist_start = cdtime.comptime(1900,1,1)
    times = np.arange(10,2100-start_year)
    for X in ["ALL","ANZDA","MADA","MXDA","NADA","OWDA"]:
        getattr(D,X).time_of_emergence(cdtime.comptime(start_year,1,1),times=times,color=colorregions(X),uncertainty="shade")
    certain=stats.norm.interval(.99)[1]
    plt.axhline(certain,lw=1,color=cm.Reds(.99),label="Virtually Certain", alpha=.5)
    plt.title("(c)")#: Time of Emergence (assuming "+str(start_year)+" start)",fontsize=8)

    ax1,ax2,ax3=plt.gcf().axes
    ax2.set_xlim(1900,2055)
    ax2.set_ylim(-3,7)
    for ax in [ax1,ax2,ax3]:
        ax.set_xlabel(ax.get_xlabel(),fontsize=8)
        ax.set_ylabel(ax.get_ylabel(),fontsize=8)
    ax3.set_xlim(start_year+10,start_year+55)
    ax2.legend(fontsize=6)
    ax3.legend(fontsize=6)
    ax1.legend(fontsize=6)

    for ax in plt.gcf().axes:
        plt.setp(ax.get_xticklabels(),fontsize=6)
        plt.setp(ax.get_yticklabels(),fontsize=6)
        plt.setp(ax.xaxis.get_label(),fontsize=6)
        plt.setp(ax.yaxis.get_label(),fontsize=6)
        
from matplotlib.patches import Polygon
def NatureRevisions_Figure6(D):
    histdata = []
    aerosoldata=[]
    recentdata=[]
    datadict={}
    recent_start = cdtime.comptime(1981,1,1)
    recent_stop=cdtime.comptime(2017,12,31)
    #historical_stop = cdtime.comptime(2005,12,31)

    aerosol_start = cdtime.comptime(1950,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)

    hist_start = cdtime.comptime(1900,1,1)
    hist_stop = cdtime.comptime(1949,12,31)
    ax1=plt.subplot(131)
    ax2=plt.subplot(132)
    ax3=plt.subplot(133)
    xc=1.5
    for X in ["ALL","ANZDA","MADA","MXDA","NADA","OWDA"]:
        hist=getattr(D,X).for_figure_4(hist_start,stop_time=hist_stop,include_trees=True,include_cru=True,include_dai=True)
        histdata+=[hist["noise"]/np.std(hist["noise"]),hist["modslopes"]/np.std(hist["noise"])]
        ax1.plot([hist["tree_rings"]/np.std(hist["noise"])],[xc],"o",color=get_dataset_color("tree"))
        ax1.plot([hist["dai"]/np.std(hist["noise"])],[xc],"o",color=get_dataset_color("dai"))
        ax1.plot([hist["cru"]/np.std(hist["noise"])],[xc],"o",color=get_dataset_color("cru"))

        aerosol=getattr(D,X).for_figure_4(aerosol_start,stop_time=aerosol_stop,include_trees=True,include_cru=True,include_dai=True)
        aerosoldata+=[aerosol["noise"]/np.std(aerosol["noise"]),aerosol["modslopes"]/np.std(aerosol["noise"])]
        ax2.plot([aerosol["tree_rings"]/np.std(aerosol["noise"])],[xc],"o",color=get_dataset_color("tree"))
        ax2.plot([aerosol["dai"]/np.std(aerosol["noise"])],[xc],"o",color=get_dataset_color("dai"))
        ax2.plot([aerosol["cru"]/np.std(aerosol["noise"])],[xc],"o",color=get_dataset_color("cru"))

        recent=getattr(D,X).for_figure_4(recent_start,stop_time=recent_stop,include_trees=False,include_cru=True,include_dai=True)
        recentdata+=[recent["noise"]/np.std(recent["noise"]),recent["modslopes"]/np.std(recent["noise"])]
        
        ax3.plot([recent["dai"]/np.std(recent["noise"])],[xc],"o",color=get_dataset_color("dai"))
        ax3.plot([recent["cru"]/np.std(recent["noise"])],[xc],"o",color=get_dataset_color("cru"))
        xc+=2
    
    #fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title('A Boxplot Example')
    #fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    boxplot_mydata(histdata,ax1)
    ax1.set_xlim(-5,5)
    region_color_axis(ax1,ticks=True)
    ax1.set_title("(a):1900-1949",fontsize=6)
    ax1.set_xlabel("S/N")
    ax1.axvline(0,ls=":",c="k")
    boxplot_mydata(aerosoldata,ax2)
    ax2.set_title("(b):1950-1975",fontsize=6)
    ax2.axvline(0,ls=":",c="k")
    ax2.set_xlim(-5,5)
    region_color_axis(ax2,ticks=False)
    ax2.set_xlabel("S/N")
    boxplot_mydata(recentdata,ax3)
    ax3.set_xlim(-5,5)
    region_color_axis(ax3,ticks=False)
    ax3.set_title("(c): 1981-2017",fontsize=6)
    ax3.set_xlabel("S/N")
    ax3.axvline(0,ls=":",c="k")
    for ax in plt.gcf().axes:
        plt.setp(ax.get_xticklabels(),fontsize=6)
        plt.setp(ax.get_yticklabels(),fontsize=6)
        plt.setp(ax.xaxis.get_label(),fontsize=6)
        plt.setp(ax.yaxis.get_label(),fontsize=6)
    
def boxplot_mydata(histdata,ax1):        
    bp = ax1.boxplot(histdata, notch=0, sym='+', vert=0, whis=1.5)
    plt.setp(bp['boxes'], color='w')
    plt.setp(bp['whiskers'], color='k')
    plt.setp(bp['fliers'], color='k', marker='+')
    plt.setp(bp['medians'], color='k')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    #ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)

    # Hide these grid behind plot objects


    # Now fill the boxes with desired colors
    boxColors = [get_dataset_color("tree_noise"), get_dataset_color("H85")]
    numBoxes = len(histdata)
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = np.column_stack([boxX, boxY])

        k = i % 2
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        # med = bp['medians'][i]
        # medianX = []
        # medianY = []
        # for j in range(2):
        #     medianX.append(med.get_xdata()[j])
        #     medianY.append(med.get_ydata()[j])
        #     ax1.plot(medianX, medianY, 'k')
        #     medians[i] = medianY[0]
        
        
def region_color_axis(ax1,ticks=True):
    testx=np.linspace(*ax1.get_xlim())
    i=0
    for X in ["ALL","ANZDA","MADA","MXDA","NADA","OWDA"]:
        ax1.fill_between(testx,np.zeros_like(testx)+.5+2*i,np.zeros_like(testx)+2.5+2*i,color=colorregions(X),alpha=0.1)
        i+=1
        if ticks:
            ax1.set_yticks(np.arange(1.5,12.5,2))
            ax1.set_yticklabels(["GDA","ANZDA","MADA","MXDA","NADA","OWDA"])
        else:
            ax1.set_yticks([])
         
    

def test_noise(D,L,marker="o"):
    early=('1400-1-1','1600-1-1')
    late = ('1650-1-1','1850-1-1')
    #earr=[]
    #larr=[]
    pvals=[]
    e=[]
    l=[]
    
    for X in ["ALL","ANZDA","MADA","MXDA","NADA","OWDA"]:
        proj=getattr(D,X).projection
        earr=b.bootstrap_slopes(proj(time=early),L).asma()
        e=np.std(earr)
        larr=b.bootstrap_slopes(proj(time=late),L).asma()
        l=np.std(larr)
        print colorregions(X)
        plt.plot([e],[l],marker=marker,c=colorregions(X),label=X)
        pvals+=[stats.ks_2samp(earr.compressed(),larr.compressed())[1]]
    return pvals
        
        
def aerosol_fingerprint(D):
    NO,YES=b.aerosol_indirect(D)
    aerosol_start = cdtime.comptime(1951,1,1)
    aerosol_stop = cdtime.comptime(1975,12,31)
    yessolver = b.Eof(YES(time=(aerosol_start,aerosol_stop)))
    nosolver = b.Eof(NO(time=(aerosol_start,aerosol_stop)))
    yesfac=da.get_orientation(yessolver)
    nofac=da.get_orientation(nosolver)

    plt.subplot(221)
    m=b.landplot(yesfac*yessolver.eofs()[0],vmin=-.1,vmax=.1)
    #m.drawcoastlines()
    m.fillcontinents(color="gray",zorder=0)
    plt.colorbar(orientation="horizontal",label="EOF loading")
    plt.title("(a): EOF1 for models with aerosol indirect effect",fontsize=8)
    plt.ylim(-60,90)

    plt.subplot(222)
    m=b.landplot(nofac*nosolver.eofs()[0],vmin=-.1,vmax=.1)
    #m.drawcoastlines()
    m.fillcontinents(color="gray",zorder=0)
    plt.ylim(-60,90)
    plt.colorbar(orientation="horizontal",label="EOF loading")
    plt.title("(b) EOF1 for models without aerosol indirect effect",fontsize=8)

    plt.subplot(223)
    Plotting.time_plot(yesfac*yessolver.pcs()[:,0],c="k")
    plt.ylabel("Temporal Amplitude")
    plt.title("(c) PC1 for models with aerosol indirect effect",fontsize=8)

    plt.subplot(224)
    Plotting.time_plot(nofac*nosolver.pcs()[:,0],c="k")
    plt.ylabel("Temporal Amplitude")
    plt.title("(d) PC1 for models without aerosol indirect effect",fontsize=8)
    for ax in plt.gcf().axes:
        plt.setp(ax.get_xticklabels(),fontsize=6)
        plt.setp(ax.get_yticklabels(),fontsize=6)
        plt.setp(ax.xaxis.get_label(),fontsize=6)
        plt.setp(ax.yaxis.get_label(),fontsize=6)
    
