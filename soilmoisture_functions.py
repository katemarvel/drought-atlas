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

### Import scipy routines for smoothing, interpolation
from scipy.interpolate import interp1d
from scipy.optimize import brentq,fminbound
import scipy.ndimage as ndimag

import CMIP5_tools as cmip5
import DA_tools as da
from Plotting import *
from seasonal_cycle_utils import mask_data

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

import BEN_PDSI as b


def soilmoisture(depth,mask=None):
    if depth == "pdsi":
        f = cdms.open("../DROUGHT_ATLAS/CMIP5/pdsi.summerseason.ensemble.hist.rcp85.nc")
        variable="pdsi"
    else:
        f = cdms.open("../DROUGHT_ATLAS/CMIP5/sm"+depth+".summerseason.ensemble.hist.rcp85.nc")
        variable = "sm"+depth
    sm = b.get_rid_of_bad(f(variable))
    sm = MV.masked_where(np.isnan(sm),sm)
    f.close()
    if mask is not None:
        sm = mask_data(sm,mask)
    return sm

def soilmoisture_fingerprints(mask,name=None,fortalk=False):
    Fingerprints={}
    depths = ["pdsi","30cm","2m"]
    if fortalk:
        letters=["","",""]
    else:
        letters = ["(a): ","(b): ","(c): "]
    pcs = []
    pclabels = []
    for depth in depths:
        i=depths.index(depth)
        if fortalk:
            plt.figure()
        else:
            plt.subplot(2,2,i+1)
        sm = soilmoisture(depth,mask=mask)
        solver = Eof(MV.average(sm,axis=0),weights='area')
        Fingerprints[depth]=solver
        fac = da.get_orientation(solver)
        if name is None:
            m=b.landplot(fac*solver.eofs()[0],vmin=-.1,vmax=.1)
            plt.colorbar(orientation='horizontal',label='EOF loading')
        else:
            m=b.plot_regional(fac*solver.eofs()[0],name,vmin=-.1,vmax=.1)
            m.drawcountries()
        m.drawcoastlines(color='gray')
        
        
        if depth is not "pdsi":
            plt.title(letters[i]+depth+" fingerprint")
        else:
            plt.title(letters[i]+" PDSI fingerprint")
        pcs+=[fac*solver.pcs()[:,0]]

    if fortalk:
        plt.figure()
    else:
        plt.subplot(2,2,4)
    for i in range(3):
        if depths[i]=="pdsi":
            label="PDSI"
        else:
            label=depths[i]
        time_plot(pcs[i],label=label,lw=3,color=cm.copper(i/2.))
    plt.legend(loc=0)
    plt.title("(d): Principal Components")
    plt.xlabel("Time")
    plt.ylabel("Temporal amplitude")
    plt.xlim(1900,2100)
    return Fingerprints
    

class SoilMoisture():
    def __init__(self,mask,data=None):
        self.mask=mask
        depths = ["30cm","2m"]
        self.depths=depths
        self.soilmoisture = {}
        self.solvers = {}
        self.mma = {}
        self.noise={}
        self.P={}
        self.OBS_PROJECTIONS={}
        self.TOE={}
        self.TOE_start='None'
        if data is not None:
            self.soilmoisture,self.mma,self.solvers = data
        else:
            for depth in depths:
                self.soilmoisture[depth]=soilmoisture(depth,mask=mask)
                self.mma[depth]=MV.average(self.soilmoisture[depth],axis=0)
                self.solvers[depth] = Eof(self.mma[depth],weights='area')
    def project_piControl_on_solver(self,depth):
        if depth in self.noise.keys():
            
            pass
        else:
            direc = "/Volumes/Marvel/PICTRL/SM"+depth+"_REGRIDDED_SUMMER/"
            files = glob.glob(direc+"*")
            npiC=len(files)

            fname=files[0]

            f=cdms.open(fname)
            piC_pdsi_regrid=f("sm"+depth+"_summer")
            piC_pdsi_regrid = MV.masked_where(np.isnan(piC_pdsi_regrid),piC_pdsi_regrid)
            mask =self.solvers[depth].eofs()[0].mask
            grid=self.soilmoisture[depth].getGrid()
            nyears = piC_pdsi_regrid.shape[0]
            
            piC_mask = mask_data(piC_pdsi_regrid,mask)
            newmask = np.prod(~piC_mask.mask,axis=0)
          
          
            solver = Eof(mask_data(self.mma[depth],newmask==0),weights='area')
            fac = da.get_orientation(solver)
            p=solver.projectField(piC_mask)[:,0]*fac
            for i in range(npiC)[1:]:
                fname=files[i]
                f=cdms.open(fname)
                piC_pdsi_regrid=f("sm"+depth+"_summer")
                piC_pdsi_regrid = MV.masked_where(np.isnan(piC_pdsi_regrid),piC_pdsi_regrid)
                nyears += piC_pdsi_regrid.shape[0]
                
                piC_mask = mask_data(piC_pdsi_regrid,mask)
                newmask = np.prod(~piC_mask.mask,axis=0)
               
               
                solver = Eof(mask_data(self.mma[depth],newmask==0),weights='area')
                fac = da.get_orientation(solver)
                f.close()
                p=MV.concatenate((p,fac*solver.projectField(piC_mask)[:,0]))
            tax=cdms.createAxis(np.arange(nyears))
            tax.designateTime()
            tax.units = 'years since 0000-7-1'
            tax.id="time"
            p.setAxis(0,tax)
            self.noise[depth] = p
    def model_projections(self,depth):
       
        if depth in self.P.keys():
            pass
        
        else:
            to_proj = mask_data(self.soilmoisture[depth],self.solvers[depth].eofs()[0].mask)
            P=MV.zeros(to_proj.shape[:2])
            for i in range(to_proj.shape[0]):
               
                tp = to_proj[i]
                mma_mask = mask_data(self.mma[depth],tp[0].mask)
                solver = Eof(mma_mask,weights='area')
                tp = mask_data(tp,solver.eofs()[0].mask)
                fac=da.get_orientation(solver)

                P[i] = solver.projectField(tp)[:,0]*fac
            P.setAxisList(to_proj.getAxisList()[:2])
            self.P[depth]=P
            self.P[depth].getTime().id="time"
    
    def sn_at_time(self,start_time,L,depth,overlapping=True):
        self.project_piControl_on_solver(depth)
       
        
       
        self.model_projections(depth)
        stop_time=start_time.add(L,cdtime.Years)
        modslopes = cmip5.get_linear_trends(self.P[depth](time=(start_time,stop_time)))
        if overlapping:
            noiseterm = b.bootstrap_slopes(self.noise[depth],L)
        else:
            noiseterm = da.get_slopes(self.noise[depth],L)/365.
        return modslopes,noiseterm

    def time_of_emergence(self,start_time,depth,times = np.arange(10,76),ax=None,**kwargs):
        
        if not (depth in self.P.keys()):
            self.model_projections(depth)
        if (not (depth in self.TOE.keys())) or (self.TOE_start != str(start_time)):
            nmod,nyears = self.P[depth].shape
            self.TOE[depth]=MV.zeros((nmod,len(times)))
            for i in range(len(times)):
                L=times[i]
                modslopes,noiseterm = self.sn_at_time(start_time,L,depth)
                sns=modslopes/np.std(noiseterm)
                self.TOE[depth][:,i]=sns
            self.TOE[depth].setAxis(0,self.P[depth].getAxis(0))
            self.TOE_start = str(start_time)
        if ax is not None:
            endyears = start_time.year+times
            ax.plot(endyears,np.ma.average(self.TOE[depth].asma(),axis=0),lw=4,label=depth+" model mean signal",**kwargs)
            ax.fill_between(endyears,np.ma.min(self.TOE[depth].asma(),axis=0),np.ma.max(self.TOE[depth].asma(),axis=0),alpha=.3,**kwargs)
            #ax.axhline(stats.norm.interval(.9)[-1],c="r",lw=3,label="90% significance threshold")
            ax.set_xlabel("Trend end year")
            ax.set_ylabel("Signal-to-noise ratio")
            plt.legend(loc=0)


    def histograms(self,start_time,stop_time,depth,overlapping=True):
        
        
        L=stop_time.year-start_time.year+1
        modslopes,noiseterm = self.sn_at_time(start_time,L,depth,overlapping=overlapping)
        ns=np.std(noiseterm)
        
        plt.hist(modslopes/ns,20,normed=True,color=cm.Oranges(.8),alpha=.5)
        lab = str(start_time.year)+"-"+str(stop_time.year)
        da.fit_normals_to_data(modslopes/ns,color=cm.Oranges(.9),lw=3,label=lab+" Model projections")

        plt.hist(noiseterm/ns,20,normed=True,color=cm.Purples(.8),alpha=.5)
        da.fit_normals_to_data(noiseterm/ns,color=cm.Purples(.9),lw=3,label="Noise")
        plt.xlabel("S/N")
        plt.ylabel("Normalized Frequency")
    def obs_SN(self,start_time,stop_time,depth,overlapping=True):
        self.project_soilmoisture("MERRA2")
        self.project_soilmoisture("GLEAM")
        L=stop_time.year-start_time.year+1
        modslopes,noiseterm = self.sn_at_time(start_time,L,depth,overlapping=overlapping)
        ns=np.std(noiseterm)
        
        plt.hist(modslopes/ns,20,normed=True,color=cm.Oranges(.8),alpha=.5)
        lab = str(start_time.year)+"-"+str(stop_time.year)
        da.fit_normals_to_data(modslopes/ns,color=cm.Oranges(.9),lw=3,label=lab+" trends in H85 projections onto fingerprint")

        plt.hist(noiseterm/ns,20,normed=True,color=cm.Purples(.8),alpha=.5)
        da.fit_normals_to_data(noiseterm/ns,color=cm.Purples(.9),lw=3,label=str(L)+"-year trends in piControl projection onto fingerprint")
        plt.xlabel("S/N")
        plt.ylabel("Normalized Frequency")

        merra=self.OBS_PROJECTIONS["MERRA2"][depth](time=(start_time,stop_time))
        gleam=self.OBS_PROJECTIONS["GLEAM"][depth](time=(start_time,stop_time))
        merrasig=cmip5.get_linear_trends(merra)/ns
        plt.axvline(merrasig,label="MERRA2",c="b",lw=3)
        gleamsig=cmip5.get_linear_trends(gleam)/ns
        plt.axvline(gleamsig,label="GLEAM",c="r",lw=3)
        plt.legend()
        
    def project_soilmoisture(self,dataset):
        mask=self.mask
        self.OBS_PROJECTIONS[dataset]={}
        surf,root=standardize_soilmoisture(dataset)
        
        surfsolver=Eof(mask_data(self.mma["30cm"],mask),weights='area') 
        surfmask=mask_data(surf,surfsolver.eofs()[0].mask)
        surfsolver=Eof(mask_data(self.mma["30cm"],surfmask[-1].mask),weights='area')
        surfanom=surfmask-MV.average(surfmask,axis=0)
        self.OBS_PROJECTIONS[dataset]["30cm"]=surfsolver.projectField(surfanom)[:,0]*da.get_orientation(surfsolver)

        rootsolver=Eof(mask_data(self.mma["2m"],mask),weights='area') 
        rootmask=mask_data(root,rootsolver.eofs()[0].mask)
        rootsolver=Eof(mask_data(self.mma["2m"],rootmask[-1].mask),weights='area')
        rootanom=rootmask-MV.average(rootmask,axis=0)
        self.OBS_PROJECTIONS[dataset]["2m"]=rootsolver.projectField(rootanom)[:,0]*da.get_orientation(rootsolver)
        

    
    


def standardize_soilmoisture(dataset):
    f = cdms.open("../DROUGHT_ATLAS/PROCESSED/ALL_ANZDA.nc")
    obs = f("pdsi")
    obs = MV.masked_where(np.isnan(obs),obs)
    obs = MV.masked_where(np.abs(obs)>90,obs)
    obs = obs(time=('1921-1-1','2000-12-31'))
    pdsi=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
    f.close()
    mu_p=MV.average(pdsi,axis=0)
    sig_p=genutil.statistics.std(pdsi,axis=0)



    f=cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/"+dataset+"_soilmoisture_summerseason.nc")
    gleam30cm=f("smsurf")
    mask=pdsi[0].mask
    gleam30cmmask = b.mask_data(gleam30cm,mask)
    mu_s = MV.average(gleam30cmmask,axis=0)
    sig_s=genutil.statistics.std(gleam30cmmask,axis=0)

    surf = (gleam30cmmask-mu_s+mu_p)*(sig_p/sig_s)


    gleam2m=f("smroot")
    mask=pdsi[0].mask
    gleam2mmask = b.mask_data(gleam2m,mask)
    mu_s2 = MV.average(gleam2mmask,axis=0)
    sig_s2=genutil.statistics.std(gleam2mmask,axis=0)
    root = (gleam2mmask-mu_s2+mu_p)*(sig_p/sig_s)
    return surf,root

    


    # surfsolver=Eof(b.mask_data(sm.mma["30cm"],mask),weights='area') 
    # surfmask=b.mask_data(surf,surfsolver.eofs()[0].mask)
    # surfsolver=Eof(b.mask_data(sm.mma["30cm"],surfmask[-1].mask),weights='area')
    # surfanom=surfmask-MV.average(surfmask,axis=0)

    # time_plot(surfsolver.projectField(surfanom)[:,0],lw=3,color=cm.copper(.2))
    # plt.title("Surface soil moisture")
    # plt.ylabel("Projection on fingerprint")

    # plt.figure()

    # rootsolver=Eof(b.mask_data(sm.mma["2m"],mask),weights='area') 
    # rootmask=b.mask_data(root,rootsolver.eofs()[0].mask)
    # rootsolver=Eof(b.mask_data(sm.mma["2m"],rootmask[-1].mask),weights='area')
    # rootanom=rootmask-MV.average(rootmask,axis=0)

    # time_plot(rootsolver.projectField(rootanom)[:,0],lw=3,color=cm.copper(.8))
    # plt.title("Root soil moisture")
    # plt.ylabel("Projection on fingerprint")
    
