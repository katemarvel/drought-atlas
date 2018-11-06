#regrid all observations to common grid and combine
#Write code to plot maps of specific regions
import CMIP5_tools as cmip5
import MV2 as MV
import glob
import cdms2 as cdms
import numpy as np
import cdutil
from seasonal_cycle_utils import mask_data
from eofs.cdms import Eof
import DA_tools as da
from seasonal_cycle_utils import mask_data
import string
import os


def djf_pdsi(X):
    """Average CMIP5 PDSI over austral summer months """
    dec = X[:,-1,:,:]
    jan = X[:,0,:,:]
    feb = X[:,1,:,:]
    djf=(dec[:-1]+jan[1:]+feb[1:])/3.
    tax = djf.getAxis(0)
    yr=int(X.getAxis(0)[0])+1
    tax.units = "years since "+str(yr)+"-1-15"
    tax.id="time"
    djf.setAxis(0,tax)
    djf.id = 'pdsi'
    return djf


def jja_pdsi(X):
    """Average CMIP5 PDSI over boreal summer months """
    jja=MV.average(X[:,5:8],axis=1)
    tax = jja.getAxis(0)
    tax.units = 'years since 0000-7-1'
    tax.id="time"
    jja.setAxis(0,tax)
    jja.id = 'pdsi'
    return jja



def summerseason_pdsi(X):
    """ Average X over the summer season: DJF in the SH, JJa in the NH"""
    nyear,nmonth,nlat,nlon = X.shape
    yax,monax,latax,lonax=X.getAxisList()
    #SS = MV.zeros((nyear-1,nlat,nlon))
    DJF = djf_pdsi(X)
    JJA = jja_pdsi(X)[1:] #get rid of 1st year since no corresponding DJF
    lon=X.getLongitude()[:]
    lat=X.getLatitude()[:]
    x,y=np.meshgrid(lon,lat)
    ny=nyear-1
    latmesh=np.repeat(y[np.newaxis],ny,axis=0)
    SS = MV.where(latmesh<0,DJF,JJA)
    SS.setAxisList([DJF.getTime(),latax,lonax])
    SS.id="pdsi"
    SS.long_name="Summer season PDSI"
    return SS

def summerseason_soilmoisture(fname):
  
    f = cdms.open(fname)
    djf_2m = f("sm2m_DJF")
    jja_2m = f("sm2m_JJA")[:,1:]

    nens,nyear,nlat,nlon=djf_2m.shape

    lon=djf_2m.getLongitude()[:]
    lat=djf_2m.getLatitude()[:]
    x,y=np.meshgrid(lon,lat)
    latmesh=np.repeat(y[np.newaxis],nyear,axis=0)
    latmesh = np.repeat(latmesh[np.newaxis],nens,axis=0)

    tax = djf_2m.getAxis(1)
    tax.units = 'years since 0000-1-1'
    tax.id="time"
    tax.designateTime()
    
    SS_2m=MV.where(latmesh<0,djf_2m,jja_2m)
    SS_2m=MV.masked_where(np.isnan(SS_2m),SS_2m)
    SS_2m.setAxisList(djf_2m.getAxisList())
    SS_2m.setAxis(1,tax)
    SS_2m.id="sm2m"


    djf_30cm = f("sm30cm_DJF")
    jja_30cm = f("sm30cm_JJA")[:,1:]
    
    SS_30cm=MV.where(latmesh<0,djf_30cm,jja_30cm)
    SS_30cm=MV.masked_where(np.isnan(SS_30cm),SS_30cm)
    SS_30cm.setAxisList(djf_30cm.getAxisList())
    SS_30cm.setAxis(1,tax)
    SS_30cm.id="sm30cm"
    f.close()

    
    return SS_2m,SS_30cm
    
    
def regrid_picontrol_summerseason(variable="pdsi"):

    direc = "/Volumes/Marvel/PICTRL/"
    cmd = "mkdir /Volumes/Marvel/PICTRL/"+string.upper(variable)+"_REGRIDDED_SUMMER/"
    os.system(cmd)
    files = glob.glob(direc+"*nc")
    npiC=len(files)
    fgrid = cdms.open("OBS/gpcp.precip.mon.mean.nc")
    grid=fgrid("precip").getGrid()

    for fname in files:   
   
        f=cdms.open(fname)
        piC_JJA=f(variable+"_JJA")[0,1:] #take first ensemble member
        piC_JJA = MV.masked_where(np.isnan(piC_JJA),piC_JJA)


        tax=cdms.createAxis(np.arange(piC_JJA.shape[0]))
        tax.designateTime()
        tax.units = 'years since 0000-7-1'
        tax.id="time"

        piC_JJA.setAxis(0,tax)

        piC_DJF=f(variable+"_DJF")[0] #take first ensemble member
        piC_DJF = MV.masked_where(np.isnan(piC_DJF),piC_DJF)
        piC_JJA.setAxis(0,tax)
        piC_DJF.setAxis(0,tax)

        piC_JJA_regrid = piC_JJA.regrid(grid,regridTool='regrid2')
        piC_DJF_regrid = piC_DJF.regrid(grid,regridTool='regrid2')   

        nyears = piC_JJA.shape[0]

        tax=cdms.createAxis(np.arange(piC_JJA.shape[0]))
        tax.designateTime()
        tax.units = 'years since 0000-7-1'
        tax.id="time"

        lon=piC_DJF_regrid.getLongitude()[:]
        lat=piC_DJF_regrid.getLatitude()[:]
        x,y=np.meshgrid(lon,lat)
        latmesh=np.repeat(y[np.newaxis],nyears,axis=0)

        SS = MV.where(latmesh<0,piC_DJF_regrid,piC_JJA_regrid)
        SS.setAxis(0,tax)
        SS.id=variable+"_summer"
        f.close()
        fgrid.close()
        fw = cdms.open(fname.replace("PICTRL/","PICTRL/"+string.upper(variable)+"_REGRIDDED_SUMMER/"),"w")
        fw.write(SS)
        fw.close()
    
    

     

def regrid_and_truncate(SS,the_grid=None):
    """ Put all CMIP5 files on the same grid"""
    if the_grid is None:
        fobs = cdms.open("OBS/gpcp.precip.mon.mean.nc")
        the_grid = fobs["precip"].getGrid()
    Xt = SS(time=('1900-1-1','2099-12-31'))
    Xr = Xt.regrid(the_grid,regridTool='regrid2')
    if the_grid is None:
        fobs.close()
    return Xr




def write_regridded_soilmoisture(label="summerseason.ensemble"):

    direc = "/Volumes/Marvel/FromBen/SCALEDSM/"
    allfiles = glob.glob(direc+"*")
    #Figure out how many ensemble members there are
    L=0
    for fname in allfiles:
        f=cdms.open(fname)
        L+=f["sm2m_DJF"].shape[0]
        f.close()
    
    
    modelnames=[]
    #Do first ensemble member, first model to set up shape for regridded arrays
    i=0
    fname=allfiles[0]
    SS2m,SS30cm=summerseason_soilmoisture(fname)
    nens=SS2m.shape[0]
    test_2m= regrid_and_truncate(SS2m[0])
    SS2m_regrid=MV.zeros((L,)+test_2m.shape)
    SS2m_regrid[i]=test_2m
    test_30cm = regrid_and_truncate(SS30cm[0])
    SS30cm_regrid=MV.zeros((L,)+test_30cm.shape)
    SS30cm_regrid[i]=test_30cm
    rip = "r1i1p1"
    modelnames+=[fname.split(".")[3]+"."+rip]
    print fname.split(".")[3]+"."+rip
    i+=1
    #if there are more ensemble members in this first model
    if nens>1:
        print "GOING INTO THIS LOOP"
        for j in range(nens)[1:]:
            SS2m_regrid[i]=regrid_and_truncate(SS2m[j])
            SS30cm_regrid[i]=regrid_and_truncate(SS30cm[j])
            rip="r"+str(j+1)+"i1p1"
            modelnames+=[fname.split(".")[3]+"."+rip]
            print fname.split(".")[3]+"."+rip
            i+=1
    #now loop over everything else
    for fname in allfiles[1:]:
        SS2m,SS30cm=summerseason_soilmoisture(fname)
        nens=SS2m.shape[0]
        for j in range(nens):
            SS2m_regrid[i]=regrid_and_truncate(SS2m[j])
            SS30cm_regrid[i]=regrid_and_truncate(SS30cm[j])
            rip="r"+str(j+1)+"i1p1"
            modelnames+=[fname.split(".")[3]+"."+rip]
            print fname.split(".")[3]+"."+rip
            i+=1
            
    #give regridded arrays the appropriate axes
    modax = cmip5.make_model_axis(modelnames)
    axlist = [modax]+test_2m.getAxisList()
    SS2m_regrid.setAxisList(axlist)
    SS2m_regrid.id="sm2m"
    SS30cm_regrid.setAxisList(axlist)
    SS30cm_regrid.id="sm30cm"
    fw2m = cdms.open("../DROUGHT_ATLAS/CMIP5/sm2m."+label+".hist.rcp85.nc","w")
    fw2m.write(SS2m_regrid)
    fw2m.close()

    fw30cm = cdms.open("../DROUGHT_ATLAS/CMIP5/sm30cm."+label+".hist.rcp85.nc","w")
    fw30cm.write(SS30cm_regrid)
    fw30cm.close()

    
    
    
        
        
            
    
    
    
    
        
        

def write_regridded_ensemble(the_grid=None,label="summerseason.ensemble"):
    """write ensemble of summer season PDSI"""
    direc = "/Volumes/Marvel/FromBen/PDSI/"
    allfiles = glob.glob(direc+"*")
    nf=len(allfiles)
    i=0
    f=cdms.open(allfiles[i])
    X=f("PDSI")
    SS=summerseason_pdsi(X)
    Xt=regrid_and_truncate(SS,the_grid=the_grid)
    
    ens = MV.zeros((nf,)+Xt.shape)+1.e20
    ens[i]=Xt
    for i in range(nf)[1:]:
        f=cdms.open(allfiles[i])
        X=f("PDSI")
        SS=summerseason_pdsi(X)
        Xt=regrid_and_truncate(SS,the_grid=the_grid)
        try:
            ens[i]=Xt
        except:
            continue
        axes = Xt.getAxisList()
        f.close()
    ens=MV.masked_where(np.abs(ens)>1.e10,ens)
    ens = MV.masked_where(np.isnan(ens),ens)
    ens.id="pdsi"
    modax = cmip5.make_model_axis(allfiles)
    ens.setAxisList([modax]+axes)
    fw = cdms.open("../DROUGHT_ATLAS/CMIP5/pdsi."+label+".hist.rcp85.nc","w")
    ens.id='pdsi'
    fw.write(ens)
    fw.close()
    return ens

def combine_observations(list_of_obs,grid='2.5',writename=None):
    """ Place observations on common grid and combine """
    #Common time axis
    starts = [cmip5.start_time(x) for x in list_of_obs]
    start = np.max(starts)
    stops = [cmip5.stop_time(x) for x in list_of_obs]
    stop = np.min(stops)
    newlist = [x(time=(start,stop)) for x in list_of_obs]
    L = newlist[0].shape[0]

     #get target grid
    if grid == '2.5':
        fgrid = cdms.open("OBS/gpcp.precip.mon.mean.nc")
        print "writing to GPCP grid"
    else:
        fgrid = cdms.open("OBS/gpcc.precip.mon.total.v7.nc")
        print "writing to gpcc grid"
    
    gpcp_grid = fgrid("precip").getGrid()
    #set up the combined observations
    combined_obs = MV.zeros((L,)+gpcp_grid.shape)+1.e20

    for X in newlist:
        Xr = X.regrid(gpcp_grid,regridTool='regrid2')
        indices=np.where(~Xr.flatten().mask)
        MV.put(combined_obs,indices,Xr.compressed())
    combined_obs = MV.masked_where(combined_obs>1.e10,combined_obs)
    combined_obs.setAxisList([newlist[0].getTime()]+gpcp_grid.getAxisList())
    combined_obs.id = "pdsi"
    
    if writename is not None:
        fw = cdms.open(writename,"w")
        fw.write(combined_obs)
        fw.close()
    fgrid.close()
    return combined_obs

def mask_models(CO,writename=None):
    """Create mask based on availability of observations at first time step and apply to CMIP5 hist+85 PDSI."""
    f = cdms.open("../DROUGHT_ATLAS/CMIP5/pdsi.summerseason.ensemble.hist.rcp85.nc")
    models = f("pdsi")
    models = MV.masked_where(np.isnan(models),models)
    masked_models = mask_data(models,CO[0].mask)
    masked_models.id='pdsi'
    if writename is not None:
        fw = cdms.open(writename,"w")
        fw.write(masked_models)
        fw.close()
    return masked_models
    
def rewrite_downloaded():
    ##NEED TO FIX TIME AXIS for all obs
    ##Transpose so axes are [time,lat,lon] and set time axis units to be years since summer 0 AD
    files = glob.glob("../DROUGHT_ATLAS/DOWNLOADED/*")
    for fil in files:
        write_transposed_and_masked(fil)
def write_transposed_and_masked(fil):
  
    atlas_name = fil.split("/")[-1].split(".")[0]
    print atlas_name
    f = cdms.open(fil)
    try:
        pdsi = f("pdsi")
        tpdsi = MV.transpose(pdsi)
    except:
        pdsi = f("PDSI")
        tpdsi = pdsi #For NADA the axes are OK
    mask_tpdsi=MV.masked_where(np.isnan(tpdsi),tpdsi)
    taxis = mask_tpdsi.getTime()
    taxis.units = "years since 0000-7-1"
    mask_tpdsi.setAxis(0,taxis)
    mask_tpdsi.id="pdsi"
    mask_tpdsi.name="pdsi"
    for key in f.attributes.keys():
        setattr(mask_tpdsi,key,f.attributes[key])
    writefname = "../DROUGHT_ATLAS/PROCESSED/"+atlas_name+".nc"
    wf = cdms.open(writefname,"w")
    wf.write(mask_tpdsi)
    wf.close()


def write_regridded():
    fnames = glob.glob("../DROUGHT_ATLAS/PROCESSED/ORIGINAL_GRID/*")
    for fil in fnames:
        atlas_name = fil.split("/")[-1].split(".")[0]
        
        f = cdms.open("../DROUGHT_ATLAS/PROCESSED/ORIGINAL_GRID/"+atlas_name+".nc")
        obs = f("pdsi")
        obs = MV.masked_where(np.isnan(obs),obs)
        obs = MV.masked_where(np.abs(obs)>90,obs)
        obs = obs(time=('1100-7-1','2020-12-31'))
        
        obs=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
        f.close()
        writefname = "../DROUGHT_ATLAS/PROCESSED/"+atlas_name+"2.5.nc"
        CO = combine_observations([obs],grid='2.5',writename=writefname)
        modelwritename = "../DROUGHT_ATLAS/CMIP5/pdsi."+atlas_name+"2.5.hist.rcp85.nc"
        mm = mask_models(CO,writename=modelwritename)
            
def write_combinations():
    """
    Make t combined NH drought atlases: ALL_OBS, which starts in 1100 and includes MADA, OWDA, and NADA, and ALL_OBS_plus_MEX, which includes the former + MXDA
    """
    list_of_obs = []
    for atlas_name in ["OWDA2.5","MADA2.5","NADA2.5"]:
        f = cdms.open("../DROUGHT_ATLAS/PROCESSED/"+atlas_name+".nc")
        obs = f("pdsi")
        obs = MV.masked_where(np.isnan(obs),obs)
        obs = MV.masked_where(np.abs(obs)>90,obs)
        obs = obs(time=('1100-7-1','2020-12-31'))
        obs=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
        f.close()
        list_of_obs+=[obs]
    writename = "../DROUGHT_ATLAS/PROCESSED/ALL_OBS.nc"
    CO = combine_observations(list_of_obs,writename=writename)
    modelwritename = "../DROUGHT_ATLAS/CMIP5/pdsi.ALL_OBS.hist.rcp85.nc"
    mm = mask_models(CO, writename = modelwritename)
    #now add in MXDA:
    atlas_name = "MXDA2.5"
    f = cdms.open("../DROUGHT_ATLAS/PROCESSED/"+atlas_name+".nc")
    obs = f("pdsi")
    obs = MV.masked_where(np.isnan(obs),obs)
    obs = MV.masked_where(np.abs(obs)>90,obs)
    obs = obs(time=('1100-7-1','2020-12-31'))
    obs=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
    f.close()
    list_of_obs+=[obs]
    writename = "../DROUGHT_ATLAS/PROCESSED/ALL_OBS_plus_MEX.nc"
    CO = combine_observations(list_of_obs,writename=writename)
    modelwritename = "../DROUGHT_ATLAS/CMIP5/pdsi.ALL_OBS_plus_MEX.hist.rcp85.nc"
    mm = mask_models(CO, writename = modelwritename)

def write_big_drought_atlas():
    list_of_obs = []
    for atlas_name in ["OWDA2.5","MADA2.5","NADA2.5","MXDA2.5","ANZDA2.5"]:
        f = cdms.open("../DROUGHT_ATLAS/PROCESSED/"+atlas_name+".nc")
        obs = f("pdsi")
        obs = MV.masked_where(np.isnan(obs),obs)
        obs = MV.masked_where(np.abs(obs)>90,obs)
        obs = obs(time=('1100-7-1','2020-12-31'))
        obs=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
        f.close()
        list_of_obs+=[obs]
    writename = "../DROUGHT_ATLAS/PROCESSED/ALL_ANZDA.nc"
    CO = combine_observations(list_of_obs,writename=writename)
    modelwritename = "../DROUGHT_ATLAS/CMIP5/pdsi.ALL_ANZDA.hist.rcp85.nc"
    mm = mask_models(CO, writename = modelwritename)
def exclude_NADA():
    """
    Make a combined NH drought atlas ALL_NONADA, which starts in 1400 and includes MADA, OWDA, ANZDA but NOT NADA
    """
    list_of_obs = []
    for atlas_name in ["OWDA2.5","MADA2.5","MXDA2.5","ANZDA2.5"]:
        f = cdms.open("../DROUGHT_ATLAS/PROCESSED/"+atlas_name+".nc")
        obs = f("pdsi")
        obs = MV.masked_where(np.isnan(obs),obs)
        obs = MV.masked_where(np.abs(obs)>90,obs)
        obs = obs(time=('1400-7-1','2020-12-31'))
        obs=mask_data(obs,obs.mask[0]) #Make all the obs have the same mask as the first datapoint
        f.close()
        list_of_obs+=[obs]
    writename = "../DROUGHT_ATLAS/PROCESSED/ALL_NONADA.nc"
    CO = combine_observations(list_of_obs,writename=writename)
    modelwritename = "../DROUGHT_ATLAS/CMIP5/pdsi.ALL_NONADA.hist.rcp85.nc"
    mm = mask_models(CO, writename = modelwritename)
     
##### MORE OBSERVATIONS TO THROW AT THIS THING #########
          
def dai_jja():
    f = cdms.open("../DROUGHT_ATLAS/pdsi.mon.mean.selfcalibrated.nc")
    dai=f("pdsi")
    cdutil.setTimeBoundsMonthly(dai)
    dai_jja = cdutil.JJA(dai)
    fgrid = cdms.open("OBS/gpcp.precip.mon.mean.nc")
    gpcp_grid = fgrid("precip").getGrid()
    fgrid.close()
    dai2=dai_jja.regrid(gpcp_grid,regridTool='regrid2')
    dai2.id="pdsi"
    for att in dai.attributes.keys():
        setattr(dai2,att,dai.attributes[att])
    fw = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/DAI_selfcalibrated.nc","w")
    fw.write(dai2)
    fw.close()
    return dai2

def cru_jja():

    f = cdms.open("../DROUGHT_ATLAS/scPDSI.cru_ts3.26early.bams2018.GLOBAL.1901.2017.nc")
    cru=f("scpdsi")
  
    cru.getTime().units = 'days since 1900-1-1'
    cdutil.setTimeBoundsMonthly(cru)
    cdutil.setTimeBoundsMonthly(cru)
    cru_jja = cdutil.JJA(cru)
    cru_jja = MV.masked_where(np.abs(cru_jja)>1000,cru_jja)

    fgrid = cdms.open("OBS/gpcp.precip.mon.mean.nc")
    gpcp_grid = fgrid("precip").getGrid()
    fgrid.close()
    cru2=cru_jja.regrid(gpcp_grid,regridTool='regrid2')
    cru2.id="pdsi"
    for att in cru.attributes.keys():
        setattr(cru2,att,cru.attributes[att])

    fw = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/CRU_selfcalibrated.nc","w")
    fw.write(cru2)
    fw.close()
    f.close()
    return cru2

#NEED TO MAKE SURE HAS THE SAME MASK AT ALL TIME STEPS


    
def soilmoisture_GLEAM():

    direc = "/Volumes/Marvel/monthly_v32a/soilmoist/"
    files = sorted(glob.glob(direc+"*nc"))
    
    L = len(files)
    i=0
    fname=files[i]
    f=cdms.open(fname)
    smroot_jja=MV.average(f("smroot")[5:8],axis=0)
    smsurf_jja=MV.average(f("smsurf")[5:8],axis=0)
    ROOT_jja = MV.zeros((L-1,)+smroot_jja.shape)
    SURF_jja = MV.zeros((L-1,)+smsurf_jja.shape)
    smroot_jf=MV.average(f("smroot")[0:2],axis=0)
    smsurf_jf=MV.average(f("smsurf")[0:2],axis=0)
    smroot_dec=f("smroot")[-1]
    smsurf_dec=f("smsurf")[-1]
    ROOT_djf = MV.zeros((L-1,)+smroot_jja.shape)
    SURF_djf = MV.zeros((L-1,)+smsurf_jja.shape)
    #ROOT[i]=smroot
    #SURF[i]=smsurf
    f.close()
    i+=1
    for fname in files[1:]:
        f=cdms.open(fname)
        smroot_jja=MV.average(f("smroot")[5:8],axis=0)
        smsurf_jja=MV.average(f("smsurf")[5:8],axis=0)
        ROOT_jja[i-1]=smroot_jja
        SURF_jja[i-1]=smsurf_jja
        smroot_jf=MV.average(f("smroot")[0:2],axis=0)
        smsurf_jf=MV.average(f("smsurf")[0:2],axis=0)
        #compute DJF average using DEC from the previous year
        ROOT_djf[i-1]=(1./3.)*(smroot_dec+2*smroot_jf)
        SURF_djf[i-1]=(1./3.)*(smsurf_dec+2*smsurf_jf)
        #update with dec of this year
        smroot_dec=f("smroot")[-1]
        smsurf_dec=f("smsurf")[-1]
        f.close()
        i+=1
        
    #set up time axis
    tax=cdms.createAxis(np.arange(L-1,dtype=np.float))
    tax.designateTime()
    tax.id='time'
    tax.units='years since 1981-7-15'
    axes = [tax,]+smroot_jja.getAxisList()

    #Get the common grid
    fgrid = cdms.open("OBS/gpcp.precip.mon.mean.nc")
    gpcp_grid = fgrid("precip").getGrid()
    fgrid.close()
    
    ROOT_jja=MV.masked_where(np.isnan(ROOT_jja),ROOT_jja)
    SURF_jja=MV.masked_where(np.isnan(SURF_jja),SURF_jja)
    ROOT_jja.id="smroot"
    ROOT_jja.setAxisList(axes)
    SURF_jja.id="smsurf"
    SURF_jja.setAxisList(axes)
   
    
    ROOT_jja_R=ROOT_jja.regrid(gpcp_grid,regridTool='regrid2')
    SURF_jja_R=SURF_jja.regrid(gpcp_grid,regridTool='regrid2')
    
    fw_jja = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_jja.nc","w")
    fw_jja.write(ROOT_jja_R)
    fw_jja.write(SURF_jja_R)
    fw_jja.close()

    ROOT_djf=MV.masked_where(np.isnan(ROOT_djf),ROOT_djf)
    SURF_djf=MV.masked_where(np.isnan(SURF_djf),SURF_djf)
    ROOT_djf.id="smroot"
    ROOT_djf.setAxisList(axes)
    SURF_djf.id="smsurf"
    SURF_djf.setAxisList(axes)
   
    
    ROOT_djf_R=ROOT_djf.regrid(gpcp_grid,regridTool='regrid2')
    SURF_djf_R=SURF_djf.regrid(gpcp_grid,regridTool='regrid2')
    
    fw_djf = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_djf.nc","w")
    fw_djf.write(ROOT_djf_R)
    fw_djf.write(SURF_djf_R)
    fw_djf.close()
        
    
def summerseason_GLEAM(jja,djf):
  
   

    nyear,nlat,nlon=djf.shape

    lon=djf.getLongitude()[:]
    lat=djf.getLatitude()[:]
    x,y=np.meshgrid(lon,lat)
    latmesh=np.repeat(y[np.newaxis],nyear,axis=0)
        
    
    SS=MV.where(latmesh<0,djf,jja)
    SS=MV.masked_where(np.isnan(SS),SS)
    SS.setAxisList(djf.getAxisList())
    SS.id=djf.id
    return SS

def write_root_and_surf_summerseason():
    f_djf = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_djf.nc")
    root_djf=f_djf("smroot")
    surf_djf=f_djf("smsurf")
    
    f_jja = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_jja.nc")
    root_jja=f_jja("smroot")
    surf_jja=f_jja("smsurf")

    fw = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_summerseason.nc","w")
    surf=summerseason_GLEAM(surf_jja,surf_djf)
    root=summerseason_GLEAM(root_jja,root_djf)
    surf.id="smsurf"
    root.id="smroot"
    fw.write(surf)
    fw.write(root)
    fw.close()

def merra2():
    fobs = cdms.open("OBS/gpcp.precip.mon.mean.nc")
    grid = fobs("precip").getGrid()
    fobs.close()
    merrafiles=sorted(glob.glob("/Volumes/Marvel/MERRA2/land/*2d_lnd_Nx*"))
    fname=merrafiles[0]
    f=cdms.open(fname)
    root_og=f("RZMC")
    root=root_og.regrid(grid,regridTool='regrid2')

    rootatts=root.getAxis(0).attributes
    surf_og=f("SFMC")
    surf=surf_og.regrid(grid,regridTool='regrid2')
    surfatts=surf.getAxis(0).attributes
    f.close()
    for fname in merrafiles[1:]:
        f=cdms.open(fname)
        root_og=f("RZMC")
        root1=root_og.regrid(grid,regridTool='regrid2')

        root=MV.concatenate((root,root1))

        surf_og=f("SFMC")
        surf1=surf_og.regrid(grid,regridTool='regrid2')
        surf=MV.concatenate((surf,surf1))
        f.close()
    tax=cdms.createAxis(np.arange(len(merrafiles)).tolist())
    tax.designateTime()
    tax.id="time"
    tax.units='months since 1980-01-01'
    surf.setAxis(0,tax)
    root.setAxis(0,tax)
    surf.id="SFMC"
    root.id="RZMC"
    return surf,root

def write_merra2(surf=None,root=None):
    if surf is None:
        surf,root = merra2()
    cdutil.setTimeBoundsMonthly(surf)
    cdutil.setTimeBoundsMonthly(root)
    fw = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/MERRA2_soilmoisture_summerseason.nc","w")
    
    djf_surf=cdutil.DJF(surf,criteriaarg=(1,None))[1:]
    jja_surf=cdutil.DJF(surf,criteriaarg=(1,None))[1:]
    ss_surf=summerseason_GLEAM(jja_surf,djf_surf)
    ss_surf.id="smsurf"
    fw.write(ss_surf)

    djf_root=cdutil.DJF(root,criteriaarg=(1,None))[1:]
    jja_root=cdutil.DJF(root,criteriaarg=(1,None))[1:]
    ss_root=summerseason_GLEAM(jja_root,djf_root)
    ss_root.id="smroot"
    fw.write(ss_root)
    fw.close()
def soilmoisture_projections():

    f = cdms.open("../DROUGHT_ATLAS/OBSERVATIONS/GLEAM_soilmoisture_summerseason.nc")
    root=f("smroot")
    surf=f("smsurf")

    ALL = b.DroughtAtlas("ALL_ANZDA")
    mask=ALL.obs[0].mask
    sm=b.SoilMoisture(mask)

    surfsolver=Eof(b.mask_data(sm.mma["30cm"],mask),weights='area') 
    surfmask=b.mask_data(surf,surfsolver.eofs()[0].mask)
    surfsolver=Eof(b.mask_data(sm.mma["30cm"],surfmask[-1].mask),weights='area')
    surfanom=surfmask-MV.average(surfmask,axis=0)

    time_plot(surfsolver.projectField(surfanom)[:,0],lw=3,color=cm.copper(.2))
    plt.title("Surface soil moisture")
    plt.ylabel("Projection on fingerprint")
    
    plt.figure()

    rootsolver=Eof(b.mask_data(sm.mma["2m"],mask),weights='area') 
    rootmask=b.mask_data(root,rootsolver.eofs()[0].mask)
    rootsolver=Eof(b.mask_data(sm.mma["2m"],rootmask[-1].mask),weights='area')
    rootanom=rootmask-MV.average(rootmask,axis=0)

    time_plot(rootsolver.projectField(rootanom)[:,0],lw=3,color=cm.copper(.8))
    plt.title("Root soil moisture")
    plt.ylabel("Projection on fingerprint")
