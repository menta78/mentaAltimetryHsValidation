import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools
from mpl_toolkits.basemap import Basemap


filterHighSsh = False
options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}


def computeSkills(obs, model, pth=99):

    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

    condition = obs>=0
    meanObs = np.nanmean(obs, where= condition)

    condition1 = obs >= pobs 
    condition2 = model >= pobs
    condition = condition1 & condition2

    nsc1 = np.nansum(np.abs(obs-model), where=condition)
    nsc2 = np.nansum(np.abs(obs-np.nanmean(obs)), where=condition)

    ssres = np.nansum((obs-model)**2, where=condition)
    sstot = np.nansum((obs-np.nanmean(obs))**2, where=condition)

    absre = np.nansum(model-obs, where=condition)
    nobs = np.nansum(obs, where=condition)

    nse  = 1 - nsc1/nsc2
    nnse = 1/(2-nse)
    r2   = 1 - ssres/sstot
    nr2 = 1/(2-r2)
    re = absre/nobs

    return nse, r2, absre, re


def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "/ERA5_schismwwm_"+strtime+"*.npy"
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    return list(itertools.chain(*fl))


def elaborateMeasuresPlot(
    startDate,
    endDate,
    hsSatAndModelDir,
    outputDir,
    filterLowHs=False,
    filterHsThreshold=0.0,
    pth = 90,
):

    def loadFile(flpth):
        print("")
        print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            cnd = np.logical_and(
                sshsat < filterSshMaximum, sshmdl < filterSshMaximum
            )
            satdts = satdts[cnd, :]
        return satdts

    fls = getFiles(hsSatAndModelDir, startDate, endDate)

    obs = np.array([])
    model = np.array([])
    Lonn = model
    Latt = model
    Indx = model

    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:,0]
        model_ = data_[:,1]
        repLon = data_[:,2]
        repLat = data_[:,3]
        repIdx = data_[:,4] 

        obs = np.concatenate([obs, obs_])
        model = np.concatenate([model, model_])
        Lonn = np.concatenate([Lonn, repLon])
        Latt = np.concatenate([Latt, repLat])
        Indx = np.concatenate([Indx, repIdx])
        
    
    uniqueIdx, jIdx = np.unique(Indx[np.isfinite(Indx)], return_index=True)

    meanModel = np.nanmean(model)
    model = model - meanModel
    obs = obs - np.nanmean(obs)
        
    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

    nseT = []
    r2T = []
    absreT = []
    reT = []
    uniqueLon = []
    uniqueLat = []    

    for i in range(len(jIdx)):
        if i == len(jIdx)-1:
            idx = jIdx[i]
            uniqueLon.append(Lonn[idx])
            uniqueLat.append(Latt[idx])
            model_ = model[idx:-1]
            obs_ = obs[idx:-1]
            nse, r2, absre, re = computeSkills(obs, model, pth)
            nseT.append(nse)
            r2T.append(r2)
            absreT.append(absre)
            reT.append(re)
            continue

        idx = jIdx[i]
        idxNext = jIdx[i+1]
        uniqueLon.append(Lonn[idx])
        uniqueLat.append(Latt[idx])
        model_ = model[idx:idxNext]
        obs_ = obs[idx:idxNext]

        nse, r2, absre, re = computeSkills(obs_, model_, pth)
        nseT.append(nse)
        r2T.append(r2)
        absreT.append(absre)
        reT.append(re)
        # print("nse  = ", nse)
        # print("r2  = ", r2)
        # print("absre  = ", absre)
        # print("re  = ", re)
        print("Running Tidal Gauge", i)

    nseT = np.array(nseT)
    r2T = np.array(r2T)
    absreT = np.array(absreT)
    reT = np.array(reT)
    uniqueLon = np.array(uniqueLon)
    uniqueLat = np.array(uniqueLat)    

    shpfile = "/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/coastline/ne_10m_coastline.shp"
    m=Basemap()

    fig, ax = plt.subplots()
    #mask = geopandas.read_file(shpfile)
    m.drawcoastlines(linewidth=0.5)
    plt.scatter(uniqueLon, uniqueLat, s=20, c=nseT, cmap="Blues")
    plt.colorbar()
    plt.clim(0, 1)
    ax.set(xlim=[-180,180], ylim=[-90,90])
    ax.set_aspect("equal", "box")
    plt.title("Nash–Sutcliffe model efficiency coefficient")
    plt.savefig('data/stats/tidalGauge_nse_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
    plt.close(fig)

    fig, ax = plt.subplots()
    #mask = geopandas.read_file(shpfile)
    m.drawcoastlines(linewidth=0.5)
    plt.scatter(uniqueLon, uniqueLat, s=20, c=r2T, cmap="Blues")
    plt.colorbar()
    plt.clim(0, 1)
    ax.set(xlim=[-180,180], ylim=[-90,90])
    ax.set_aspect("equal", "box")
    plt.title("R2")
    plt.savefig('data/stats/tidalGauge_r2_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
    plt.close(fig)

    fig, ax = plt.subplots()
    #mask = geopandas.read_file(shpfile)
    m.drawcoastlines(linewidth=0.5)
    plt.scatter(uniqueLon, uniqueLat, s=20, c=absreT, cmap="Blues")
    plt.colorbar()
    plt.clim(-50, 50)
    ax.set(xlim=[-180,180], ylim=[-90,90])
    ax.set_aspect("equal", "box")
    plt.title("Absolute Bias")
    plt.savefig('data/stats/tidalGauge_absre_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
    plt.close(fig)

    fig, ax = plt.subplots()
    #mask = geopandas.read_file(shpfile)
    m.drawcoastlines(linewidth=0.5)
    plt.scatter(uniqueLon, uniqueLat, s=20, c=reT, cmap="Blues")
    plt.colorbar()
    plt.clim(-1, 1)
    ax.set(xlim=[-180,180], ylim=[-90,90])
    ax.set_aspect("equal", "box")
    plt.title("Relative Bias")
    plt.savefig('data/stats/tidalGauge_re_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
    plt.close(fig)

