import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools


filterHighSsh = False

def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "/ERA5_schismwwm_"+strtime+"*.npy"
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    return list(itertools.chain(*fl))


def elaborateMeasures(
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

   r2lst = []
   #looping on tidal gauge files
    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:,0]
        obs_ = obs_ - np.nanmean(obs_)
        model_ = data_[:,1]
        model_ = model_ - np.nanmean(model_)

       #obs = np.concatenate([obs, obs_])
       #model = np.concatenate([model, model_])

       # computing r2 (and other measures) gauge by gauge
        pmodel_ = np.nanpercentile(model_, pth)
        pobs_ = np.nanpercentile(obs_, pth)

        condition1_ = obs_ >= pobs_ 
        condition2_ = model_ >= pobs_
        condition_ = condition1_ & condition2_

        ssres_ = np.nansum((obs_-model_)**2, where=condition_)
        sstot = np.nansum((obs_-np.nanmean(obs_))**2, where=condition_)

        r2_ = 1 - ssres_/sstot_
        r2lst.append(r2)

   r2 = np.mean(np.array(r2lst))

   #meanModel = np.nanmean(model)
   #model = model - meanModel
   #    
    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

   #condition = obs>=0

   #meanObs = np.nanmean(obs, where= condition)

    condition1 = obs >= pobs 
    condition2 = model >= pobs
    condition = condition1 & condition2

    nsc1 = np.nansum(np.abs(obs-model), where=condition)
    nsc2 = np.nansum(np.abs(obs-np.nanmean(obs)), where=condition)
    N = np.sum(condition)

    ssres = np.nansum((obs-model)**2, where=condition)
    sstot = np.nansum((obs-np.nanmean(obs))**2, where=condition)

    absre_ = np.nansum(model-obs, where=condition)
    nobs = np.nansum(obs, where=condition)

    absre = absre_/N
    nse  = 1 - nsc1/nsc2
    nnse = 1/(2-nse)
    r2   = 1 - ssres/sstot
    nr2 = 1/(2-r2)
    re = absre_/nobs
    rmse = np.sqrt(ssres/N)


    print("ssres", ssres)
    print("sstot", sstot)
    print("N", N)

    print("nse = ",nse, "nnse = ", nnse, "r2 = ", r2, "nr2 = ", nr2)
    print("rmse = ", rmse)
    print("abs bias = ", absre, "relative bias = ", re)

    with open('data/stats/tidalGauge_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".txt", 'w') as f:
        f.write("nse = " + str(nse)+"\n")
        f.write("nnse = " + str(nnse)+"\n")
        f.write("r2 = " + str(r2)+"\n")
        f.write("nr2 = " + str(nr2)+"\n")
        f.write("rmse = " + str(rmse)+"\n")
        f.write("abs bias = " + str(absre)+"\n")
        f.write("rel bias = " + str(re)+"\n")
