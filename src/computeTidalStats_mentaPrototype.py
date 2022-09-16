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
        #print("    loading file " + flpth)
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
    nselst = []
    ablst = []
    rblst = []
    rmselst = []

   #looping on tidal gauge files
    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:,0]
        obs_ = obs_ - np.nanmean(obs_)
        model_ = data_[:,1]
        #model_ = model_ - np.nanmean(model_)

       #obs = np.concatenate([obs, obs_])
       #model = np.concatenate([model, model_])

       # computing r2 (and other measures) gauge by gauge
        pmodel_ = np.nanpercentile(model_, pth)
        pobs_ = np.nanpercentile(obs_, pth)

        if np.abs(pobs_ - pmodel_) > 0.1:
            print(pobs_, pmodel_)
            continue

        condition1_ = obs_ >= pobs_ 
        condition2_ = model_ >= model_
        condition_ = condition1_ & condition2_

        N = np.nansum(condition_)

        ssres_ = np.nansum((obs_-model_)**2, where=condition_)
        sstot_ = np.nansum((obs_-np.nanmean(obs_))**2, where=condition_)
        nsc1   = np.nansum(np.abs(obs_-model_), where=condition_)
        nsc2   = np.nansum(np.abs(obs_-np.nanmean(obs_)), where=condition_)

        absre_ = np.nansum(model_ - obs_, where=condition_)
        nobs = np.nansum(obs_, where=condition_)

        r2_ = 1 - ssres_/sstot_
        nse_ = 1 - nsc1/nsc2
        ab_ = absre_/N
        rb_ = absre_/nobs * 100
        rmse_ = np.sqrt(ssres_/N)


        r2lst.append(r2_)
        nselst.append(nse_)
        ablst.append(ab_)
        rmselst.append(rmse_)
        rblst.append(rb_)

    print(len(r2lst))
    r2 = np.mean(np.array(r2lst))
    nse = np.mean(np.array(nselst))
    ab = np.mean(np.array(ablst))
    rb = np.mean(np.array(rblst))
    rmse = np.mean(np.array(rmselst))

    nnse = 1/(2-nse)
    nr2 = 1/(2-r2)

    print("nse = ",nse, "nnse = ", nnse, "r2 = ", r2, "nr2 = ", nr2)
    print("rmse = ", rmse)
    print("abs bias = ", ab, "relative bias = ", rb)

    with open('data/stats/tidalGauge_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".txt", 'w') as f:
        f.write("nse = " + str(nse)+"\n")
        f.write("nnse = " + str(nnse)+"\n")
        f.write("r2 = " + str(r2)+"\n")
        f.write("nr2 = " + str(nr2)+"\n")
        f.write("rmse = " + str(rmse)+"\n")
        f.write("abs bias = " + str(ab)+"\n")
        f.write("rel bias = " + str(rb)+"\n")
