import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools

import src.utils as utils


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
    pearsonlst = []

   #looping on tidal gauge files
    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:,0]
        obs_ = obs_ - np.nanmean(obs_)
        model_ = data_[:,1]
        model_ = model_ - np.nanmean(model_)
        
        stats = utils.computeStats(obs_, model_, pth)

    r2lst.append(stats["r2"])
    nselst.append(stats["NSE"])
    ablst.append(stats["Absolute_Bias"])
    rmselst.append(stats["RMSE"])
    rblst.append(stats["Relative_Bias"])
    pearsonlst.append(stats["Pearson"])

    r2 = np.mean(np.array(r2lst))
    nse = np.mean(np.array(nselst))
    ab = np.mean(np.array(ablst))
    rb = np.mean(np.array(rblst))
    rmse = np.mean(np.array(rmselst))
    pearson = np.mean(np.array(pearsonlst))

    nnse = 1/(2-nse)
    nr2 = 1/(2-r2)

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("N = ", stats["N"])
    print("nse = ",nse)
    print("nnse = ", nnse)
    print("r2 = ", r2)
    print("nr2 = ", nr2)
    print("rmse = ", rmse)
    print("abs bias = ", ab)
    print("relative bias = ", rb)
    print("Pearson  Correlation =", pearson)
    print("=========================================")


    with open('data/stats/tidalGauge_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".txt", 'w') as f:
        f.write("nse = " + str(nse)+"\n")
        f.write("nnse = " + str(nnse)+"\n")
        f.write("r2 = " + str(r2)+"\n")
        f.write("nr2 = " + str(nr2)+"\n")
        f.write("rmse = " + str(rmse)+"\n")
        f.write("abs bias = " + str(ab)+"\n")
        f.write("rel bias = " + str(rb)+"\n")
        f.write("pearson = " + str(pearson)+"\n")
