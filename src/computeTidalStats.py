import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools

import src.utils_test as utils

filterHighSsh = False


def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "/ERA5_schismwwm_" + strtime + "*.npy"
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
    pth=90,
):
    def loadFile(flpth):
        # print("    loading file " + flpth)
        satdts = np.load(flpth)
        # if filterHighSsh:
        #     sshsat = satdts[:, 0]
        #     sshmdl = satdts[:, 1]
        #     cnd = np.logical_and(
        #         sshsat < filterSshMaximum, sshmdl < filterSshMaximum
        #     )
        #     satdts = satdts[cnd, :]
        return satdts

    fls = getFiles(hsSatAndModelDir, startDate, endDate)

    obs = np.array([])
    model = np.array([])

    r2lst = []
    nselst = []
    ablst = []
    rblst = []
    rmselst = []
    nrmselst = []
    pearsonlst = []

    # looping on tidal gauge files
    obs = np.array([])
    model = np.array([])

    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:, 0]
        model_ = data_[:, 1]

        model = np.concatenate((model, model_), axis=0)
        obs = np.concatenate((obs, obs_), axis=0)

        if len(obs) < 1:
            continue
        r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_ = utils.computeStats(obs, model, pth)

        r2lst.append(r2_)
        nselst.append(nse_)
        ablst.append(ab_)
        rmselst.append(rmse_)
        nrmselst.append(nrmse_)
        rblst.append(rb_)
        pearsonlst.append(pearson_)


    r2array = np.array(r2lst)
    nsearray = np.array(nselst)
    abarray = np.array(ablst)
    rbarray = np.array(rblst)
    rmsearray = np.array(rmselst)
    nrmsearray = np.array(nrmselst)
    pearsonarray = np.array(pearsonlst)
    r2 = np.mean(r2array)
    nse = np.mean(nsearray)
    ab = np.mean(abarray)
    rb = np.mean(rbarray)
    rmse = np.mean(rmsearray)
    nrmse = np.mean(nrmsearray)
    pearson = np.mean(pearsonarray)

    nnse = 1 / (2 - nse)
    nr2 = 1 / (2 - r2)

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("nse = ", nse)
    print("nnse = ", nnse)
    print("r2 = ", r2)
    print("nr2 = ", nr2)
    print("rmse = ", rmse)
    print("nrmse = ", nrmse)
    print("abs bias = ", ab)
    print("relative bias = ", rb)
    print("Pearson  Correlation =", pearson)
    print("=========================================")

    with open(
        "data/stats/tidalGauge_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt",
        "w",
    ) as f:
        f.write("nse = " + str(nse) + "\n")
        f.write("nnse = " + str(nnse) + "\n")
        f.write("r2 = " + str(r2) + "\n")
        f.write("nr2 = " + str(nr2) + "\n")
        f.write("rmse = " + str(rmse) + "\n")
        f.write("nrmse = " + str(nrmse) + "\n")
        f.write("abs bias = " + str(ab) + "\n")
        f.write("rel bias = " + str(rb) + "\n")
        f.write("pearson = " + str(pearson) + "\n")
