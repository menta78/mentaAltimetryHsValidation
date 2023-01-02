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

    nbi, absBias, nrmse, hh = utils.computeStatsHs(obs, model, pth)

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("nbi = ", nbi)
    print("absolute bias = ", absBias)
    print("nrmse = ", nrmse)
    print("hh = ", hh)
    print("=========================================")

    with open(
        "data/stats/buoys_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt",
        "w",
    ) as f:
        f.write("nbi = " + str(nbi) + "\n")
        f.write("absBias = " + str(absBias) + "\n")
        f.write("nrmse = " + str(nrmse) + "\n")
        f.write("hh = " + str(hh) + "\n")

