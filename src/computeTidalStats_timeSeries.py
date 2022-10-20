import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools

import src.utils as utils

def elaborateMeasures(
    target,
    timeSerieModel,
    timeSerieTidal,
    startDate,
    endDate,
    outputDir,
    pth = 90,
):

    r2lst = []
    nselst = []
    ablst = []
    rblst = []
    rmselst = []
    pearsonlst = []

    stats = utils.computeStats(timeSerieTidal, timeSerieModel, pth)

    r2lst.append(stats["r2"])
    nselst.append(stats["NSE"])
    ablst.append(stats["Absolute_Bias"])
    rmselst.append(stats["RMSE"])
    rblst.append(stats["Relative_Bias"])
    pearsonlst.append(stats["Pearson"])

    print(len(r2lst))
    r2 = np.mean(np.array(r2lst))
    nse = np.mean(np.array(nselst))
    ab = np.mean(np.array(ablst))
    rb = np.mean(np.array(rblst))
    rmse = np.mean(np.array(rmselst))
    pearson = np.mean(np.array(pearsonlst))

    nnse = 1/(2-nse)
    nr2 = 1/(2-r2)

    N = stats["N"]

    print("=========================================")
    print("N = ", N)
    print("nse = ",nse)
    print("nnse = ", nnse)
    print("r2 = ", r2)
    print("nr2 = ", nr2)
    print("rmse = ", rmse)
    print("abs bias = ", ab)
    print("relative bias = ", rb)
    print("Pearson  Correlation =", pearson)
    print("=========================================")

    fileStats = os.path.join(outputDir,'tidalGauge_pth_'+str(pth)+"_"+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".txt")

    with open(fileStats, 'w') as f:
        f.write("nse = " + str(nse)+"\n")
        f.write("nnse = " + str(nnse)+"\n")
        f.write("r2 = " + str(r2)+"\n")
        f.write("nr2 = " + str(nr2)+"\n")
        f.write("rmse = " + str(rmse)+"\n")
        f.write("abs bias = " + str(ab)+"\n")
        f.write("rel bias = " + str(rb)+"\n")
        f.write("rel bias = " + str(pearson)+"\n")
