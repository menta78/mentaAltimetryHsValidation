import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools


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

    obs_ = timeSerieTidal
    obs_ = obs_ - np.nanmean(obs_)
    model_ = timeSerieModel
    model_ = model_ - np.nanmean(model_)

    # computing r2 (and other measures) gauge by gauge
    pmodel_ = np.nanpercentile(model_, pth)
    pobs_ = np.nanpercentile(obs_, pth)

    condition1_ = obs_ >= pobs_ 
    condition2_ = model_ >= model_
    condition_ = condition1_ & condition2_

    N = np.nansum(condition_)

    ssres_ = np.nansum((obs_-model_)**2, initial=0, where=condition_)
    sstot_ = np.nansum((obs_-np.nanmean(obs_))**2,  initial=0, where=condition_)
    nsc1   = np.nansum(np.abs(obs_-model_), initial = 0, where=condition_)
    nsc2   = np.nansum(np.abs(obs_-np.nanmean(obs_)), initial=0, where=condition_)

    absre_ = np.nansum(model_ - obs_, initial=0, where=condition_)
    nobs = np.nansum(obs_, initial=0, where=condition_)

    sigmaObs = np.sqrt(np.nansum((obs_-np.nanmean(obs_))**2,  initial=0, where=condition_))
    sigmaModel = np.sqrt(np.nansum((model_-np.nanmean(model_))**2,  initial=0, where=condition_))
    cov_ = np.nansum((obs_-np.nanmean(obs_))*(model_-np.nanmean(model_)), initial=0, where=condition_)

    r2_ = 1 - ssres_/sstot_
    nse_ = 1 - nsc1/nsc2
    ab_ = absre_/N
    rb_ = absre_/nobs * 100
    rmse_ = np.sqrt(ssres_/N)
    pearson_ = cov_/(sigmaModel*sigmaObs)

    r2lst.append(r2_)
    nselst.append(nse_)
    ablst.append(ab_)
    rmselst.append(rmse_)
    rblst.append(rb_)
    pearsonlst.append(pearson_)

    print(len(r2lst))
    r2 = np.mean(np.array(r2lst))
    nse = np.mean(np.array(nselst))
    ab = np.mean(np.array(ablst))
    rb = np.mean(np.array(rblst))
    rmse = np.mean(np.array(rmselst))
    pearson = np.mean(np.array(pearsonlst))

    nnse = 1/(2-nse)
    nr2 = 1/(2-r2)

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
