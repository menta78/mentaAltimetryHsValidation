import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib.colorbar import Colorbar
import matplotlib.gridspec as gridspec
import itertools
from mpl_toolkits.basemap import Basemap

import src.utils as utils


filterHighSsh = False
options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
title_font = {
    "size": "12",
    "color": "black",
    "weight": "normal",
    "verticalalignment": "bottom",
}


# def computeSkills(obs, model, meanTidal, meanModel, pth=99):

#     obs = obs - meanModel

#     pmodel = np.nanpercentile(model, pth)
#     pobs = np.nanpercentile(obs, pth)

#     condition = obs>=0
#     meanObs = np.nanmean(obs, where= condition)

#     condition1 = obs >= pobs
#     condition2 = model >= pmodel
#     condition = condition1 & condition2
#     #condition = condition1 | condition2

#     nsc1 = np.nansum(np.abs(obs-model), where=condition)
#     nsc2 = np.nansum(np.abs(obs-np.nanmean(obs)), where=condition)
#     N = np.sum(condition)

#     ssres = np.nansum((obs-model)**2, where=condition)
#     sstot = np.nansum((obs-np.nanmean(obs))**2, where=condition)

#     absre_ = np.nansum(model-obs, where=condition)
#     obsTot = np.nansum(obs, where=condition)

#     sigmaObs = np.sqrt(np.nansum((obs-np.nanmean(obs))**2,  initial=0, where=condition_))
#     sigmaModel = np.sqrt(np.nansum((model-np.nanmean(model))**2,  initial=0, where=condition_))
#     cov_ = np.nansum((obs_-np.nanmean(obs))*(model-np.nanmean(model)), initial=0, where=condition_)

#     absre = absre_/N
#     nse  = 1 - nsc1/nsc2
#     nnse = 1/(2-nse)
#     r2   = 1 - ssres/sstot
#     nr2 = 1/(2-r2)
#     rebias = absre_/obsTot * 100
#     pearson = cov_/(sigmaModel*sigmaObs)

#     return nnse, nr2, absre, rebias, pearson


def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "/ERA5_schismwwm_" + strtime + "*.npy"
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
    pth=90,
):
    def loadFile(flpth):
        # print("")
        # print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            # cnd = np.logical_and(
            #     sshsat < filterSshMaximum, sshmdl < filterSshMaximum
            # )
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
        obs_ = data_[:, 0]
        model_ = data_[:, 1]
        repLon = data_[:, 2]
        repLat = data_[:, 3]
        repIdx = data_[:, 4]

        obs = np.concatenate([obs, obs_])
        model = np.concatenate([model, model_])
        Lonn = np.concatenate([Lonn, repLon])
        Latt = np.concatenate([Latt, repLat])
        Indx = np.concatenate([Indx, repIdx])

    uniqueIdx, jIdx = np.unique(Indx[np.isfinite(Indx)], return_index=True)

    r2lst = []
    nselst = []
    ablst = []
    rblst = []
    rmselst = []
    nrmselst = []
    pearsonlst = []
    uniqueLon = []
    uniqueLat = []

    # computing skills for each tidal gauge considering the 95th percentile of each whole
    # time series
    for i in range(len(jIdx)):
        if i == len(jIdx) - 1:
            idx = jIdx[i]
            uniqueLon.append(Lonn[idx])
            uniqueLat.append(Latt[idx])
            model_ = model[idx:-1] - np.nanmean(model[idx:-1])
            obs_ = obs[idx:-1] - np.nanmean(obs[idx:-1])
            # nse, r2, absre, re = computeSkills(obs, model, meanTidals, meanModels, pth)
            if len(obs_) < 1:
                continue

            r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_ = utils.computeStats(
                obs_, model_, pth
            )

            r2lst.append(r2_)
            nselst.append(nse_)
            ablst.append(ab_)
            rmselst.append(rmse_)
            nrmselst.append(nrmse_)
            rblst.append(rb_)
            pearsonlst.append(pearson_)
            continue

        idx = jIdx[i]
        idxNext = jIdx[i + 1]
        uniqueLon.append(Lonn[idx])
        uniqueLat.append(Latt[idx])
        model_ = model[idx:idxNext] - np.nanmean(model[idx:idxNext])
        obs_ = obs[idx:idxNext] - np.nanmean(obs[idx:-1])

        if len(obs_) < 10:
            continue

        r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_ = utils.computeStats(
            obs_, model_, pth
        )

        r2lst.append(r2_)
        nselst.append(nse_)
        ablst.append(ab_)
        rmselst.append(rmse_)
        nrmselst.append(nrmse_)
        rblst.append(rb_)
        pearsonlst.append(pearson_)

        print("Running Tidal Gauge", i)

    nseT = np.array(nselst)
    r2T = np.array(r2lst)
    absreT = np.array(ablst)
    reT = np.array(rblst)
    rmseT = np.array(rmselst)
    nrmseT = np.array(nrmselst)
    pearsonT = np.array(pearsonlst)

    uniqueLon = np.array(uniqueLon)
    uniqueLat = np.array(uniqueLat)

    # shpfile = "/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/coastline/ne_10m_coastline.shp"
    m = Basemap()

    # ==========================================================================================
    # NSE
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt1 = axMap.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=nseT,
        cmap="RdBu",
        edgecolors="black",
        linewidths=0.5,
        vmin=0,
        vmax=1,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title(
        "Normalized Nash–Sutcliffe model efficiency coefficient", **title_font
    )

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt1, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_nse_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)

    # ==========================================================================================
    # R2
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=(1, 0.05))

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt2 = plt.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=r2T,
        cmap="RdBu",
        edgecolors="black",
        linewidths=0.5,
        vmin=0,
        vmax=1,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("Normalized R2", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt2, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_r2_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.show()
    plt.close(fig)

    # ==========================================================================================
    # Absolute Bias
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=(1, 0.05))

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt3 = plt.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=absreT,
        cmap="RdBu",
        edgecolors="black",
        linewidths=0.5,
        vmin=-0.5,
        vmax=0.5,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("Absolute Bias", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt3, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_absre_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)

    # ==========================================================================================
    # Relative Bias
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=(1, 0.05))

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt4 = plt.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=reT,
        cmap="RdBu",
        edgecolors="black",
        linewidths=0.5,
        vmin=-100,
        vmax=100,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("Relative Bias (%)", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt4, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_re_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)

    # ==========================================================================================
    # Pearson
    # ==========================================================================================
    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=(1, 0.05))

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt5 = plt.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=pearsonT,
        cmap="RdBu",
        edgecolors="black",
        linewidths=0.5,
        vmin=0,
        vmax=1,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("Pearson Correlation Coefficient", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt5, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_pearson_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)

    # ==========================================================================================
    # NRMSE
    # ==========================================================================================
    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=(1, 0.05))

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])
    m.drawcoastlines(linewidth=0.5)
    plt6 = plt.scatter(
        uniqueLon,
        uniqueLat,
        s=20,
        c=nrmselst,
        cmap="Reds",
        edgecolors="black",
        linewidths=0.5,
        vmin=0,
        vmax=1,
    )
    axMap.set(xlim=[-180, 180], ylim=[-90, 90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("NRMSE", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(
        ax=axCb, mappable=plt6, orientation="vertical"
    )  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/tidalGauge_nrmse_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)
