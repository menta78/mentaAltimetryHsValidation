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

    hhlst = []
    nbilst = []
    ablst = []
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

            nbi_, absBias_, nrmse_, hh_ = utils.computeStatsHs(
                obs_, model_, pth
            )

            nbilst.append(nbi_)
            hhlst.append(hh_)
            ablst.append(absBias_)
            nrmselst.append(nrmse_)
            continue

        idx = jIdx[i]
        idxNext = jIdx[i + 1]
        uniqueLon.append(Lonn[idx])
        uniqueLat.append(Latt[idx])
        model_ = model[idx:idxNext] - np.nanmean(model[idx:idxNext])
        obs_ = obs[idx:idxNext] - np.nanmean(obs[idx:-1])

        if len(obs_) < 10:
            continue

        nbi_, absBias_, nrmse_, hh_ = utils.computeStatsHs(
                obs_, model_, pth
            )

        nbilst.append(nbi_)
        hhlst.append(hh_)
        ablst.append(absBias_)
        nrmselst.append(nrmse_)

        print("Running Tidal Gauge", i)

    nbiT = np.array(nbilst)
    hhT = np.array(hhlst)
    absreT = np.array(ablst)
    nrmseT = np.array(nrmselst)

    uniqueLon = np.array(uniqueLon)
    uniqueLat = np.array(uniqueLat)

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
