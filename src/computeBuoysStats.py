import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib.colorbar import Colorbar
import matplotlib.gridspec as gridspec
import itertools
from mpl_toolkits.basemap import Basemap
import itertools

import src.utils as utils

filterHighSsh = True
options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
title_font = {
    "size": "12",
    "color": "black",
    "weight": "normal",
    "verticalalignment": "bottom",
}


def getFiles(rootdir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = rootdir + f"/{startDate.year}/Buoys_*" + strtime + "*.npy"
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    #flat_list = [item for sublist in fl for item in sublist]

    return list(itertools.chain(*fl))

def organize_strings(strings_list):
    organized_list = {}
    for string in strings_list:
        match = re.search(r'station_(\d+)', string)
        if match:
            station_number = int(match.group(1))
            if station_number not in organized_list:
                organized_list[station_number] = []
            organized_list[station_number].append(string)
    return organized_list


def elaborateMeasures(
    startDate,
    endDate,
    hsSatAndModelDir,
    outputDir,
    filterLowHs=False,
    filterHsThreshold=0.0,
    pth=90,
    nminobs=1,
):
    def loadFile(flpth):
        # print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            cnd = np.logical_and(
                sshsat > filterHsThreshold, sshmdl > filterHsThreshold
            )
            satdts = satdts[cnd, :]
        return satdts

    fls = getFiles(hsSatAndModelDir, startDate, endDate)
    fls_organized = organize_strings(fls)

    obs = np.array([])
    model = np.array([])
    time = np.array([])
    # Lonn = np.array([])
    # Latt = np.array([])
    # Indx = np.array([])
    # Node = np.array([])

    nbi, absBias, nrmse, hh, lon, lat = [], [], [], [], [], []

    for key in fls_organized:
        fls_station = fls_organized[key]
        # computing skills in each station after combining the records from different years
        for f in fls_station:
            data = loadFile(f)
            obs_   = data[:, 0]
            model_ = data[:, 1]
            if len(model_) > 0:
                _lon, _lat = data[0,2], data[0,3]
            # repLon = data[:, 2]
            # repLat = data[:, 3]
            # repIdx = data[:, 4]
            # repNode = data[:,5]
            time_ = data[:,6]


            obs = np.concatenate([obs, obs_])
            model = np.concatenate([model, model_])
            # Lonn = np.concatenate([Lonn, repLon])
            # Latt = np.concatenate([Latt, repLat])
            # Indx = np.concatenate([Indx, repIdx])
            # Node = np.concatenate([Node, repNode])
            time = np.concatenate([time, time_])

            # _obs = np.zeros([obs.shape[0], 1])
            # _model = np.zeros([obs.shape[0], 1]) 
            # # _lon = np.zeros([obs.shape[0], 1])
            # # _lat = np.zeros([obs.shape[0], 1]) 
            # # _idx = np.zeros([obs.shape[0], 1])
            # # _node = np.zeros([obs.shape[0], 1]) 
            # _time = np.zeros([obs.shape[0], 1]) 

            # _obs[:, 0] = obs
            # _model[:, 0] = model
            # _lon[:, 0] = Lonn
            # _lat[:, 0] = Latt
            # _idx[:, 0] = Indx
            # _node[:, 0] = Node
            # _time[:, 0] = time

            if len(model) < nminobs:
                continue
            _nbi, _absBias, _nrmse, _hh = utils.computeStatsHs(obs, model, pth)
            nbi.append(_nbi)
            absBias.append(_absBias)
            nrmse.append(_nrmse)
            hh.append(_hh)
            lon.append(_lon)
            lat.append(_lat)    


    print(nbi)
    print(absBias)

    nbiMean = np.mean(np.array(nbi))
    hhMean = np.mean(np.array(hh))
    absBiasMean = np.mean(np.array(absBias))
    nrmseMean = np.mean(np.array(nrmse))

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("nbi = ", nbiMean)
    print("absolute bias = ", absBiasMean)
    print("nrmse = ", nrmseMean)
    print("hh = ", hhMean)
    print("=========================================")

    with open(
        "data/stats/buoys_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt",
        "w",
    ) as f:
        f.write("nbi = " + str(nbiMean) + "\n")
        f.write("absBias = " + str(absBiasMean) + "\n")
        f.write("nrmse = " + str(nrmseMean) + "\n")
        f.write("hh = " + str(hhMean) + "\n")


    m = Basemap()

    # ==========================================================================================
    # NBI
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])

    plt1 = axMap.scatter(
        lon,
        lat,
        s=20,
        c=nbi,
        cmap="RdYlBu",
        ##edgecolors="black",
        vmin=-1, 
        vmax=1,
    )

    m.drawcoastlines(linewidth=0.5)
    #m.fillcontinents(color="gray")
    # axMap.set(xlim=[-180,180], ylim=[-90,90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("NBI", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/buoys-hs-nbi_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)

    # ==========================================================================================
    # HH
    # ==========================================================================================

    fig, ax = plt.subplots(figsize=(9, 4))
    grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

    # Map and scatter plot
    axMap = plt.subplot(grd[0, 0])

    plt1 = axMap.scatter(
        lon,
        lat,
        s=20,
        c=hh,
        cmap="rainbow",
        #edgecolors="black",
        vmin=0, 
        vmax=1,
    )
    m.drawcoastlines(linewidth=0.5)
    #m.fillcontinents(color="gray")
    # axMap.set(xlim=[-180,180], ylim=[-90,90])
    axMap.set_aspect("equal", "box")
    axMap.set_title("HH", **title_font)

    # Colorbar
    axCb = plt.subplot(grd[0, 1])
    cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

    plt.savefig(
        "data/stats/buoys-hs-hh_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
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
        lon,
        lat,
        s=20,
        c=absBias,
        cmap="RdYlBu",
        #edgecolors="black",
        linewidths=0.5,
        vmin=-1,
        vmax=1,
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
        "data/stats/buoys-hs-absre_"
        + "pth_"
        + str(pth)
        + "_"
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
        lon,
        lat,
        s=20,
        c=nrmse,
        cmap="rainbow",
        #edgecolors="black",
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
        "data/stats/buoys-hs-nrmse_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig
    )
    plt.close(fig)
