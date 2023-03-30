import os, re, glob, csv
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
        pthfile = (
            rootdir
            + f"/{startDate.year}/TidalGauge_GESLA_station_*"
            + strtime
            + "*.npy"
        )
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    # flat_list = [item for sublist in fl for item in sublist]

    return list(itertools.chain(*fl))


def organize_strings(strings_list):
    organized_list = {}
    for string in strings_list:
        match = re.search(r"station_(\d+)", string)
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
    time_window = 1,
):
    def loadFile(flpth):
        # print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            cnd = np.logical_and(sshsat > -100, sshmdl > -100)
            satdts = satdts[cnd, :]
            cnd = np.logical_and(satdts[:, 0] < 100, satdts[:, 1] < 100)
            satdts = satdts[cnd, :]
        return satdts

    fls = getFiles(hsSatAndModelDir, startDate, endDate)
    fls_organized = organize_strings(fls)


    r2, nse, ab, rb, rmse, nrmse, pearson, lon, lat = [], [], [], [], [], [], [], [], []

    for key in fls_organized:
        fls_station = fls_organized[key]
        print(f"running skills for station {key}")
        # computing skills in each station after combining the records from different years

        #initialize time, obs and model
        obs = np.array([])
        model = np.array([])
        time = np.array([])

        for f in fls_station:
            print(f"componing files: ", f)
            data = loadFile(f)
            obs_ = data[:, 0]
            model_ = data[:, 1]
            if len(model_) > 0:
                _lon, _lat = data[0, 2], data[0, 3]
            time_ = data[:, 6]

            obs = np.concatenate([obs, obs_])
            model = np.concatenate([model, model_])
            time = np.concatenate([time, time_])

        if len(model) < nminobs:
            continue

        r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_, ref1, ref2 = utils.computeStats(
            obs, model, pth, time_window=time_window
        )
        #r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_, obs1, model1 = utils.computeStats(
        #    obs, model, pth, time_window=time_window
        #)

        print("pearson = ", pearson_)
        """
        if pearson_ < 0.5:
            fig, ax = plt.subplots()
            ax.plot(obs1, label="observation")
            ax.plot(model1, label="model", alpha=0.5)
            ax.legend()
            plt.show()
        """

        # Skip adding the result if r2_ is None (indicating empty obs or model)
        if r2_ is None:
            continue

        r2.append(r2_)
        nse.append(nse_)
        ab.append(ab_)
        rb.append(rb_)
        rmse.append(rmse_)
        nrmse.append(nrmse_)
        pearson.append(pearson_)
        lon.append(_lon)
        lat.append(_lat)

    r2array = np.array(r2)
    nsearray = np.array(nse)
    abarray = np.array(ab)
    rbarray = np.array(rb)
    rmsearray = np.array(rmse)
    nrmsearray = np.array(nrmse)
    pearsonarray = np.array(pearson)

    # Here we save the skills for each tidal gauge
    header = ["lon", "lat", "absBias", "rmse", "nrmse", "pearson"]

    outFile = os.path.join(
        outputDir,
        "TidalGauge_skills_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

    with open(outFile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for i in range(len(lon)):
            writer.writerow(
                [
                    lon[i],
                    lat[i],
                    abarray[i],
                    rmsearray[i],
                    nrmsearray[i],
                    pearsonarray[i],
                ]
            )

    r2Mean = np.mean(r2array)
    nseMean = np.mean(nsearray)
    abMean = np.mean(abarray)
    rbMean = np.mean(rbarray)
    rmseMean = np.mean(rmsearray)
    nrmseMean = np.mean(nrmsearray)
    pearsonMean = np.mean(pearsonarray)

    nnseMean = 1 / (2 - nseMean)
    nr2Mean = 1 / (2 - r2Mean)

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("nse = ", nseMean)
    print("nnse = ", nnseMean)
    print("r2 = ", r2Mean)
    print("nr2 = ", nr2Mean)
    print("rmse = ", rmseMean)
    print("nrmse = ", nrmseMean)
    print("abs bias = ", abMean)
    print("relative bias = ", rbMean)
    print("Pearson  Correlation =", pearsonMean)
    print("=========================================")

    with open(
        "data/stats/tidalGauge_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt",
        "w",
    ) as f:
        f.write("nse = " + str(nseMean) + "\n")
        f.write("nnse = " + str(nnseMean) + "\n")
        f.write("r2 = " + str(r2Mean) + "\n")
        f.write("nr2 = " + str(nr2Mean) + "\n")
        f.write("rmse = " + str(rmseMean) + "\n")
        f.write("nrmse = " + str(nrmseMean) + "\n")
        f.write("abs bias = " + str(abMean) + "\n")
        f.write("rel bias = " + str(rbMean) + "\n")
        f.write("pearson = " + str(pearsonMean) + "\n")

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
        lon,
        lat,
        s=20,
        c=nse,
        cmap="rainbow",
        # edgecolors="black",
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
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
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
        lon,
        lat,
        s=20,
        c=r2,
        cmap="rainbow",
        # edgecolors="black",
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
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
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
        c=ab,
        cmap="RdYlBu",
        # edgecolors="black",
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
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
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
        lon,
        lat,
        s=20,
        c=rb,
        cmap="RdYlBu",
        # edgecolors="black",
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
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
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
        lon,
        lat,
        s=20,
        c=pearson,
        cmap="rainbow",
        # edgecolors="black",
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
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
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
        # edgecolors="black",
        linewidths=0.5,
        vmin=0,
        vmax=100,
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
        "data/stats/tidalGauge_nrmse_newdef_"
        + "pth_"
        + str(pth)
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".png",
        **options_savefig,
    )
    plt.close(fig)