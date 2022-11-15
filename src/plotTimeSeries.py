import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib.colorbar import Colorbar
import matplotlib.gridspec as gridspec
import itertools
from mpl_toolkits.basemap import Basemap


filterHighSsh = False
options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
title_font = {
    "size": "12",
    "color": "black",
    "weight": "normal",
    "verticalalignment": "bottom",
}


def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return array[idx]


def find_by_index(array, value):
    value_ = find_nearest(array, value)
    index = np.where(array == value_)
    return index


def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "/ERA5_schismwwm_" + strtime + "*.npy"
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    return list(itertools.chain(*fl))


def timeSeriesPlot(
    startDate,
    endDate,
    hsSatAndModelDir,
    outputDir,
    target,
    meanFileTidals,
    meanModelFile,
    pth=90,
):
    def loadFile(flpth):
        print("")
        print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            cnd = np.logical_and(sshsat < filterSshMaximum, sshmdl < filterSshMaximum)
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

    # load mean of each tidal
    print(meanFileTidals)
    meanTidals = np.load(meanFileTidals)

    # load mean of each node in the model
    meanModels = 0

    uniqueIdx, jIdx = np.unique(Indx[np.isfinite(Indx)], return_index=True)

    # meanModel = np.nanmean(model)
    # model = model - meanModel

    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

    nseT = []
    r2T = []
    absreT = []
    reT = []
    uniqueLon = []
    uniqueLat = []

    lon, lat = target[0], target[1]
    indexLon = find_by_index(Lonn, lon)
    indexLat = find_by_index(Latt, lat)
    index = np.intersect1d(indexLon, indexLat)
    print(index)
    obsTarget = obs[index]
    modelTarget = model[index]

    fig, ax = plt.subplots()
    # It's missing time array
    ax.plot(obsTarget, label="observation")
    ax.plot(modelTarget, label="model", alpha=0.5)
    ax.legend()
    plt.savefig(
        "/home/vousdmi/Desktop/tidalVSmodel_lon="
        + str(lon)
        + "_lat="
        + str(lat)
        + ".png",
        **options_savefig
    )
    plt.close(fig)
