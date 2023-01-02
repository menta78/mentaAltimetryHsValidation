import os, re
from datetime import datetime
import numpy as np
import h5py

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM

# from src.computeTidalStats import elaborateMeasures
from src.computeTidalStats import elaborateMeasures
from src.plotStatsTidals import elaborateMeasuresPlot

import src.utils as utils


def get_serie_gesla(fileName):
    f = h5py.File(fileName, "r")
    data = f.get("GESELD")

    npoints = data["residual"].shape[0]  # number of stations

    res = []
    time = []
    lon = []
    lat = []

    for i in range(npoints):
        ref = f["GESELD"]["longitude"][i][0]
        lon.append(f[ref][0][0])
        ref = f["GESELD"]["latitude"][i][0]
        lat.append(f[ref][0][0])
        ref = f["GESELD"]["residual"][i][0]
        res.append(np.array(f[ref][0]))
        ref = f["GESELD"]["time"][i][0]
        time.append(np.array(f[ref][0]))
    return lon, lat, res, time


def getMeanEachTidal(fileName):
    _, lat, res, _ = get_serie_gesla(fileName)
    meanEachTidal = np.zeros(len(lat))
    i = 0
    for resTidal in res:
        meanEachTidal[i] = np.nanmean(resTidal)
        i += 1

    return meanEachTidal


rootDir = os.path.dirname(os.path.realpath(__file__))

(
    tidalGaugeDataDir,
    crsSatDataDir,
    crsWaveSatDataDir,
    modelNcFilesDir,
    hsModelAndSatObsSshDir,
    hsModelAndSatObsHsDir,
    hsModelAndSatObsTidalDir,
    statsDir,
) = utils.load_paths(rootDir)


# time interval
startDate, endDate = datetime(2002, 3, 22), datetime(2009, 12, 30)
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 31)
startDate, endDate = datetime(1995, 1, 1), datetime(2009, 12, 30)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 32

# Percentile
pth = 95

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

r2Compute = False
if r2Compute:
    elaborateMeasures(
        startDate,
        endDate,
        hsModelAndSatObsTidalDir,
        statsDir,
        pth=pth,
    )

r2ComputePlot = True

if r2ComputePlot:
    elaborateMeasuresPlot(
        startDate,
        endDate,
        hsModelAndSatObsTidalDir,
        statsDir,
        pth=pth,
        nminobs = 250,
    )
