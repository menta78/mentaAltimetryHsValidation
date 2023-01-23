import os, re
from src.coarsenCmemsSshSatData import coarsenCmemsSshSatData
from src.interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
    interpolateModelToCoarsenedSatData_schismWWM,
)

import src.computeSshStats as computeSshStats
from datetime import datetime

import src.utils as utils


rootDir = os.path.dirname(os.path.realpath(__file__))

(
    tidalGaugeDataDir,
    buoysDir,
    crsSatDataDir,
    crsWaveSatDataDir,
    modelNcFilesDir,
    hsModelAndSatObsSshDir,
    hsModelAndSatObsHsDir,
    hsModelAndSatObsTidalDir,
    hsModelAndSatObsBuoysDir,
    statsDir,
) = utils.load_paths(rootDir)


# time interval
startDate, endDate = datetime(1995, 1, 1), datetime(2000, 1, 1)

overwriteExisting = True

pth = 95

nminobs = 250

# number of processes to be used for the interpolation
nParWorker = 8

# threshold above which hs should be considered
filterSshMaximum = 250
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

# computing the statistics
dx, dy = 1.0, 1.0
computeSshStats.maskPointsCloseToTheCoast = False

latlims = [-90, 90]
computeSshStats.elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsSshDir,
    statsDir,
    dx=dx,
    dy=dy,
    filterHighSsh=True,
    filterSshMaximum=filterSshMaximum,
    latlims=latlims,
    pth=pth,
    nminobs = nminobs
)

