import os, re
from src.coarsenCmemsSshSatData import coarsenCmemsSshSatData
from src.interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
    interpolateModelToCoarsenedSatData_schismWWM,
)

import src.computeHsStats_test as computeHsStats
from datetime import datetime

import src.utils as utils


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
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
startDate, endDate = datetime(2012, 1, 1), datetime(2019, 12, 31)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
startDate, endDate = datetime(2003, 1, 1), datetime(2009, 12, 30)

overwriteExisting = True

pth = 0

# number of processes to be used for the interpolation
nParWorker = 8

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

# computing the statistics
dx, dy = 2.0, 2.0
computeHsStats.maskPointsCloseToTheCoast = False

latlims = [-90, 90]
computeHsStats.elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsHsDir,
    statsDir,
    dx=dx,
    dy=dy,
    filterHighSsh=True,
    filterSshMaximum=filterSshMaximum,
    latlims=latlims,
    pth = pth,
)
