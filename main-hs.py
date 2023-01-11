import os, re
from src.coarsenSatData import coarsenSatData
from src.interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
    interpolateModelToCoarsenedSatData_schismWWM,
)
import src.computeSshStats as computeSshStats
import src.computeStats as computeStats
from datetime import datetime

import src.utils as utils
# Directory where the raw globwave files are located
#rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"

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
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
startDate, endDate = datetime(2012, 1, 1), datetime(2019, 12, 31)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
overwriteExisting = True

# number of processes to be used for the interpolation
nParWorker = 8

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

doCoarsenSatData = True
if doCoarsenSatData:
    # coarsening the sat data.
    # This must be done because single alt observation are noisy and too numerous.
    # The data are averaged on a latitudinal tract with size latdelta.
    # You can do this once, then you can disable the coarseining, unless you want to change latdelta, or the time extent of the sat data
    rawSatDataDir = os.path.join(rootDir, "data/rawDataSWH")
    startDate, endDate = datetime(2003, 1, 1), datetime(2003, 12, 31)
    assert os.path.exists(rawSatDataDir) == True
    latdelta = 0.5
    coarsenSatData(
        rawSatDataDir, crsSatDataDir, startDate, endDate, latdelta
    )
fjrijir


doInterpolateModelToSat = True
if doInterpolateModelToSat:
    # interpolating the model hs along the sat tracks
    interpolateModelToCoarsenedSatData_schismWWM(
        crsSatDataDir,
        modelNcFilesDir,
        hsModelAndSatObsDir,
        boundaries,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )


# computing the statistics
dx, dy = 1.0, 1.0
computeSshStats.maskPointsCloseToTheCoast = False

latlims = [-90, 90]
computeSshStats.elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsDir,
    statsDir,
    dx=dx,
    dy=dy,
    filterHighSsh=True,
    filterSshMaximum=filterSshMaximum,
    latlims=latlims,
    pth = 95,
)
