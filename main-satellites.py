from src.coarsenCmemsSshSatData import coarsenCmemsSshSatData
from src.interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
    interpolateModelToCoarsenedSatData_schismWWM,
)
import src.computeSshStats as computeSshStats
import src.computeStats as computeStats
from datetime import datetime

# Directory where the raw globwave files are located
#rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"

rootDir = "/home/guille/Projects/mentaAltimetryHsValidation/"

rawSatDataDir = rootDir + "data/rawData"

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
crsSatDataDir = rootDir + "data/crsSatData"

# Directory where the model nc files are located
modelNcFilesDir = rootDir + "data/schismwwm"

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = (
    rootDir + "data/satModelPairs/"
)

# Directory where the stats are generated
statsDir = rootDir + "data/stats/"

# time interval
startDate, endDate = datetime(2000, 2, 1), datetime(2000, 3, 30)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 4

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

doCoarsenSatData = False
if doCoarsenSatData:
    # coarsening the sat data.
    # This must be done because single alt observation are noisy and too numerous.
    # The data are averaged on a latitudinal tract with size latdelta.
    # You can do this once, then you can disable the coarseining, unless you want to change latdelta, or the time extent of the sat data
    latdelta = 0.5
    coarsenCmemsSshSatData(
        rawSatDataDir, crsSatDataDir, startDate, endDate, latdelta
    )

doInterpolateModelToSat = False
if doInterpolateModelToSat:
    # interpolating the model hs along the sat tracks
    interpolateModelTocoarsenCmemsSshSatData_schismWWM(
#    interpolateModelToCoarsenedSatData_schismWWM(
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
