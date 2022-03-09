import sys

sys.path.append("../altSshValidation/")
from coarsenCmemsSshSatData import coarsenCmemsSshSatData
from interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
)
import computeSshStats
from datetime import datetime

# Directory where the raw globwave files are located
#rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"
rawSatDataDir = "/media/ggarcia/TOSHIBA EXT/SATELLITES/SEALEVEL_GLO"

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
crsSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/crsSatData"

# Directory where the model nc files are located
modelNcFilesDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/outputs_final/"

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = (
    "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/satModelPairs/"
)

# Directory where the stats are generated
statsDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/stats/"

# time interval
startDate, endDate = datetime(2000, 3, 1), datetime(2000, 3, 30)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 4

# threshold above which hs should be considered
filterHsThreshold = 0.5
filterLowHs = filterHsThreshold > 0

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
    filterLowHs=filterLowHs,
    filterHsThreshold=filterHsThreshold,
    latlims=latlims,
)
