import sys
sys.path.append('../altHsValidation/')
from coarsenSatData import coarsenSatData
from interpolateModelToCoarsenedSatData import interpolateModelToCoarsenedSatData_schismWWM
import computeStats
from datetime import datetime

# Directory where the raw globwave files are located
rawSatDataDir = '/mnt/d/DATA/CMEMS-swh'

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
crsSatDataDir = '/mnt/d/DATA/CMEMS-swh/crsSatData'

# Directory where the model nc files are located
modelNcFilesDir = '/work/opa/lm09621/work/WWMV/schismWwmIce/run_experiment_spart/outputs_final/'

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = '/work/opa/lm09621/data/waves/satData/models/schismWwm_spart/202108/satModelPairs/'

# Directory where the stats are generated
statsDir = '/work/opa/lm09621/data/waves/satData/models/schismWwm_spart/202108/stats/'

# time interval 
startDate, endDate = datetime(2000, 2, 1), datetime(2025, 3, 31)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 4

# threshold above which hs should be considered
filterHsThreshold = .5
filterLowHs = filterHsThreshold > 0

# set this if you need to limit your analysis to a subdomain
boundaries = None

doCoarsenSatData = True
if doCoarsenSatData:
 # coarsening the sat data.
 # This must be done because single alt observation are noisy and too numerous.
 # The data are averaged on a latitudinal tract with size latdelta.
 # You can do this once, then you can disable the coarseining, unless you want to change latdelta, or the time extent of the sat data
  latdelta = 0.4
  coarsenSatData(rawSatDataDir, crsSatDataDir, startDate, endDate, latdelta)

exit()

doInterpolateModelToSat = True
if doInterpolateModelToSat:
 # interpolating the model hs along the sat tracks
  interpolateModelToCoarsenedSatData_schismWWM(crsSatDataDir, modelNcFilesDir, hsModelAndSatObsDir, boundaries, startDate, endDate, overwriteExisting=overwriteExisting, nParWorker=nParWorker)

# computing the statistics
dx, dy = 1., 1.
computeStats.maskPointsCloseToTheCoast = False
latlims = [-90, 90]
computeStats.elaborateMeasures(startDate, endDate, hsModelAndSatObsDir, statsDir, 
  dx=dx, dy=dy, filterLowHs=filterLowHs, filterHsThreshold=filterHsThreshold, latlims=latlims)


