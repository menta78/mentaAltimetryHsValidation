import os, re
from datetime import datetime
import numpy as np
import h5py

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM
from src.computeTidalStats import elaborateMeasures

def get_serie_gesla(fileName):
    f = h5py.File(fileName,'r')
    data = f.get("GESELD")

    npoints = data["residual"].shape[0] # number of stations

    res   = []
    time  = []
    lon   = []
    lat   = []

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


# Directory where the raw globwave files are located
# rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"

rootDir = "/mnt/c/Users/ggarc/OneDrive/Documents/Projects/mentaAltimetryHsValidation"

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = os.path.join(rootDir, "data/satModelPairs/")

# Directory where the stats are generated
#statsDir = os.path.join(rootDir, "data/stats/")
statsDir = "../../experiments/jrc/stats/"

# time interval
startDate, endDate = datetime(2000, 3, 28), datetime(2000, 3, 30)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 4

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None


filterTidalGaugeData = False

if filterTidalGaugeData:
    pathname = os.path.join(tidalGaugeDataDir, "GESLAv1_withResiduals.mat")
    lonTidal, latTidal, resTidal, timeTidal = get_serie_gesla(pathname)


doInterpolateModelToSat = False
if doInterpolateModelToSat:
    varsTidal = [lonTidal, latTidal, resTidal, timeTidal]
    # interpolating the model ssh along tidal gauges
    interpolateModelToTidalGauge_schismWWM(
        varsTidal,
        modelNcFilesDir,
        statsDir,
        boundaries,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )


r2Compute = True
if r2Compute:
    elaborateMeasures(
        startDate,
        endDate,
        statsDir,
        pth = 95,
    ) 
