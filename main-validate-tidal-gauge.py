import os, re
from datetime import datetime
import numpy as np
import h5py

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM
# from src.computeTidalStats import elaborateMeasures
from src.computeTidalStats_mentaPrototype import elaborateMeasures
from src.plotStatsTidals import elaborateMeasuresPlot


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

def getMeanEachTidal(fileName):
    _, lat, res, _ = get_serie_gesla(fileName)
    meanEachTidal = np.zeros(len(lat))
    i = 0
    for resTidal in res:
        meanEachTidal[i] = np.nanmean(resTidal)
        i += 1

    return meanEachTidal


# Directory where the raw globwave files are located
# rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"

rootDir = "/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation"

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = os.path.join(rootDir, "data/TidalModelPairs/")
hsModelAndSatObsDir = os.path.join(rootDir, "data/tidalInterpolations/")
assert os.path.exists(hsModelAndSatObsDir) == True

# Directory where the stats are generated
statsDir = os.path.join(rootDir, "data/stats")
#statsDir = "../../experiments/jrc/stats/"

# Directory to Tidal Gauge
pathname = os.path.join(tidalGaugeDataDir, "GESLAv1_withResiduals.mat")

meanFileTidals = os.path.join(tidalGaugeDataDir, "meanResiduals.npy")
meanModelFile = os.path.join(rootDir, "data/elev/elevmean.nc")


# time interval
startDate, endDate = datetime(2002, 3, 22), datetime(2009, 12, 30)
startDate, endDate = datetime(2006, 12, 20), datetime(2007, 12, 29)
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
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


doGetMeanTidal = False
if doGetMeanTidal:
    meanEachTidal = getMeanEachTidal(pathname)
    print(meanEachTidal.shape[:])
    np.save(meanFileTidals, meanEachTidal)


doInterpolateModelToSat = 0
if doInterpolateModelToSat:
    lonTidal, latTidal, resTidal, timeTidal = get_serie_gesla(pathname)

    varsTidal = [lonTidal, latTidal, resTidal, timeTidal]
    # interpolating the model ssh along tidal gauges
    interpolateModelToTidalGauge_schismWWM(
        varsTidal,
        modelNcFilesDir,
        hsModelAndSatObsDir,
        meanModelFile,
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
        hsModelAndSatObsDir,
        statsDir,
        pth = pth,
    ) 



r2ComputePlot = False

if r2ComputePlot:
    elaborateMeasuresPlot(
        startDate,
        endDate,
        hsModelAndSatObsDir,
        statsDir,
        pth = pth,
) 
