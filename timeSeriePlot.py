import math
import os
from datetime import datetime, timedelta

import h5py
import numpy as np
import scipy.io
from dotenv import load_dotenv
from matplotlib import pyplot as plt

import src.utils as utils

dotenv_path = os.path.join(os.path.dirname(__file__), ".env")
load_dotenv(dotenv_path)

rootDir = os.environ.get("ROOT_DIR")

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")

hsModelAndSatObsDir = os.path.join(rootDir, "data/TidalModelPairs/")

# time interval
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
overwriteExisting = False

# Get ssh from model NetCDFs and plot time serie
extension = ".nc"
fltPattern = "ERA5_schismwwm_"
flsPath = utils.getFiles(modelNcFilesDir, startDate, endDate, fltPattern, extension)

timeModel, lonModel, latModel, elev = utils.getModelVariables(flsPath)


# Get Tomas data
timSerieFl = "/mnt/c/Users/ggarc/OneDrive/Documents/experiments/GLOBAL_REANALISIS_timeSeries/WL2012_2014.mat"
f = scipy.io.loadmat(timSerieFl)
wl = f["WL"]
time = f["time"]
points = f["point2findGWL"]
npoint = 11
target = points[npoint, :]
print(target)

node = utils.find_closest_node(lonModel, latModel, target)
print(lonModel[node], latModel[node])




# Elevation Time Serie
elevTimeSerie = elev[:, node]
print(elevTimeSerie.shape[:])

plot = True
if plot:
    pathname = "data/tidalGaugeData/GESLAv1_withResiduals.mat"
    timSerieFl = "data/WL2012_2014.mat"
    
    plt.plot(elevTimeSerie)
    plt.plot(wl[npoint, :], alpha=0.4)
    plt.show()
