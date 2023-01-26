import os, re
from datetime import datetime
import numpy as np
import numpy.ma as ma
import gzip
import glob
import netCDF4
from datetime import datetime, timedelta

from matplotlib import pyplot as plt
import geopandas

from src.interpolateModelToBouy import interpolateModelToTidalGauge_schismWWM
from src.computeBuoysStats import elaborateMeasures
from src.plotStatsBuoys import elaborateMeasuresPlot
import itertools

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

startDate, endDate = datetime(1985, 1, 1), datetime(2015, 1, 1)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 32

nminobs = 1

# Percentile
pth = 0

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

hsModelAndSatObsBuoysDir = os.path.join(rootDir, "data/buoyModelPairs_before2015")
assert os.path.exists(hsModelAndSatObsBuoysDir) == True



elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsBuoysDir,
    statsDir,
    pth = pth,
    nminobs = nminobs
) 


# r2ComputePlot = True
# if r2ComputePlot:
#     elaborateMeasuresPlot(
#         startDate,
#         endDate,
#         hsModelAndSatObsBuoysDir,
#         statsDir,
#         None,
#         None,
#         pth = pth,
#         nminobs = 250,
# ) 
