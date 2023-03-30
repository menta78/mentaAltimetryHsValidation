import os, re
from datetime import datetime
import numpy as np
import h5py

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM

# from src.computeTidalStats import elaborateMeasures
from src.computeTidalStats import elaborateMeasures
from src.plotStatsTidals import elaborateMeasuresPlot

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
#startDate, endDate = datetime(2005, 1, 1), datetime(2015, 1, 1)
startDate, endDate = datetime(1995, 1, 1), datetime(2015, 1, 1)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 32

# Percentile
pth = 0

nminobs = 1

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None

elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsTidalDir,
    statsDir,
    pth=pth,
    nminobs = nminobs
)

# r2ComputePlot = True

# if r2ComputePlot:
#     elaborateMeasuresPlot(
#         startDate,
#         endDate,
#         hsModelAndSatObsTidalDir,
#         statsDir,
#         pth=pth,
#         nminobs = 250,
#     )
