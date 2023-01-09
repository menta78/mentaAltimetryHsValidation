import os, re, glob, itertools
from datetime import datetime
import numpy as np
import h5py

from src.interpolatePredictionToScatterData import interpolateModelToScatterData

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

# number of processes to be used for the interpolation
nParWorker = 8

doInterpolateModelTidalGauge = False
if doInterpolateModelTidalGauge:
    startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)

    overwriteExisting = False
    scatterDataType = "TidalGauge_GESLA"

    pathname = os.path.join(tidalGaugeDataDir, "GESLAv1_withResiduals.mat")

    lonTidal, latTidal, resTidal, timeTidal = utils.get_serie_gesla(pathname)

    varsTidal = [lonTidal, latTidal, resTidal, timeTidal]
    varsModel = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]

    # interpolating the model ssh along tidal gauges
    interpolateModelToScatterData(
        varsTidal,
        varsModel,
        modelNcFilesDir,
        hsModelAndSatObsTidalDir,
        scatterDataType,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )


doInterpolateModelBuoys = True
if doInterpolateModelBuoys:
    startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)

    overwriteExisting = False
    scatterDataType = "Buoys_CMEMS"

    # get buoy files
    fl = []
    pthfile = os.path.join(buoysDir, "*.nc")
    fileFound = glob.glob(pthfile)
    fl.append(fileFound)
    lstBuoys = list(itertools.chain(*fl))

    lonBuoys, latBuoys, timeBuoys, hsBuoys = utils.get_hs_buoy(lstBuoys)
    print(len(lonBuoys))
    assert len(lonBuoys) == len(latBuoys) == len(timeBuoys) == len(timeBuoys)

    varsBuoys = [lonBuoys, latBuoys, hsBuoys, timeBuoys]

    varsModel = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "WWM_1", "time"]

    # interpolating the model ssh along tidal gauges
    interpolateModelToScatterData(
        varsBuoys,
        varsModel,
        modelNcFilesDir,
        hsModelAndSatObsBuoysDir,
        scatterDataType,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )

fjrifji

r2Compute = True
if r2Compute:
    elaborateMeasures(
        startDate,
        endDate,
        hsModelAndSatObsDir,
        statsDir,
        pth=pth,
    )

r2ComputePlot = False

if r2ComputePlot:
    elaborateMeasuresPlot(
        startDate,
        endDate,
        hsModelAndSatObsDir,
        statsDir,
        meanFileTidals,
        meanFileModel,
        pth=pth,
    )
