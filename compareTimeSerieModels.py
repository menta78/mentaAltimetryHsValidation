import csv
import math
import os
from datetime import datetime, timedelta

import h5py
import numpy as np
import scipy.io
from dotenv import load_dotenv
from matplotlib import pyplot as plt
<<<<<<< HEAD
=======
from src.readModelFast import readNcSchism
from scipy.interpolate import interp1d
>>>>>>> d15e7daae56a682931ccf688ffe618b4396fda9d

import src.utils as utils

from src.computeTidalStats_timeSeries import elaborateMeasures


dotenv_path = os.path.join(os.path.dirname(__file__), ".env")
load_dotenv(dotenv_path)

rootDir = os.environ.get("ROOT_DIR")

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")

hsModelAndSatObsDir = os.path.join(rootDir, "data/TidalModelPairs/")

statsDir = os.path.join(rootDir, "data/stats")

timSerieFl = os.path.join(
    rootDir, "data/GLOBAL_REANALISIS_timeSeries_Tomas/WL2012_2014.mat"
)

timeSeriefile = os.path.join(rootDir, "data/timeSeriesLorenzoModel.csv")


assert os.path.exists(tidalGaugeDataDir) == True
assert os.path.exists(hsModelAndSatObsDir) == True
assert os.path.exists(timSerieFl) == True

# time interval
#startDate, endDate = datetime(2013, 1, 2), datetime(2014, 12, 30)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
overwriteExisting = False

storeData=False
if storeData:
    # Get ssh from model NetCDFs and plot time serie
    extension = ".nc"
    fltPattern = "ERA5_schismwwm_"
    flsPath = utils.getFiles(modelNcFilesDir, startDate, endDate, fltPattern, extension)

    timeModel, lonModel, latModel, elev = utils.getModelVariables(flsPath)

    # Get Tomas data
    f = scipy.io.loadmat(timSerieFl)
    wl = f["WL"]
    time = f["time"]
    points = f["point2findGWL"]
    # npoint = 1
    # target = points[npoint, :]
    # print("== Target Tomas' model ==", target)

    #sample
    node = utils.find_closest_node(lonModel, latModel, points[0, :])
    sample = elev[:, node]

<<<<<<< HEAD
npoints = points.shape[0]
ntime = sample.shape[0] # all the time series will have the same time dimension but it doesnt
# coincide with the time dimension of Tomas'. It will match only if we are running the same period
allTimeSeries = np.zeros([npoints, ntime])
=======
    npoints = points.shape[0]
    ntime = sample.shape[0] # all the time series will have the same time dimension but it doesnt
    # coincide with the time dimension of Tomas'. It will match only if we are running the same period
    npoints = 3
    allTimeSeries = np.zeros([npoints, ntime])
>>>>>>> d15e7daae56a682931ccf688ffe618b4396fda9d


    for i in range(npoints):
        target = points[i, :]
        node = utils.find_closest_node(lonModel, latModel, target)
        elevTimeSerie = elev[:, node]
        allTimeSeries[i,:] = elevTimeSerie
        print(i)


    np.savetxt(timeSeriefile, allTimeSeries, delimiter=",")

    print("===== Time series stored =====")


interpolate = True

if interpolate:
    pathname = os.path.join(tidalGaugeDataDir, "GESLAv1_withResiduals.mat")
    filePoints = os.path.join(rootDir, "data/point2findGWL.csv")
    fileTimeSeries = os.path.join(rootDir, "data/timeSeriesLorenzoModel.csv")

    point2analyse = 59

    extension = ".nc"
    fltPattern = "ERA5_schismwwm_"
    flsPath = utils.getFiles(modelNcFilesDir, startDate, endDate, fltPattern, extension)

    timeModel, lonModel, latModel, elev = utils.getModelVariables(flsPath)

    points = np.genfromtxt(filePoints, delimiter=',')
    timeSeriesModel = np.genfromtxt(fileTimeSeries, delimiter=',')
    
    #target = points[point2analyse, :]
    lonTidal, latTidal, resTidal, timeTidal = utils.get_serie_gesla(pathname)

    for i in range(points.shape[0]):
        target = points[i, :]
        node = utils.find_closest_node(lonTidal, latTidal, target)
        tmstmpTidal_ = timeTidal[node][:]
        # Get timestamps recorded in the stations that belongs to the array of timestamps of the schism model
        tmstmpTidal_, isSubset, idxMin, idxMax = utils.get_subsest_list(tmstmpTidal_, min(timeModel), max(timeModel))
        if isSubset and (len(tmstmpTidal_) > 100):
            print("target = ", target)
            print(idxMin, idxMax)
            tmstmpTidal = np.asarray(tmstmpTidal_)
            break

    nodeModel = utils.find_closest_node(lonModel, latModel, target)
    timeSerieModel = elev[:, nodeModel]

    intpltr = interp1d(timeModel, timeSerieModel, bounds_error=False)

    # Model elevation interpolated at tidal time
    intp = intpltr(tmstmpTidal)

    res = np.array(resTidal[node][idxMin:idxMax])
    timeSerieTidal = np.asarray(resTidal[node][idxMin:idxMax])

    print(tmstmpTidal.shape[:], timeSerieTidal.shape[:], intp.shape[:])

    finalCSV = np.zeros((intp.shape[0],3,))
    finalCSV[:, 0] = tmstmpTidal
    finalCSV[:, 1] = timeSerieTidal
    finalCSV[:, 2] = intp

    #np.savetxt(rf"data/timeSerieInterpolated_point_{target}.csv", finalCSV, delimiter=",")

    elaborateMeasures(target, intp, timeSerieTidal, startDate, endDate, statsDir, pth = 90)


plot = False
if plot:
    options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
    title_font = {
        "size": "12",
        "color": "black",
        "weight": "normal",
        "verticalalignment": "bottom",
    }

    fig, ax = plt.subplots()
    # It's missing time array
    ax.plot(timeModel, elevTimeSerie, label="Lorenzo's Model")
    ax.plot(time, wl[npoint, :], label="Tomas' Model", alpha=0.5)
    ax.legend()
    plt.show()
    plt.savefig(
        "/home/vousdmi/Desktop/tomasVSlorenzo_lon="
        + str(target[0])
        + "_lat="
        + str(target[1])
        + ".png",
        **options_savefig
    )
    plt.close(fig)
