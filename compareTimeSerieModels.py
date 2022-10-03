import csv
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

timSerieFl = os.path.join(
    rootDir, "data/GLOBAL_REANALISIS_timeSeries_Tomas/WL2012_2014.mat"
)


assert os.path.exists(tidalGaugeDataDir) == True
assert os.path.exists(hsModelAndSatObsDir) == True
assert os.path.exists(timSerieFl) == True

# time interval
startDate, endDate = datetime(2013, 1, 2), datetime(2014, 12, 30)
startDate, endDate = datetime(2013, 1, 2), datetime(2013, 1, 10)
overwriteExisting = False

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


for i in range(points.shape[0]):
    target = points[i, :]
    node = utils.find_closest_node(lonModel, latModel, target)
    elevTimeSerie = elev[:, node]
    with open("data/timeSeriesLorenzoModel.csv", "a", encoding="UTF8") as f:
        # create the csv writer
        writer = csv.writer(f)
        # write a row to the csv file
        writer.writerow(elevTimeSerie)

exit


# Elevation Time Serie
elevTimeSerie = elev[:, node]

plot = True
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
