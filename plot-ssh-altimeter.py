import os, re
from src.coarsenCmemsSshSatData import coarsenCmemsSshSatData
from src.interpolateModelToCoarsenedSatData import (
    interpolateModelTocoarsenCmemsSshSatData_schismWWM,
    interpolateModelToCoarsenedSatData_schismWWM,
)

import src.computeHsStats_test as computeHsStats
from datetime import datetime

import src.utils as utils
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colorbar import Colorbar
import matplotlib.gridspec as gridspec
import itertools
from mpl_toolkits.basemap import Basemap


rootDir = os.path.dirname(os.path.realpath(__file__))

(
    tidalGaugeDataDir,
    crsSatDataDir,
    crsWaveSatDataDir,
    modelNcFilesDir,
    hsModelAndSatObsSshDir,
    hsModelAndSatObsHsDir,
    hsModelAndSatObsTidalDir,
    statsDir,
) = utils.load_paths(rootDir)

options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
title_font = {
    "size": "12",
    "color": "black",
    "weight": "normal",
    "verticalalignment": "bottom",
}

# time interval
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
startDate, endDate = datetime(2012, 1, 1), datetime(2019, 12, 31)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
startDate, endDate = datetime(2003, 1, 1), datetime(2009, 12, 30)

pth= 95


lonFile = os.path.join(statsDir, "lons-ssh.csv")
latFile = os.path.join(statsDir, "lats-ssh.csv")

pearsonFile = os.path.join(
        statsDir,
        "SSH-pearson-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

r2File = os.path.join(
        statsDir,
        "SSH-r2-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

rmseFile = os.path.join(
        statsDir,
        "SSH-rmse-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

biasFile = os.path.join(
        statsDir,
        "SSH-bias-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

# check that data exists
assert os.path.exists(lonFile) == True
assert os.path.exists(latFile) == True
assert os.path.exists(rmseFile) == True
assert os.path.exists(r2File) == True
assert os.path.exists(pearsonFile) == True
assert os.path.exists(biasFile) == True

lon = np.genfromtxt(lonFile)
lat = np.genfromtxt(latFile)
rmse = np.genfromtxt(rmseFile)
pearson = np.genfromtxt(pearsonFile)
r2 = np.genfromtxt(r2File)
bias = np.genfromtxt(biasFile)

Latt, Lonn = np.meshgrid(lat, lon)

mask = bias == 99999
bias = np.ma.masked_array(bias, mask)
rmse = np.ma.masked_array(rmse, mask)
r2 = np.ma.masked_array(r2, mask)
pearson = np.ma.masked_array(pearson, mask)


m = Basemap()

# ==========================================================================================
# NRMSE
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    rmse, cmap="summer", origin="lower", extent=[-180, 180, -90, 90], vmin=-5, vmax=5
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("RMSE", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/SSH-rmse_"
    + "pth_"
    + str(pth)
    + "_"
    + startDate.strftime("%Y%m%d")
    + "_"
    + endDate.strftime("%Y%m%d")
    + ".png",
    **options_savefig
)
plt.close(fig)

# ==========================================================================================
# R2
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    r2, cmap="RdBu", origin="lower", extent=[-180, 180, -90, 90], vmin=0, vmax=1
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("r2", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/SSH-r2_"
    + "pth_"
    + str(pth)
    + "_"
    + startDate.strftime("%Y%m%d")
    + "_"
    + endDate.strftime("%Y%m%d")
    + ".png",
    **options_savefig
)
plt.close(fig)


# ==========================================================================================
# PEARSON
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    pearson, cmap="summer", origin="lower", extent=[-180, 180, -90, 90], vmin=0, vmax=1
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("Pearson Coefficient", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/SSH-pearson_"
    + "pth_"
    + str(pth)
    + "_"
    + startDate.strftime("%Y%m%d")
    + "_"
    + endDate.strftime("%Y%m%d")
    + ".png",
    **options_savefig
)
plt.close(fig)

# ==========================================================================================
# BIAS
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    bias, cmap="RdBu", origin="lower", extent=[-180, 180, -90, 90], vmin=-1, vmax=1
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("BIAS", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/SSH-bias_"
    + "pth_"
    + str(pth)
    + "_"
    + startDate.strftime("%Y%m%d")
    + "_"
    + endDate.strftime("%Y%m%d")
    + ".png",
    **options_savefig
)
plt.close(fig)
