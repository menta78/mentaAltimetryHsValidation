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


options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
title_font = {
    "size": "12",
    "color": "black",
    "weight": "normal",
    "verticalalignment": "bottom",
}

# time interval
startDate, endDate = datetime(2002, 1, 1), datetime(2005, 1, 1)

pth = 95


lonFile = os.path.join(statsDir, "lons.csv")
latFile = os.path.join(statsDir, "lats.csv")

nrmseFile = os.path.join(
        statsDir,
        "HS-nrmse-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

hhFile = os.path.join(
        statsDir,
        "HS-hh-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

nbiFile = os.path.join(
        statsDir,
        "HS-nbi-altimeter_"
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
        "HS-bias-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

# check that data exists
print(nrmseFile)
assert os.path.exists(lonFile) == True
assert os.path.exists(latFile) == True
assert os.path.exists(nrmseFile) == True
assert os.path.exists(hhFile) == True
assert os.path.exists(nbiFile) == True
assert os.path.exists(biasFile) == True

lon = np.genfromtxt(lonFile)
lat = np.genfromtxt(latFile)
nrmse = np.genfromtxt(nrmseFile)
nbi = np.genfromtxt(nbiFile)
hh = np.genfromtxt(hhFile)
bias = np.genfromtxt(biasFile)

Latt, Lonn = np.meshgrid(lat, lon)

mask = bias == 99999
bias = np.ma.masked_array(bias, mask)
nrmse = np.ma.masked_array(nrmse, mask)
nbi = np.ma.masked_array(nbi, mask)
hh = np.ma.masked_array(hh, mask)


m = Basemap()

# ==========================================================================================
# NRMSE
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    nrmse, cmap="rainbow", origin="lower", extent=[-180, 180, -90, 90], vmin=0, vmax=0.5
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("NRMSE", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/HS-nrmse_"
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
# NBI
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    nbi, cmap="RdYlBu", origin="lower", extent=[-180, 180, -90, 90], vmin=-0.5, vmax=0.5
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("NBI", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/HS-nbi_"
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
# HH
# ==========================================================================================

fig, ax = plt.subplots(figsize=(9, 4))
grd = gridspec.GridSpec(1, 2, wspace=0.025, width_ratios=[1, 0.05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(
    hh, cmap="rainbow", origin="lower", extent=[-180, 180, -90, 90], vmin=0, vmax=1
)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color="gray")
# axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("HH", **title_font)

# Colorbar
axCb = plt.subplot(grd[0, 1])
cb = Colorbar(ax=axCb, mappable=plt1, orientation="vertical")  # , ticklocation = 'top')

plt.savefig(
    "data/stats/HS-hh_"
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
    bias, cmap="RdYlBu", origin="lower", extent=[-180, 180, -90, 90], vmin=-0.5, vmax=0.5
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
    "data/stats/HS-bias_"
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
