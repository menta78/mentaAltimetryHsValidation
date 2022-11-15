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

options_savefig = {'dpi': 150, "bbox_inches": "tight", "transparent": False}
title_font = {'size':'12', 'color':'black', 'weight':'normal',
                      'verticalalignment':'bottom'}

# time interval
startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
startDate, endDate = datetime(2012, 1, 1), datetime(2019, 12, 31)
startDate, endDate = datetime(2003, 1, 1), datetime(2009, 12, 30)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)




lonFile = os.path.join(statsDir, "lons.csv")
latFile = os.path.join(statsDir, "lats.csv")

nrmseFile = os.path.join(
        statsDir,
        "HS-nrmse-altimeter_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

hhFile = os.path.join(
        statsDir,
        "HS-hh-altimeter_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

nbiFile = os.path.join(
        statsDir,
        "HS-nbi-altimeter_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

biasFile = os.path.join(
        statsDir,
        "HS-bias-altimeter_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".csv",
    )

# check that data exists
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

m=Basemap()

#==========================================================================================
# NRMSE
#==========================================================================================

fig, ax = plt.subplots(figsize=(9,4))
grd = gridspec.GridSpec(1, 2, wspace=.025, width_ratios=[1, .05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(nrmse, 
        cmap='Reds',
        origin='lower',
        extent=[-180,180, -90,90],
        vmin=0, vmax=1)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='gray')
#axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("NRMSE", **title_font)

# Colorbar
axCb = plt.subplot(grd[0,1])
cb = Colorbar(ax = axCb, mappable = plt1, orientation = 'vertical') #, ticklocation = 'top')

plt.savefig('data/stats/HS-nrmse_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
plt.close(fig)

#==========================================================================================
# NBI
#==========================================================================================

fig, ax = plt.subplots(figsize=(9,4))
grd = gridspec.GridSpec(1, 2, wspace=.025, width_ratios=[1, .05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(nbi, 
        cmap='Reds',
        origin='lower',
        extent=[-180,180, -90,90],
        vmin=0, vmax=1)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='gray')
#axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("NBI", **title_font)

# Colorbar
axCb = plt.subplot(grd[0,1])
cb = Colorbar(ax = axCb, mappable = plt1, orientation = 'vertical') #, ticklocation = 'top')

plt.savefig('data/stats/HS-nbi_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
plt.close(fig)


#==========================================================================================
# HH
#==========================================================================================

fig, ax = plt.subplots(figsize=(9,4))
grd = gridspec.GridSpec(1, 2, wspace=.025, width_ratios=[1, .05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(hh, 
        cmap='Reds',
        origin='lower',
        extent=[-180,180, -90,90],
        vmin=0, vmax=1)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='gray')
#axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("HH", **title_font)

# Colorbar
axCb = plt.subplot(grd[0,1])
cb = Colorbar(ax = axCb, mappable = plt1, orientation = 'vertical') #, ticklocation = 'top')

plt.savefig('data/stats/HS-hh_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
plt.close(fig)

#==========================================================================================
# BIAS
#==========================================================================================

fig, ax = plt.subplots(figsize=(9,4))
grd = gridspec.GridSpec(1, 2, wspace=.025, width_ratios=[1, .05])

# Map and scatter plot
axMap = plt.subplot(grd[0, 0])

plt1 = axMap.imshow(bias, 
        cmap='RdBu',
        origin='lower',
        extent=[-180,180, -90,90],
        vmin=-1, vmax=1)
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='gray')
#axMap.set(xlim=[-180,180], ylim=[-90,90])
axMap.set_aspect("equal", "box")
axMap.set_title("BIAS", **title_font)

# Colorbar
axCb = plt.subplot(grd[0,1])
cb = Colorbar(ax = axCb, mappable = plt1, orientation = 'vertical') #, ticklocation = 'top')

plt.savefig('data/stats/HS-bias_'+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d") + ".png", **options_savefig)
plt.close(fig)