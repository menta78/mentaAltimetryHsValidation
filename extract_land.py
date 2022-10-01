import os, re
import netCDF4
from datetime import datetime
import numpy as np
import h5py

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM
from src.computeTidalStats import elaborateMeasures
from src.plotStatsTidals import elaborateMeasuresPlot
from src.plotTimeSeries import timeSeriesPlot
from alphaBetaLab.alphaBetaLab.abTriangularMesh import loadFromGr3File


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

rootDir = "/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation"

gridFile = os.path.join(rootDir, "data/hgrid.gr3")
ncFile = os.path.join(rootDir, "data/schismwwm/ERA5_schismwwm_20030303.nc")

varNames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]

# tidal position
lonInterest = -8.4006 #-3.9657
latInterest = 43.3674 # 48.7184
target = [lonInterest, latInterest]


# extract land nodes
meshObj = loadFromGr3File(gridFile)
land = meshObj.landBoundaryNodes 

# extract model elev data
ds = netCDF4.Dataset(ncFile)
lon = ds.variables[varNames[0]][:]
lat = ds.variables[varNames[1]][:]
hs = ds.variables[varNames[2]][0, :]


idxLon, lonNode = find_nearest(lon, lonInterest)
idxLat, latNode = find_nearest(lat, latInterest)

print(idxLon, lonNode)
print(idxLat, latNode)


if idxLon and idxLat in land:
    print("touch land")
    print(land)
