import os, re, glob, itertools
from datetime import datetime, timedelta

import netCDF4
import numpy as np
import h5py

def get_serie_gesla(fileName):
    f = h5py.File(fileName,'r')
    data = f.get("GESELD")

    npoints = data["residual"].shape[0] # number of stations

    res   = []
    time  = []
    lon   = []
    lat   = []

    for i in range(npoints):
        ref = f["GESELD"]["longitude"][i][0]
        lon.append(f[ref][0][0])
        ref = f["GESELD"]["latitude"][i][0]
        lat.append(f[ref][0][0])
        ref = f["GESELD"]["residual"][i][0]
        res.append(np.array(f[ref][0]))
        ref = f["GESELD"]["time"][i][0]
        time.append(np.array(f[ref][0]))
    return lon, lat, res, time


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def find_closest_node(lonM, latM, target):
    dsts  = []
    nodes = []

    i = 0
    for lon, lat in zip(lonM, latM):
        coordinate = np.array([lon, lat])
        dsts.append(np.linalg.norm(coordinate - target))
        nodes.append(i)
        i += 1

    idX = np.argmin(np.array(dsts))
    node = nodes[idX]
    return node


def get_subsest_list(array, a, b):
    if a >= b:
        print("first value must be smaller than second")
        return [], False, 0, 0

    if a < min(array) or b > max(array):
        return [], False, 0, 0

    i = np.argmin(np.abs(np.array(array)-a))
    j = np.argmin(np.abs(np.array(array)-b))

    return array[i:j], True, i, j

def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac

def toJulian(tarray_):
    i = 0
    tarray = tarray_
    # Convert to julian days
    for tt in tarray_:
        tarray[i] = datetime2matlabdn(tt)
        i += 1
    return tarray

def getFiles(hsSatAndModelDir, startDate, endDate, fltPattern, extension):
    fl = []
    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = os.path.join(hsSatAndModelDir, fltPattern + strtime + "*" + extension)
        fileFound = glob.glob(pthfile)
        fl.append(fileFound)
        startDate += timedelta(days=1)

    return list(itertools.chain(*fl))


def getNcCalendar(tmnc):
    try:
        return tmnc.calendar
    except:
        return "standard"


def getModelVariables(flsPath, varNames=None):
    varNames = (
        varNames if not varNames is None else ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]
    )

    nfiles = len(flsPath)

    ds = netCDF4.Dataset(flsPath[0])
    timeVarName = varNames[3]
    try:
        tmnc = ds.variables[timeVarName]
        clndr = getNcCalendar(tmnc)
    except:
        pass

    tmmdl = toJulian(netCDF4.num2date(
        tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
    ))

    lon = ds.variables[varNames[0]][:]
    lat = ds.variables[varNames[1]][:]
    var = ds.variables[varNames[2]][:,:]

    for i in range(1, nfiles):
        print(flsPath[i])
        ds = netCDF4.Dataset(flsPath[i])
        timeVarName = varNames[3]
        try:
            tmnc = ds.variables[timeVarName]
        except:
            print("something wrong in file " + fpth)
            outFlPath = "none"
            return outFlPath
    
        tmmdl_ = toJulian(netCDF4.num2date(
            tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
        ))

        var_ = ds.variables[varNames[2]][:,:]

        tmmdl = np.concatenate([tmmdl, tmmdl_], axis=0)
        var = np.concatenate([var, var_], axis=0)

    return tmmdl, lon, lat, var




def computeStats(obs, model, pth):

    # computing r2 (and other measures) gauge by gauge
    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

    condition1_ = obs >= pobs 
    condition2_ = model >= model
    condition_ = condition1_ & condition2_

    N = np.nansum(condition_)

    ssres_ = np.nansum((obs-model)**2, initial=0, where=condition_)
    sstot_ = np.nansum((obs-np.nanmean(obs))**2,  initial=0, where=condition_)
    nsc1   = np.nansum(np.abs(obs-model), initial = 0, where=condition_)
    nsc2   = np.nansum(np.abs(obs-np.nanmean(obs)), initial=0, where=condition_)

    absre_ = np.nansum(model - obs, initial=0, where=condition_)
    nobs = np.nansum(obs, initial=0, where=condition_)

    sigmaObs = np.sqrt(np.nansum((obs-np.nanmean(obs))**2,  initial=0, where=condition_))
    sigmaModel = np.sqrt(np.nansum((model-np.nanmean(model))**2,  initial=0, where=condition_))
    cov_ = np.nansum((obs-np.nanmean(obs))*(model-np.nanmean(model)), initial=0, where=condition_)

    r2 = 1 - ssres_/sstot_
    nse = 1 - nsc1/nsc2
    ab = absre_/N
    rb = absre_/nobs * 100
    rmse = np.sqrt(ssres_/N)
    pearson = cov_/(sigmaModel*sigmaObs)

    print(r2, nse)
    stats = {}
    stats["r2"] = r2
    stats["NSE"] = nse
    # stats["Normalized_NSE"] = nnse
    # stats["Normalized_r2"] = nr2
    stats["Absolute_Bias"] = ab
    stats["Relative_Bias"] = rb
    stats["RMSE"] = rmse
    stats["Pearson"] = pearson
    stats["N"] = N

    return stats

