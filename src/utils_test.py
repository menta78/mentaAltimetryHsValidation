import glob
import itertools
import os
import re
from datetime import datetime, timedelta

from scipy import signal, stats

import h5py
import netCDF4
import numpy as np


def get_serie_gesla(fileName):
    f = h5py.File(fileName, "r")
    data = f.get("GESELD")

    npoints = data["residual"].shape[0]  # number of stations

    res = []
    time = []
    lon = []
    lat = []

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
    dsts = []
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

    i = np.argmin(np.abs(np.array(array) - a))
    j = np.argmin(np.abs(np.array(array) - b))

    return array[i:j], True, i, j


def datetime2matlabdn(dt):
    ord = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt - datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / (
        24.0 * 60.0 * 60.0
    )
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
        varNames
        if not varNames is None
        else ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]
    )

    nfiles = len(flsPath)

    ds = netCDF4.Dataset(flsPath[0])
    timeVarName = varNames[3]
    try:
        tmnc = ds.variables[timeVarName]
        clndr = getNcCalendar(tmnc)
    except:
        pass

    tmmdl = toJulian(
        netCDF4.num2date(tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False)
    )

    lon = ds.variables[varNames[0]][:]
    lat = ds.variables[varNames[1]][:]
    var = ds.variables[varNames[2]][:, :]

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

        tmmdl_ = toJulian(
            netCDF4.num2date(
                tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
            )
        )

        var_ = ds.variables[varNames[2]][:, :]

        tmmdl = np.concatenate([tmmdl, tmmdl_], axis=0)
        var = np.concatenate([var, var_], axis=0)

    return tmmdl, lon, lat, var


def computeStatsHs(obs_, model_, pth):

    # computing r2 (and other measures) gauge by gauge
    pmodel = np.nanpercentile(model_, pth)
    pobs = np.nanpercentile(obs_, pth)

    # print(pmodel, pobs)
    # considering the data above the 95th percentile of the observation
    cnd = np.logical_and(obs_ >= pobs, model_ >= model_)
    model = model_[cnd]
    obs = obs_[cnd]

    N = len(obs)

    # for i in range(N):
    #     print(obs[i], model[i])

    devSum = np.nansum(model - obs)
    sqDevSum = np.nansum((model - obs) ** 2)
    obsSum = np.nansum(obs)
    # print(obsSum)

    cov = np.nansum(obs * model)

    nbi = devSum / obsSum
    absBias = devSum / N
    nrmse = np.sqrt(sqDevSum / np.nansum(obs**2))
    hh = np.sqrt(sqDevSum / cov)

    return nbi, absBias, nrmse, hh


def computeStats(obs_, model_, pth):

    not_nan_ind = ~np.isnan(obs_)
    obs__ = signal.detrend(obs_[not_nan_ind]) 

    model__ = model_[not_nan_ind]

    meanobs = np.nanmean(obs__)
    meanmodel = np.nanmean(model__)

    model_ = model__ - meanmodel
    obs_ = obs__ - meanobs

    # computing r2 (and other measures) gauge by gauge
    pmodel = np.nanpercentile(model_, pth)
    pobs = np.nanpercentile(obs_, pth)

    # print(pmodel, pobs)
    # considering the data above the 95th percentile of the observation
    cnd = np.logical_and(obs_ >= pobs, model_ >= model_)
    model = model_[cnd]
    obs = obs_[cnd]

    N = len(obs)

    # for i in range(N):
    #     print("===")
    #     print(model[i], obs[i])
    #     print(np.nanmean(model_), np.nanmean(obs_))
    #     print("===")

    # for i in range(N):
    #     print(obs[i], model[i])

    ssres_ = np.nansum((obs - model) ** 2)
    sstot_ = np.nansum((obs) ** 2)
    #sstot_ = np.nansum((obs - np.nanmean(obs_)) ** 2)
    nsc1 = np.nansum(np.abs(obs - model))
    #nsc2 = np.nansum(np.abs(obs - np.nanmean(obs_)))
    nsc2 = np.nansum(np.abs(obs))

    absre_ = np.nansum(model - obs)
    nobs = np.nansum(obs)

    #sigmaObs = np.sqrt(np.nansum((obs - np.nanmean(obs_)) ** 2))
    sigmaObs = np.sqrt(np.nansum((obs) ** 2))

    #sigmaModel = np.sqrt(np.nansum((model - np.nanmean(model_)) ** 2))

    sigmaModel = np.sqrt(np.nansum((model) ** 2))

    #cov_ = np.nansum((obs - np.nanmean(obs_)) * (model - np.nanmean(model_)))

    cov_ = np.nansum((obs) * (model))

    r2 = 1 - ssres_ / sstot_
    nse = 1 - nsc1 / nsc2
    ab = absre_ / N
    rb = absre_ / nobs * 100
    nrmse = np.sqrt(ssres_ / N)
    rmse = np.sqrt(ssres_)
    pearson = cov_ / (sigmaModel * sigmaObs)


    return r2, nse, ab, rb, rmse, nrmse, pearson


def load_paths(rootDir):
    # Directory where the coarsened satellite data are located
    # the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
    tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")
    crsSatDataDir = os.path.join(rootDir, "data/crsSatData/")
    crsWaveSatDataDir = os.path.join(rootDir, "data/crsSatDataWaves/")
    assert os.path.exists(tidalGaugeDataDir) == True
    assert os.path.exists(crsSatDataDir) == True
    assert os.path.exists(crsWaveSatDataDir) == True

    # Directory where the model nc files are located
    modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")
    assert os.path.exists(modelNcFilesDir) == True

    # Directory where the pairs observation/model are to be generated
    hsModelAndSatObsSshDir = os.path.join(rootDir, "data/satModelPairs/")
    hsModelAndSatObsHsDir = os.path.join(rootDir, "data/satWaveModelPairs/")
    hsModelAndSatObsTidalDir = os.path.join(rootDir, "data/tidalInterpolations/")
    assert os.path.exists(hsModelAndSatObsSshDir) == True
    assert os.path.exists(hsModelAndSatObsHsDir) == True
    assert os.path.exists(hsModelAndSatObsTidalDir) == True

    # Directory where the stats are generated
    statsDir = os.path.join(rootDir, "data/stats")
    assert os.path.exists(statsDir) == True

    # # Directory to Tidal Gauge
    # pathname = os.path.join(tidalGaugeDataDir, "GESLAv1_withResiduals.mat")

    return (
        tidalGaugeDataDir,
        crsSatDataDir,
        crsWaveSatDataDir,
        modelNcFilesDir,
        hsModelAndSatObsSshDir,
        hsModelAndSatObsHsDir,
        hsModelAndSatObsTidalDir,
        statsDir,
    )
