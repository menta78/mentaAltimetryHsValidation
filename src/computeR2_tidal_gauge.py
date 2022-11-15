from http.client import NotConnected
import multiprocessing
import os
import re
from datetime import datetime
from xmlrpc.client import INVALID_ENCODING_CHAR
from GeslaDataset.gesla import GeslaDataset

import netCDF4
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation
from scipy.interpolate import RegularGridInterpolator, interp1d


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def check_times(minTime, maxTime, datetimeArray):

    for time in datetimeArray:
        if minTime <= time <= maxTime:
            return True
        else:
            pass

    return False


def get_tidal(filename):
    """return tidal gauge data

    Info:
        output["sea_level"]
        output["longitude"]
        output["latitude"]
        output["date_time"]

    Args:
        g3 (_type_): _description_
        filename (_type_): _description_

    Returns:
        xarray: tidal gauge data
    """
    return _g3.files_to_xarray([filename])


def convert_datetime_nctime(nctime):
    """_summary_

    Args:
        nctime (xarray): _description_

    Returns:
        float: _description_
    """

    def getNcCalendar(nctime):
        try:
            return nctime.calendar
        except:
            return "standard"

    clndr = getNcCalendar(nctime)

    date = netCDF4.num2date(nctime, units=nctime.units, calendar=clndr)
    return date


def convert_timestamp_nc_calendar(nctime, timestamp):
    """_summary_

    Args:
        nctime (xarray): _description_
        timestamp (str): _description_

    Returns:
        float: _description_
    """

    def getNcCalendar(nctime):
        try:
            return nctime.calendar
        except:
            return "standard"

    clndr = getNcCalendar(nctime)

    if isinstance(timestamp, datetime):
        ttimestamp = timestamp
    else:
        ttimestamp = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")

    numdate = netCDF4.date2num(ttimestamp, units=nctime.units, calendar=clndr)
    return numdate


def get_datetime_model(modelDir, filenames):
    """_summary_

    Args:
        modelDir (str): _description_
        filename (array): _description_

    Returns:
        array: _description_
    """

    i = 0

    for filename in filenames:
        if i == 0:
            print("    ... loading model nc")
            fpth = os.path.join(modelDir, filename)
            ds = netCDF4.Dataset(fpth)
            # timeVarName = varNames[3]
            timeVarName = "time"
            try:
                auxtmnc = ds.variables[timeVarName]
            except:
                print("something wrong in file " + fpth)
                outFlPath = "none"
                return outFlPath
            i += 1
        else:
            print("    ... loading model nc")
            fpth = os.path.join(modelDir, filename)
            ds = netCDF4.Dataset(fpth)
            # timeVarName = varNames[3]
            timeVarName = "time"
            try:
                tmnc = ds.variables[timeVarName]
            except:
                print("something wrong in file " + fpth)
                outFlPath = "none"
                return outFlPath
            auxtmnc = np.concatenate(auxtmnc, tmnc)

    return auxtmnc


def get_all_ncfiles(modelDir, ncExpression):
    flpattern = "(.*)\.nc$"
    mdlfl0 = [f for f in os.listdir(modelDir) if re.match(flpattern, f)]
    dts = [int(re.match(ncExpression, fn).groups(0)[0]) for fn in mdlfl0]
    iii = np.argsort(dts)
    mdlfl = np.array(mdlfl0)[iii]
    return mdlfl


def filter_ncfiles(modelDir, modelList, startDate, endDate):
    """_summary_

    Args:
        modelList (_type_): _description_
        startDate (_type_): _description_
        endDate (_type_): _description_
    """

    idx2delete = []
    for i in range(len(modelList)):
        filename = modelList[i]
        fpth = os.path.join(modelDir, filename)
        ds = netCDF4.Dataset(fpth)
        # timeVarName = varNames[3]
        timeVarName = "time"
        try:
            tmnc = ds.variables[timeVarName]
        except:
            print("something wrong in file " + fpth)
            outFlPath = "none"
            return outFlPath
        numStart = convert_timestamp_nc_calendar(tmnc, startDate)
        numEnd = convert_timestamp_nc_calendar(tmnc, endDate)
        if not check_times(numStart, numEnd, tmnc):
            idx2delete.append(i)

    return np.delete(modelList, idx2delete)


def get_ssh_time_model(modelDir, modelListFiltered):
    i = 0
    for filename in modelListFiltered:
        if i == 0:
            fpth = os.path.join(modelDir, filename)
            ds = netCDF4.Dataset(fpth)
            # timeVarName = varNames[3]
            sshVarName = "hvel"
            timeVarName = "time"
            try:
                ssh = ds.variables[sshVarName][:, :, 0, 0]
                # ssh = ds.variables[sshVarName][:,:]
                time = ds.variables[timeVarName]
                i += 1
            except:
                print("something wrong in file " + fpth)
                outFlPath = "none"
                return outFlPath
        else:
            fpth = os.path.join(modelDir, filename)
            ds = netCDF4.Dataset(fpth)
            try:
                ssh_ = ds.variables[sshVarName][:, :, 0, 0]
                # ssh_ = ds.variables[sshVarName][:,:]
                time_ = ds.variables[timeVarName]
            except:
                print("something wrong in file " + fpth)
                outFlPath = "none"
                return outFlPath

            ssh = np.concatenate([ssh, ssh_], axis=0)
            time = np.concatenate([time, time_], axis=0)
    return ssh, time


def get_latlon_model(modelDir, filename):
    fpth = os.path.join(modelDir, filename)
    ds = netCDF4.Dataset(fpth)
    # timeVarName = varNames[3]
    lonVarName = "SCHISM_hgrid_node_x"
    latVarName = "SCHISM_hgrid_node_y"
    try:
        lon = ds.variables[lonVarName][:]
        lat = ds.variables[latVarName][:]
    except:
        print("something wrong in file " + fpth)
        outFlPath = "none"
        return outFlPath
    return lat, lon


def interpolate_ssh_tidal_points(latModel, lonModel, sshModel, nctime, tidalGauge):
    # import matplotlib.pyplot as plt
    lonTidal = tidalGauge["longitude"].to_numpy()
    latTidal = tidalGauge["latitude"].to_numpy()

    # Prepare times
    timeTidal = tidalGauge["date_time"].to_numpy()
    timeTidalTimestamp_ = [t.astype(datetime) for t in timeTidal]
    timeTidalTimestamp = [t / 1e9 for t in timeTidalTimestamp_]
    timeModel_ = convert_datetime_nctime(nctime)
    timeModel = [
        datetime.strptime(
            t.strftime("%Y-%m-%d %H:%M:%S"), "%Y-%m-%d %H:%M:%S"
        ).timestamp()
        for t in timeModel_
    ]

    minIdx = np.where(np.in1d(timeTidalTimestamp, [timeModel[0]]))[0]
    maxIdx = np.where(np.in1d(timeTidalTimestamp, [timeModel[-1]]))[0]
    timeTidalTimestamp = timeTidalTimestamp[minIdx[0] : maxIdx[0]]

    # Interpolating in space
    nobs = len(timeTidalTimestamp)
    ntime = len(timeModel)
    intp0 = np.zeros([ntime, nobs]) * np.nan  # number of times, number of tidal gauges
    triObj = Triangulation(latModel, lonModel)
    # plt.triplot(triObj)
    # plt.savefig("mesh.png")

    intp = np.zeros([ntime, 1]) * np.nan

    for itm in range(ntime):
        print("interpolating ", itm)
        sshii = sshModel[itm, :]
        try:
            sshii = sshii.filled(np.nan)
        except:
            pass
        intpltr = LinearTriInterpolator(triObj, sshii)
        intp0[itm, :] = intpltr(latTidal, lonTidal)

        # intpltr = interp1d(timeTidalTimestamp, intp0[itm, :])
        # intp[itm] = intpltr(timeModel[itm])

    #    print(intp[:])
    #    return

    intp = np.zeros([nobs, 1]) * np.nan
    for iobs in range(nobs):
        intpltr = interp1d(timeModel, intp0[:, iobs])
        intp[iobs] = intpltr(timeTidalTimestamp[iobs])

    sl = tidalGauge["sea_level"].to_numpy()
    sl = sl[minIdx[0] : maxIdx[0]]
    for iobs in range(nobs):
        print(datetime.fromtimestamp(timeTidalTimestamp[iobs]), intp[iobs], sl[iobs])


def get_tidal_variables(tidalGauge, timeModel):
    lonTidal = tidalGauge["longitude"].to_numpy()  # a number
    latTidal = tidalGauge["latitude"].to_numpy()
    timeTidal = tidalGauge["date_time"].to_numpy()
    timeTidalConverted_ = [datetime.utcfromtimestamp(int(t) / 1e9) for t in timeTidal]
    timeTidalConverted = [
        convert_timestamp_nc_calendar(timeModel, t) for t in timeTidalConverted_
    ]
    return lonTidal, latTidal, timeTidalConverted


""" def compute_r2(g3, modelDir, ncExpression=r"schout_([0-9]*)\_compressed.nc", varNames=None, startDate, endDate):
    varNames = (
        varNames if not varNames is None else ["longitude", "latitude", "hs", "time"]
    )
    global _g3, _modelDir, _ncExpression, _varNames, _startDate, _endDate
    _g3 = g3
    _modelDir = modelDir
    _ncExpression = ncExpression
    _varNames = varNames
    _startDate = startDate
    _endDate = endDate

    # Retrieve lat lon array from model
    #latModel, lonModel = get_latlon_model() """
