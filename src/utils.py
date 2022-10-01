import os, re, glob, itertools
from datetime import datetime, timedelta

import netCDF4
import numpy as np


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
        print("something wrong in file " + fpth)
        outFlPath = "none"
        return outFlPath

    tmmdl = toJulian(netCDF4.num2date(
        tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
    ))

    lon = ds.variables[varNames[0]][:]
    lat = ds.variables[varNames[1]][:]
    var = ds.variables[varNames[2]][:,:]

    for i in range(1, nfiles):
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


