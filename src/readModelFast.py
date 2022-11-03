import multiprocessing
import os
import re
from datetime import datetime, timedelta
#from GeslaDataset.gesla import GeslaDataset

import netCDF4
import numpy as np
import src.utils as utils



def __elabFile(mdlF, mdlFPrev=None, mdlFNext=None, varNames=None):
    # xr_tidalGauge corresponds to one tidal gauge instrument
    varNames = (
        varNames if not varNames is None else ["longitude", "latitude", "hs", "time"]
    )
    print("  elaborating file " + mdlF)

    print("    ... loading model nc")
    fpth = os.path.join(_modelNcFileDir, mdlF)
    ds = netCDF4.Dataset(fpth)
    timeVarName = varNames[3]
    try:
        tmnc = ds.variables[timeVarName]
    except:
        print("something wrong in file " + fpth)
        outFlPath = "none"
        return outFlPath

    def getNcCalendar(tmnc):
        try:
            return tmnc.calendar
        except:
            return "standard"

    clndr = getNcCalendar(tmnc)
    try:
        tmstrt = netCDF4.num2date(
            tmnc[0], tmnc.units, clndr, only_use_cftime_datetimes=False
        )
        print("First time record schism file ", fpth, " is ", tmstrt)
    except:
        print("something wrong in file " + fpth)
        outFlPath = "none"
        return "case " + str(mdlF)


    outFlName = mdlF.replace(".nc", "_" + "tidalGauge")
    outFlPath = os.path.join(_destDir, outFlName)
    if (not _overwriteExisting) and os.path.isfile(outFlPath + ".npy"):
        print("      file " + outFlPath + " already exists. Skipping")
        ds.close()
        return outFlPath

    if not _startDate is None and not _endDate is None:
        tmend = netCDF4.num2date(
            tmnc[-1], tmnc.units, clndr, only_use_cftime_datetimes=False
        )
        if (_endDate < tmstrt) or (_startDate > tmend):
            print("  file out of date range. Skipping")
            ds.close()
            return "case " + str(mdlF)

    tmmdl = toJulian(netCDF4.num2date(
        tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
    ))
  
    _lock.acquire()
    try:
        lon = ds.variables[varNames[0]][:]
        lat = ds.variables[varNames[1]][:]
        hs = ds.variables[varNames[2]][:,:]

        usePrev = not mdlFPrev is None
        if usePrev:
            fpth_ = os.path.join(_modelNcFileDir, mdlFPrev)
            try:
                ds_ = netCDF4.Dataset(fpth_)
                tmnc_ = ds_.variables[timeVarName]
                tmmdl_ = toJulian(netCDF4.num2date(
                    tmnc_[-1],
                    tmnc_.units,
                    getNcCalendar(tmnc_),
                    only_use_cftime_datetimes=False,
                ))
                hs_ = ds_.variables[varNames[2]][-1, :]
                hs_ = hs_.reshape([1, len(hs_)])
                hs = np.concatenate([hs_, hs], axis=0)
                tmmdl = np.insert(tmmdl, 0, tmmdl_)
                ds_.close()
            except:
                pass

        useNext = not mdlFNext is None
        if useNext:
            fpth_ = os.path.join(_modelNcFileDir, mdlFNext)
            try:
                ds_ = netCDF4.Dataset(fpth_)
                tmnc_ = ds_.variables[timeVarName]
                tmmdl_ = toJulian(netCDF4.num2date(
                    tmnc_[0],
                    tmnc_.units,
                    getNcCalendar(tmnc_),
                    only_use_cftime_datetimes=False,
                ))
                hs_ = ds_.variables[varNames[2]][0, :]
                hs_ = hs_.reshape([1, len(hs_)])
                hs = np.concatenate([hs, hs_], axis=0)
                tmmdl = np.append(tmmdl, tmmdl_)
                ds_.close()
            except:
                pass

        ds.close()
    finally:
        _lock.release()

    return tmmdl, lon, lat, hs


def __elabFile_schismWWM(args):
    mdlFPrev = args[0]
    mdlF = args[1]

    mdlFNext = args[2]
    return __elabFile(mdlF, mdlFPrev=mdlFPrev, mdlFNext=mdlFNext, varNames=_varnames)


def readNcSchism(
    modelNcFileDir,
    destDir,
    startDate=None,
    endDate=None,
    flpattern="(.*)\.nc$",
    overwriteExisting=True,
    nParWorker=4,
):
    print("interpolating model data to tidal gauge")
    mdlfl0 = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    dts = [
        int(re.match("ERA5_schismwwm_([0-9]*)\.nc", fn).groups(0)[0])
        for fn in mdlfl0
    ]
    iii = np.argsort(dts)
    mdlfl = np.array(mdlfl0)[iii]
    mdlflpre = mdlfl[:-1].tolist()
    mdlflpre.insert(0, None)
    mdlflnext = mdlfl[1:].tolist()
    mdlflnext.append(None)

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _startDate, _endDate, dtmngr, _varnames
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _varnames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]


    # dtmngr = coarsenSatData.satDataManager(crsSatDataDir, boundaries)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    fls = zip(mdlflpre, mdlfl, mdlflnext)
    print(mdlflpre, mdlfl, mdlflnext)
    flItr = map_(__elabFile_schismWWM, fls)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)

    return fls