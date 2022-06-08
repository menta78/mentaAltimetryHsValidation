import multiprocessing
import os
import re
from datetime import datetime
from GeslaDataset.gesla import GeslaDataset

import netCDF4
import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation
from scipy.interpolate import RegularGridInterpolator, interp1d

GRID_TYPE_REGULAR = "regular"
GRID_TYPE_UNSTRUCT = "unstructured"


def __elabFile(mdlF, filenameTidalGauge, mdlFPrev=None, mdlFNext=None, varNames=None):
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
    except:
        print("something wrong in file " + fpth)
        outFlPath = "none"
        return outFlPath

    #outFlName = mdlF.replace(".nc", "_sshModelAndTidalObs_" + tmstrt.strftime("%Y%m%d"))
    outFlName = mdlF.replace(".nc", "_sshModelAndTidalObs_" + filenameTidalGauge)
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
            return outFlPath

    tmmdl = netCDF4.num2date(
        tmnc[:], tmnc.units, clndr, only_use_cftime_datetimes=False
    )

    _lock.acquire()
    try:
        lon = ds.variables[varNames[0]][:]
        lat = ds.variables[varNames[1]][:]
        hs = ds.variables[varNames[2]][:]

        usePrev = not mdlFPrev is None
        if usePrev:
            fpth_ = os.path.join(_modelNcFileDir, mdlFPrev)
            try:
                ds_ = netCDF4.Dataset(fpth_)
                tmnc_ = ds_.variables[timeVarName]
                tmmdl_ = netCDF4.num2date(
                    tmnc_[-1],
                    tmnc_.units,
                    getNcCalendar(tmnc_),
                    only_use_cftime_datetimes=False,
                )
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
                tmmdl_ = netCDF4.num2date(
                    tmnc_[0],
                    tmnc_.units,
                    getNcCalendar(tmnc_),
                    only_use_cftime_datetimes=False,
                )
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

    if _gridType == GRID_TYPE_UNSTRUCT:
        triObj = Triangulation(lat, lon)
    elif _gridType == GRID_TYPE_REGULAR:
        grdPoints = (lat.filled(), lon.filled())
    else:
        raise Exception("unsupported grid type: " + _gridType)

    lonSat = xr_tidalGauge["longitude"] # a number
    latSat = xr_tidalGauge["latitude"]
    timeSat = xr_tidalGauge["date_time"]

    print("    ... interpolating to sat point for each time step")
    ntmmdl = len(tmmdl)
    nobs = timeSat.shape[0]
    intp0 = np.zeros([ntmmdl, nobs]) * np.nan
    for tm, itm in zip(tmmdl, range(ntmmdl)):
        print("      processing time " + str(tm))
        hsii = hs[itm, :]
        try:
            hsii = hsii.filled(np.nan)
        except:
            pass
        if _gridType == GRID_TYPE_UNSTRUCT:
            intpltr = LinearTriInterpolator(triObj, hsii)
            # TODO:
            # latSat y lonSat en tidal gauge data son un valor entonces
            # tengo que repetirlo tantas veces como observaciones haya
            # NOTA. Eso se hará arriba
            intp0[itm, :] = intpltr(latSat, lonSat)
        elif _gridType == GRID_TYPE_REGULAR:
            intpltr = RegularGridInterpolator(
                grdPoints, hsii, bounds_error=False, fill_value=np.nan
            )
            intp0[itm, :] = intpltr((latSat, lonSat))

    print("    ... interpolating on time")
    tmstmpmdl = [datetime.timestamp(t) for t in tmmdl]
    tmstmpsat = timeSat
    intp = np.zeros([nobs, 1]) * np.nan
    for iobs in range(nobs):
        intpltr = interp1d(tmstmpmdl, intp0[:, iobs])
        intp[iobs] = intpltr(tmstmpsat[iobs])

    dtSatAndMod = np.concatenate([xr_tidalGauge["sea_level"], intp], 1)

    print("    ... saving output file")
    np.save(outFlPath, dtSatAndMod)
    return outFlPath


def __elabFile_schismWWM(args):
    mdlFPrev = args[0]
    mdlF = args[1]
    mdlFNext = args[2]
    return __elabFile(mdlF, mdlFPrev=mdlFPrev, mdlFNext=mdlFNext, varNames=_varnames)


def interpolateModelToTidalGauge_schismWWM(
    crsSatDataDir,
    modelNcFileDir,
    destDir,
    boundaries=None,
    startDate=None,
    endDate=None,
    flpattern="(.*)\.nc$",
    overwriteExisting=True,
    nParWorker=4,
):
    print("interpolating model data to tidal gauge")
    mdlfl0 = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    dts = [
        int(re.match("schout_([0-9]*)\_compressed.nc", fn).groups(0)[0])
        for fn in mdlfl0
    ]
    iii = np.argsort(dts)
    mdlfl = np.array(mdlfl0)[iii]
    mdlflpre = mdlfl[:-1].tolist()
    mdlflpre.insert(0, None)
    mdlflnext = mdlfl[1:].tolist()
    mdlflnext.append(None)

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _startDate, _endDate, dtmngr, _gridType, _varnames
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _gridType = GRID_TYPE_UNSTRUCT
    _varnames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "hvel", "time"]
    # dtmngr = coarsenSatData.satDataManager(crsSatDataDir, boundaries)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    fls = zip(mdlflpre, mdlfl, mdlflnext)
    flItr = map_(__elabFile_schismWWM, fls)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)

    return fls
