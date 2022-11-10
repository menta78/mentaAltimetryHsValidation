import os, re
from datetime import datetime

import numpy as np
import netCDF4
from matplotlib.tri import Triangulation, LinearTriInterpolator
from scipy.interpolate import interp1d, RegularGridInterpolator
import multiprocessing

import src.coarsenSatData as coarsenSatData
import src.coarsenCmemsSshSatData as coarsenCmemsSshSatData


GRID_TYPE_REGULAR = "regular"
GRID_TYPE_UNSTRUCT = "unstructured"


def __elabFile(mdlF, mdlFPrev=None, mdlFNext=None, varNames=None):
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

    outFlName = mdlF.replace(".nc", "_hsModelAndSatObs_" + tmstrt.strftime("%Y%m%d"))
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

    dtSatAll = dtmngr.getAllSatDataBtwDt(tmmdl[0], tmmdl[-1])
    dtSatAll = dtmngr.stackData(dtSatAll)
    if (dtSatAll is None) or (dtSatAll.shape[0] == 0):
        print("      no sat data found for " + mdlF)
        ds.close()
        return outFlPath

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

    # Substract the model's mean of each node
    if _meanModelFile:
        ds = netCDF4.Dataset(_meanModelFile)
        try:
            meanelev = ds.variables["elev"][0,:]
        except:
            print("something wrong in file " + meanModelFile)
            outFlPath = "none"
            return outFlPath

        for i in range(hs.shape[-1]):
            hs[:,i] = hs[:,i] - meanelev[i]


    lonSat = dtSatAll[:, 1]
    latSat = dtSatAll[:, 2]

    print("    ... interpolating to sat point for each time step")
    ntmmdl = len(tmmdl)
    nobs = dtSatAll.shape[0]
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
            intp0[itm, :] = intpltr(latSat, lonSat)
        elif _gridType == GRID_TYPE_REGULAR:
            intpltr = RegularGridInterpolator(
                grdPoints, hsii, bounds_error=False, fill_value=np.nan
            )
            intp0[itm, :] = intpltr((latSat, lonSat))

    print("    ... interpolating on time")
    tmstmpmdl = [datetime.timestamp(t) for t in tmmdl]
    tmstmpsat = dtSatAll[:, 0]
    intp = np.zeros([nobs, 1]) * np.nan
    for iobs in range(nobs):
        intpltr = interp1d(tmstmpmdl, intp0[:, iobs])
        intp[iobs] = intpltr(tmstmpsat[iobs])

    dtSatAndMod = np.concatenate([dtSatAll, intp], 1)

    print("    ... saving output file")
    np.save(outFlPath, dtSatAndMod)
    return outFlPath


def interpolateModelToCoarsenedSatData_WW3(
    crsSatDataDir,
    modelNcFileDir,
    destDir,
    boundaries=None,
    startDate=None,
    endDate=None,
    gridType=GRID_TYPE_UNSTRUCT,
    flpattern="(.*)\.nc$",
    overwriteExisting=True,
    nParWorker=12,
):
    print("interpolating model data to sat")
    mdlfl = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    mdlfl.sort()
    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _startDate, _endDate, dtmngr, _gridType
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _gridType = gridType
    dtmngr = coarsenSatData.satDataManager(crsSatDataDir, boundaries)
    print(dtmngr)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    flItr = map_(__elabFile, mdlfl)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)


interpolateModelToCoarsenedSatData = interpolateModelToCoarsenedSatData_WW3


def __elabFile_schismWWM(args):
    mdlFPrev = args[0]
    mdlF = args[1]
    mdlFNext = args[2]
    return __elabFile(mdlF, mdlFPrev=mdlFPrev, mdlFNext=mdlFNext, varNames=_varnames)


def interpolateModelTocoarsenCmemsSshSatData_schismWWM(
    crsSatDataDir,
    modelNcFileDir,
    destDir,
    meanModelFile,
    boundaries=None,
    startDate=None,
    endDate=None,
    flpattern="(.*)\.nc$",
    overwriteExisting=True,
    nParWorker=4,
):
    print("interpolating model data to sat")
    mdlfl0 = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    dts = [
        int(re.match("ERA5_schismwwm_([0-9]*)\.nc", fn).groups(0)[0])
        for fn in mdlfl0
    ]
    print(dts)
    iii = np.argsort(dts)
    mdlfl = np.array(mdlfl0)[iii]
    mdlflpre = mdlfl[:-1].tolist()
    mdlflpre.insert(0, None)
    mdlflnext = mdlfl[1:].tolist()
    mdlflnext.append(None)

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _meanModelFile, _startDate, _endDate, dtmngr, _gridType, _varnames
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _meanModelFile = meanModelFile
    _startDate, _endDate = startDate, endDate
    _gridType = GRID_TYPE_UNSTRUCT
    _varnames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]
    dtmngr = coarsenCmemsSshSatData.satDataManager(crsSatDataDir, boundaries)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    fls = zip(mdlflpre, mdlfl, mdlflnext)
    flItr = map_(__elabFile_schismWWM, fls)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)


def interpolateModelToCoarsenedSatData_schismWWM(
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
    print("interpolating model data to sat")
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

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _startDate, _endDate, dtmngr, _gridType, _varnames
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _gridType = GRID_TYPE_UNSTRUCT
    _varnames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "WWM_1", "time"]
    dtmngr = coarsenSatData.satDataManager(crsSatDataDir, boundaries)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    fls = zip(mdlflpre, mdlfl, mdlflnext)
    flItr = map_(__elabFile_schismWWM, fls)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)


def interpolateModelToCoarsenedSatData_WWMExperimental(
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
    print("interpolating model data to sat")
    mdlfl0 = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    dts = [int(re.match("WWM_out_([0-9]*)\.nc", fn).groups(0)[0]) for fn in mdlfl0]
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
    _varnames = ["lon", "lat", "HS", "ocean_time"]
    dtmngr = coarsenSatData.satDataManager(crsSatDataDir, boundaries)
    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map
    fls = zip(mdlflpre, mdlfl, mdlflnext)
    flItr = map_(__elabFile_schismWWM, fls)
    for mdlF in flItr:
        print("  file successfully saved: " + mdlF)


if __name__ == "__main__":
    crsSatDataDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION/redSeaAndPersicGulf/coarsenedSatData"

    # RED SEA 10 years
    modelNcFileDir = "/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_red_sea_unst/results/"
    destDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst/hsModelAndSatObs"
    # modelNcFileDir = '/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_red_sea_unst_noobst/results/'
    # destDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/hsModelAndSatObs'
    boundaries = [[30, 32], [30, 12.7], [43.8, 12.7], [43.8, 32]]
    startDate, endDate = datetime(2000, 1, 1), datetime(2010, 1, 1)

    # PERSIC GULF 10 years
    modelNcFileDir = "/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_persic_gulf_unst/results/"
    destDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst/hsModelAndSatObs"
    modelNcFileDir = "/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_persic_gulf_unst_noobst/results/"
    destDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst_noobst/hsModelAndSatObs"
    boundaries = [[47, 23.5], [47, 31], [56.5, 31], [56.5, 23.5]]
    startDate, endDate = datetime(2000, 1, 1), datetime(2010, 1, 1)

    # RED SEA 40 YEARS
    modelNcFileDir = "/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_red_sea_unst/results/"
    destDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION/hindcast_red_sea_unst/hsModelAndSatObs"
    # modelNcFileDir = '/media/lmentaschi/TOSHIBA EXT/ClimateRuns/WW3_UNST/hindcast_red_sea_unst_noobst/results/'
    # destDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/hsModelAndSatObs'
    boundaries = [[30, 32], [30, 12.7], [43.8, 12.7], [43.8, 32]]
    startDate, endDate = datetime(1977, 1, 1), datetime(2021, 1, 1)

    interpolateModelToCoarsenedSatData(
        crsSatDataDir, modelNcFileDir, destDir, boundaries, startDate, endDate
    )
