import multiprocessing
import os
import re
from datetime import datetime, timedelta
#from GeslaDataset.gesla import GeslaDataset

import netCDF4
import numpy as np
import src.utils as utils
from matplotlib.tri import LinearTriInterpolator, Triangulation
from scipy.interpolate import RegularGridInterpolator, interp1d

GRID_TYPE_REGULAR = "regular"
GRID_TYPE_UNSTRUCT = "unstructured"


def read_file(line, fileName = "data/tidalGaugeData/stations-tidalGauge.txt"):
    # open the sample file used
    file = open(fileName)
    
    # read the content of the file opened
    content = file.readlines()

    return content[line].rstrip('\n')


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

def get_subsest_list(array, a, b):
    if a >= b:
        print("first value must be smaller than second")
        return [], False, 0, 0

    if a < min(array) or b > max(array):
        return [], False, 0, 0

    i = np.argmin(np.abs(np.array(array)-a))
    j = np.argmin(np.abs(np.array(array)-b))

    return array[i:j], True, i, j


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


    if _meanModelFile:
        # Substract the model's mean of each node
        ds = netCDF4.Dataset(_meanModelFile)
        try:
            meanelev = ds.variables["elev"][0,:]
        except:
            print("something wrong in file " + meanModelFile)
            outFlPath = "none"
            return outFlPath

        for i in range(hs.shape[-1]):
            hs[:,i] = hs[:,i] - meanelev[i]

    print("    ... Interpolant schism model for each time step")

    # Generate an interpolant for each time step of the schism model
    nobs = len(_lonTidal)
    ntmmdl = len(tmmdl)
  
    resT  = []
    intpT = []
    lonT = []
    latT = []
    indxT = []

    j = 0

    for i in range(nobs):
        print("    ... interpolating on time", flush=True)
        
        tmstmpTidal_ = _timeTidal[i][:]

        # Get timestamps recorded in the stations that belongs to the array of timestamps of the schism model
        tmstmpTidal, isSubset, idxMin, idxMax = get_subsest_list(tmstmpTidal_, min(tmmdl), max(tmmdl))
        
        if not isSubset:
            print("Missing data station " + read_file(i))
            continue
        
        ntTidal = len(tmstmpTidal)
        intp = np.zeros([ntTidal, 1]) * np.nan

        node = utils.find_closest_node(lon, lat, [_lonTidal[i], _latTidal[i]])
        print("------ node = ", node, flush=True)
        intp0 = hs[:, node]

        intpltr = interp1d(tmmdl, intp0, bounds_error=False)

        intp = intpltr(tmstmpTidal)
        print(intp0, intp)
        jfrijfri
        res = np.array(_resTidal[i][idxMin:idxMax])

        lonn = np.array(_lonTidal[i])
        latt = np.array(_latTidal[i])

        aux_ = res*0 + 1

        Lonn = lonn*aux_
        Latt = latt*aux_
        Indx = j*aux_

        resT.append(res.tolist())
        intpT.append(intp.tolist())
        lonT.append(Lonn.tolist())
        latT.append(Latt.tolist())
        indxT.append(Indx.tolist())

        j += 1

    resT = [item for sublist in resT for item in sublist]
    intpT = [item for sublist in intpT for item in sublist]
    lonT = [item for sublist in lonT for item in sublist]
    latT = [item for sublist in latT for item in sublist]
    indxT = [item for sublist in indxT for item in sublist]

    resTT = np.zeros([len(resT), 1])
    intpTT = np.zeros([len(intpT), 1])
    lonTT = np.zeros([len(lonT), 1])
    latTT = np.zeros([len(latT), 1])
    indxTT = np.zeros([len(indxT), 1])


    resTT[:,0] = np.array(resT)
    intpTT[:,0] = np.array(intpT)
    lonTT[:,0] = np.array(lonT)
    latTT[:,0] = np.array(latT)
    indxTT[:,0] = np.array(indxT)


    dtSatAndMod = np.concatenate([resTT, intpTT, lonTT, latTT, indxTT], axis=1)
    print("    ... saving output file in ", _destDir)
    #outFlName = mdlF.replace(".nc", "_" + read_file(i) + "_tidalGauge_" + tmstrt.strftime("%Y%m%d") )
    outFlName = mdlF.replace(".nc", "_" + "tidalGauge")
    outFlPath = os.path.join(_destDir, outFlName)
    print(outFlName)
    np.save(outFlPath, dtSatAndMod)
    
    return outFlPath


def __elabFile_schismWWM(args):
    mdlFPrev = args[0]
    mdlF = args[1]

    mdlFNext = args[2]
    return __elabFile(mdlF, mdlFPrev=mdlFPrev, mdlFNext=mdlFNext, varNames=_varnames)


def interpolateModelToTidalGauge_schismWWM(
    varsTidal,
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

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _meanModelFile, _startDate, _endDate, dtmngr, _gridType, _varnames, _lonTidal, _latTidal, _resTidal, _timeTidal
    _overwriteExisting = overwriteExisting
    _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _meanModelFile = meanModelFile
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _gridType = GRID_TYPE_UNSTRUCT
    _varnames = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "elev", "time"]
    _lonTidal = varsTidal[0]
    _latTidal = varsTidal[1]
    _resTidal = varsTidal[2]
    _timeTidal= varsTidal[3]

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
