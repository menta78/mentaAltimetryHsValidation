import multiprocessing
import os
import re
from datetime import datetime, timedelta

# from GeslaDataset.gesla import GeslaDataset

import netCDF4
import numpy as np
import src.utils as utils
from matplotlib.tri import LinearTriInterpolator, Triangulation
from scipy.interpolate import RegularGridInterpolator, interp1d

GRID_TYPE_REGULAR = "regular"
GRID_TYPE_UNSTRUCT = "unstructured"


def read_file(line, fileName="data/tidalGaugeData/stations-tidalGauge.txt"):
    # open the sample file used
    file = open(fileName)

    # read the content of the file opened
    content = file.readlines()

    return content[line].rstrip("\n")


def find_closest_node(target):
    dsts = []
    nodes = []

    i = 0
    for lon, lat in zip(_lonGrid, _latGrid):
        coordinate = np.array([lon, lat])
        dsts.append(np.linalg.norm(coordinate - target))
        nodes.append(i)
        i += 1

    idX = np.argmin(np.array(dsts))
    node = nodes[idX]
    print(node)

    return node


def load_model_variables(mdlfl):
    nfiles = len(mdlfl)

    ds = netCDF4.Dataset(os.path.join(_modelNcFileDir, mdlfl[0]))
    timeVarName = _varsModel[-1]

    tmnc = ds.variables[timeVarName]
    tmmdl = utils.toJulian(
        netCDF4.num2date(
            tmnc[:], tmnc.units, _timeCalendar, only_use_cftime_datetimes=False
        )
    )

    var = ds.variables[_varsModel[2]][:, :]

    for i in range(1, nfiles):
        flsPath = os.path.join(_modelNcFileDir, mdlfl[i])
        ds = netCDF4.Dataset(flsPath)
        timeVarName = _varsModel[3]
        try:
            tmnc = ds.variables[timeVarName]
        except:
            print("something wrong in file " + fpth)
            outFlPath = "none"
            return outFlPath

        tmmdl_ = utils.toJulian(
            netCDF4.num2date(
                tmnc[:], tmnc.units, _timeCalendar, only_use_cftime_datetimes=False
            )
        )

        var_ = ds.variables[_varsModel[2]][:, :]

        tmmdl = np.concatenate([tmmdl, tmmdl_], axis=0)
        var = np.concatenate([var, var_], axis=0)

    return [tmmdl, var]


def interpolate_model_scatter_observations(
    nodeGrid, scatterObservation, timeModel, gridVar
):
    timeScatterData, isSubset, idxMin, idxMax = utils.get_subsest_list(
        _timeScatterData[scatterObservation], min(timeModel), max(timeModel)
    )

    if not isSubset:
        print(f"Missing data station {scatterObservation}")
        return None

    ntObs = len(timeScatterData)
    intp = np.zeros([ntObs, 1]) * np.nan

    print("------ node = ", nodeGrid, flush=True)
    intp0 = gridVar[:, nodeGrid]

    intpltr = interp1d(timeModel, intp0, bounds_error=False)

    intp = intpltr(timeScatterData)

    scatterValue = np.array(_obsScatterData[scatterObservation][idxMin:idxMax])

    # for i in range(len(intp)):
    #     print(intp[i], scatterValue[i])

    scatterValueTT = np.zeros([ntObs, 1])
    intpTT = np.zeros([ntObs, 1])
    lonTT = np.zeros([ntObs, 1])
    latTT = np.zeros([ntObs, 1])
    indxTT = np.zeros([ntObs, 1])
    nodeTT = np.zeros([ntObs, 1])

    aux_ = intp * 0 + 1

    lon = _lonScatterData[scatterObservation] * aux_
    lat = _latScatterData[scatterObservation] * aux_
    indx = scatterObservation * aux_
    node = nodeGrid * aux_

    scatterValueTT[:, 0] = np.array(scatterValue)
    intpTT[:, 0] = np.array(intp)
    lonTT[:, 0] = np.array(lon)
    latTT[:, 0] = np.array(lat)
    indxTT[:, 0] = np.array(indx)
    nodeTT[:, 0] = np.array(node)

    # print(intp.shape[:])
    # print(scatterValue.shape[:])
    # print(intp.shape[:])
    # print(lon.shape[:])
    # print(indx.shape[:])
    # print(node.shape[:])

    dtSatAndMod = np.concatenate(
        [scatterValueTT, intpTT, lonTT, latTT, indxTT, nodeTT], axis=1
    )
    print("    ... saving output file in ", _destDir)
    # outFlName = mdlF.replace(".nc", "_" + read_file(i) + "_tidalGauge_" + tmstrt.strftime("%Y%m%d") )
    outFlName = f"{_scatterDataType}_station_{scatterObservation}"
    outFlPath = os.path.join(_destDir, outFlName)
    print(outFlName)
    np.save(outFlPath, dtSatAndMod)


def _interpModel(args):
    nodeGrid = args[0]
    scatterObservation = args[1]

    return interpolate_model_scatter_observations(
        nodeGrid, scatterObservation, _timeModel, _gridVar
    )


def interpolateModelToScatterData(
    varsScatterObs,
    varsModel,
    modelNcFileDir,
    destDir,
    scatterDataType,
    startDate=None,
    endDate=None,
    flpattern="(.*)\.nc$",
    overwriteExisting=True,
    nParWorker=4,
):
    print("interpolating prediction in grid data to observation as scatter data")
    mdlfl0 = [f for f in os.listdir(modelNcFileDir) if re.match(flpattern, f)]
    dts = [
        int(re.match("ERA5_schismwwm_([0-9]*)\.nc", fn).groups(0)[0]) for fn in mdlfl0
    ]
    iii = np.argsort(dts)
    mdlfl = np.array(mdlfl0)[iii]
    # mdlflpre = mdlfl[:-1].tolist()
    # mdlflpre.insert(0, None)
    # mdlflnext = mdlfl[1:].tolist()
    # mdlflnext.append(None)

    global _overwriteExisting, _lock, _modelNcFileDir, _destDir, _startDate, _endDate, _lonScatterData, _latScatterData, _obsScatterData, _timeScatterData, _lonGrid, _latGrid, _varsModel, _timeCalendar, _timeModel, _gridVar
    _overwriteExisting = overwriteExisting
    # _lock = multiprocessing.Lock()
    _modelNcFileDir = modelNcFileDir
    _destDir = destDir
    _startDate, _endDate = startDate, endDate
    _lonScatterData = varsScatterObs[0]
    _latScatterData = varsScatterObs[1]
    _obsScatterData = varsScatterObs[2]
    _timeScatterData = varsScatterObs[3]
    _varsModel = varsModel

    ds = netCDF4.Dataset(os.path.join(modelNcFileDir, mdlfl[0]))
    _lonGrid = ds.variables[varsModel[0]][:]
    _latGrid = ds.variables[varsModel[1]][:]
    tmnc = ds.variables[varsModel[-1]]

    try:
        _timeCalendar = tmnc.calendar
    except:
        _timeCalendar = "standard"

    ds.close()

    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map

    # localize the grid nodes associated to each location of scatter observations
    positionsScatterData = zip(_lonScatterData, _latScatterData)
    nodesGrid = map_(find_closest_node, positionsScatterData)
    if nParWorker > 1:
        p.terminate()  # it's necessary to terminate before declaring new global variables and open new pools
    # store nodes in a file
    nodesFile = os.path.join(destDir, f"{scatterDataType}_nodes_{startDate.strftime('%Y%m%d')}_{endDate.strftime('%Y%m%d')}.npy")
    np.save(nodesFile, nodesGrid)


    # load grid variables from model to interpolate
    modelVariables = load_model_variables(mdlfl)
    global _timeModel, _gridVar, _scatterDataType
    _timeModel = modelVariables[0]
    _gridVar = modelVariables[1]
    _scatterDataType = scatterDataType

    if nParWorker > 1:
        p = multiprocessing.Pool(nParWorker)
        map_ = p.map
    else:
        map_ = map

    # interpolate selected variable for each scatter observation and time record
    scatterObservations = list(range(len(nodesGrid)))
    observations = zip(nodesGrid, scatterObservations)
    map_(_interpModel, observations)

    if nParWorker > 1:
        p.terminate()

    return fls
