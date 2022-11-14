import numpy as np
import xarray as xa
import src.geodiccaCoastLine


def mapByLonLatCumm(
    mp, times, lons, lats, msrs, mods, mapdx, mapdy, lonlims=[-180, 180], latlims=[-90, 90]
):
    # import pdb; pdb.set_trace()
    maplons = np.arange(lonlims[0], lonlims[1], mapdx)
    maplats = np.arange(latlims[0], latlims[1], mapdy)

    #mp = {}

    for tm, lon, lat, msr, mod in zip(times, lons, lats, msrs, mods):
        if (mod > 20) or (msr > 20):
            continue
        if (not (lonlims[0] <= lon) and (lon <= lonlims[1])) or (
            not (latlims[0] <= lat) and (lat <= latlims[1])
        ):
            continue
        ix = int((lon - lonlims[0]) // mapdx)
        iy = int((lat - latlims[0]) // mapdy)

        lst = mp.get((ix, iy), [[], [], [], [], []]) # return empty if value doesn't exists
        mp[(ix, iy)] = lst

        _tms = lst[0]
        _lons = lst[1]
        _lats = lst[2]
        _msrs = lst[3]
        _mods = lst[4]

        _tms.append(tm)
        _lons.append(lon)
        _lats.append(lat)
        _msrs.append(msr)
        _mods.append(mod)

    return maplons, maplats, mp

def mapByLonLat(
    times, lons, lats, msrs, mods, mapdx, mapdy, lonlims=[-180, 180], latlims=[-90, 90]
):
    # import pdb; pdb.set_trace()
    maplons = np.arange(lonlims[0], lonlims[1], mapdx)
    maplats = np.arange(latlims[0], latlims[1], mapdy)

    mp = {}

    for tm, lon, lat, msr, mod in zip(times, lons, lats, msrs, mods):
        if (mod > 20) or (msr > 20):
            continue
        if (not (lonlims[0] <= lon) and (lon <= lonlims[1])) or (
            not (latlims[0] <= lat) and (lat <= latlims[1])
        ):
            continue
        ix = int((lon - lonlims[0]) // mapdx)
        iy = int((lat - latlims[0]) // mapdy)

        lst = mp.get((ix, iy), [[], [], [], [], []]) # return empty if value doesn't exists
        mp[(ix, iy)] = lst

        _tms = lst[0]
        _lons = lst[1]
        _lats = lst[2]
        _msrs = lst[3]
        _mods = lst[4]

        _tms.append(tm)
        _lons.append(lon)
        _lats.append(lat)
        _msrs.append(msr)
        _mods.append(mod)

    return maplons, maplats, mp


def createCoastlinePointsMask(maplons, maplats, resl="h", minDistFromCoast=20000):
    print("  creating mask of points close to the coast")
    cstln = geodiccaCoastLine.coastLine(
        resl, llcrnr=[min(maplons), min(maplats)], urcrnr=[max(maplons), max(maplats)]
    )
    mask = np.ones([len(maplats), len(maplons)])
    npt = len(maplons) * len(maplats)
    dx = maplons[1] - maplons[0]
    dy = maplats[1] - maplats[0]
    ipt = 0
    for ix in range(len(maplons)):
        for iy in range(len(maplats)):
            lon = maplons[ix] + dx / 2
            lat = maplats[iy] - dy / 2
            xdst, ydst = cstln.getDistanceVectorFromCoast(lon, lat)
            if np.sqrt(xdst**2.0 + ydst**2.0) < minDistFromCoast:
                mask[iy, ix] = 0
            ipt += 1
            if ipt % 100 == 0:
                print("  done " + str(ipt / npt * 100) + "%")
    return mask


def createBathyMask(maplons, maplats, minSeaDepth, bathyFile):
    btmds = xa.open_dataset(bathyFile)
    intpBtmDs = btmds.interp(x=maplons, y=maplats)
    intpBtm = np.array(intpBtmDs.variables["z"])
    mask = intpBtm < minSeaDepth
    return mask


nminobs = 0


def computeBiasNrmse(lons, lats, mapdata):
    # import pdb; pdb.set_trace()
    bias = np.ones((len(lats), len(lons))) * 99999
    nrmse = np.ones((len(lats), len(lons))) * 99999
    dtcount = np.ones((len(lats), len(lons))) * 99999
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue
            _bias = (np.mean(mods) - np.mean(msrs)) / np.mean(msrs)
            _nrmse = (sum((mods - msrs) ** 2.0) / sum(msrs**2.0)) ** 0.5
            bias[iy, ix] = _bias
            nrmse[iy, ix] = _nrmse
            dtcount[iy, ix] = len(mods)
    mask = bias == 99999
    bias = np.ma.masked_array(bias, mask)
    nrmse = np.ma.masked_array(nrmse, mask)
    dtcount = np.ma.masked_array(dtcount, mask)
    return bias, nrmse, dtcount


def computeCumDeviations(lons, lats, mapdata):
    obsSum = np.ones((len(lats), len(lons))) * np.nan
    sqObsSum = np.ones((len(lats), len(lons))) * np.nan
    sqModSum = np.ones((len(lats), len(lons))) * np.nan
    devSum = np.ones((len(lats), len(lons))) * np.nan
    sqDevSum = np.ones((len(lats), len(lons))) * np.nan
    mdlByObsSum = np.ones((len(lats), len(lons))) * np.nan
    dtcount = np.ones((len(lats), len(lons))) * np.nan
    _consideredCells = 0
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue
            _obsSum = np.sum(msrs)
            _sqObsSum = np.sum(msrs**2.0)
            _sqModSum = np.sum(mods**2.0)
            _devSum = np.sum(mods - msrs)
            _sqDevSum = np.sum((mods - msrs) ** 2.0)
            _mdlByObsSum = np.sum(msrs * mods)
            obsSum[iy, ix] = _obsSum
            sqObsSum[iy, ix] = _sqObsSum
            sqModSum[iy, ix] = _sqModSum
            devSum[iy, ix] = _devSum
            sqDevSum[iy, ix] = _sqDevSum
            mdlByObsSum[iy, ix] = _mdlByObsSum
            dtcount[iy, ix] = len(mods)
#            print(iy, ix)
            _consideredCells += 1
    print("        considered cells: " + str(_consideredCells))
    return obsSum, sqObsSum, sqModSum, devSum, sqDevSum, mdlByObsSum, dtcount
    

def computeMean_cell(lons, lats, mapdata, mapDataAll):

    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            tm = np.array(data[0])
            lon = np.array(data[1])
            lat = np.array(data[2])
            msrs = np.array(data[3])
            mods = np.array(data[4])

            if len(mods) < nminobs:
                continue
            
            # tm_mean_ =  np.nanmean(tm)
            # lon_mean_ = np.nanmean(lon)
            # lat_mean_ = np.nanmean(lat)
            # sat_mean_ = np.nanmean(msrs)
            # mod_mean_ = np.nanmean(mods)
            tm_mean_ =  tm
            lon_mean_ = np.nanmean(lon)
            lat_mean_ = np.nanmean(lat)
            sat_mean_ = msrs
            mod_mean_ = mods

            lst = mapDataAll.get((ix, iy), [[], [], [], [], []]) # return empty if value doesn't exists
            lst[0].append(tm_mean_)
            lst[1].append(lon_mean_)
            lst[2].append(lat_mean_)
            lst[3].append(sat_mean_)
            lst[4].append(mod_mean_)

            mapDataAll[(ix, iy)] = lst

    return mapDataAll

def computeMean(lons, lats, mapdata):
    _mod_mean = 0
    _sat_mean = 0
    _consideredCells = 0
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue
            _mod_mean += np.sum(mods)
            _sat_mean += np.sum(msrs)
            _consideredCells += 1
    print("        considered cells: " + str(_consideredCells))
    sat_mean = _sat_mean/_consideredCells
    mod_mean = _mod_mean/_consideredCells
    return sat_mean, mod_mean

def computeSkillsSsh(lons, lats, mapdata):
    sqObsSum = np.ones((len(lats), len(lons))) * np.nan
    dtcount = np.ones((len(lats), len(lons))) * np.nan
    _consideredCells = 0
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue
            _sqDevSum = np.sum((mods - msrs) ** 2.0)
            sqDevSum[iy, ix] = _sqDevSum
            dtcount[iy, ix] = len(mods)
#            print(iy, ix)
            _consideredCells += 1
    print("        considered cells: " + str(_consideredCells))
    return sqDevSum, dtcount


def computeMaxima(lons, lats, mapdata):
    mxMsrs = np.ones((len(lats), len(lons))) * np.nan
    mxMdl = np.ones((len(lats), len(lons))) * np.nan
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue

            mxMsrs[iy, ix] = np.nanmax(msrs)
            mxMdl[iy, ix] = np.nanmax(mods)
    return mxMsrs, mxMdl


def computeMaxima(lons, lats, mapdata):
    mxMsrs = np.ones((len(lats), len(lons))) * np.nan
    mxMdl = np.ones((len(lats), len(lons))) * np.nan
    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            mods = np.array(data[4])
            if len(mods) < nminobs:
                continue

            mxMsrs[iy, ix] = np.nanmax(msrs)
            mxMdl[iy, ix] = np.nanmax(mods)
    return mxMsrs, mxMdl
