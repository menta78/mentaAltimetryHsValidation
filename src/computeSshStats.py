import os, re, sys
import glob
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt

import src.coarsenSatData


import src.mapByLonLat as mll


maskPointsCloseToTheCoast = True
minDistFromCoast = 20000
minObsForValidation = 0
maskShallowPoints = False
minSeaDepth = -5
bathyFile = "/home/lmentaschi/usr/WaveWatchIII/gridgen1.1/reference_data/etopo2.nc"


def elaborateMeasures(
    startDate,
    endDate,
    hsSatAndModelDir,
    outputDir,
    latlims=[-63, 63],
    lonlims=[-180, 180],
    filterLowHs=False,
    filterHsThreshold=0.0,
    dx=0.2,
    dy=0.2,
    pth = 90, 
):

    years = []


    def iterateByYear():
        fls_ = [f for f in os.listdir(hsSatAndModelDir) if re.match("(.*)\.npy", f)]
        pttrn = "(.*)_hsModelAndSatObs_([0-9]{8}).npy"
        dtstrs = [re.match(pttrn, f).groups(0)[1] for f in fls_]
        dts = [datetime.strptime(dtstr, "%Y%m%d") for dtstr in dtstrs]
        iii = np.argsort(dts)
        dts = np.array(dts)[iii]
        fls_ = np.array(fls_)[iii]
        fls = [
            fls_[i]
            for i in range(len(fls_))
            if (startDate <= dts[i] and endDate >= dts[i])
        ]

        loopYear = None

        def loadFile(flpth):
            print("")
            print("    loading file " + flpth)
            satdts = np.load(flpth)
            if filterLowHs:
                hssat = satdts[:, 3]
                hsmdl = satdts[:, 4]
                cnd = np.logical_and(
                    hssat > filterHsThreshold, hsmdl > filterHsThreshold
                )
                satdts = satdts[cnd, :]
            return satdts

        def _elab_ssh(dts, lons, lats, saths, modhs):
            maplons, maplats, mapdata = mll.mapByLonLat(
                dts, lons, lats, saths, modhs, dx, dy, lonlims=lonlims, latlims=latlims
            )
            return (
                maplons,
                maplats,
                mapdata,
            )

        def _elab(dts, lons, lats, saths, modhs):
            maplons, maplats, mapdata = mll.mapByLonLat(
                dts, lons, lats, saths, modhs, dx, dy, lonlims=lonlims, latlims=latlims
            )
            (
                obsSum,
                sqObsSum,
                devSum,
                sqDevSum,
                mdlByObsSum,
                dtcount,
            ) = mll.computeCumDeviations(maplons, maplats, mapdata)
            obsMax, mdlMax = mll.computeMaxima(maplons, maplats, mapdata)
            return (
                maplons,
                maplats,
                obsSum,
                sqObsSum,
                devSum,
                sqDevSum,
                mdlByObsSum,
                dtcount,
                obsMax,
                mdlMax,
            )

        for f in fls:
            fpth = os.path.join(hsSatAndModelDir, f)
            satdts_ = loadFile(fpth)
            dts = satdts_[:, 0]
            currYr = datetime.fromtimestamp(dts[0]).year
            if loopYear is None:
                loopYear = currYr
                years.append(loopYear)
                satdts = None
                print("  loading year " + str(loopYear))
            if loopYear == currYr:
                satdts = (
                    satdts_ if satdts is None else np.concatenate([satdts, satdts_], 0)
                )
            else:
                # loopYear != currYr
                # yielding the current year
                dts = satdts[:, 0]
                lons = satdts[:, 1]
                lats = satdts[:, 2]
                saths = satdts[:, 3]
                modhs = satdts[:, 4]
#                yield _elab_ssh(dts, lons, lats, saths, modhs)
                yield _elab_ssh(dts, lons, lats, saths, modhs)
                fnrinfri
                satdts = satdts_
                loopYear = currYr
                years.append(loopYear)
                print("  loading year " + str(loopYear))
        # yielding the last year
        dts = satdts[:, 0]
        lons = satdts[:, 1]
        lats = satdts[:, 2]
        saths = satdts[:, 3]
        modhs = satdts[:, 4]
        yield _elab_ssh(dts, lons, lats, saths, modhs)

    #(
    #    obsSum,
    #    sqObsSum,
    #    devSum,
    #    sqDevSum,
    #    dtcount,
    #    obsMaxSum,
    #    mdlMaxSum,
    #    obsTotMax,
    #    mdlTotMax,
    #) = (None, None, None, None, None, None, None, None, None)
    for blob in iterateByYear():
        (
            maplons,
            maplats,
            mp,
           # _sqObsSum,
           # _devSum,
           # _sqDevSum,
           # _mdlByObsSum,
           # _dtcount,
           # _obsMax,
           # _mdlMax,
        ) = blob
      #  print(blob)
        #if _obsSum is None:
        #    continue
        #if obsSum is None:
        #    (
        #        obsSum,
        #        sqObsSum,
        #        devSum,
        #        sqDevSum,
        #        mdlByObsSum,
        #        dtcount,
        #        obsMaxSum,
        #        mdlMaxSum,
        #        obsTotMax,
        #        mdlTotMax,
        #    ) = (
        #        _obsSum,
        #        _sqObsSum,
        #        _devSum,
        #        _sqDevSum,
        #        _mdlByObsSum,
        #        _dtcount,
        #        _obsMax,
        #        _mdlMax,
        #        _obsMax,
        #        _mdlMax,
        #    )
        #else:
        #    obsSum = np.nansum([obsSum, _obsSum], 0)
        #    sqObsSum = np.nansum([sqObsSum, _sqObsSum], 0)
        #    devSum = np.nansum([devSum, _devSum], 0)
        #    sqDevSum = np.nansum([sqDevSum, _sqDevSum], 0)
        #    mdlByObsSum = np.nansum([mdlByObsSum, _mdlByObsSum], 0)
        #    dtcount = np.nansum([dtcount, _dtcount], 0)
        #    obsMaxSum = np.nansum([obsMaxSum, _obsMax], 0)
        #    mdlMaxSum = np.nansum([mdlMaxSum, _mdlMax], 0)
        #    obsTotMax = np.nanmax([obsTotMax, _obsMax], 0)
        #    mdlTotMax = np.nanmax([mdlTotMax, _mdlMax], 0)

    if not latlims is None:
        cnd = np.logical_and(maplats >= latlims[0], maplats <= latlims[1])
        maplats = maplats[cnd]
        #obsSum = obsSum[cnd, :]
        #sqObsSum = sqObsSum[cnd, :]
        #devSum = devSum[cnd, :]
        #sqDevSum = sqDevSum[cnd, :]
        #mdlByObsSum = mdlByObsSum[cnd, :]
        #dtcount = dtcount[cnd, :]
        #obsMaxSum = obsMaxSum[cnd, :]
        #mdlMaxSum = mdlMaxSum[cnd, :]
        #obsTotMax = obsTotMax[cnd, :]
        #mdlTotMax = mdlTotMax[cnd, :]

    if maskPointsCloseToTheCoast:
        msk = mll.createCoastlinePointsMask(
            maplons, maplats, resl="h", minDistFromCoast=minDistFromCoast
        )
        cnd = msk == 0
        obsSum[cnd] = np.nan
        sqObsSum[cnd] = np.nan
        devSum[cnd] = np.nan
        sqDevSum[cnd] = np.nan
        mdlByObsSum[cnd] = np.nan
        dtcount[cnd] = np.nan
        obsMaxSum[cnd] = np.nan
        mdlMaxSum[cnd] = np.nan
        obsTotMax[cnd] = np.nan
        mdlTotMax[cnd] = np.nan

    if maskShallowPoints:
        msk = mll.createBathyMask(maplons, maplats, minSeaDepth, bathyFile)
        cnd = msk == 0
        obsSum[cnd] = np.nan
        sqObsSum[cnd] = np.nan
        devSum[cnd] = np.nan
        sqDevSum[cnd] = np.nan
        mdlByObsSum[cnd] = np.nan
        dtcount[cnd] = np.nan
        obsMaxSum[cnd] = np.nan
        mdlMaxSum[cnd] = np.nan
        obsTotMax[cnd] = np.nan
        mdlTotMax[cnd] = np.nan


    # Compute rmse and extract mean from data
    deviation_ = np.ones((len(maplats), len(maplons))) * np.nan
    _consideredCells = 0
    data_sat_ = []
    data_sat = []
    data_mod_ = []
    data_mod = []
    for ix in range(len(maplons)):
        for iy in range(len(maplats)):
            data = mp.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3])
            data_sat_.append(msrs)
            mods = np.array(data[4])
            data_mod_.append(mods)

            msrs_sum = np.sum(np.array(data[3]))
            mods_sum = np.sum(np.array(data[4]))

            deviation_[iy, ix] = np.sum((mods - msrs) ** 2.0)
#            print(deviation_[iy, ix])
#            print(iy, ix)
            _consideredCells += 1

    deviation = deviation_/(_consideredCells)
    sat_mean = msrs_sum / _consideredCells
    mod_mean = mods_sum / _consideredCells

    print("considered cells: ", _consideredCells)

    # flatting array
    for i in range(len(data_sat_)):
      data_sat.extend(data_sat_[i])
    for i in range(len(data_mod_)):
      data_mod.extend(data_mod_[i])


    pth_data_mod = np.percentile(data_mod, pth)
    pth_data_sat = np.percentile(data_sat, pth)

    filter_data_sat = []
    filter_data_mod = []

    # compute NSE and Coefficient of determination
    nse_ = 0
    ssres = 0
    sstot = 0
    _consideredCellsFilter = 0
    for ix in range(len(maplons)):
        for iy in range(len(maplats)):
            data = mp.get((ix, iy))
            if not data:
                continue
            msrs = np.array(data[3]) - sat_mean
            mods = np.array(data[4]) - mod_mean
            for i in range(len(msrs)):
                if msrs[i] >= pth_data_sat and mods[i] >= pth_data_sat:
                    msrs_ = np.sum(msrs[i])
                    mods_ = np.sum(mods[i])
                    msrs_sum = np.sum(msrs_)
                    mods_sum = np.sum(mods_)

                    # in this case, msrs_sum mean is zero? or is the mean before substract the mean?
                    nse_ = np.sum((mods - msrs)**2)/np.sum((msrs_sum - sat_mean)**2)
                    ssres = np.sum(mods - msrs)
                    sstot = np.sum(msrs - sat_mean)
                    #_consideredCellsFilter += 1
                else:
                    continue

#    print("considered cells after percentile filtering: ", _consideredCellsFilter)
    nse = 1 - nse_
    r2 = 1 - ssres/sstot
    print("NSE = ", nse)
    print("R2 = ", r2)

    #cnd = dtcount < minObsForValidation
    #obsSum[cnd] = np.nan
    #sqObsSum[cnd] = np.nan
    #devSum[cnd] = np.nan
    #sqDevSum[cnd] = np.nan
    #mdlByObsSum[cnd] = np.nan
    #dtcount[cnd] = np.nan
    #obsMaxSum[cnd] = np.nan
    #mdlMaxSum[cnd] = np.nan
    #obsTotMax[cnd] = np.nan
    #mdlTotMax[cnd] = np.nan


    rmseTot = np.sqrt(np.nansum(deviation))
    totIndStr = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOTAL ERROR INDICATORS:
rmse: {rmseTot:2.5f}
NSE: {nse:2.5f}
R2: {r2:2.5f}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
    totIndStr = totIndStr.format(
        rmseTot=rmseTot,
        nse=nse,
        r2=r2,
    )
    print("")
    print(totIndStr)
    print("")
    totIndsFilePath = os.path.join(outputDir, "totalIndicators.txt")
    with open(totIndsFilePath, "w") as f:
        f.write(totIndStr)
        f.close()

    rmse = np.sqrt(deviation)

    # saving to files
    try:
        os.makedirs(outputDir)
    except:
        pass
    np.savetxt(os.path.join(outputDir, "lons.csv"), maplons)
    np.savetxt(os.path.join(outputDir, "lats.csv"), maplats)
    np.savetxt(os.path.join(outputDir, "rmse.csv"), rmse)

    print("output dir: " + outputDir)
