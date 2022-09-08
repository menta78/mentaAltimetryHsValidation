import os, re
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt

import src.coarsenSatData


import src.mapByLonLat as mll


maskPointsCloseToTheCoast = False
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
    filterHighSsh=False,
    filterSshMaximum=0.0,
    dx=0.2,
    dy=0.2,
    pth = 90,
):

    years = []
    satssh = []
    modssh = []

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
            if filterHighSsh:
                sshsat = satdts[:, 3]
                sshmdl = satdts[:, 4]
                cnd = np.logical_and(
                    sshsat < filterSshMaximum, sshmdl < filterSshMaximum
                )
                satdts = satdts[cnd, :]
            return satdts

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
                satssh_ = satdts[:, 3]
                satssh.append(satssh_)
                modssh_ = satdts[:, 4]
                modssh.append(modssh_)

                yield _elab(dts, lons, lats, satssh_, modssh_)
                satdts = satdts_
                loopYear = currYr
                years.append(loopYear)
                print("  loading year " + str(loopYear))
        # yielding the last year
        dts = satdts[:, 0]
        lons = satdts[:, 1]
        lats = satdts[:, 2]
        satssh_ = satdts[:, 3]
        satssh.append(satssh_)
        modssh_ = satdts[:, 4]
        modssh.append(modssh_)
        yield _elab(dts, lons, lats, satssh_, modssh_)

    (
        obsSum,
        sqObsSum,
        devSum,
        sqDevSum,
        dtcount,
        obsMaxSum,
        mdlMaxSum,
        obsTotMax,
        mdlTotMax,
    ) = (None, None, None, None, None, None, None, None, None)
    for blob in iterateByYear():
        (
            maplons,
            maplats,
            _obsSum,
            _sqObsSum,
            _devSum,
            _sqDevSum,
            _mdlByObsSum,
            _dtcount,
            _obsMax,
            _mdlMax,
        ) = blob
        if _obsSum is None:
            continue
        if obsSum is None:
            (
                obsSum,
                sqObsSum,
                devSum,
                sqDevSum,
                mdlByObsSum,
                dtcount,
                obsMaxSum,
                mdlMaxSum,
                obsTotMax,
                mdlTotMax,
            ) = (
                _obsSum,
                _sqObsSum,
                _devSum,
                _sqDevSum,
                _mdlByObsSum,
                _dtcount,
                _obsMax,
                _mdlMax,
                _obsMax,
                _mdlMax,
            )
        else:
            obsSum = np.nansum([obsSum, _obsSum], 0)
            sqObsSum = np.nansum([sqObsSum, _sqObsSum], 0)
            devSum = np.nansum([devSum, _devSum], 0)
            sqDevSum = np.nansum([sqDevSum, _sqDevSum], 0)
            mdlByObsSum = np.nansum([mdlByObsSum, _mdlByObsSum], 0)
            dtcount = np.nansum([dtcount, _dtcount], 0)
            obsMaxSum = np.nansum([obsMaxSum, _obsMax], 0)
            mdlMaxSum = np.nansum([mdlMaxSum, _mdlMax], 0)
            obsTotMax = np.nanmax([obsTotMax, _obsMax], 0)
            mdlTotMax = np.nanmax([mdlTotMax, _mdlMax], 0)

    if not latlims is None:
        cnd = np.logical_and(maplats >= latlims[0], maplats <= latlims[1])
        maplats = maplats[cnd]
        obsSum = obsSum[cnd, :]
        sqObsSum = sqObsSum[cnd, :]
        devSum = devSum[cnd, :]
        sqDevSum = sqDevSum[cnd, :]
        mdlByObsSum = mdlByObsSum[cnd, :]
        dtcount = dtcount[cnd, :]
        obsMaxSum = obsMaxSum[cnd, :]
        mdlMaxSum = mdlMaxSum[cnd, :]
        obsTotMax = obsTotMax[cnd, :]
        mdlTotMax = mdlTotMax[cnd, :]

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


    # flatting array
    for i in range(len(satssh)):
      satssh.extend(satssh[i])
      modssh.extend(modssh[i])

    satssh_ = np.array(satssh)[0]
    modssh_ = np.array(modssh)[0]

    modssh_mean = np.nanmean(modssh_)
    satssh_mean = np.nanmean(satssh_)

    satssh = satssh_  - satssh_mean
    modssh = modssh_  - modssh_mean

    # compute percentiles.
    pobs = np.nanpercentile(satssh, pth)
    pmod = np.nanpercentile(modssh, pth)
    print("pth observations = ", pobs)
    print("pth model = ", pmod)

    condition = satssh>=0
    meanObs = np.nanmean(satssh, where= condition)

    condition1 = satssh >= pobs 
    condition2 = modssh >= pmod
    conditionPth = condition1 & condition2
    print(condition1, condition2, conditionPth)

    cnd = dtcount < minObsForValidation
    conditionNumberSamples = condition 

    sqDevSum[cnd] = np.nan
    dtcount[cnd] = np.nan

    N = np.nansum(conditionNumberSamples)
    rmseTot = np.sqrt(np.nansum(sqDevSum) / N, where = conditionPth)

    # computing skills
    """
    filt_satssh = []
    filt_modssh = []
    for i in range(len(satssh)):
        if satssh[i] >= pth_satssh and modssh[i] >= pth_modssh:
            print(satssh[i],  modssh[i])
            filt_satssh.append(satssh[i])
            filt_modssh.append(modssh[i])

    filt_data_sat = np.array(filt_satssh)
    filt_data_mod = np.array(filt_modssh)

#    mean_filt_data_sat = np.nanmean(filt_data_sat)
#    print(mean_filt_data_sat)
    mean_filt_data_sat = 0
    N = len(filt_data_sat)
    nse1 = 0
    nse2 = 0
    ssres = 0
    sstot = 0

    for i in range(N):
        nse1 += (filt_data_mod[i] - filt_data_sat[i])**2
        nse2 += (filt_data_sat[i] - mean_filt_data_sat)**2
        ssres += (filt_data_mod[i] - filt_data_sat[i])**2
        sstot += (filt_data_sat[i] - mean_filt_data_sat)**2
    """

    nsc1 = np.nansum(np.abs(satssh-modssh), where=conditionPth)
    nsc2 = np.nansum(np.abs(satssh-np.nanmean(satssh)), where=conditionPth)

    ssres = np.nansum((satssh-modssh)**2, where=conditionPth)
    sstot = np.nansum((satssh-np.nanmean(satssh))**2, where=conditionPth)

    absre_ = np.nansum(modssh-satssh, where=conditionPth)
    nobs = np.nansum(satssh, where=conditionPth)

    print(absre_, nobs)

    absre = absre_/N
    nse  = 1 - nsc1/nsc2
    nnse = 1/(2-nse)
    r2   = 1 - ssres/sstot
    nr2 = 1/(2-r2)
    reb = absre_/nobs * 100
    rmseTot = np.sqrt(ssres/N)
    totIndStr = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOTAL ERROR INDICATORS:
rmse: {rmseTot:2.5f}
NSE: {nse:2.5f}
R2: {r2:2.5f}
NNSE: {nnse:2.5f}
NR2: {nr2:2.5f}
Bias: {absre:2.5f}
RelBias: {reb:2.5f}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
    totIndStr = totIndStr.format(
        rmseTot = rmseTot,
        nse     = nse,
        r2      = r2,
        nnse     = nnse,
        nr2      = nr2,
        absre = absre,
        reb = reb,
    )
    print("")
    print(totIndStr)
    print("")
    totIndsFilePath = os.path.join(outputDir, "totalIndicators.txt")
    with open(totIndsFilePath, "w") as f:
        f.write(totIndStr)
        f.close()

    rmse = np.sqrt(sqDevSum / N)

    # saving to files
    try:
        os.makedirs(outputDir)
    except:
        pass
    np.savetxt(os.path.join(outputDir, "lons.csv"), maplons)
    np.savetxt(os.path.join(outputDir, "lats.csv"), maplats)
    np.savetxt(os.path.join(outputDir, "rmse.csv"), rmse)
    np.savetxt(os.path.join(outputDir, "dtcount.csv"), dtcount)

    print("output dir: " + outputDir)


if __name__ == "__main__":
    startYear = 2000
    endYear = 2009

    # RED SEA 10 YEARS
    lonlims = [30, 43.8]
    latlims = [13.4, 32]
    # lonlims = [36.5, 43.8]
    # latlims = [13.4, 21]
    hsSatAndModelDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst/hsModelAndSatObs/"
    outputDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst/stats/"
    hsSatAndModelDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/hsModelAndSatObs/"
    outputDir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/stats/"

    # elaborateMeasures(startYear, endYear, hsSatAndModelDir, outputDir, lonlims=lonlims, latlims=latlims)

    """
  # PERSIC GULF 10 YEARS
  lonlims = [47, 56.5]
  latlims = [23.5, 31]
 #lonlims = [36.5, 43.8]
 #latlims = [13.4, 21]
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst/stats/'
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst_noobst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst_noobst/stats/'
  """

    elaborateMeasures(
        startYear,
        endYear,
        hsSatAndModelDir,
        outputDir,
        lonlims=lonlims,
        latlims=latlims,
    )
