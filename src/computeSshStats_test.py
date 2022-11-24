import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools

import src.coarsenSatData


import src.mapByLonLat as mll

import src.utils as utils


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
    pth=90,
    nminobs = 2,
):

    years = []
    satssh = []
    modssh = []
    dts = []
    Lons = []
    Lats = []
    Dts = []
    Maplons = []
    Maplats = []
    Mapdata = []

    # def iterateByYear():
    #     fls_ = [f for f in os.listdir(hsSatAndModelDir) if re.match("(.*)\.npy", f)]
    #     pttrn = "(.*)_hsModelAndSatObs_([0-9]{8}).npy"
    #     dtstrs = [re.match(pttrn, f).groups(0)[1] for f in fls_]
    #     dts = [datetime.strptime(dtstr, "%Y%m%d") for dtstr in dtstrs]
    #     iii = np.argsort(dts)
    #     dts = np.array(dts)[iii]
    #     fls_ = np.array(fls_)[iii]
    #     fls = [
    #         fls_[i]
    #         for i in range(len(fls_))
    #         if (startDate <= dts[i] and endDate >= dts[i])
    #     ]

    #     loopYear = None

    #     def loadFile(flpth):
    #         print("")
    #         print("    loading file " + flpth)
    #         satdts = np.load(flpth)
    #         if filterHighSsh:
    #             sshsat = satdts[:, 3]
    #             sshmdl = satdts[:, 4]
    #             cnd = np.logical_and(
    #                 sshsat < filterSshMaximum, sshmdl < filterSshMaximum
    #             )
    #             satdts = satdts[cnd, :]
    #         return satdts

    #     def _elab(dts, lons, lats, saths, modhs):
    #         maplons, maplats, mapdata = mll.mapByLonLat(dts, lons, lats, saths, modhs, dx, dy, lonlims=lonlims, latlims=latlims)
    #         Maplons.append(maplons)
    #         Maplats.append(maplats)
    #         Mapdata.append(mapdata)
    #         #return maplons, maplats, mapdata
    #         # (
    #         #     obsSum,
    #         #     sqObsSum,
    #         #     sqModSum,
    #         #     devSum,
    #         #     sqDevSum,
    #         #     mdlByObsSum,
    #         #     dtcount,
    #         # ) = mll.computeCumDeviations(maplons, maplats, mapdata)
    #         # obsMax, mdlMax = mll.computeMaxima(maplons, maplats, mapdata)
    #         # return (
    #         #     maplons,
    #         #     maplats,
    #         #     obsSum,
    #         #     sqObsSum,
    #         #     devSum,
    #         #     sqDevSum,
    #         #     mdlByObsSum,
    #         #     dtcount,
    #         #     obsMax,
    #         #     mdlMax,
    #         # )

    #     for f in fls:
    #         fpth = os.path.join(hsSatAndModelDir, f)
    #         satdts_ = loadFile(fpth)
    #         dts = satdts_[:, 0]
    #         currYr = datetime.fromtimestamp(dts[0]).year
    #         if loopYear is None:
    #             loopYear = currYr
    #             years.append(loopYear)
    #             satdts = None
    #             print("  loading year " + str(loopYear))
    #         if loopYear == currYr:
    #             satdts = (
    #                 satdts_ if satdts is None else np.concatenate([satdts, satdts_], 0)
    #             )
    #         else:
    #             # loopYear != currYr
    #             # yielding the current year
    #             dts_ = satdts[:, 0]
    #             lons_ = satdts[:, 1]
    #             lats_ = satdts[:, 2]
    #             satssh_ = satdts[:, 3]
    #             modssh_ = satdts[:, 4]

    #             yield _elab(dts_, lons_, lats_, satssh_, modssh_)
    #             satdts = satdts_
    #             loopYear = currYr
    #             years.append(loopYear)
    #             print("  loading year " + str(loopYear))
    #     # yielding the last year
    #     dts = satdts[:, 0]
    #     lons = satdts[:, 1]
    #     lats = satdts[:, 2]
    #     satssh_ = satdts[:, 3]
    #     modssh_ = satdts[:, 4]

    #     satssh.append(satssh_)
    #     modssh.append(modssh_)
    #     Lons.append(lons)
    #     Lats.append(lats)
    #     Dts.append(dts)

    #     yield _elab(dts, lons, lats, satssh_, modssh_)

    # (
    #     obsSum,
    #     sqObsSum,
    #     devSum,
    #     sqDevSum,
    #     dtcount,
    #     obsMaxSum,
    #     mdlMaxSum,
    #     obsTotMax,
    #     mdlTotMax,
    # ) = (None, None, None, None, None, None, None, None, None)
    # for blob in iterateByYear():
    #     pass
    # #     (
    # #         maplons,
    # #         maplats,
    # #         _obsSum,
    # #         _sqObsSum,
    # #         _devSum,
    # #         _sqDevSum,
    # #         _mdlByObsSum,
    # #         _dtcount,
    # #         _obsMax,
    # #         _mdlMax,
    # #     ) = blob
    # #     if _obsSum is None:
    # #         continue
    # #     if obsSum is None:
    # #         (
    # #             obsSum,
    # #             sqObsSum,
    # #             devSum,
    # #             sqDevSum,
    # #             mdlByObsSum,
    # #             dtcount,
    # #             obsMaxSum,
    # #             mdlMaxSum,
    # #             obsTotMax,
    # #             mdlTotMax,
    # #         ) = (
    # #             _obsSum,
    # #             _sqObsSum,
    # #             _devSum,
    # #             _sqDevSum,
    # #             _mdlByObsSum,
    # #             _dtcount,
    # #             _obsMax,
    # #             _mdlMax,
    # #             _obsMax,
    # #             _mdlMax,
    # #         )
    # #     else:
    # #         obsSum = np.nansum([obsSum, _obsSum], 0)
    # #         sqObsSum = np.nansum([sqObsSum, _sqObsSum], 0)
    # #         devSum = np.nansum([devSum, _devSum], 0)
    # #         sqDevSum = np.nansum([sqDevSum, _sqDevSum], 0)
    # #         mdlByObsSum = np.nansum([mdlByObsSum, _mdlByObsSum], 0)
    # #         dtcount = np.nansum([dtcount, _dtcount], 0)
    # #         obsMaxSum = np.nansum([obsMaxSum, _obsMax], 0)
    # #         mdlMaxSum = np.nansum([mdlMaxSum, _mdlMax], 0)
    # #         obsTotMax = np.nanmax([obsTotMax, _obsMax], 0)
    # #         mdlTotMax = np.nanmax([mdlTotMax, _mdlMax], 0)

    # # if not latlims is None:
    # #     cnd = np.logical_and(maplats >= latlims[0], maplats <= latlims[1])
    # #     maplats = maplats[cnd]
    # #     obsSum = obsSum[cnd, :]
    # #     sqObsSum = sqObsSum[cnd, :]
    # #     devSum = devSum[cnd, :]
    # #     sqDevSum = sqDevSum[cnd, :]
    # #     mdlByObsSum = mdlByObsSum[cnd, :]
    # #     dtcount = dtcount[cnd, :]
    # #     obsMaxSum = obsMaxSum[cnd, :]
    # #     mdlMaxSum = mdlMaxSum[cnd, :]
    # #     obsTotMax = obsTotMax[cnd, :]
    # #     mdlTotMax = mdlTotMax[cnd, :]

    # # if maskPointsCloseToTheCoast:
    # #     msk = mll.createCoastlinePointsMask(
    # #         maplons, maplats, resl="h", minDistFromCoast=minDistFromCoast
    # #     )
    # #     cnd = msk == 0
    # #     obsSum[cnd] = np.nan
    # #     sqObsSum[cnd] = np.nan
    # #     devSum[cnd] = np.nan
    # #     sqDevSum[cnd] = np.nan
    # #     mdlByObsSum[cnd] = np.nan
    # #     dtcount[cnd] = np.nan
    # #     obsMaxSum[cnd] = np.nan
    # #     mdlMaxSum[cnd] = np.nan
    # #     obsTotMax[cnd] = np.nan
    # #     mdlTotMax[cnd] = np.nan

    # # if maskShallowPoints:
    # #     msk = mll.createBathyMask(maplons, maplats, minSeaDepth, bathyFile)
    # #     cnd = msk == 0
    # #     obsSum[cnd] = np.nan
    # #     sqObsSum[cnd] = np.nan
    # #     devSum[cnd] = np.nan
    # #     sqDevSum[cnd] = np.nan
    # #     mdlByObsSum[cnd] = np.nan
    # #     dtcount[cnd] = np.nan
    # #     obsMaxSum[cnd] = np.nan
    # #     mdlMaxSum[cnd] = np.nan
    # #     obsTotMax[cnd] = np.nan
    # #     mdlTotMax[cnd] = np.nan

    def getFiles(hsSatAndModelDir, startDate, endDate):
        fl = []

        while startDate <= endDate:
            strtime = startDate.strftime("%Y%m%d")
            pthfile = hsSatAndModelDir + "/ERA5_schismwwm_" + strtime + "*.npy"
            fileFound = glob.glob(pthfile)
            fl.append(fileFound)
            startDate += timedelta(days=1)

        return list(itertools.chain(*fl))

    def loadFile(flpth):
        print("    loading file " + flpth)
        satdts = np.load(flpth)
        return satdts

    fls = getFiles(hsSatAndModelDir, startDate, endDate)

    obs = np.array([])
    model = np.array([])

    r2lst = []
    nselst = []
    ablst = []
    rblst = []
    rmselst = []
    nrmselst = []
    pearsonlst = []

    # looping on tidal gauge files
    obs = np.array([])
    model = np.array([])
    dts = np.array([])
    mapdata = {}

    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:, 3]
        model_ = data_[:, 4]
        model = np.concatenate((model, model_), axis=0)
        obs = np.concatenate((obs, obs_), axis=0)

    fltr = np.logical_and(obs >= -10, obs <= 10)
    obs=obs[fltr]
    model=model[fltr]
    meanobs = np.nanmean(obs)
    meanmodel = np.nanmean(model)

    for f in fls:
        data_ = loadFile(f)
        dts_ = data_[:, 0]
        lons = data_[:, 1]
        lats = data_[:, 2]
        obs_ = data_[:, 3]-meanobs
        model_ = data_[:, 4]-meanmodel

        maplons, maplats, mapdata = mll.mapByLonLatCumm(
            mapdata,
            dts_,
            lons,
            lats,
            obs_,
            model_,
            dx,
            dy,
            lonlims=lonlims,
            latlims=latlims,
        )
        # mapdata = mll.computeMean_cell(lons, lats, mapdata, mapdata)

        # model = np.concatenate((model, model_), axis=0)
        # obs = np.concatenate((obs, obs_), axis=0)
        # dts = np.concatenate((dts, dts_), axis=0)

    print(meanobs, meanmodel)

    fig, ax = plt.subplots()
    ix = 16
    iy = 17
    data = mapdata.get((ix, iy))
    obs = np.array(data[3])
    model = np.array(data[4])
    print(np.mean(mapdata.get((ix, iy))[1]), np.mean(mapdata.get((ix, iy))[2]))
    r2_, nse_, ab_, rb_, rmse_, nrmse_, pearson_ = utils.computeStats(obs, model, pth)
    print("r2 = ", r2_)
    print("pearson = ", pearson_)
    print("bias = ", ab_)
    print("rmse = ", rmse_)
    ax.plot(obs, label="observation")
    ax.plot(model, label="model", alpha=0.5)
    ax.legend()
    plt.savefig("testtt.png")


    r2 = np.ones((len(maplats), len(maplons))) * 99999
    nse = np.ones((len(maplats), len(maplons))) * 99999
    nrmse = np.ones((len(maplats), len(maplons))) * 99999
    bias = np.ones((len(maplats), len(maplons))) * 99999
    rmse = np.ones((len(maplats), len(maplons))) * 99999
    pearson = np.ones((len(maplats), len(maplons))) * 99999
    relbias = np.ones((len(maplats), len(maplons))) * 99999

    for ix in range(len(lons)):
        for iy in range(len(lats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            obs = np.array(data[3])
            model = np.array(data[4])

            if len(model) <= nminobs:
                continue

            r2_, nse_, absBias_, rb_, rmse_, nrmse_, pearson_ = utils.computeStats(
                obs, model, pth
            )
            """
            if (r2_ > 0.7) & (pearson_ > 0.9):
                print("ix = ", ix, "iy = ", iy)
                fjrijfri
            """

            r2[iy, ix] = r2_
            nse[iy, ix] = nse_
            bias[iy, ix] = absBias_
            relbias[iy, ix] = rb_
            rmse[iy, ix] = rmse_
            nrmse[iy, ix] = nrmse_
            pearson[iy, ix] = pearson_

    mask = bias == 99999
    bias = np.ma.masked_array(bias, mask)
    nrmse = np.ma.masked_array(nrmse, mask)
    relbias = np.ma.masked_array(relbias, mask)
    r2 = np.ma.masked_array(r2, mask)
    pearson = np.ma.masked_array(pearson, mask)
    nse = np.ma.masked_array(nse, mask)

    # modssh_mean = np.nanmean(modssh_)
    # satssh_mean = np.nanmean(satssh_)

    # satssh = satssh_ - satssh_mean
    # modssh = modssh_  # - modssh_mean

    # r2, nse, ab, rb, rmse, nrmse, pearson = utils.computeStats(satssh, modssh, pth)
    nnse = 1 / (2 - nse)
    nr2 = 1 / (2 - r2)

    abTot = np.nanmean(bias)
    relabTot = np.nanmean(relbias)
    rmseTot = np.nanmean(rmse)
    nrmseTot = np.nanmean(nrmse)
    pearsonTot = np.nanmean(pearson)
    nseTot = np.nanmean(nse)
    nnseTot = np.nanmean(nnse)
    r2Tot = np.nanmean(r2)
    nr2Tot = np.nanmean(nr2)

    totIndStr = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOTAL ERROR INDICATORS:
rmse: {rmseTot:2.5f}
NSE: {nseTot:2.5f}
R2: {r2Tot:2.5f}
NNSE: {nnseTot:2.5f}
NR2: {nr2Tot:2.5f}
Bias: {abTot:2.5f}
RelBias: {relabTot:2.5f}
Pearson: {pearsonTot:2.5f}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
    totIndStr = totIndStr.format(
        rmseTot=rmseTot,
        nseTot=nseTot,
        r2Tot=r2Tot,
        nnseTot=nnseTot,
        nr2Tot=nr2Tot,
        abTot=abTot,
        relabTot=relabTot,
        pearsonTot=pearsonTot,
    )
    print("")
    print(totIndStr)
    print("")

    statFile = (
        "ssh-altimeter_"
        + "pth_"
        + str(pth)
        + "_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt"
    )
    totIndsFilePath = os.path.join(outputDir, statFile)
    with open(totIndsFilePath, "w") as f:
        f.write(totIndStr)
        f.close()

    # rmse = np.sqrt(sqDevSum)
    # bias = devSum / obsSum
    # absBias = devSum / dtcount
    # nrmse = np.sqrt(sqDevSum / sqObsSum)
    # hh = np.sqrt(sqDevSum / mdlByObsSum)

    # saving to files
    try:
        os.makedirs(outputDir)
    except:
        pass

    # lons = Maplons[0]
    # lats = Maplats[0]
    # nminobs = 0

    # ModData = np.ones((len(lats), len(lons))) * np.nan
    # ObsData = ModData

    # modd = []
    # obss = []

    # # initalize new map
    # mp = {}
    # # for ix in range(len(lons)):
    # #     for iy in range(len(lats)):
    # #         lst = mp.get((ix, iy), [[], []])
    # #         mp[(ix, iy)] = lst

    # # store all data´
    # _obs = []
    # _mod = []
    # for mapdata in Mapdata:
    #     for ix in range(len(lons)):
    #         for iy in range(len(lats)):
    #             data = mapdata.get((ix, iy))
    #             if not data:
    #                 continue
    #             msrs = data[3]
    #             mods = data[4]
    #             if len(mods) < nminobs:
    #                 continue

    #             assert len(mods) == len(msrs)

    #             lst = mp.get((ix, iy), [[], []])
    #             _obs.append(msrs)
    #             _mod.append(mods)
    #             lst[0] = _obs
    #             lst[1] = _mod
    #             mp[(ix, iy)] = lst

    # Mod = np.ones((len(lats), len(lons))) * np.nan
    # Obs = np.ones((len(lats), len(lons))) * np.nan

    # for ix in range(len(lons)):
    #     for iy in range(len(lats)):

    #         data = mp.get((ix, iy))
    #         if not data:
    #             continue

    #         obs = [item for sublist in data[0] for item in sublist]
    #         mod = [item for sublist in data[1] for item in sublist]

    #         Mod[iy, ix] = np.array(mod)[0]
    #         Obs[iy, ix] = np.array(obs)[0]

    # print(Mod[:,0])
    # fjrifir

    np.savetxt(os.path.join(outputDir, "lons-ssh.csv"),maplons)
    np.savetxt(os.path.join(outputDir, "lats-ssh.csv"), maplats)
    # np.savetxt(os.path.join(outputDir, "rmse-ssh-altimeter_"+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".csv"), rmse)
    # np.savetxt(os.path.join(outputDir, "bias-ssh-altimeter_"+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".csv"), bias)

    np.savetxt(
        os.path.join(
            outputDir,
            "SSH-rmse-altimeter_"
            + "pth_"
            + str(pth)
            + "_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        rmse,
    )

    np.savetxt(
        os.path.join(
            outputDir,
            "SSH-pearson-altimeter_"
            + "pth_"
            + str(pth)
            + "_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        pearson,
    )

    np.savetxt(
        os.path.join(
            outputDir,
            "SSH-r2-altimeter_"
            + "pth_"
            + str(pth)
            + "_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        r2,
    )

    np.savetxt(
        os.path.join(
            outputDir,
            "SSH-bias-altimeter_"
            + "pth_"
            + str(pth)
            + "_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        bias,
    )

    # np.savetxt(os.path.join(outputDir, "dtcount.csv"), dtcount)

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
