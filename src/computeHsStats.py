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
    nminobs = 1,
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

    ablst = []
    nrmselst = []
    hhlst = []
    nbilst = []

    # looping on tidal gauge files
    obs = np.array([])
    model = np.array([])
    dts = np.array([])
    mapdata = {}

    print(fls)
    for f in fls:
        data_ = loadFile(f)
        dts_ = data_[:, 0]
        lons = data_[:, 1]
        lats = data_[:, 2]
        obs_ = data_[:, 3]
        model_ = data_[:, 4]

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


    r2lst = []
    nselst = []
    ablst = []
    rblst = []
    rmselst = []
    nrmselst = []
    pearsonlst = []

    

    mpabBias = {}
    mpnrmse = {}
    mpnbi = {}
    mphh = {}

    bias = np.ones((len(maplats), len(maplons))) * 99999
    nrmse = np.ones((len(maplats), len(maplons))) * 99999
    nbi = np.ones((len(maplats), len(maplons))) * 99999
    hh = np.ones((len(maplats), len(maplons))) * 99999
    for ix in range(len(maplons)):
        for iy in range(len(maplats)):
            data = mapdata.get((ix, iy))
            if not data:
                continue
            obs = np.array(data[3])
            model = np.array(data[4])

            if len(model) <= nminobs:
                continue

            _nbi, _absBias, _nrmse, _hh = utils.computeStatsHs(obs, model, pth)

            bias[iy, ix] = _absBias
            nbi[iy, ix] = _nbi
            hh[iy, ix] = _hh
            nrmse[iy, ix] = _nrmse

            # lst = mpabBias.get((ix, iy), []) # return empty if value doesn't exists
            # lst = nbi
            # mpabBias[(ix, iy)] = lst

            # lst = mpnbi.get((ix, iy), []) # return empty if value doesn't exists
            # lst = nbi
            # mpnbi[(ix, iy)] = lst

            # lst = mphh.get((ix, iy), []) # return empty if value doesn't exists
            # lst = hh
            # mphh[(ix, iy)] = lst

            # lst = mpnrmse.get((ix, iy), []) # return empty if value doesn't exists
            # lst = nrmse
            # mpnrmse[(ix, iy)] = lst

            ablst.append(_absBias)
            nrmselst.append(_nrmse)
            nbilst.append(_nbi)
            hhlst.append(_hh)

    abArray = np.array(ablst)
    nrmseArray = np.array(nrmselst)
    nbiArray = np.array(nbilst)
    hhArray = np.array(hhlst)

    abTot = np.nanmean(abArray)
    hhTot = np.nanmean(hhArray)
    nbiTot = np.nanmean(nbiArray)
    nrmseTot = np.nanmean(nrmseArray)

    # modssh_mean = np.nanmean(modssh_)
    # satssh_mean = np.nanmean(satssh_)

    # satssh = satssh_ - satssh_mean
    # modssh = modssh_  # - modssh_mean

    print("=========================================")
    print("Start date ", startDate, "   :::::::   End date", endDate)
    print("Percentile ", pth)
    print("nrmse = ", nrmseTot)
    print("abs bias = ", abTot)
    print("hh = ", hhTot)
    print("nbi =", nbiTot)
    print("=========================================")

    statFile = (
        "hs-altimeter_"
        + startDate.strftime("%Y%m%d")
        + "_"
        + endDate.strftime("%Y%m%d")
        + ".txt"
    )
    totIndsFilePath = os.path.join(outputDir, statFile)
    with open(totIndsFilePath, "w") as f:
        f.write("nrmse = " + str(nrmseTot) + "\n")
        f.write("abs bias = " + str(abTot) + "\n")
        f.write("nbi = " + str(nbiTot) + "\n")
        f.write("hh = " + str(hhTot) + "\n")

    # saving to files
    try:
        os.makedirs(outputDir)
    except:
        pass

    np.savetxt(os.path.join(outputDir, "lons.csv"), maplons)
    np.savetxt(os.path.join(outputDir, "lats.csv"), maplats)
    # np.savetxt(os.path.join(outputDir, "rmse-ssh-altimeter_"+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".csv"), rmse)
    # np.savetxt(os.path.join(outputDir, "bias-ssh-altimeter_"+startDate.strftime("%Y%m%d")+"_"+endDate.strftime("%Y%m%d")+".csv"), bias)
    np.savetxt(
        os.path.join(
            outputDir,
            "HS-nrmse-altimeter_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        nrmse,
    )
    np.savetxt(
        os.path.join(
            outputDir,
            "HS-hh-altimeter_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        hh,
    )
    np.savetxt(
        os.path.join(
            outputDir,
            "HS-nbi-altimeter_"
            + startDate.strftime("%Y%m%d")
            + "_"
            + endDate.strftime("%Y%m%d")
            + ".csv",
        ),
        nbi,
    )
    np.savetxt(
        os.path.join(
            outputDir,
            "HS-bias-altimeter_"
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
