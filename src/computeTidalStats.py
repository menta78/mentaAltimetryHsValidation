import os, re, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import itertools


filterHighSsh = False

def getFiles(hsSatAndModelDir, startDate, endDate):
    fl = []

    while startDate <= endDate:
        strtime = startDate.strftime("%Y%m%d")
        pthfile = hsSatAndModelDir + "*"+strtime+".npy"
        print(pthfile)
        fl.append(glob.glob(pthfile))
        startDate += timedelta(days=1)

    return list(itertools.chain(*fl))


def elaborateMeasures(
    startDate,
    endDate,
    outputDir,
    filterLowHs=False,
    filterHsThreshold=0.0,
    pth = 90,
):

    def loadFile(flpth):
        print("")
        print("    loading file " + flpth)
        satdts = np.load(flpth)
        if filterHighSsh:
            sshsat = satdts[:, 0]
            sshmdl = satdts[:, 1]
            cnd = np.logical_and(
                sshsat < filterSshMaximum, sshmdl < filterSshMaximum
            )
            satdts = satdts[cnd, :]
        return satdts

    fls = getFiles(outputDir, startDate, endDate)

    obs = np.array([])
    model = np.array([])

    for f in fls:
        data_ = loadFile(f)
        obs_ = data_[:,0]
        model_ = data_[:,1]

        obs = np.concatenate([obs, obs_])
        model = np.concatenate([model, model_])

    meanModel = np.nanmean(model)
    model = model - meanModel
        
    pmodel = np.nanpercentile(model, pth)
    pobs = np.nanpercentile(obs, pth)

    nsc1 = 0
    ssres = 0
    nsc2 = 0
    sstot = 0
    meanObs = 0
    nobs = 0
    valuesObs = []

    for (obsi, modeli) in zip(obs, model):
        if obsi >= 0:
            valuesObs.append(obsi)
            nobs += 1
        else:
            continue

        meanObs = sum(valuesObs)/nobs

    print(meanObs)

    # meanObs = np.nanmean(obs)
    # print(meanObs)

    for (obsi, modeli) in zip(obs, model):
        if obsi >= pobs and modeli >= pobs:
        #if obsi >= 0 and modeli >= 0 and modeli == modeli:
            #print(obsi, modeli)
            nsc1 += np.abs(obsi-modeli)
            nsc2 += np.abs(obsi-meanObs)

            ssres += (obsi-modeli)**2
            sstot += (obsi-meanObs)**2


    # nsc1 = np.nansum(np.abs(obs-model))
    # nsc2 = np.nansum(np.abs(obs-np.nanmean(obs)))

    # ssres = np.nansum((obs-model)**2)
    # sstot = np.nansum((obs-np.nanmean(obs))**2)

    nse  = 1 - nsc1/nsc2
    nnse = 1/(2-nse)
    r2   = 1 - ssres/sstot
    nr2 = 1/(2-r2)

    print("nse = ",nse, "nnse = ", nnse, "r2 = ", r2, "nr2 = ", nr2)
