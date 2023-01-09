import csv
import math
import os
from datetime import datetime, timedelta

import h5py
import numpy as np
import scipy.io
import glob
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import itertools
import src.utils as utils
import netCDF4

from src.computeTidalStats_timeSeries import elaborateMeasures

lon = []
lat = []
time = []
hsBuoy = []

hsVars = ["VHM0"]

def datetime2matlabdn(dt):
    ord = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt - datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / (
        24.0 * 60.0 * 60.0
    )
    return mdn.toordinal() + frac


def toJulian(tarray_):
    i = 0
    tarray = tarray_
    # Convert to julian days
    for tt in tarray_:
        tarray[i] = datetime2matlabdn(tt)
        i += 1
    return tarray

def getFiles(hsModelAndSatObsBuoysDir):
    fl = []

    pthfile = os.path.join(hsModelAndSatObsBuoysDir, "*.nc")
    fileFound = glob.glob(pthfile)
    fl.append(fileFound)

    return list(itertools.chain(*fl))


def get_hs_buoy(lstBuoys):
    i = 0
    for buoy in lstBuoys:
        print(buoy)
        
        try:
            ds = netCDF4.Dataset(buoy)
            for hsvar in hsVars:
                if hsvar in ds.variables:
                    i += 1
                    varHs = hsvar
            if i == 0:
                print("ignoring ...")
                continue
            else:
                i = 0

            depth_qflag = ds["DEPH_QC"][:][0]
            position_qflag = ds["POSITION_QC"][:][0]
            time_qflag = ds["TIME_QC"][:][0]
            hs_qflag = ds[varHs+"_QC"][:][0][0]
            
            if (position_qflag == 1) and (time_qflag == 1) and (hs_qflag == 1):
                tmnc_ = ds["TIME"]
                hs_ = ds[varHs][:]
                hs_ = hs_.flatten()
                hs__ = hs_[~hs_.mask]
                tmnc__ = tmnc_[:]
                tmnc = tmnc__[~hs_.mask]
                tmnc = tmnc.flatten()
                tmmdl_ = toJulian(netCDF4.num2date(tmnc[:], tmnc_.units, tmnc_.calendar, only_use_cftime_datetimes=False))
                tmmdl = tmmdl_.data # convert masked array to "usual" array

                hs__ = ma.getdata(hs_)
                hs = hs__.flatten()

                hsBuoy.append(hs.tolist())
                time.append(tmmdl.tolist())
                lon.append(ds["LONGITUDE"][0])
                lat.append(ds["LATITUDE"][0])
                # TODO
                # en utils crear una función para que los datos tengan la misma estructura que en tidalGauges
                # no puede haber fill values. Si algo es nan, debe eliminarse y también en el array de tiempos
                # no puede quedarse un res vacio. Tiene que tener siempre valor
            else:
                print("ignoring...")
        except:
            pass

    return lon, lat, time, hsBuoy


rootDir = os.path.dirname(os.path.realpath(__file__))

(
    tidalGaugeDataDir,
    buoysDir,
    crsSatDataDir,
    crsWaveSatDataDir,
    modelNcFilesDir,
    hsModelAndSatObsSshDir,
    hsModelAndSatObsHsDir,
    hsModelAndSatObsTidalDir,
    hsModelAndSatObsBuoysDir,
    statsDir,
) = utils.load_paths(rootDir)


startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
overwriteExisting = False

interpolate = True

if interpolate:

    target = [-131.339644, 49.335392]
    lstBuoys = getFiles(buoysDir)
    lonBuoys, latBuoys, timeBuoys, hsBuoys = get_hs_buoy(lstBuoys)

    print(lonBuoys)
    fjrijfri

    extension = ".nc"
    fltPattern = "ERA5_schismwwm_"
    flsPath = utils.getFiles(modelNcFilesDir, startDate, endDate, fltPattern, extension)

    timeModel, lonModel, latModel, elev = utils.getModelVariables(flsPath)

    for i in range(10):
        node = utils.find_closest_node(lonBuoys, latBuoys, target)
        tmstmpTidal_ = timeBuoys[node][:]
        # Get timestamps recorded in the stations that belongs to the array of timestamps of the schism model
        tmstmpTidal_, isSubset, idxMin, idxMax = utils.get_subsest_list(
            tmstmpTidal_, min(timeModel), max(timeModel)
        )
        if isSubset and (len(tmstmpTidal_) > 10):
            print("target = ", target)
            print(idxMin, idxMax)
            tmstmpTidal = np.asarray(tmstmpTidal_)
            break

    nodeModel = utils.find_closest_node(lonModel, latModel, target)
    timeSerieModel = elev[:, nodeModel]

    intpltr = interp1d(timeModel, timeSerieModel, bounds_error=False)

    # Model elevation interpolated at tidal time
    intp = intpltr(tmstmpTidal)

    res = np.array(hsBuoys[node][idxMin:idxMax])
    timeSerieTidal = np.asarray(hsBuoys[node][idxMin:idxMax])

    print(tmstmpTidal.shape[:], timeSerieTidal.shape[:], intp.shape[:])



plot = False
if plot:
    options_savefig = {"dpi": 150, "bbox_inches": "tight", "transparent": False}
    title_font = {
        "size": "12",
        "color": "black",
        "weight": "normal",
        "verticalalignment": "bottom",
    }

    fig, ax = plt.subplots()
    # It's missing time array
    ax.plot(timeModel, elevTimeSerie, label="Lorenzo's Model")
    ax.plot(time, wl[npoint, :], label="Tomas' Model", alpha=0.5)
    ax.legend()
    plt.show()
    plt.savefig(
        "/home/vousdmi/Desktop/tomasVSlorenzo_lon="
        + str(target[0])
        + "_lat="
        + str(target[1])
        + ".png",
        **options_savefig
    )
    plt.close(fig)
