import os, re
from datetime import datetime
import numpy as np
import numpy.ma as ma
import gzip
import glob
import netCDF4
from datetime import datetime, timedelta

from matplotlib import pyplot as plt
import geopandas

from src.interpolateModelToBouy import interpolateModelToTidalGauge_schismWWM
from src.computeBuoysStats import elaborateMeasures
from src.plotStatsBuoys import elaborateMeasuresPlot
import itertools

import src.utils as utils


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

def getFiles(hsSatAndModelDir):
    fl = []

    pthfile = os.path.join(hsSatAndModelDir, "*.nc")
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


# time interval
startDate, endDate = datetime(2002, 3, 22), datetime(2009, 12, 30)
startDate, endDate = datetime(2006, 12, 20), datetime(2007, 12, 29)
startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 29)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 1

# Percentile
pth = 95

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None


doInterpolateModelToSat = False
if doInterpolateModelToSat:
    lstBuoys = getFiles(buoysDir)
    lonBuoys, latBuoys, timeBuoys, hsBuoys = get_hs_buoy(lstBuoys)
    print(len(lonBuoys), len(latBuoys), len(timeBuoys), len(timeBuoys))
    assert len(lonBuoys) == len(latBuoys) == len(timeBuoys) == len(timeBuoys)
    
    varsBuoys = [lonBuoys, latBuoys, hsBuoys, timeBuoys]
    # interpolating the model ssh along tidal gauges
    interpolateModelToTidalGauge_schismWWM(
        varsBuoys,
        modelNcFilesDir,
        hsModelAndSatObsBuoysDir,
        None,
        boundaries,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )


r2Compute = False
if r2Compute:
    elaborateMeasures(
        startDate,
        endDate,
        hsModelAndSatObsBuoysDir,
        statsDir,
        pth = pth,
    ) 

r2ComputePlot = True

if r2ComputePlot:
    elaborateMeasuresPlot(
        startDate,
        endDate,
        hsModelAndSatObsBuoysDir,
        statsDir,
        None,
        None,
        pth = pth,
        nminobs = 10,
) 
