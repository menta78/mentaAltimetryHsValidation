import os, re
from datetime import datetime
from GeslaDataset.gesla import GeslaDataset

from src.interpolateModelToTidelGauge import interpolateModelToTidalGauge_schismWWM

# Directory where the raw globwave files are located
# rawSatDataDir = "/home/ggarcia/Projects/mentaAltimetryHsValidation/satData/rawData"

rootDir = "/mnt/c/Users/ggarc/OneDrive/Documents/Projects/mentaAltimetryHsValidation"

# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelNcFilesDir = os.path.join(rootDir, "data/schismwwm")

# Directory where the pairs observation/model are to be generated
hsModelAndSatObsDir = os.path.join(rootDir, "data/satModelPairs/")

# Directory where the stats are generated
statsDir = os.path.join(rootDir, "data/stats/")

# time interval
startDate, endDate = datetime(2000, 2, 1), datetime(2000, 3, 30)
overwriteExisting = False

# number of processes to be used for the interpolation
nParWorker = 4

# threshold above which hs should be considered
filterSshMaximum = 100
filterHighSsh = True

# set this if you need to limit your analysis to a subdomain
boundaries = None


filterTidalGaugeData = False

if filterTidalGaugeData:
    south_lat = 20
    north_lat = 60
    west_lon = 20
    east_lon = 60

    meta_file = os.path.join(tidalGaugeDataDir, "GESLA3_ALL.csv")
    data_path = os.path.join(tidalGaugeDataDir, "GESLA3.0_ALL/")

    g3 = GeslaDataset(meta_file=meta_file, data_path=data_path)

    data = g3.load_lat_lon_range(
        south_lat=south_lat,
        north_lat=north_lat,
        west_lon=west_lon,
        east_lon=east_lon,
    )
    print(data.site_name.values)
else:
    meta_file = os.path.join(tidalGaugeDataDir, "GESLA3_ALL.csv")
    data_path = os.path.join(tidalGaugeDataDir, "GESLA3.0_ALL/")

    g3 = GeslaDataset(meta_file=meta_file, data_path=data_path)

    filenames = [
        'johnston-052a-usa-uhslc',
        'malakal-007b-plw-uhslc',
        'alboran-alb-esp-cmems',
    ]

    xr_dataset = g3.files_to_xarray(filenames)

print(xr_dataset["sea_level"].shape[:])
kjoeifr

doInterpolateModelToSat = True
if doInterpolateModelToSat:
    # interpolating the model ssh along tidal gauges
    interpolateModelToTidalGauge_schismWWM(
        xr_dataset,
        modelNcFilesDir,
        boundaries,
        startDate,
        endDate,
        overwriteExisting=overwriteExisting,
        nParWorker=nParWorker,
    )


# computing the statistics
""" dx, dy = 1.0, 1.0
computeSshStats.maskPointsCloseToTheCoast = False

latlims = [-90, 90]
computeSshStats.elaborateMeasures(
    startDate,
    endDate,
    hsModelAndSatObsDir,
    statsDir,
    dx=dx,
    dy=dy,
    filterHighSsh=True,
    filterSshMaximum=filterSshMaximum,
    latlims=latlims,
    pth = 95,
) """
