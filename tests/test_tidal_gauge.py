import os
from datetime import datetime
import netCDF4

import pytest

from GeslaDataset.gesla import GeslaDataset
import src.computeR2_tidal_gauge as tg

rootDir = "/mnt/c/Users/ggarc/OneDrive/Documents/Projects/mentaAltimetryHsValidation"
# Directory where the coarsened satellite data are located
# the coarsening is performed by the coarsenSatData function. If you already performed this operation you don't need to repeat it
tidalGaugeDataDir = os.path.join(rootDir, "data/tidalGaugeData")

# Directory where the model nc files are located
modelDir = os.path.join(rootDir, "data/schismwwm")

ncExpression="schout_([0-9]*)\_compressed.nc"


def test_check_times():
    datetimeTidal = [
        datetime(2000, 1, 1),
        datetime(2000, 1, 2),
        datetime(2000, 1, 3),
        datetime(2000, 1, 4),
    ]
    datetimeModel = [datetime(2000, 1, 2), datetime(2000, 1, 3)]
    assert tg.check_times(datetimeTidal[0], datetimeTidal[-1], datetimeModel) == True

    datetimeTidal = [
        datetime(2000, 1, 1),
        datetime(2000, 1, 2),
        datetime(2000, 1, 3),
        datetime(2000, 1, 4),
    ]
    datetimeModel = [datetime(2000, 2, 2), datetime(2000, 2, 3)]
    assert tg.check_times(datetimeTidal[0], datetimeTidal[-1], datetimeModel) == False

    datetimeTidal = [
        datetime(2000, 1, 1),
        datetime(2000, 1, 2),
        datetime(2000, 1, 3),
        datetime(2000, 1, 4),
    ]
    datetimeModel = [datetime(1999, 12, 31), datetime(2000, 1, 1)]
    assert tg.check_times(datetimeTidal[0], datetimeTidal[-1], datetimeModel) == True

    datetimeTidal = [
        datetime(2000, 1, 1),
        datetime(2000, 1, 2),
        datetime(2000, 1, 3),
        datetime(2000, 1, 4),
    ]
    datetimeModel = [datetime(2000, 1, 4), datetime(2000, 1, 5)]
    assert tg.check_times(datetimeTidal[0], datetimeTidal[-1], datetimeModel) == True


def test_filter_ncfiles():

    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-28 03:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    assert tg.filter_ncfiles(modelDir, mdlList, startDate, endDate) == ['schout_88_compressed.nc']

    
    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-30 03:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    actual =  tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)
    expected = ['schout_88_compressed.nc','schout_89_compressed.nc', 'schout_90_compressed.nc']
    assert all([a == b for a, b in zip(actual, expected)])

    startDate = "2000-01-28 00:00:00"
    endDate = "2000-01-30 03:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    assert len(tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)) == 0


def test_combine_ssh():
    tol = 1e-5
    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-28 12:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    modelListFiltered = tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)
    ssh, time = tg.get_ssh_time_model(modelDir, modelListFiltered)
    assert ssh[0,0] == pytest.approx(-0.001303, rel=tol, abs=tol)
    assert ssh[-1,0] == pytest.approx(-0.003422, rel=tol, abs=tol)
    assert ssh.shape[0] == 8
    assert time[0] == 7527600.000000
    assert time[-1] == 7603200.000000

    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-29 12:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    modelListFiltered = tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)
    ssh, time = tg.get_ssh_time_model(modelDir, modelListFiltered)
    assert ssh[0,0] == pytest.approx(-0.001303, rel=tol, abs=tol)
    assert ssh[-1,0] == pytest.approx(-0.005541, rel=tol, abs=tol)
    assert ssh.shape[0] == 16
    assert time[0] == 7527600.000000
    assert time[-1] == 7689600


@pytest.mark.skip(reason="not used")
def test_get_tidal_variables():
    meta_file = os.path.join(tidalGaugeDataDir, "GESLA3_ALL.csv")
    data_path = os.path.join(tidalGaugeDataDir, "GESLA3.0_ALL/")
    g3 = GeslaDataset(meta_file=meta_file, data_path=data_path)

    filenames = [
#        'johnston-052a-usa-uhslc',
        'st_helena-sth-gbr-noc',
#        'malakal-007b-plw-uhslc',
#        'alboran-alb-esp-cmems',
    ]
    tidalGauge = g3.files_to_xarray(filenames)

    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-28 12:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    modelListFiltered = tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)
    fpth = os.path.join(modelDir, modelListFiltered[0])
    ds = netCDF4.Dataset(fpth)
    timeModel = ds["time"]
    #lonTidal, latTidal, timeTidal = tg.get_tidal_variables(tidalGauge, timeModel)
    #assert timeTidal[211103] == 7528838


def test_interpolate_ssh():
    startDate = "2000-03-28 00:00:00"
    endDate = "2000-03-28 12:00:00"
    mdlList = tg.get_all_ncfiles(modelDir, ncExpression)
    modelListFiltered = tg.filter_ncfiles(modelDir, mdlList, startDate, endDate)


    ssh, time = tg.get_ssh_time_model(modelDir, modelListFiltered)
    lon, lat = tg.get_latlon_model(modelDir, modelListFiltered[0]) 

    meta_file = os.path.join(tidalGaugeDataDir, "GESLA3_ALL.csv")
    data_path = os.path.join(tidalGaugeDataDir, "GESLA3.0_ALL/")
    g3 = GeslaDataset(meta_file=meta_file, data_path=data_path)

    filenames = [
#        'johnston-052a-usa-uhslc',
#        'st_helena-sth-gbr-noc',
          'palmeira-235a-cpv-uhslc',
#        'alboran-alb-esp-cmems',
    ]

    tidalGauge = g3.files_to_xarray(filenames)

    sshMean = ssh.mean()
    ssh = ssh - sshMean
    tg.interpolate_ssh_tidal_points(lat, lon, ssh, time, tidalGauge)

def test_r2():
    pass

