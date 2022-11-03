import numpy as np
import h5py
from datetime import datetime, timedelta
import math

def get_serie_gesla(fileName):
    f = h5py.File(fileName,'r')
    data = f.get("GESELD")

    npoints = data["residual"].shape[0] # number of stations

    res   = []
    time  = []
    lon   = []
    lat   = []
    stations = []

    #print(np.array(f["GESELD"]["station_name"]).item().decode('utf-8'))
    for i in range(4):
        ref = f["GESELD"]["station_name"][i][0]
        station = np.array(f[ref][0][0].item())
        stations.append(station)
        ref = f["GESELD"]["longitude"][i][0]
        lon.append(f[ref][0][0])
        ref = f["GESELD"]["latitude"][i][0]
        lat.append(f[ref][0][0])
        ref = f["GESELD"]["residual"][i][0]
        res.append(np.array(f[ref][0]))
        ref = f["GESELD"]["time"][i][0]
        time.append(np.array(f[ref][0]))
    return lon, lat, res, time, stations


def get_timeSerie(fileName):
    f = h5py.File(fileName,'r')
    data = f.get("GESELD")

    npoints = data["residual"].shape[0] # number of stations

    res   = []
    time  = []
    lon   = []
    lat   = []
    stations = []

    #print(np.array(f["GESELD"]["station_name"]).item().decode('utf-8'))
    for i in range(4):
        ref = f["GESELD"]["station_name"][i][0]
        station = np.array(f[ref][0][0].item())
        stations.append(station)
        ref = f["GESELD"]["longitude"][i][0]
        lon.append(f[ref][0][0])
        ref = f["GESELD"]["latitude"][i][0]
        lat.append(f[ref][0][0])
        ref = f["GESELD"]["residual"][i][0]
        res.append(np.array(f[ref][0]))
        ref = f["GESELD"]["time"][i][0]
        time.append(np.array(f[ref][0]))
    return lon, lat, res, time, stations

#pathname = "data/tidalGaugeData/GESLAv1_withResiduals.mat"
timSerieFl = "data/WL2012_2014.mat"

on, lat, res, time, stationss = get_serie_gesla(pathname)
#print(res[:][:])
#print(len(lon), len(res))
print(stationss)


# def get_subsest_list(array, a, b):
#     if a >= b:
#         print("first value must be smaller than second")
#         return None
#     if a <= max(array) and a>= min(array) and b <= max(array) and b>= min(array):
#         i = np.argmin(np.abs(np.array(array)-a))
#         j = np.argmin(np.abs(np.array(array)-b))
#         return array[i:j]
#     elif a <= max(array) and a>= min(array) and b >= max(array):
#         i = np.argmin(np.abs(np.array(array)-a))
#         return array[i:-1]
#     elif a <= min(array) and b <= max(array) and b>= min(array):
#         j = np.argmin(np.abs(np.array(array)-b))
#         return array[0:j]    


# dt = datetime(1993,1,1,0,0,0)

# ord = dt.toordinal()
# mdn = dt + timedelta(days = 366)
# frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)

# print(mdn.toordinal() + frac)



# kk = np.zeros([3,1])
# kk2 = np.zeros([3,1])

# kk[:,0] = np.array([1,2,3])

# kk2[:,0] = np.array([4,5,6])

# print(np.concatenate([kk,kk2], 1))


def find_nearest(a, value1, value2):
    return np.argmin(np.abs(np.array(a)-value1)), np.argmin(np.abs(np.array(a)-value2))

k = [10.1,11.1,12.1, 13.1, 14.1, 15.1]

print(find_nearest(k, 11, 15.0))

print(find_nearest(k, 10, 15.0))

print(find_nearest(k, 10, 16.0))

print(find_nearest(k, 12, 16.0))



def nash_sutcliffe_coefficient(data):
    obs   = data[::6,0]
    model = data[::6,1]

    nsc1 = np.sum(np.abs(obs-model))
    nsc2 = np.sum(np.abs(obs-np.nanmean(obs)))

    ssres = np.sum((obs-model)**2)
    sstot = np.sum((obs-np.nanmean(obs))**2)

    nse = 1 - nsc1/nsc2
    nnse = 1/(2-nse)
    r2 = 1 - ssres/sstot

    return nse, nnse, r2

def read_file(line, fileName = "data/tidalGaugeData/stations-tidalGauge.txt"):
    # open the sample file used
    file = open(fileName)
    
    # read the content of the file opened
    content = file.readlines()

    return content[line].rstrip('\n')

print(read_file(1))
