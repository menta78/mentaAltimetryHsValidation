import numpy as np
import h5py
from datetime import datetime, timedelta
import math
import matplotlib.pyplot as plt

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
    for i in range(npoints):
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
        time.append(f[ref][0])
    return lon, lat, res, time, stations

def flatten(l):
    return [item for sublist in l for item in sublist]

pathname = "data/tidalGaugeData/GESLAv1_withResiduals.mat"
on, lat, res, time, stationss = get_serie_gesla(pathname)

timme = np.array(flatten(time))

intervalos = range(int(min(timme)), int(max(timme)) + 2) #calculamos los extremos de los intervalos

print(intervalos)


plt.hist(x=timme, bins=intervalos, color='#F2AB6D', rwidth=0.85)
plt.title('Histogram tidal gauge records')
plt.xlabel('time')
plt.ylabel('Frequency')

plt.show() #dibujamos el histograma