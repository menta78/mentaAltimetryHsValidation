import matplotlib.pyplot as plt
import geopandas
import numpy as np
import h5py


def get_serie_gesla(fileName):
    f = h5py.File(fileName, "r")
    data = f.get("GESELD")

    npoints = data["residual"].shape[0]  # number of stations

    res = []
    time = []
    lon = []
    lat = []
    stations = []

    # print(np.array(f["GESELD"]["station_name"]).item().decode('utf-8'))
    for i in range(npoints):
        ref = f["GESELD"]["longitude"][i][0]
        lon.append(f[ref][0][0])
        ref = f["GESELD"]["latitude"][i][0]
        lat.append(f[ref][0][0])
    return lon, lat


pathname = "data/tidalGaugeData/GESLAv1_withResiduals.mat"
lon, lat = get_serie_gesla(pathname)

shpfile = "/mnt/c/Users/ggarc/OneDrive/Documents/external_data/JRC/VALIDATION/GESLA/ne_10m_coastline/ne_10m_coastline.shp"


fig, ax = plt.subplots()
mask = geopandas.read_file(shpfile)
mask.plot(facecolor=None, edgecolor="grey", ax=ax)

plt.scatter(lon, lat, s=20)
ax.set(xlim=[-180, 180], ylim=[-90, 90])
plt.show()
