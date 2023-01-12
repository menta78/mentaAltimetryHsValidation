from matplotlib import pyplot as plt
import numpy as np


flpth = "data/tidalModelPairsNew/1995/TidalGauge_GESLA_station_78_19950101_19960101.npy"
flpth = "data/buoyModelPairs/1998/Buoys_CMEMS_station_625_19980101_19990101.npy"

satdts = np.load(flpth)
sshsat = satdts[:, 0]
sshmdl = satdts[:, 1]
if "TidalGauge" in flpth:
    sshsat = sshsat - np.nanmean(satdts[:, 0])
    sshmdl = sshmdl - np.nanmean(satdts[:, 1])
    cnd = np.logical_and(sshsat > -100, sshmdl > -100)
    satdts = satdts[cnd, :]
    cnd = np.logical_and(satdts[:, 0] < 100, satdts[:, 1] < 100)
    satdts = satdts[cnd, :]
else:
    filterHsThreshold = 0
    cnd = np.logical_and(sshsat > filterHsThreshold, sshmdl > filterHsThreshold)
    satdts = satdts[cnd, :]

obsTarget = satdts[:,0]
modelTarget = satdts[:,1]

fig, ax = plt.subplots()
ax.plot(obsTarget, label="observation")
ax.plot(modelTarget, label="model", alpha=0.5)
ax.legend()
plt.show()
