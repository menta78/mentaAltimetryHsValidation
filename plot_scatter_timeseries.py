from matplotlib import pyplot as plt
import numpy as np


flpth = "data/buoyModelPairs/Buoys_CMEMS_station_97.npy"

satdts = np.load(flpth)
sshsat = satdts[:, 0]
sshmdl = satdts[:, 1]
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
