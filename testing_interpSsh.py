import os, re

from datetime import datetime

import src.utils as utils
import numpy as np

from matplotlib import pyplot as plt

rootDir = os.path.dirname(os.path.realpath(__file__))

(
    tidalGaugeDataDir,
    crsSatDataDir,
    crsWaveSatDataDir,
    modelNcFilesDir,
    hsModelAndSatObsSshDir,
    hsModelAndSatObsHsDir,
    hsModelAndSatObsTidalDir,
    statsDir,
) = utils.load_paths(rootDir)


# time interval
# startDate, endDate = datetime(1995, 1, 1), datetime(1999, 12, 30)
# startDate, endDate = datetime(2012, 1, 1), datetime(2019, 12, 31)
# startDate, endDate = datetime(2003, 12, 20), datetime(2003, 12, 28)
# startDate, endDate = datetime(2003, 1, 1), datetime(2009, 12, 30)

startDate = datetime(2003, 1, 1)


def loadFile(flpth):
    satdts = np.load(flpth)
    return satdts


dateTime="20031123"
fileNPY=f"ERA5_schismwwm_{dateTime}_hsModelAndSatObs_{dateTime}.npy"


testFile = os.path.join(hsModelAndSatObsSshDir, fileNPY)
data = loadFile(testFile)
print(data.shape[:])

obs_ = data[:,3]
model_ = data[:,4]

filter_arr = np.logical_and(obs_ >= -10, obs_ <= 10)

obs = obs_[filter_arr]
model = model_[filter_arr]

obs = obs-np.nanmean(obs)
model = model-np.nanmean(model)

pobs = np.nanpercentile(obs, 0)
flr = np.logical_and(obs>=pobs, model>=model)

obs = obs[flr]
model =model[flr]


print(obs.shape)

fig, ax = plt.subplots()
ax.plot(obs, label="observation")
ax.plot(model, label="model", alpha=0.5)
ax.legend()
plt.show()
