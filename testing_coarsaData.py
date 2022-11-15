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


testFile = os.path.join(crsSatDataDir, "SSH_GLO_PHY_L3_2003.npy")
data = loadFile(testFile)

obs_ = data[:,3]

filter_arr = np.logical_and(obs_ >= -10, obs_ <= 10)

obs = obs_[filter_arr]

print(obs.shape[::100])

plt.plot(obs[0:2000], '.')
plt.show()