import os, re
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt

import coarsenSatData


import mapByLonLat as mll


maskPointsCloseToTheCoast = True
minDistFromCoast = 20000
minObsForValidation = 0
maskShallowPoints = False
minSeaDepth = -5
bathyFile = '/home/lmentaschi/usr/WaveWatchIII/gridgen1.1/reference_data/etopo2.nc'


def elaborateMeasures(startDate, endDate, hsSatAndModelDir, outputDir, latlims=[-63, 63], lonlims=[-180, 180], filterLowHs=False, filterHsThreshold=.0, dx=.2, dy=.2):

  years = []
 
  def iterateByYear():
    fls_ = [f for f in os.listdir(hsSatAndModelDir) if re.match('(.*)\.npy', f)]
    pttrn = '(.*)_hsModelAndSatObs_([0-9]{8}).npy'
    dtstrs = [re.match(pttrn, f).groups(0)[1] for f in fls_]
    dts = [datetime.strptime(dtstr, '%Y%m%d') for dtstr in dtstrs]
    iii = np.argsort(dts)
    dts = np.array(dts)[iii]
    fls_ = np.array(fls_)[iii]
    fls = [fls_[i] for i in range(len(fls_)) if (startDate <= dts[i] and endDate >= dts[i])]

    loopYear = None

    def loadFile(flpth):
      print('')
      print('    loading file ' + flpth)
      satdts = np.load(flpth)
      if filterLowHs:
        hssat = satdts[:,3]
        hsmdl = satdts[:,4]
        cnd = np.logical_and(hssat > filterHsThreshold, hsmdl > filterHsThreshold)
        satdts = satdts[cnd, :]
      return satdts
 
    def _elab(dts, lons, lats, saths, modhs):
      maplons, maplats, mapdata = mll.mapByLonLat(dts, lons, lats, saths, modhs, dx, dy, lonlims=lonlims, latlims=latlims)
      obsSum, sqObsSum, devSum, sqDevSum, mdlByObsSum, dtcount =\
                         mll.computeCumDeviations(maplons, maplats, mapdata)
      obsMax, mdlMax = mll.computeMaxima(maplons, maplats, mapdata)
      return maplons, maplats, obsSum, sqObsSum, devSum, sqDevSum, mdlByObsSum, dtcount, obsMax, mdlMax

    for f in fls:
      fpth = os.path.join(hsSatAndModelDir, f)
      satdts_ = loadFile(fpth)
      dts = satdts_[:,0]
      currYr = datetime.fromtimestamp(dts[0]).year
      if loopYear is None:
        loopYear = currYr
        years.append(loopYear)
        satdts = None
        print('  loading year ' + str(loopYear))
      if loopYear == currYr:
        satdts = satdts_ if satdts is None else np.concatenate([satdts, satdts_], 0)
      else:
       # loopYear != currYr
       # yielding the current year
        dts = satdts[:,0]
        lons = satdts[:,1]
        lats = satdts[:,2]
        saths = satdts[:,3]
        modhs = satdts[:,4]
        yield _elab(dts, lons, lats, saths, modhs)
        satdts = satdts_
        loopYear = currYr
        years.append(loopYear)
        print('  loading year ' + str(loopYear))
   # yielding the last year
    dts = satdts[:,0]
    lons = satdts[:,1]
    lats = satdts[:,2]
    saths = satdts[:,3]
    modhs = satdts[:,4]
    yield _elab(dts, lons, lats, saths, modhs)
    


  obsSum, sqObsSum, devSum, sqDevSum, dtcount, obsMaxSum, mdlMaxSum, obsTotMax, mdlTotMax = None, None, None, None, None, None, None, None, None
  for blob in iterateByYear():
    maplons, maplats, _obsSum, _sqObsSum, _devSum, _sqDevSum, _mdlByObsSum, _dtcount, _obsMax, _mdlMax = blob
    if _obsSum is None:
      continue
    if obsSum is None:
      obsSum, sqObsSum, devSum, sqDevSum, mdlByObsSum, dtcount, obsMaxSum, mdlMaxSum, obsTotMax, mdlTotMax =\
                 _obsSum, _sqObsSum, _devSum, _sqDevSum, _mdlByObsSum,_dtcount, _obsMax, _mdlMax, _obsMax, _mdlMax
    else:
      obsSum = np.nansum([obsSum, _obsSum], 0)
      sqObsSum = np.nansum([sqObsSum, _sqObsSum], 0)
      devSum = np.nansum([devSum, _devSum], 0)
      sqDevSum = np.nansum([sqDevSum, _sqDevSum], 0)
      mdlByObsSum = np.nansum([mdlByObsSum, _mdlByObsSum], 0)
      dtcount = np.nansum([dtcount, _dtcount], 0)
      obsMaxSum = np.nansum([obsMaxSum, _obsMax], 0)
      mdlMaxSum = np.nansum([mdlMaxSum, _mdlMax], 0)
      obsTotMax = np.nanmax([obsTotMax, _obsMax], 0)
      mdlTotMax = np.nanmax([mdlTotMax, _mdlMax], 0)

  if not latlims is None:
    cnd = np.logical_and(maplats >= latlims[0], maplats <= latlims[1])
    maplats = maplats[cnd]
    obsSum = obsSum[cnd, :]
    sqObsSum = sqObsSum[cnd, :]
    devSum = devSum[cnd, :]
    sqDevSum = sqDevSum[cnd, :]
    mdlByObsSum = mdlByObsSum[cnd, :]
    dtcount = dtcount[cnd, :]
    obsMaxSum = obsMaxSum[cnd, :]
    mdlMaxSum = mdlMaxSum[cnd, :]
    obsTotMax = obsTotMax[cnd, :]
    mdlTotMax = mdlTotMax[cnd, :]

  if maskPointsCloseToTheCoast:
    msk = mll.createCoastlinePointsMask(maplons, maplats, resl='h', minDistFromCoast=minDistFromCoast)
    cnd = msk == 0
    obsSum[cnd] = np.nan
    sqObsSum[cnd] = np.nan
    devSum[cnd] = np.nan
    sqDevSum[cnd] = np.nan
    mdlByObsSum[cnd] = np.nan
    dtcount[cnd] = np.nan
    obsMaxSum[cnd] = np.nan
    mdlMaxSum[cnd] = np.nan
    obsTotMax[cnd] = np.nan
    mdlTotMax[cnd] = np.nan
  
  if maskShallowPoints:
    msk = mll.createBathyMask(maplons, maplats, minSeaDepth, bathyFile)
    cnd = msk == 0
    obsSum[cnd] = np.nan
    sqObsSum[cnd] = np.nan
    devSum[cnd] = np.nan
    sqDevSum[cnd] = np.nan
    mdlByObsSum[cnd] = np.nan
    dtcount[cnd] = np.nan
    obsMaxSum[cnd] = np.nan
    mdlMaxSum[cnd] = np.nan
    obsTotMax[cnd] = np.nan
    mdlTotMax[cnd] = np.nan

  cnd = dtcount < minObsForValidation
  obsSum[cnd] = np.nan
  sqObsSum[cnd] = np.nan
  devSum[cnd] = np.nan
  sqDevSum[cnd] = np.nan
  mdlByObsSum[cnd] = np.nan
  dtcount[cnd] = np.nan
  obsMaxSum[cnd] = np.nan
  mdlMaxSum[cnd] = np.nan
  obsTotMax[cnd] = np.nan
  mdlTotMax[cnd] = np.nan
    

  nbiTot = np.nansum(devSum) / np.nansum(obsSum)
  absBiasTot = np.nansum(devSum) / np.nansum(dtcount)
  nrmseTot = np.sqrt( np.nansum(sqDevSum) / np.nansum(sqObsSum) )
  hhTot = np.sqrt( np.nansum(sqDevSum) / np.nansum(mdlByObsSum) )
  nbiYMaxTot = np.nansum(mdlMaxSum - obsMaxSum) / np.nansum(obsMaxSum)
  nbiTotMaxTot = np.nansum(mdlTotMax - obsTotMax) / np.nansum(obsTotMax)
  totIndStr = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOTAL ERROR INDICATORS:
nbi:   {nbiTot:2.5f}
absbi: {absbiTot:2.5f}
nrmse: {nrmseTot:2.5f}
hh: {hhTot:2.5f}
nbiYMax: {nbiYMaxTot:2.5f}
nbiTotMax: {nbiTotMaxTot:2.5f}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
  totIndStr = totIndStr.format(nbiTot=nbiTot, nrmseTot=nrmseTot, hhTot=hhTot, nbiYMaxTot=nbiYMaxTot, nbiTotMaxTot=nbiTotMaxTot, absbiTot=absBiasTot)
  print('')
  print(totIndStr)
  print('')
  totIndsFilePath = os.path.join(outputDir, 'totalIndicators.txt')
  with open(totIndsFilePath, 'w') as f:
    f.write(totIndStr)
    f.close()

  bias = devSum / obsSum
  absBias = devSum / dtcount
  nrmse = np.sqrt(sqDevSum / sqObsSum)
  hh = np.sqrt(sqDevSum / mdlByObsSum)

  obsYearMax = obsMaxSum/len(years)
  mdlYearMax = mdlMaxSum/len(years)

  
  #saving to files
  try:
    os.makedirs(outputDir)
  except:
    pass
  np.savetxt(os.path.join(outputDir, 'lons.csv'), maplons)
  np.savetxt(os.path.join(outputDir, 'lats.csv'), maplats)
  np.savetxt(os.path.join(outputDir, 'bias.csv'), bias)
  np.savetxt(os.path.join(outputDir, 'absbias.csv'), absBias)
  np.savetxt(os.path.join(outputDir, 'nrmse.csv'), nrmse)
  np.savetxt(os.path.join(outputDir, 'hh.csv'), hh)
  np.savetxt(os.path.join(outputDir, 'dtcount.csv'), dtcount)
  np.savetxt(os.path.join(outputDir, 'obsYearMaxMean.csv'), obsYearMax)
  np.savetxt(os.path.join(outputDir, 'modelYearMaxMean.csv'), mdlYearMax)
  np.savetxt(os.path.join(outputDir, 'obsTotMax.csv'), obsTotMax)
  np.savetxt(os.path.join(outputDir, 'modelTotMax.csv'), mdlTotMax)

  print("output dir: " + outputDir)
 

if __name__ == '__main__':
  startYear = 2000
  endYear = 2009

  # RED SEA 10 YEARS
  lonlims = [30, 43.8]
  latlims = [13.4, 32]
 #lonlims = [36.5, 43.8]
 #latlims = [13.4, 21]
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst/stats/'
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_red_sea_unst_noobst/stats/'

 #elaborateMeasures(startYear, endYear, hsSatAndModelDir, outputDir, lonlims=lonlims, latlims=latlims)

  """
  # PERSIC GULF 10 YEARS
  lonlims = [47, 56.5]
  latlims = [23.5, 31]
 #lonlims = [36.5, 43.8]
 #latlims = [13.4, 21]
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst/stats/'
  hsSatAndModelDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst_noobst/hsModelAndSatObs/'
  outputDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION_10yrs/hindcast_persic_gulf_unst_noobst/stats/'
  """

  elaborateMeasures(startYear, endYear, hsSatAndModelDir, outputDir, lonlims=lonlims, latlims=latlims)


