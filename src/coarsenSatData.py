import os
import re
from datetime import datetime
import time

import numpy as np
import netCDF4
from matplotlib import path

import src.geodiccaCircUtils as circularUtils


def coarsenSatData(rootdir, outputdir, startdate, enddate, latdelta, areaRectangles=[]):
    # areaRectangles is a list of rectangles [minlon minlat, maxlon, maxlat]

    drs = list([d for d in os.listdir(rootdir) if re.match("[0-9]{4}", d)])
    drs.sort()

    os.chdir(rootdir)

    satelliteMap = {
        1: "CMEMS_SWH",
    }
    satIds = list(satelliteMap.keys())
    satIds.sort()

    tic = time.time()

    try:
        os.makedirs(outputdir)
    except:
        pass

    for dr in drs:
        print("elaborating directory " + dr)
        drpth = os.path.join(rootdir, dr)
        monthdirs = [d for d in os.listdir(drpth) if re.match("[0-9]{2}", d)]
        monthdirs.sort()

        year = dr
        dataBySatellite = {}

        for mdr in monthdirs:
            actDate = datetime(int(dr), int(mdr), 1)
            if (startdate > actDate) or (enddate <= actDate):
                print(str(actDate) + " out of time range. Skipping")
                continue

            mdrpath = os.path.join(drpth, mdr)
            files = os.listdir(mdrpath)
            files.sort()
            for f in files:
                print("elaborating file " + f)
                fpath = os.path.join(mdrpath, f)
                ds = netCDF4.Dataset(fpath)
                #sats = np.array(ds.variables["satellite"])
                timevar = ds.variables["time"]
                times = np.array(timevar)
                # lons = np.array(ds.variables["lon"])
                # lats = np.array(ds.variables["lat"])
                # hss = np.array(ds.variables["swh"])
                lons = np.array(ds.variables["longitude"])
                lats = np.array(ds.variables["latitude"])
                hss = np.array(ds.variables["VAVH"])
                for areaRectangle in areaRectangles:
                    minlon = areaRectangle[0]
                    minlat = areaRectangle[1]
                    maxlon = areaRectangle[2]
                    maxlat = areaRectangle[3]
                    cndlon = np.logical_and(minlon <= lons, lons < maxlon)
                    cndlat = np.logical_and(minlat <= lats, lats < maxlat)
                    cnd = np.logical_and(cndlon, cndlat)
                    lons = lons[cnd]
                    lats = lats[cnd]
                    times = times[cnd]
                    hss = hss[cnd]
                    #sats = sats[cnd]

                for satid in satIds:
                    # stms = times[sats == satid]
                    # slons = lons[sats == satid]
                    # slats = lats[sats == satid]
                    # shss = hss[sats == satid]

                    stms = times
                    slons = lons
                    slats = lats
                    shss = hss

                    # average data are by time, lon, lat, hs
                    avgDt = dataBySatellite.get(satid, [[], [], [], []])
                    dataBySatellite[satid] = avgDt

                    icoord = 0
                    ilat = 0
                    indexes = []
                    while icoord < len(slats):
                        # assiming a polar orbit and slicing by latitude
                        lat0 = int(np.floor(slats[ilat] / latdelta))
                        icoord = ilat
                        for lat in slats[ilat:]:
                            if int(np.floor(lat / latdelta)) == lat0:
                                indexes.append(icoord)
                                icoord += 1
                            else:
                                ltms = stms[indexes]
                                llons = slons[indexes]
                                llats = slats[indexes]
                                lhss = shss[indexes]

                                avgtm = np.mean(ltms)
                                avgtm = netCDF4.num2date(
                                    avgtm,
                                    timevar.units,
                                    only_use_cftime_datetimes=False,
                                )
                                avglon = circularUtils.angleMeanDeg(llons)
                                avglat = circularUtils.angleMeanDeg(llats)
                                avghs = np.mean(lhss)

                                avgDt[0].append(datetime.timestamp(avgtm))
                                avgDt[1].append(avglon)
                                avgDt[2].append(avglat)
                                avgDt[3].append(avghs)

                                ilat = icoord
                                indexes = []
                                break

                ds.close()

        print("  ... year completed. Saving output ...")
        for satid in dataBySatellite.keys():
            data = dataBySatellite[satid]
            if not len(data[0]):
                continue
            satname = satelliteMap[satid]
            outputfile = os.path.join(outputdir, satname + "_" + year)
            print("saving file for satellite " + satname)
            np.save(outputfile, np.array(data).transpose())

    toc = time.time()
    print("time elapsed: {t} s".format(t=toc - tic))


class satDataManager:
    def __init__(self, dataDir, boundaries=None):
        self.dataDir = dataDir
        if not boundaries is None:
            self.boundaries = path.Path(boundaries)
        else:
            self.boundaries = None

    def getSatNames(self):
        satNames = [
            re.sub("(.*)_(.*)(\.npy)", "\\1", f)
            for f in os.listdir(self.dataDir)
            if re.match("(.*)_(.*)\.npy", f)
        ]
        return set(satNames)

    def getObsYears(self):
        yrs = [
            int(re.sub("(.*)_(.*)(\.npy)", "\\2", f))
            for f in os.listdir(self.dataDir)
            if re.match("(.*)_(.*)\.npy", f)
        ]
        yrs = np.unique(yrs)
        return yrs

    def getSatData(self, satname, year):
        yearstr = str(year)
        flname = satname + "_" + yearstr + ".npy"
        flpath = os.path.join(self.dataDir, flname)
        if not os.path.isfile(flpath):
            return None
        print("  ... loading file " + flname)
        dtnp = np.load(flpath)
        if not self.boundaries is None:
            lon = dtnp[:, 1]
            lat = dtnp[:, 2]
            pts = np.array([lon, lat]).transpose()
            cnd = self.boundaries.contains_points(pts)
            dtnp = dtnp[cnd, :]
        return dtnp

    def getAllSatDataOfYr(self, year):
        result = {}
        allNone = True
        satNames = self.getSatNames()
        for satname in satNames:
            dt = self.getSatData(satname, year)
            result[satname] = dt
            if not dt is None:
                allNone = False
        if allNone:
            return None
        else:
            return result

    def getAllSatDataBtwDt(self, startDate, endDate):
        result = {}
        allNone = True
        satNames = self.getSatNames()
        startYear = startDate.year
        endYear = endDate.year
        startTms = datetime.timestamp(startDate)
        endTms = datetime.timestamp(endDate)
        yrs = self.getObsYears()
        for satname in satNames:
            for year in range(startYear, endYear + 1):
                dty = self.getSatData(satname, year)
                if dty is None:
                    continue
                if dty.ndim == 1:
                    dty = np.expand_dims(dty, 0)
                tmstmp = dty[:, 0]
                cnd = np.logical_and(startTms <= tmstmp, tmstmp < endTms)
                dty = dty[cnd, :]
                dt0 = result.get(satname, None)
                dt = dty if dt0 is None else np.concatenate([dt0, dty], 0)
                result[satname] = dt
        return result

    def stackData(self, data):
        rslt = None
        dts = data.values()
        for dt in dts:
            if rslt is None:
                rslt = dt
            else:
                rslt = np.concatenate([rslt, dt], 0)
        return rslt

    def getAllSatData(self):
        result = {}
        allNone = True
        satNames = self.getSatNames()
        yrs = self.getObsYears()
        for satname in satNames:
            for year in yrs:
                dty = self.getSatData(satname, year)
                if dty is None:
                    continue
                if dty.ndim == 1:
                    dty = np.expand_dims(dty, 0)
                dt0 = result.get(satname, None)
                dt = dty if dt0 is None else np.concatenate([dt0, dty], 0)
                result[satname] = dt
        return result

    def getMergedArray(self, dataBySat):
        result = None
        for dt in dataBySat.values():
            result = dt if result is None else np.concatenate([result, dt], 0)
            if result.ndim == 1:
                result = np.expand_dims(result, 0)
        return result


def plotCoarsenedDataStats(dataDir, dlon, dlat, vmin, vmax):
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    dataMngr = satDataManager(dataDir)
    data = dataMngr.getAllSatData()
    data = dataMngr.getMergedArray(data)
    tmstmp = data[:, 0]
    lon = data[:, 1]
    lat = data[:, 2]
    hs = data[:, 3]

    minlon, maxlon = min(lon), max(lon)
    lonmap = np.arange(minlon, maxlon, dlon)
    lonii = np.digitize(lon, lonmap)
    lonbin = lonmap[lonii - 1]

    minlat, maxlat = min(lat), max(lat)
    latmap = np.arange(minlat, maxlat, dlat)
    latii = np.digitize(lat, latmap)
    latbin = latmap[latii - 1]

    cnt = np.zeros([latmap.shape[0], lonmap.shape[0]])
    for ilon, ilat in zip(lonii - 1, latii - 1):
        cnt[ilat, ilon] += 1
    cnt[cnt == 0] = np.nan

    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent([minlon - 1, maxlon + 1, minlat - 1, maxlat + 1], ccrs.PlateCarree())
    plt.pcolormesh(
        lonmap, latmap, cnt, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax
    )
    ax.coastlines(resolution="10m")
    ax.gridlines(draw_labels=True)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    import pdb

    pdb.set_trace()

    rootdir = "/media/lmentaschi/TOSHIBA EXT/Globwave/"
    outputdir = "/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION/redSeaAndPersicGulf/coarsenedSatData"
    # startdate = datetime(2000, 1, 1)
    # enddate = datetime(2010, 1, 1)
    startdate = datetime(1979, 1, 1)
    enddate = datetime(2021, 1, 1)
    latdelta = 0.05
    areaRectangles = [(30, 12, 57, 31)]
    coarsenSatData(rootdir, outputdir, startdate, enddate, latdelta, areaRectangles)

# crsDataDir = '/media/lmentaschi/TOSHIBA EXT/analysis/WW3_UNST_VALIDATION/redSeaAndPersicGulf/coarsenedSatData'
# dlon, dlat = .05, .05
# vmin, vmax = 0, 200
# plotCoarsenedDataStats(crsDataDir, dlon, dlat, vmin, vmax)
