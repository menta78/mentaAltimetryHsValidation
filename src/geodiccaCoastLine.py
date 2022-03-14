from src.geodiccaUtility import *

# from alphaBetaLab import abFixBasemap
from mpl_toolkits import basemap
import math


class coastLine:
    def __init__(self, resolution="c", llcrnr=None, urcrnr=None):
        """
        resolution: resolution of the coastline.
          Can be ``c`` (crude), ``l`` (low), ``i`` (intermediate),
          ``h`` (high), ``f`` (full) or None.
        """
        if llcrnr is None:
            self.basemap = basemap.Basemap(resolution=resolution)
        else:
            self.basemap = basemap.Basemap(
                #resolution=resolution,
                resolution=resolution,
                llcrnrlon=llcrnr[0],
                llcrnrlat=llcrnr[1],
                urcrnrlon=urcrnr[0],
                urcrnrlat=urcrnr[1],
            )

        self.coastsegs = self.basemap.coastsegs

    def filterRect(self, minlon, minlat, maxlon, maxlat):
        # building the list of segs containing points in the desired area
        segs = []
        for cs in self.coastsegs:
            for p in cs:
                if (minlon <= p[0] <= maxlon) and (minlat <= p[1] <= maxlat):
                    segs.append(cs)
                    break

        # building segs containing only points in the desired area
        cut_segs = []
        for cs in segs:
            cut_seg = []
            cut_segs.append(cut_seg)
            for p in cs:
                if (minlon <= p[0] <= maxlon) and (minlat <= p[1] <= maxlat):
                    cut_seg.append(p)
        self.coastsegs = cut_segs
        return cut_segs

    def getClosestPointTripletteAndSeg(self, lon, lat):
        mindst = 100000000.0
        closestPtTrpl = None
        closest_seg = None
        for coastseg in self.coastsegs:
            lngth = len(coastseg)
            for i in range(1, lngth):
                pt_1 = coastseg[i - 1]
                pt0 = coastseg[i]
                pt1 = coastseg[(i + 1) % lngth]
                lon0 = pt0[0]
                lat0 = pt0[1]
                dst = ((lon - lon0) ** 2.0 + (lat - lat0) ** 2.0) ** 0.5
                if dst < mindst:
                    mindst = dst
                    closestPtTrpl = (pt_1, pt0, pt1)
                    closest_seg = coastseg
        return closestPtTrpl, closest_seg

    def getDistanceVectorFromCoast(self, lon, lat):
        """
        getDistanceVectorFromCoast: returns the distance of point (lon, lat)
        from the coast. The output is a vector with the x (longitude) and y (latitude)
        components of the distance. The measure unit of the output distance is meters.
        THIS ROUTINE USES AN EUCLIDEAN APPROXIMATION OF THE EARTH SURFACE,
        VALID ONLY FOR SHORT DISTANCES.
        """

        xyintp = [0.0, 0.0]

        def pointInFrontOfCoastSegment(coastpt0, coastpt1):
            result, xyintp[0], xyintp[1] = normalIntersectSegment(
                (lon, lat), coastpt0, coastpt1
            )
            return result

        def computeDistanceDeg():
            return (lon - xyintp[0], lat - xyintp[1])

        closestPtTrpl, _ = self.getClosestPointTripletteAndSeg(lon, lat)
        pt0, pt1, pt2 = closestPtTrpl[0], closestPtTrpl[1], closestPtTrpl[2]
        if pointInFrontOfCoastSegment(pt0, pt1):
            distanceDeg = computeDistanceDeg()
        elif pointInFrontOfCoastSegment(pt1, pt2):
            distanceDeg = computeDistanceDeg()
        else:
            # the nearest part of coast is the cusp
            distanceDeg = (lon - pt1[0], lat - pt1[1])

        xdist = (distanceDeg[0] * cos(lat * deToRa)) * earthCircle / 360.0
        ydist = distanceDeg[1] * earthCircle / 360.0
        return (xdist, ydist)

    def getMathDirNormalToCoast(self, lon, lat):
        dist = self.getDistanceVectorFromCoast(lon, lat)
        dr = math.atan2(-dist[1], -dist[0])
        return dr

    def getDegNDirNormalToCoast(self, lon, lat):
        mathdir = self.getMathDirNormalToCoast(lon, lat)
        return mathDirToDegNPropDir(mathdir)
