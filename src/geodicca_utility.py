from math import *
import numpy as np

earthCircle = 40000000.0
earthRadius = earthCircle / (2 * pi)
deToRa = pi / 180.0
raToDe = 180.0 / pi

seaMileToKm = 1.852
kmToSeaMile = 1.0 / seaMileToKm


def mathDirToDegNPropDir(mathdir):
    return ((-mathdir + pi / 2.0) * raToDe) % 360.0


def mathDirToDegNOrigDir(mathdir):
    return -mathDirToDegNPropDir(mathdir) % 360.0


def greatCircleAngleRad(lon1, lat1, lon2, lat2):
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    delta = abs(
        2
        * asin(
            (sin(dlat / 2) ** 2.0 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2.0) ** 0.5
        )
    )
    return delta


def greatCircleAngleDeg(lon1, lat1, lon2, lat2):
    radlon1 = lon1 * deToRa
    radlat1 = lat1 * deToRa
    radlon2 = lon2 * deToRa
    radlat2 = lat2 * deToRa
    anglerad = greatCircleAngleRad(radlon1, radlat1, radlon2, radlat2)
    return anglerad * raToDe


def greatCircleDistanceRad(lon1, lat1, lon2, lat2):
    angle = greatCircleAngleRad(lon1, lat1, lon2, lat2)
    return angle * earthRadius


def greatCircleDistanceDeg(lon1, lat1, lon2, lat2):
    angle = greatCircleAngleDeg(lon1, lat1, lon2, lat2)
    return angle * deToRa * earthRadius


def pointInsidePolygon(x, y, poly):
    """
    Determine if a point is inside a given polygon or not
    Polygon is a list of (x,y) pairs.

    thanks to Patrick Jordan
    http://www.ariel.com.au/a/python-point-int-poly.html
    """

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
                        p1x, p1y = p2x, p2y

    return inside


def getNormalVersor(lon0, lat0, lon1, lat1):
    dy = lat1 - lat0
    avglat = (lat0 + lat1) / 2.0 * deToRa
    assert (
        avglat != 0
    ), "getNormalVersor: this function is affordable only far from the poles"
    dx = (lon1 - lon0) * cos(avglat)
    vx = -dy / cos(avglat)
    vy = dx
    l = (vx**2.0 + vy**2.0) ** 0.5
    versor = (vx / l, vy / l)
    return versor


def normalIntersectSegment(normalPt, segpt0, segpt1):
    """
    normalIntersectSegment: returns True if the recta passing through normalPt,
    normal to the segment with endpoints segpt0 and segpt1, intersects the segment
    """
    avgLat = (normalPt[1] + segpt0[1] + segpt1[1]) / 3.0 * deToRa
    assert abs(avgLat) < pi / 2.0 * (
        89.0 / 90.0
    ), "normalIntersectSegment: funcion invalid at the poles"

    xpt0 = (segpt0[0] - normalPt[0]) * cos(avgLat)
    ypt0 = segpt0[1] - normalPt[1]

    xpt1 = (segpt1[0] - normalPt[0]) * cos(avgLat)
    ypt1 = segpt1[1] - normalPt[1]

    if (xpt0 != xpt1) and (ypt0 != ypt1):
        m = (ypt1 - ypt0) / (xpt1 - xpt0)
        l = ypt0 - m * xpt0

        mn = -1 / m

        xpintp = l / (mn - m)
        ypintp = mn * xpintp
    elif xpt0 == xpt1:
        xpintp = xpt0
        ypintp = 0.0
    elif ypt0 == ypt1:
        ypintp = ypt0
        xpintp = 0.0

    minx = min(xpt0, xpt1)
    maxx = max(xpt0, xpt1)
    miny = min(ypt0, ypt1)
    maxy = max(ypt0, ypt1)
    result = (minx <= xpintp <= maxx) and (miny <= ypintp <= maxy)

    ypintp = ypintp + normalPt[1]
    xpintp = xpintp / cos(avgLat) + normalPt[0]

    return result, xpintp, ypintp


def rectanglesIntersect(minx0, maxx0, miny0, maxy0, minx1, maxx1, miny1, maxy1):
    segspoly = [
        (minx1, miny1),
        (minx1, maxy1),
        (maxx1, maxy1),
        (maxx1, miny1),
        (minx1, miny1),
    ]

    oneAngleContained = (
        pointInsidePolygon(minx0, miny0, segspoly)
        or pointInsidePolygon(minx0, maxy0, segspoly)
        or pointInsidePolygon(maxx0, maxy0, segspoly)
        or pointInsidePolygon(maxx0, miny0, segspoly)
    )
    domainContainsSegs = (
        (minx0 <= minx1) and (maxx0 >= maxx1) and (miny0 <= miny1) and (maxy0 >= maxy1)
    )
    segsContainsDomain = (
        (minx0 >= minx1) and (maxx0 <= maxx1) and (miny0 >= miny1) and (maxy0 <= maxy1)
    )
    crossDomains1 = (
        (minx0 <= minx1) and (maxx0 >= maxx1) and (miny0 >= miny1) and (maxy0 <= maxy1)
    )
    crossDomains2 = (
        (minx0 >= minx1) and (maxx0 <= maxx1) and (miny0 <= miny1) and (maxy0 >= maxy1)
    )
    return (
        oneAngleContained
        or domainContainsSegs
        or segsContainsDomain
        or crossDomains1
        or crossDomains2
    )


def polygonSurface(poly):
    """
    Computation of the surface enclosed in a polygon applying Green's theorem
    """

    def segments(p):
        return zip(p, p[1:] + [p[0]])

    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for ((x0, y0), (x1, y1)) in segments(poly)))


def getIntersection(seg00, seg01, seg10, seg11):
    x00 = seg00[0]
    y00 = seg00[1]
    x01 = seg01[0]
    y01 = seg01[1]
    x10 = seg10[0]
    y10 = seg10[1]
    x11 = seg11[0]
    y11 = seg11[1]
    m0 = (y01 - y00) / (x01 - x00) if x01 != x00 else None
    m1 = (y11 - y10) / (x11 - x10) if x11 != x10 else None
    if x01 == x00:
        x = x00
        if m1 != None:
            y = y10 + m1 * (x - x10)
        else:
            return
    elif x11 == x10:
        x = x10
        if m0 != None:
            y = y00 + m0 * (x - x00)
        else:
            return
    elif m0 != m1:
        x = (y10 - y00 - (m1 * x10 - m0 * x00)) / (m0 - m1)
        y = y00 + m0 * (x - x00)
    else:
        return
    if (
        (min(x00, x01) <= x <= max(x00, x01))
        and (min(x10, x11) <= x <= max(x10, x11))
        and (min(y00, y01) <= y <= max(y00, y01))
        and (min(y10, y11) <= y <= max(y10, y11))
    ):
        return x, y


def almostEqual(a, b, sigdigits=7):
    return (a == b) or (
        int(round(a * 10**sigdigits)) == int(round(b * 10**sigdigits))
    )


def getPtInsidePolygon(poly):
    """
    Tries to find a point contained in polygon poly.
    """

    def addPtsOnSectionToPoly(minx, maxx, miny, maxy):
        result = [p for p in poly]
        if result[0] != result[-1]:
            result.append(result[0])
        segs = [(sg[0], sg[1]) for sg in zip(result[:-1], result[1:])]
        bsg1 = ((minx, miny), (minx, maxy))
        bsg2 = ((minx, maxy), (maxx, maxy))
        bsg3 = ((maxx, miny), (maxx, maxy))
        bsg4 = ((minx, miny), (maxx, miny))
        for sg in segs:
            sgp1 = sg[0]
            sgp2 = sg[1]
            if sgp2[0] >= sgp1[0]:
                if sgp2[1] >= sgp1[1]:
                    blst = [bsg1, bsg4, bsg2, bsg3]
                else:
                    blst = [bsg1, bsg2, bsg3, bsg4]
            else:
                if sgp2[1] >= sgp1[1]:
                    blst = [bsg3, bsg4, bsg1, bsg2]
                else:
                    blst = [bsg2, bsg3, bsg1, bsg4]
            for bsg in blst:
                intersectPt = getIntersection(bsg[0], bsg[1], sg[0], sg[1])
                if (
                    intersectPt
                    and (not intersectPt in result)
                    and not (
                        (
                            almostEqual(intersectPt[0], sgp1[0])
                            and almostEqual(intersectPt[1], sgp1[1])
                        )
                        or (
                            almostEqual(intersectPt[0], sgp2[0])
                            and almostEqual(intersectPt[1], sgp2[1])
                        )
                    )
                ):
                    ptIndx = result.index(sg[0]) + 1
                    result.insert(ptIndx, intersectPt)
        return result

    def getSubPolygonsInSect(minx, maxx, miny, maxy):
        poly_ = addPtsOnSectionToPoly(minx, maxx, miny, maxy)
        subplgs = []
        currentsubpl = None
        for p in poly_:
            insect = (minx <= p[0] <= maxx) and (miny <= p[1] <= maxy)
            if insect:
                if not currentsubpl:
                    currentsubpl = []
                    subplgs.append(currentsubpl)
                currentsubpl.append(p)
            else:
                currentsubpl = None
        return subplgs

    def getLargestSubPolygonInSect(minx, maxx, miny, maxy):
        result = None
        maxsfc = 0.0
        subplgs = getSubPolygonsInSect(minx, maxx, miny, maxy)
        for sp in subplgs:
            actsfc = polygonSurface(sp)
            if actsfc > maxsfc:
                maxsfc = actsfc
                result = sp
        return result

    def getPtInsideSect(minx, maxx, miny, maxy):
        poly = addPtsOnSectionToPoly(minx, maxx, miny, maxy)
        xs = [p[0] for p in poly if ((minx <= p[0] <= maxx) and (miny <= p[1] <= maxy))]
        ys = [p[1] for p in poly if ((minx <= p[0] <= maxx) and (miny <= p[1] <= maxy))]
        lngth = float(len(xs))
        x = sum(xs) / lngth
        y = sum(ys) / lngth
        return x, y

    def getPolySurfaceInSect(minx, maxx, miny, maxy):
        poly = addPtsOnSectionToPoly(minx, maxx, miny, maxy)
        p_ = [
            (p[0], p[1])
            for p in poly
            if ((minx <= p[0] <= maxx) and (miny <= p[1] <= maxy))
        ]
        if not p_:
            if (
                pointInsidePolygon(minx, miny, poly)
                or pointInsidePolygon(minx, maxy, poly)
                or pointInsidePolygon(maxx, maxy, poly)
                or pointInsidePolygon(maxx, miny, poly)
            ):
                p_ = [(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)]
        if p_:
            return polygonSurface(p_)
        else:
            return 0.0

    if not poly:
        return None, None
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    minx = min(xs)
    maxx = max(xs)
    miny = min(ys)
    maxy = max(ys)
    cnt = 0
    maxcnt = len(poly) * 100
    dividex = True
    while True:
        x, y = getPtInsideSect(minx, maxx, miny, maxy)
        if pointInsidePolygon(x, y, poly):
            return x, y

        if dividex:
            halfx = float(maxx + minx) / 2.0
            polySfcInSect1 = getPolySurfaceInSect(minx, halfx, miny, maxy)
            polySfcInSect2 = getPolySurfaceInSect(halfx, maxx, miny, maxy)
            if polySfcInSect1 >= polySfcInSect2:
                maxx = halfx
            else:
                minx = halfx
        else:
            halfy = float(maxy + miny) / 2.0
            polySfcInSect1 = getPolySurfaceInSect(minx, maxx, miny, halfy)
            polySfcInSect2 = getPolySurfaceInSect(minx, maxx, halfy, maxy)
            if polySfcInSect1 >= polySfcInSect2:
                maxy = halfy
            else:
                miny = halfy

        dividex = not dividex

        if cnt == maxcnt:
            raise Exception("getPtInsidePolygon: unable to converge")
        cnt += 1


def getPtOutsidePolygon(poly, distance):
    def getMostDistantPolyPoint(ptInside):
        rslt = None
        dst2min = None
        xi, yi = ptInside
        for x, y in poly:
            dst2 = (x - xi) ** 2.0 + (y - yi) ** 2.0
            if (not dst2min) or (dst2min > dst2):
                dst2min = dst2
                rslt = x, y
        return rslt

    ptInside = getPtInsidePolygon(poly)
    polyPt = getMostDistantPolyPoint(ptInside)
    xi, yi = ptInside
    x1, y1 = polyPt
    if not almostEqual(xi, x1):
        m = (y1 - yi) / (x1 - xi)
        distxmod = (distance**2.0 / (1 + m**2.0)) ** 0.5
        distx = distxmod if x1 >= xi else -distxmod
        disty = m * distx
        xout = x1 + distx
        yout = y1 + disty
    else:
        xout = x1
        yout = y1 + distance
    return xout, yout


def segmentLength(pt0, pt1):
    if not (almostEqual(pt0[0], pt1[0], 10) and almostEqual(pt0[1], pt1[1], 10)):
        return ((pt0[0] - pt1[0]) ** 2.0 + (pt0[1] - pt1[1]) ** 2.0) ** 0.5
    else:
        return 0


def angleBetweenPoints(ptv, pt1, pt2):
    """
    angleBetweenPoints:
    returns the angle < pi between points ptv (xv, yv),  pt1 (x1, y1) and pt2 (x2, y2).
    ptv is the vertex.
    """
    lv1 = segmentLength(ptv, pt1)
    lv2 = segmentLength(ptv, pt2)
    l12 = segmentLength(pt1, pt2)
    if lv1 and lv2:
        angle = acos((lv1**2.0 + lv2**2.0 - l12**2.0) / (2.0 * lv1 * lv2))
        return angle
    else:
        return float("NaN")


def getPointOnLine(pt0, pt1, distanceFromPt0):
    """
    getPointOnLine:
    return a point of straight line passing through points pt0 and pt1, at a distance
    of distanceFromPt0 from point pt0.
    A positive sign of distanceFromPt0 involves that the returned point is in the
    same direction as pt1 with respect to pt0.
    """
    x0, y0 = pt0
    x1, y1 = pt1
    if not almostEqual(x0, x1):
        m = (y1 - y0) / (x1 - x0)
        dx = (distanceFromPt0**2.0 / (1.0 + m**2.0)) ** 0.5
        dx = dx if x1 > x0 else -dx
        dx = dx if distanceFromPt0 > 0 else -dx
        dy = m * dx
    else:
        dx = 0
        dy = distanceFromPt0 if y1 > y0 else -distanceFromPt0
    return x0 + dx, y0 + dy


def smoothSteepAnglesOfPolygon(points, steepangle, maxCutLength=0):
    """
    smoothSteepAnglesOfPolygon:
    useful when working with unstructured grids. If a coastline polygon
    has too steep angles, involving a "fractal" generation of triangles,
    this function smoothes them adding intermediate points.
    steepangle: minumum accepted angle in radiants.
    """
    points_ = [points[0]]
    for pt0, pt1, pt2 in zip(points[:-2], points[1:-1], points[2:]):
        angle = angleBetweenPoints(pt1, pt0, pt2)
        if not isnan(angle) and (angle < steepangle):
            l01 = segmentLength(pt0, pt1)
            l12 = segmentLength(pt1, pt2)
            cutLength = min(l01, l12) / 2.0
            if maxCutLength:
                cutLength = min(cutLength, maxCutLength)
            pt0_ = getPointOnLine(pt1, pt0, cutLength)
            pt1_ = getPointOnLine(pt1, pt2, cutLength)
            points_.append(pt0_)
            points_.append(pt1_)
        else:
            points_.append(pt1)
    points_.append(points[-1])
    return points_


def bias(srs, refSrs):
    """
    bias: computes bias betwee srs and refSrs
    """
    result = 0.0
    n = float(len(srs))
    for s, r in zip(srs, refSrs):
        result += (s - r) / n
    return result


def nbias(srs, refSrs):
    """
    nbias: computes normalized bias between srs and refSrs
    """
    if not len(srs):
        return 0.0

    bi = bias(srs, refSrs)
    avgref = sum(refSrs) / len(refSrs)
    if avgref:
        nbi = bi / avgref
        if abs(nbi) > 1.0e10:
            # avgref is almost 0.
            nbi = 0.0
        return nbi
    else:
        return 0.0


def rmse(srs, refSrs):
    """
    rmse: root mean square error
    """
    if not len(srs):
        return 0.0
    sqsum = 0.0
    n = float(len(srs))
    for s, r in zip(srs, refSrs):
        sqsum += (s - r) ** 2.0
    sqavg = sqsum / n
    return sqavg**0.5


def nrmse(srs, refSrs):
    """
    rmse: normalized root mean square error
    """
    sqsum = 0.0
    sqrefsum = 0.0
    for s, r in zip(srs, refSrs):
        sqsum += (s - r) ** 2.0
        sqrefsum += r**2.0
    if sqrefsum:
        return (sqsum / sqrefsum) ** 0.5
    else:
        return 0.0


def hh(srs, refSrs):
    """
    hh: simmetrically normalized root mean square error:
    hh = [ sum(srs[i] - refSrs[i])**2. / sum(srs[i] * refSrs[i])  ]**0.5
    """
    if not len(srs):
        return 0
    sqsum = 0.0
    srsbyrefsum = 0.0
    n = float(len(srs))
    for s, r in zip(srs, refSrs):
        sqsum += (s - r) ** 2.0
        srsbyrefsum += s * r
    if srsbyrefsum:
        h = (sqsum / srsbyrefsum) ** 0.5
        if abs(h) > 1.0e10:
            # srsbyrefsum is almost 0.
            h = 0.0
        return h
    else:
        return 0.0
