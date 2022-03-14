from math import *
import numpy as np
from src.geodiccaUtility import *


def angleDiff(angle1, angle2):
    diff = angle1 - angle2
    while diff > pi:
        diff -= 2 * pi
    while diff < -pi:
        diff += 2 * pi
    return diff


def angleDiffDeg(angle1, angle2):
    angle1 = angle1 * deToRa
    angle2 = angle2 * deToRa
    diffRad = angleDiff(angle1, angle2)
    return diffRad * raToDe


def angleMean(angles):
    l = len(angles)
    if not l:
        return 9999.0

    x = [cos(a) for a in angles]
    y = [sin(a) for a in angles]
    meanx = sum(x) / l
    meany = sum(y) / l

    return atan2(meany, meanx)


def angleMeanDeg(angles):
    radangles = [a * deToRa for a in angles]
    radmean = angleMean(radangles)
    return radmean * raToDe


def angleBiasDeg(angles1, angles2):
    # assert len(angles1) == len(angles2), 'angle arrays must be of the same length'

    # return angleDiffDeg(angleMeanDeg(angles1), angleMeanDeg(angles2))
    diffs = [angleDiffDeg(a1, a2) for a1, a2 in zip(angles1, angles2)]
    return sum(diffs) / len(diffs)


def angleRMSEDeg(angles1, angles2):
    assert len(angles1) == len(angles2), "angle arrays must be of the same length"

    l = len(angles1)
    if not l:
        return 0.0

    diffs = [angleDiffDeg(ang1, ang2) for ang1, ang2 in zip(angles1, angles2)]
    sqdiffs = [d**2.0 for d in diffs]
    return (sum(sqdiffs) / l) ** 0.5


def _adjustLonLatPoly(poly, polyWest, polyEast):
    if polyWest <= polyEast:
        _poly = poly
    else:
        _poly = []
        for p in poly:
            (px, py) = p
            while px < polyWest:
                px += 360.0
            _poly.append((px, py))
    return _poly


def adjustLons(lons, west, east):
    if west <= east:
        _lons = lons
    else:
        _lons = []
        for l in lons:
            while l < west:
                l += 360.0
            _lons.append(l)
    return _lons


def pointInsidePolygonLonLatDeg(x, y, poly, polyWest, polyEast):
    """
    Determine if a point is inside a given polygon or not.
    parameters polyWest and polyEast are necessary to let the
    system know what are the westmost and eastmost points of the
    domain. Without those parameters there are ambiguous situations,
    e.g. poly = ((-90, ymin), (-90, ymax), (90, ymax), (90, ymin), (-90, ymin))
    """
    _poly = _adjustLonLatPoly(poly, polyWest, polyEast)
    if polyWest > polyEast:
        while x < polyWest:
            x += 360.0
    return pointInsidePolygon(x, y, _poly)


def getPtInsidePolygonLonLatDeg(poly, polyWest, polyEast):
    """
    Gets point inside polygon.
    parameters polyWest and polyEast are necessary to let the
    system know what are the westmost and eastmost points of the
    domain. Without those parameters there are ambiguous situations,
    e.g. poly = ((-90, ymin), (-90, ymax), (90, ymax), (90, ymin), (-90, ymin))
    """
    _poly = _adjustLonLatPoly(poly, polyWest, polyEast)
    return getPtInsidePolygon(_poly)


def rawConvertMetersToLonLat(x, y, lon0):
    lat = y / (deToRa * earthRadius)
    lon = x / (deToRa * np.cos(lat * deToRa) * earthRadius) + lon0
    return lon, lat


def rawConvertLonLatToMeters(lon, lat, lon0):
    x = (lon - lon0) * deToRa * np.cos(lat * deToRa) * earthRadius
    y = lat * deToRa * earthRadius
    return x, y


def rawConvertLonLatPointsToMeters(lons, lats, minlon, maxlon):
    lons = adjustLons(lons, minlon, maxlon)
    xs = []
    ys = []
    for lon, lat in zip(lons, lats):
        x, y = rawConvertLonLatToMeters(lon, lat)
        xs.append(x)
        ys.append(y)
    return xs, ys


def rawConvertLonLatPolyToMeters(poly, minlon, maxlon):
    lons = [p[0] for p in poly]
    lats = [p[1] for p in poly]
    xs, ys = rawConvertLonLatPointsToMeters(lons, lats, minlon, maxlon)
    polym = [(x, y) for (x, y) in zip(xs, ys)]
    return polym


def getPtOutsidePolygonLonLatDeg(poly, polyWest, polyEast, distanceInMeters):
    """
    Gets point outside polygon at a distance of distanceInMeters.
    parameters polyWest and polyEast are necessary to let the
    system know what are the westmost and eastmost points of the
    domain. Without those parameters there are ambiguous situations,
    e.g. poly = ((-90, ymin), (-90, ymax), (90, ymax), (90, ymin), (-90, ymin))
    """
    polym = rawConvertLonLatPolyToMeters(poly, polyWest, polyEast)
    xm, ym = getPtOutsidePolygon(polym, distanceInMeters)
    return rawConvertMetersToLonLat(xm, ym)


def rectanglesIntersectLonLatDeg(
    minx0, maxx0, miny0, maxy0, minx1, maxx1, miny1, maxy1
):
    while maxx0 < minx0:
        maxx0 += 360.0
    while maxx1 < minx1:
        maxx1 += 360.0
    return rectanglesIntersect(minx0, maxx0, miny0, maxy0, minx1, maxx1, miny1, maxy1)


def polygonSurfaceLonLatDeg(poly, polyWest, polyEast):
    """
    Computation of the surface enclosed in a polygon applying Green's theorem.
    Polygon vertex are in deg, returned surface is in squared meters.
    A euclidean geometry is used, hence use only for small surfaces.
    """
    polym = rawConvertLonLatPolyToMeters(poly, polyWest, polyEast)
    return polygonSurface(polym)


def liesBetweenLonDeg(west, x, east):
    if hasattr(x, "__iter__"):
        return np.array([liesBetweenLonDeg(west, xi, east) for xi in x])
    else:
        while east < west:
            east += 360.0
            while x < west:
                x += 360.0
        return west <= x <= east
