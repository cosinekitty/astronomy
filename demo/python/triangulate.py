#!/usr/bin/env python3
import sys
from astronomy import *

UsageText = r'''
USAGE:  triangulate.py  lat1 lon1 elv1 az1 alt1  lat2 lon2 elv2 az2 alt2

Calculate the best-fit location of a point as observed
from two different locations on or near the Earth's surface.

lat1, lat2 = Geographic latitudes in degrees north of the equator.
lon1, lon2 = Geographic longitudes in degrees east of the prime meridian.
elv1, elv2 = Elevations above sea level in meters.
az1,  az2  = Azimuths toward observed object in degrees clockwise from north.
alt1, alt2 = Altitude angles toward observed object in degrees above horizon.

This program extrapolates lines in the given directions from the two
geographic locations and finds the location in space where they
come closest to intersecting. It then prints out the coordinates
of that triangulation point, along with the error radius in meters.
'''

def DotProduct(a, b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def Intersect(pos1, dir1, pos2, dir2):
    F = DotProduct(dir1, dir2)
    amb = Vector(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z, pos1.t)
    E = DotProduct(dir1, amb)
    G = DotProduct(dir2, amb)
    # 0 = E + u  - Fv
    # 0 = G + Fu - v
    denom = 1.0 - F*F
    if denom == 0.0:
        print('Cannot solve because directions are parallel.')
    else:
        u = (F*G - E) / denom
        v = G + F*u
        print('u={}, v={}'.format(u, v))

if __name__ == '__main__':
    # Validate and parse command line arguments.

    if len(sys.argv) != 11:
        print(UsageText)
        sys.exit(1)

    lat1, lon1, elv1, az1, alt1, lat2, lon2, elv2, az2, alt2 = [float(x) for x in sys.argv[1:]]
    obs1 = Observer(lat1, lon1, elv1)
    obs2 = Observer(lat2, lon2, elv2)

    # Convert geographic coordinates into 3D vectors.
    # We can use an arbitrary time because we don't care about the Earth's rotation
    # with respect to extraterrestrial bodies.

    time = Time(0.0)

    pos1 = ObserverVector(time, obs1, True)
    print('Observer #1 = {}'.format(obs1))
    print('Position #1 = ({:0.6e}, {:0.6e}, {:0.6e})'.format(pos1.x, pos1.y, pos1.z))

    pos2 = ObserverVector(time, obs2, True)
    print('Observer #2 = {}'.format(obs2))
    print('Position #2 = ({:0.6e}, {:0.6e}, {:0.6e})'.format(pos2.x, pos2.y, pos2.z))

    # Convert horizontal coordinates into unit direction vectors.

    hor1 = Spherical(alt1, az1, 1.0)
    hor2 = Spherical(alt2, az2, 1.0)
    dir1 = VectorFromHorizon(hor1, time, Refraction.Airless)
    dir2 = VectorFromHorizon(hor2, time, Refraction.Airless)
    print('Direction #1 = ({:0.6e}, {:0.6e}, {:0.6e})'.format(dir1.x, dir1.y, dir1.z))
    print('Direction #2 = ({:0.6e}, {:0.6e}, {:0.6e})'.format(dir2.x, dir2.y, dir2.z))

    # Solve for the target point.
    Intersect(pos1, dir1, pos2, dir2)

    # Solve again with the inputs reversed, as a sanity check.
    Intersect(pos2, dir2, pos1, dir1)

    sys.exit(0)
