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

Verbose = False

def Debug(text):
    if Verbose:
        print(text)

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
        Debug('Cannot solve because directions are parallel.')
        return None, None

    u = (F*G - E) / denom
    v = G + F*u
    Debug('u={}, v={}'.format(u, v))
    if u < 0 or v < 0:
        Debug('Lines of sight do not converge.')
        return None, None

    a = Vector(pos1.x + u*dir1.x, pos1.y + u*dir1.y, pos1.z + u*dir1.z, pos1.t)
    b = Vector(pos2.x + v*dir2.x, pos2.y + v*dir2.y, pos2.z + v*dir2.z, pos2.t)
    c = Vector((a.x + b.x)/2, (a.y + b.y)/2, (a.z + b.z)/2, a.t)
    Debug('c = {}'.format(c))
    # Calculate the error distance in meters between the two skew lines.
    dist = (KM_PER_AU * 1000) * math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2)
    # Convert vector back to geographic coordinates
    obs = VectorObserver(c, True)
    return obs, dist

def DirectionVector(time, observer, altitude, azimuth):
    # Convert horizontal angles to a horizontal unit vector.
    hor = Spherical(altitude, azimuth, 1.0)
    hvec = VectorFromHorizon(hor, time, Refraction.Airless)    
    # Find the rotation matrix that converts horizontal vectors to equatorial vectors.
    rot = Rotation_HOR_EQD(time, observer)
    # Rotate the horizontal (HOR) vector to an equator-of-date (EQD) vector.
    evec = RotateVector(rot, hvec)
    return evec

if __name__ == '__main__':
    # Validate and parse command line arguments.

    # The '-v' option enables debug prints.
    args = sys.argv[1:]
    if len(args) > 0 and args[0] == '-v':
        Verbose = True
        args = args[1:]

    if len(args) != 10:
        print(UsageText)
        sys.exit(1)

    lat1, lon1, elv1, az1, alt1, lat2, lon2, elv2, az2, alt2 = [float(x) for x in args]
    obs1 = Observer(lat1, lon1, elv1)
    obs2 = Observer(lat2, lon2, elv2)

    # Convert geographic coordinates into 3D vectors.
    # We can use an arbitrary time because we don't care about the Earth's rotation
    # with respect to extraterrestrial bodies.

    time = Time(0.0)

    pos1 = ObserverVector(time, obs1, True)
    Debug('Observer #1 = {}'.format(obs1))
    Debug('Position #1 = ({:0.6f}, {:0.6f}, {:0.6f})'.format(pos1.x * KM_PER_AU, pos1.y * KM_PER_AU, pos1.z * KM_PER_AU))

    pos2 = ObserverVector(time, obs2, True)
    Debug('Observer #2 = {}'.format(obs2))
    Debug('Position #2 = ({:0.6f}, {:0.6f}, {:0.6f})'.format(pos2.x * KM_PER_AU, pos2.y * KM_PER_AU, pos2.z * KM_PER_AU))

    # Convert horizontal coordinates into unit direction vectors in EQD orientation.
    dir1 = DirectionVector(time, obs1, alt1, az1)
    dir2 = DirectionVector(time, obs2, alt2, az2)

    Debug('Direction #1 = ({:0.6f}, {:0.6f}, {:0.6f})'.format(dir1.x, dir1.y, dir1.z))
    Debug('Direction #2 = ({:0.6f}, {:0.6f}, {:0.6f})'.format(dir2.x, dir2.y, dir2.z))

    # Solve for the target point.
    obs, error = Intersect(pos1, dir1, pos2, dir2)
    if obs is None:
        print('ERROR: Could not find intersection.')
        sys.exit(1)

    print('Solution #1 = {}, err = {:0.3f} meters'.format(obs, error))

    # Solve again with the inputs reversed, as a sanity check.
    check_obs, check_error = Intersect(pos2, dir2, pos1, dir1)
    if check_obs is None:
        print('INTERNAL ERROR: inconsistent solution.')
        sys.exit(1)

    print('Solution #2 = {}, err = {:0.3f} meters'.format(check_obs, check_error))

    sys.exit(0)
