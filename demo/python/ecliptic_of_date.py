#!/usr/bin/env python3
#
#   ecliptic_of_date.py  -  by Don Cross  -  2023-10-02
#
#   Example program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   This program shows how to calculate true ecliptic coordinates of date (ECT)
#   for the Sun, Moon, and planets for an observer at a given location
#   on the Earth.
#
import sys
from astronomy import Body, Equator, Rotation_EQD_ECT, RotateVector, SphereFromVector
from astro_demo_common import ParseArgs

if __name__ == '__main__':
    observer, time = ParseArgs(sys.argv)
    bodyList = [
        Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
        Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
    ]

    # Create a rotation matrix that converts equator-of-date (EQD) vectors
    # to true ecliptic of date (ECT) vectors.
    rot = Rotation_EQD_ECT(time)

    print('{:9s} {:>8s} {:>8s} {:>11s}'.format('body', 'lon', 'lat', 'dist'))
    for body in bodyList:
        # Calculate equatorial coordinates of date = EQD (ofdate parameter is True).
        eqd = Equator(body, time, observer, True, True)

        # The returned object `eqd` is of type `Equatorial`, which contains a vector.
        # Use the rotation matrix to convert that vector from (EQD) to ecliptic of date (ECT).
        ect = RotateVector(rot, eqd.vec)

        # Convert vector to spherical coordinate angles: ecliptic longitude, ecliptic latitude.
        sphere = SphereFromVector(ect)

        # Print the output for this body
        print('{:9s} {:8.3f} {:8.3f} {:11.6f}'.format(body.name, sphere.lon, sphere.lat, sphere.dist))

    sys.exit(0)
