#!/usr/bin/env python3
#
#   galactic.py  -  by Don Cross - 2021-06-14
#
#   Example Python program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   Demo of converting galactic coordinates to horizontal coordinates.
#

import sys
import math
from astronomy import *

UsageText = r'''
USAGE: galactic olat olon glat glon [yyyy-mm-ddThh:mm:ssZ]

where

    olat = observer's latitude on the Earth
    olon = observer's longitude on the Earth
    glat = IAU 1958 galatic latitude of the target
    glon = IAU 1958 galatic longitude of the target
    yyyy-mm-ddThh:mm:ssZ = optional UTC date/time

Given the galactic coordinates of a point source in the sky,
this program calculates horizontal aiming coordinates for an
observer on or near the Earth's surface.

If the date/time is given on the command line, it is used.
Otherwise, the computer's current date/time is used.
'''


def GalacticToHorizontal(time, observer, glat, glon):
    # Calculate a matrix that converts galactic coordinates
    # to J2000 equatorial coordinates.
    rot = Rotation_GAL_EQJ()

    # Adjust the rotation matrix to convert galactic to horizontal.
    adjust_rot = Rotation_EQJ_HOR(time, observer)
    rot = CombineRotation(rot, adjust_rot)

    # Convert the galactic coordinates from angles to a unit vector.
    gsphere = Spherical(glat, glon, 1.0)
    gvec = VectorFromSphere(gsphere, time)

    # Use the rotation matrix to convert the galactic vector to a horizontal vector.
    hvec = RotateVector(rot, gvec)

    # Convert the horizontal vector back to angular coordinates.
    # Assume this is a radio source (not optical), do not correct for refraction.
    hsphere = HorizonFromVector(hvec, Refraction.Airless)
    return hsphere.lat, hsphere.lon


if __name__ == '__main__':
    if len(sys.argv) not in [5, 6]:
        print(UsageText)
        sys.exit(1)

    olat = float(sys.argv[1])
    olon = float(sys.argv[2])
    observer = Observer(olat, olon)

    glat = float(sys.argv[3])
    glon = float(sys.argv[4])

    if len(sys.argv) > 5:
        time = Time.Parse(sys.argv[5])
    else:
        time = Time.Now()

    altitude, azimuth = GalacticToHorizontal(time, observer, glat, glon)
    print('altitude = {:10.3f}, azimuth = {:10.3f}'.format(altitude, azimuth))
    sys.exit(0)
