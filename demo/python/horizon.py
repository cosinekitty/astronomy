#!/usr/bin/env python3
#
#    horizon.py  -  by Don Cross - 2019-12-18
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    This is a more advanced example. It shows how to use coordinate
#    transforms and a binary search to find the two azimuths where the
#    ecliptic intersects with an observer's horizon at a given date and time.
#
#    To execute, run the command:
#
#    python3 horizon.py latitude longitude [yyyy-mm-ddThh:mm:ssZ]
#
import sys
from astronomy import *
from astro_demo_common import ParseArgs
from typing import Tuple

NUM_SAMPLES = 4

def ECLIPLON(i: int) -> float:
    return (360.0 * i) / NUM_SAMPLES


def HorizontalCoords(ecliptic_longitude: float, time: Time, rot_ecl_hor: RotationMatrix) -> Spherical:
    eclip = Spherical(
        0.0,                    # being "on the ecliptic plane" means ecliptic latitude is zero.
        ecliptic_longitude,
        1.0                     # any positive distance value will work fine.
    )

    # Convert ecliptic angular coordinates to ecliptic vector.
    ecl_vec = VectorFromSphere(eclip, time)

    # Use the rotation matrix to convert ecliptic vector to horizontal vector.
    hor_vec = RotateVector(rot_ecl_hor, ecl_vec)

    # Find horizontal angular coordinates, correcting for atmospheric refraction.
    return HorizonFromVector(hor_vec, Refraction.Normal)


def SearchCrossing(time: Time, rot_ecl_hor: RotationMatrix, e1: float, e2: float) -> Tuple[float, Spherical]:
    tolerance = 1.0e-6      # one-millionth of a degree is close enough!
    # Binary search: find the ecliptic longitude such that the horizontal altitude
    # ascends through a zero value. The caller must pass e1, e2 such that the altitudes
    # bound zero in ascending order.
    while True:
        e3 = (e1 + e2) / 2.0
        h3 = HorizontalCoords(e3, time, rot_ecl_hor)
        if abs(e2-e1) < tolerance:
            return (e3, h3)
        if h3.lat < 0.0:
            e1 = e3
        else:
            e2 = e3


def FindEclipticCrossings(observer: Observer, time: Time) -> int:
    # The ecliptic is a celestial circle that describes the mean plane of
    # the Earth's orbit around the Sun. We use J2000 ecliptic coordinates,
    # meaning the x-axis is defined to where the plane of the Earth's
    # equator on January 1, 2000 at noon UTC intersects the ecliptic plane.
    # The positive x-axis points toward the March equinox.
    # Calculate a rotation matrix that converts J2000 ecliptic vectors
    # to horizontal vectors for this observer and time.
    rot = Rotation_ECL_HOR(time, observer)

    # Sample several points around the ecliptic.
    # Remember the horizontal coordinates for each sample.
    hor = [HorizontalCoords(ECLIPLON(i), time, rot) for i in range(NUM_SAMPLES)]

    for i in range(NUM_SAMPLES):
        a1 = hor[i].lat
        a2 = hor[(i+1) % NUM_SAMPLES].lat
        e1 = ECLIPLON(i)
        e2 = ECLIPLON(i+1)
        if a1 * a2 <= 0.0:
            if a2 > a1:
                (ex, h) = SearchCrossing(time, rot, e1, e2)
            else:
                (ex, h) = SearchCrossing(time, rot, e2, e1)

            if h.lon > 0.0 and h.lon < 180.0:
                direction = 'ascends'
            else:
                direction = 'descends'
            print('Ecliptic longitude {:0.4f} {} through horizon az {:0.4f}, alt {:0.5g}'.format(ex, direction, h.lon, h.lat))

    return 0

if __name__ == '__main__':
    observer, time = ParseArgs(sys.argv)
    sys.exit(FindEclipticCrossings(observer, time))
