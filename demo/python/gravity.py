#!/usr/bin/env python3
import sys
from astronomy import ObserverGravity

UsageText = r'''

    USAGE:

    gravity.py latitude height

    Calculates the gravitational acceleration experienced
    by an observer on the surface of the Earth at the specified
    latitude (degrees north of the equator) and height
    (meters above sea level).
    The output is the gravitational acceleration in m/s^2.

'''

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(UsageText)
        sys.exit(1)

    latitude = float(sys.argv[1])
    if latitude < -90.0 or latitude > +90.0:
        print("ERROR: Invalid latitude '{}'. Must be a number between -90 and +90.".format(sys.argv[1]))
        sys.exit(1)

    height = float(sys.argv[2])
    MAX_HEIGHT_METERS = 100000.0
    if height < 0.0 or height > MAX_HEIGHT_METERS:
        print("ERROR: Invalid height '{}'. Must be a number between 0 and {}.".format(sys.argv[1], MAX_HEIGHT_METERS))
        sys.exit(1)

    gravity = ObserverGravity(latitude, height)
    print("latitude = {:8.4f},  height = {:6.0f},  gravity = {:8.6f}".format(latitude, height, gravity))
    sys.exit(0)
