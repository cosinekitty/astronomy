#!/usr/bin/env python3
#
# earth_shadow.py  -  Don Cross  -  2020-05-11
#
# Example program for Astronomy Engine:
# https://github.com/cosinekitty/astronomy
#
# Given a date and time, figure out how far the Moon's center
# is from the Earth's shadow axis (the ray passing from
# the center of the Sun through the center of the Earth).
# This is an experiment to help me figure out a good
# lunar eclipse algorithm.

import sys
import math
from astronomy import Time, Body, HelioVector, GeoMoon

def main(text):
    time = Time.Parse(text)

    # Find heliocentric vector for Earth and Moon at this time.
    e = HelioVector(Body.Earth, time)
    print('Earth', e.x, e.y, e.z)

    m = GeoMoon(time)
    print('Moon', m.x, m.y, m.z)

    # Find the parametric value u where the Earth's shadow ray
    # passes closest to the center of the Moon.
    u = (e.x*m.x + e.y*m.y + e.z*m.z) / (e.x*e.x + e.y*e.y + e.z*e.z)
    print('u=', u)

    # Find the distance between the shadow ray and the Moon.
    r = math.sqrt((u*e.x - m.x)**2 + (u*e.y - m.y)**2 + (u*e.z - m.z)**2)
    r *= 1.4959787069098932e+8      # convert AU to km
    print('r=', r, 'km')

    # TODO: Calculate the umbral and penumbral radii at the parametric location u.
    return 0

if __name__ == '__main__':
    if len(sys.argv) == 2:
        text = sys.argv[1]
    else:
        text = '2021-05-26T11:19:00Z'
    sys.exit(main(text))
