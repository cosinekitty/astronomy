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

KM_PER_AU = 1.4959787069098932e+8
SUN_RADIUS_AU  = 4.6505e-3
SUN_RADIUS_KM = KM_PER_AU * SUN_RADIUS_AU
ATMOSPHERE_KM = 65.4
EARTH_RADIUS_KM = 6371.0 + ATMOSPHERE_KM
MOON_RADIUS_KM = 1737.4

def main(text):
    time = Time.Parse(text)
    print('time=', time, 'ut=', time.ut, 'tt=', time.tt, 'dt=', 86400*(time.tt-time.ut))

    # Find heliocentric vector for Earth and Moon at this time.
    e = HelioVector(Body.Earth, time)
    print('Earth', e.x, e.y, e.z)

    m = GeoMoon(time)
    print('Moon', m.x, m.y, m.z)

    # Not needed... just checking distance ratios.
    edist = math.sqrt(e.x**2 + e.y**2 + e.z**2)
    mdist = math.sqrt((e.x+m.x)**2 + (e.y+m.y)**2 + (e.z+m.z)**2)
    print('Earth distance=', edist, '  Moon distance=', mdist, '  ratio=', mdist/edist)
    print('Sun radius km=', SUN_RADIUS_KM)

    # Find the parametric value u where the Earth's shadow ray
    # passes closest to the center of the Moon.
    u = (e.x*m.x + e.y*m.y + e.z*m.z) / (e.x*e.x + e.y*e.y + e.z*e.z)
    print('u=', u)

    # Find the distance between the shadow ray and the center of the Moon.
    r = KM_PER_AU * math.sqrt((u*e.x - m.x)**2 + (u*e.y - m.y)**2 + (u*e.z - m.z)**2)
    print('r=', r, 'km  r-M=', r-MOON_RADIUS_KM, '  r+M=', r+MOON_RADIUS_KM)

    # Penumbra radius
    pr = -SUN_RADIUS_KM + (1+u)*(SUN_RADIUS_KM + EARTH_RADIUS_KM)
    print('penumbra radius =', pr)

    # Umbra radius
    ur = SUN_RADIUS_KM - (1+u)*(SUN_RADIUS_KM - EARTH_RADIUS_KM)
    print('umbra radius =', ur, 'DELTA=', r+MOON_RADIUS_KM-ur)

    # Is the Moon totally inside the umbra?
    umbra1 = -ur <= (r-MOON_RADIUS_KM) <= +ur
    umbra2 = -ur <= (r+MOON_RADIUS_KM) <= +ur
    pen1 = -pr <= (r-MOON_RADIUS_KM) <= +pr
    pen2 = -pr <= (r+MOON_RADIUS_KM) <= +pr
    if umbra1 and umbra2:
        print('Total eclipse')
    elif umbra1 or umbra2:
        print('Partial eclipse')
    elif pen1 or pen2:
        print('Penumbral eclipse')
    else:
        print('No eclipse')
    return 0

def FindPeak(text1, text2):
    t1 = Time.Parse(text1)
    t2 = Time.Parse(text2)
    dt = 1.0 / (24.0 * 3600.0)
    t = t1
    tmin = None
    errmin = None
    errprev = None
    translist = []
    while t.ut <= t2.ut:
        e = HelioVector(Body.Earth, t)
        m = GeoMoon(t)
        u = (e.x*m.x + e.y*m.y + e.z*m.z) / (e.x*e.x + e.y*e.y + e.z*e.z)
        ur = SUN_RADIUS_KM - (1+u)*(SUN_RADIUS_KM - EARTH_RADIUS_KM)
        r = KM_PER_AU * math.sqrt((u*e.x - m.x)**2 + (u*e.y - m.y)**2 + (u*e.z - m.z)**2)
        err = r + MOON_RADIUS_KM - ur
        if (errmin is None) or (err < errmin):
            tmin = t
            errmin = err
        if (errprev is not None) and (errprev * err) <= 0.0:
            print('{} transition'.format(t))
            translist.append(t)
        t = t.AddDays(dt)
        errprev = err

    if errmin is None:
        print('NOT FOUND')
    else:
        print('{} PEAK @ {:0.2f} km'.format(tmin, errmin))

    if len(translist) == 2:
        totalmins = (translist[1].tt - translist[0].tt) * (24.0 * 60.0)
        print('Total eclipse lasted {:0.2f} minutes.'.format(totalmins))

if __name__ == '__main__':
    if len(sys.argv) == 3:
        sys.exit(FindPeak(sys.argv[1], sys.argv[2]))

    if len(sys.argv) == 2:
        text = sys.argv[1]
    else:
        text = '2021-05-26T11:18:40Z'
    sys.exit(main(text))
