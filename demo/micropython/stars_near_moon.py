#!/usr/bin/env python3
#
#   stars_near_moon.py  -  by Don Cross - 2021-11-08
#
#   Example program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   Given an observer's geographic latitude and longitude,
#   and an optional date and time, this program displays a
#   a list of bright stars that appear near the Moon in the sky.
#
import sys
import csv
from astronomy import GeoMoon, ObserverVector, Spherical, VectorFromSphere, AngleBetween
from astro_demo_common import ParseArgs

if __name__ == '__main__':
    MAG_LIMIT = 4.0     # dimmest star we will consider
    SEP_LIMIT = 1.0     # max angular separation (degrees)
    observer, time = ParseArgs(["stars_near_moon.py", "30","-80", "2021-11-08T23:00:00Z"])
    print('UTC date = {}'.format(time))
    print()

    # Get the Moon's equatorial coordinates as seen by
    # the given observer at the specified time.
    moon_vec = GeoMoon(time) - ObserverVector(time, observer, False)

    with open('../../generate/hyg_v36_1.csv') as starfile:
        reader = csv.DictReader(starfile)
        lnum = 0
        for row in reader:
            lnum += 1
            if lnum > 1:        # skip the Sun, which does not have fixed RA/DEC
                mag = float(row['mag'])
                if mag <= MAG_LIMIT:
                    star_ra = float(row['ra'])
                    star_dec = float(row['dec'])
                    # Convert equatorial angular coordinates to vector
                    star_vec = VectorFromSphere(Spherical(star_dec, star_ra*15.0, 1.0), time)
                    angle = AngleBetween(moon_vec, star_vec)
                    if angle <= SEP_LIMIT:
                        name = row['bf']
                        if row['proper']:
                            name += ' (' + row['proper'] + ')'
                        print('{} mag={}, angle={:0.3f} deg'.format(name, mag, angle))
    sys.exit(0)
