#!/usr/bin/env python3
#
#   usr_star.py  -  by Emard - 2026-01-20
#
#   Example program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   for user defined star ra/dec, calculate az/alt position
#   at some date and geographic lat/lon
#
import sys
from astronomy import Body, Body_name, Refraction, Equator, Horizon, DefineStar
from astro_demo_common import ParseArgs

def run():
    observer, time = ParseArgs(["pos.py", "+46", "+16", "2026-01-29T23:00:00.0Z"])
    print('UTC date = {}'.format(time))
    print()
    # star with RA and Dec (e.g., Sirius)
    star_ra  =    6+45/60+ 9/3600  # [h] right ascension hour angle
    star_dec = -(16+42/60+58/3600) # [°] declination degrees
    DefineStar(body=Body.Star1,ra=star_ra,dec=star_dec,distanceLightYears=8.6)
    # ra=101.2870833° dec=-16.7161111° Sirius
    #DefineStar(body=Body.Star1,ra=101.2870833/15,dec=-16.7161111,distanceLightYears=8.6)
    # ra=101.58051° dec=-16.7541° Sirius, stellarim on date
    DefineStar(body=Body.Star2,ra=101.58051/15,dec=-16.7541,distanceLightYears=8.6)
    # ra= 88.7929167° dec=  7.4069444° Betelgeuse
    DefineStar(body=Body.Star3,ra=88.7929167/15,dec=7.4069444,distanceLightYears=900)
    print('BODY           RA°     DEC°      AZ°     ALT°')
    body_list = [
        Body.Star1, Body.Star2, Body.Star3,
    ]
    body_comment = { Body.Star1:"Sirius J2000", Body.Star2:"Sirius on date by stellarium (wrong)", Body.Star3:"Betelgeuse J2000" }
    for body in body_list:
        equ_2000 = Equator(body, time, observer, ofdate=False, aberration=True)
        equ_ofdate = Equator(body, time, observer, ofdate=True, aberration=True)
        hor = Horizon(time, observer, equ_ofdate.ra, equ_ofdate.dec, Refraction.Normal)
        print('{:<8} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:s}'.format(Body_name[body], equ_2000.ra*15, equ_2000.dec, hor.azimuth, hor.altitude, body_comment[body]))
run()
