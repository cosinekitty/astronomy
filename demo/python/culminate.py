#!/usr/bin/env python3
#
# culminate.py  -  Don Cross  -  2019-08-10
#
# Example program for Astronomy Engine:
# https://github.com/cosinekitty/astronomy
#
# This example program shows how to calculate the time
# the Sun, Moon, and planets will next reach their highest point in the sky
# as seen by an observer at a given location on the Earth.
# This is called culmination, and is found by finding when
# each body's "hour angle" is 0.
#
# Having an hour angle of 0 is another way of saying that the body is
# crossing the meridian, the imaginary semicircle in the sky that passes
# from due north on the horizon, through the zenith (straight up),
# toward due south on the horizon. At this moment the body appears to
# have an azimuth of either 180 degrees (due south) or 0 (due north).

import sys
from astronomy import Body, SearchHourAngle
from astro_demo_common import ParseArgs

if __name__ == '__main__':
    observer, time = ParseArgs(sys.argv)
    body_list = [
        Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
        Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
    ]
    print('search   :', time)
    for body in body_list:
        evt = SearchHourAngle(body, observer, 0.0, time)
        print('{:<8s} : {}  altitude={:6.2f}  azimuth={:7.2f}'.format(body.name, evt.time, evt.hor.altitude, evt.hor.azimuth))
    sys.exit(0)
