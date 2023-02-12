#!/usr/bin/env python3
#
#   solar_time.py  -  by Don Cross - 2023-02-12
#
#   Example program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   Given an observer's geographic latitude and longitude,
#   and an optional date and time, this program displays
#   true solar time for that observer and time.
#
import sys
from math import fmod
from astronomy import Body, HourAngle
from astro_demo_common import ParseArgs

if __name__ == '__main__':
    observer, time = ParseArgs(sys.argv)
    ha = HourAngle(Body.Sun, time, observer)
    solarTimeHours = fmod(ha + 12.0, 24.0)
    milli = int(round(solarTimeHours * 3.6e+6))
    second = milli // 1000
    milli %= 1000
    minute = second // 60
    second %= 60
    hour = minute // 60
    minute %= 60
    hour %= 24
    print('True solar time = {:0.4f} hours ({:02d}:{:02d}:{:02d}.{:03d})'.format(solarTimeHours, hour, minute, second, milli))
    sys.exit(0)
