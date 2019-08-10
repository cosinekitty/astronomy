#!/usr/bin/env python3
#
#   seasons.py  -  Don Cross  -  2019-08-10
#
#   Example C program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   This program demonstrates how to calculate the
#   equinoxes and solstices for a year.
#
import sys
from astronomy import Seasons

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: seasons.py year')
        sys.exit(1)
    year = int(sys.argv[1])
    seasons = Seasons(year)
    print('March equinox     :', seasons.mar_equinox)
    print('June solstice     :', seasons.jun_solstice)
    print('September equinox :', seasons.sep_equinox)
    print('December solstice :', seasons.dec_solstice)
    sys.exit(0)
