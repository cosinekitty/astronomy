#!/usr/bin/env python3
#
#    riseset.py  -  by Don Cross - 2019-07-28
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    This program calculates sunrise, sunset, moonrise, and moonset
#    times for an observer at a given latitude and longitude.
#
#    To execute, run the command:

#    python3 riseset.py latitude longitude [yyyy-mm-ddThh:mm:ssZ]
#
import sys
from astronomy import Body, Direction, SearchRiseSet
from astro_demo_common import ParseArgs

def PrintEvent(name, time):
    if time is None:
        raise Exception('Failure to calculate ' + name)
    print('{:<8s} : {}'.format(name, time))

def main(args):
    observer, time = ParseArgs(args)
    sunrise  = SearchRiseSet(Body.Sun,  observer, Direction.Rise, time, 300)
    sunset   = SearchRiseSet(Body.Sun,  observer, Direction.Set,  time, 300)
    moonrise = SearchRiseSet(Body.Moon, observer, Direction.Rise, time, 300)
    moonset  = SearchRiseSet(Body.Moon, observer, Direction.Set,  time, 300)
    PrintEvent('search',   time)
    PrintEvent('sunrise',  sunrise)
    PrintEvent('sunset',   sunset)
    PrintEvent('moonrise', moonrise)
    PrintEvent('moonset',  moonset)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
