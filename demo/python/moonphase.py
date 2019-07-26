#!/usr/bin/env python3
#
#    moonphase.py  -  by Don Cross - 2019-07-26
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    This program calculates the Moon's phase for a given date and time,
#    or for the computer's current date and time if none is given.
#    It also finds the dates and times of the subsequent 10 quarter phase changes.
#
#    To execute, run the command:
#    python3 moonphase.py [date]
#
import sys
import astronomy
from astro_demo_common import ParseTime

def QuarterName(quarter):
    return [
        'New Moon',
        'First Quarter',
        'Full Moon',
        'Third Quarter'
    ][quarter]

def main(args):
    if len(args) == 1:
        time = astronomy.Time.Now()
    elif len(args) == 2:
        time = ParseTime(args[1])
    else:
        print('USAGE: {} [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        return 1
    phase = astronomy.MoonPhase(time)
    print("{} : Moon's phase angle = {:0.6f} degrees.".format(time, phase))
    print()
    print('The next 10 lunar quarters are:')
    for i in range(10):
        if i == 0:
            mq = astronomy.SearchMoonQuarter(time)
        else:
            mq = astronomy.NextMoonQuarter(mq)
        print('{} : {}'.format(mq.time, QuarterName(mq.quarter)))
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
