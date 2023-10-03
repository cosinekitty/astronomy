#!/usr/bin/env python3
#
#    moonphase.py  -  by Don Cross - 2019-07-26
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    This program calculates the Moon's ecliptic phase and illumination percentage
#    for a given date and time, or for the computer's current date and time if none is given.
#    It also finds the dates and times of the subsequent 10 quarter phase changes.
#
#    To execute, run the command:
#    python3 moonphase.py [date]
#
import sys
import astronomy
from typing import List

def QuarterName(quarter: int) -> str:
    return [
        'New Moon',
        'First Quarter',
        'Full Moon',
        'Third Quarter'
    ][quarter]

def main(args: List[str]) -> int:
    if len(args) == 1:
        time = astronomy.Time.Now()
    elif len(args) == 2:
        time = astronomy.Time.Parse(args[1])
    else:
        print('USAGE: {} [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        return 1
    # Calculate the Moon's ecliptic phase angle,
    # which ranges from 0 to 360 degrees.
    #   0 degrees = new moon,
    #  90 degrees = first quarter,
    # 180 degrees = full moon,
    # 270 degrees = third quarter.
    phase = astronomy.MoonPhase(time)
    print("{} : Moon's ecliptic phase angle = {:0.3f} degrees.".format(time, phase))

    # Calculate the fraction of the Moon's disc
    # that appears illuminated, as seen from the Earth.
    illum = astronomy.Illumination(astronomy.Body.Moon, time)
    print("{} : Moon's illuminated fraction = {:0.2f}%.".format(time, 100.0 * illum.phase_fraction))
    print()

    # Predict when the next 10 lunar quarter phases will happen.
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
