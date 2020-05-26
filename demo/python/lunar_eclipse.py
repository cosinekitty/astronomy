#!/usr/bin/env python3
#
#    lunar_eclipse.py  -  by Don Cross - 2020-05-16
#
#    Example Python program for Astronomy Engine:
#    https://github.com/cosinekitty/astronomy
#
#    Searches for the next 10 partial/total lunar eclipses after
#    the current date, or a date specified on the command line.
#
#    To execute, run the command:
#    python3 lunar_eclipse.py [date]
#
import sys
from astronomy import Time, SearchLunarEclipse, NextLunarEclipse, EclipseKind


def PrintEclipse(e):
    # Calculate beginning/ending of different phases
    # of an eclipse by subtracting/adding the peak time
    # with the number of minutes indicated by the "semi-duration"
    # fields sd_partial and sd_total.
    MINUTES_PER_DAY = 24 * 60
    p1 = e.peak.AddDays(-e.sd_partial / MINUTES_PER_DAY)
    print('{}  Partial eclipse begins.'.format(p1))
    if e.sd_total > 0.0:
        t1 = e.peak.AddDays(-e.sd_total / MINUTES_PER_DAY)
        print('{}  Total eclipse begins.'.format(t1))
    print('{}  Peak of {} eclipse.'.format(e.peak, e.kind.name.lower()))
    if e.sd_total > 0.0:
        t2 = e.peak.AddDays(+e.sd_total / MINUTES_PER_DAY)
        print('{}  Total eclipse ends.'.format(t2))
    p2 = e.peak.AddDays(+e.sd_partial / MINUTES_PER_DAY)
    print('{}  Partial eclipse ends.'.format(p2))
    print()


def main(args):
    if len(args) == 1:
        time = Time.Now()
    elif len(args) == 2:
        time = Time.Parse(args[1])
    else:
        print('USAGE: {} [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        return 1

    count = 0
    e = SearchLunarEclipse(time)
    while count < 10:
        if e.kind != EclipseKind.Penumbral:
            count += 1
            PrintEclipse(e)
        e = NextLunarEclipse(e.peak)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
