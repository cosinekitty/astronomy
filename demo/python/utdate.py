#!/usr/bin/env python3
#
#   utdate.py  -  Don Cross  -  2020-01-03
#
#   Convert J2000 UT day number to UTC date/time.
#   Also displays Terrestrial Time (TT).
#
import sys
import astronomy

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('USAGE: utdate.py ut')
        sys.exit(1)
    ut = float(sys.argv[1])
    time = astronomy.Time(ut)
    print(time)
    print('TT = {:0.16f}'.format(time.tt))
    sys.exit(0)
