#!/usr/bin/env python3
#
#   lunarday.py  -  Example of calculating lunar days.
#
import sys
import math
from astronomy import Time, MoonPhase, SearchMoonPhase

def main(args):
    if len(args) == 1:
        time = Time.Now()
    elif len(args) == 2:
        time = Time.Parse(args[1])
    else:
        print('USAGE: {} [yyyy-mm-ddThh:mm:ssZ]'.format(args[0]))
        return 1

    # Get the Moon phase at the starting time, so we know what
    # multiple of 12 degrees to search for next.
    phase = MoonPhase(time)
    print("{} : Moon's phase angle = {:0.6f} degrees.".format(time, phase))
    print()
    print('The next 30 lunar days are:')

    target = 12.0 * math.ceil(phase / 12.0)
    for _ in range(30):
        time = SearchMoonPhase(target, time, 2.0)
        if time is None:
            print('SEARCH FAILURE.')
            return 0
        d = 1 + int(target / 12)
        print('{} : Lunar Day {}'.format(time, d))
        target = math.fmod(target + 12.0, 360.0)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
