#!/usr/bin/env python3
#
#   solar_eclipse.py  -  by Don Cross  -  2023-10-03
#
#   Example program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   Predicts the next 10 solar eclipses visible from a given location on the Earth.
#
import sys
from astronomy import LocalSolarEclipseInfo, SearchLocalSolarEclipse, NextLocalSolarEclipse
from astro_demo_common import ParseArgs


def PrintEclipse(count: int, eclipse: LocalSolarEclipseInfo) -> None:
    print('Eclipse #{} : {}, obscuration = {:0.3f}'.format(count, eclipse.kind.name, eclipse.obscuration))
    print('{}  Partial begins at altitude = {:0.2f}'.format(eclipse.partial_begin.time, eclipse.partial_begin.altitude))
    if eclipse.total_begin is not None:
        print('{}  Total begins at altitude = {:0.2f}'.format(eclipse.total_begin.time, eclipse.total_begin.altitude))
    print('{}  Peak at altitude = {:0.2f}'.format(eclipse.peak.time, eclipse.peak.altitude))
    if eclipse.total_end is not None:
        print('{}  Total ends at altitude = {:0.2f}'.format(eclipse.total_end.time, eclipse.total_end.altitude))
    print('{}  Partial ends at altitude = {:0.2f}'.format(eclipse.partial_end.time, eclipse.partial_end.altitude))
    print()

if __name__ == '__main__':
    observer, startTime = ParseArgs(sys.argv)
    count = 0
    while count < 10:
        if count == 0:
            eclipse = SearchLocalSolarEclipse(startTime, observer)
        else:
            eclipse = NextLocalSolarEclipse(eclipse.peak.time, observer)
        # We ignore eclipses that occur while the center of the
        # Sun is below the horizon at the peak of the eclipse.
        if eclipse.peak.altitude > 0.0:
            count += 1
            PrintEclipse(count, eclipse)
    sys.exit(0)
