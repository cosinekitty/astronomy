#!/usr/bin/env python3
#
#   pairlon.py  -  Don Cross  -  2021-04-26
#
#   Graphs the relative longitudes of two bodies
#   as seen from the Earth, over the specified time span.
#

import sys
import matplotlib.pyplot as plt
from astronomy import Time, Body, PairLongitude

def GraphPairLongitude(body1, body2, t1, t2, npoints):
    xlist = []
    ylist = []
    for i in range(npoints+1):
        ut = t1.ut + ((i / npoints) * (t2.ut - t1.ut))
        time = Time(ut)
        lon = PairLongitude(body1, body2, time)
        xlist.append(ut)
        ylist.append(lon)
    plt.plot(xlist, ylist, 'b.')
    plt.show()
    plt.close('all')

if __name__ == '__main__':
    if not (5 <= len(sys.argv) <= 6):
        print('USAGE:  pairlon.py body1 body2 yyyy-mm-ddThh:mm:ssZ yyyy-mm-ddThh:mm:ssZ [npoints]')
        sys.exit(1)
    body1 = Body[sys.argv[1]]
    body2 = Body[sys.argv[2]]
    time1 = Time.Parse(sys.argv[3])
    time2 = Time.Parse(sys.argv[4])
    if len(sys.argv) > 5:
        npoints = int(sys.argv[5])
    else:
        npoints = 1200
    GraphPairLongitude(body1, body2, time1, time2, npoints)
    sys.exit(0)
