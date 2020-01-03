#!/usr/bin/env python3
#
#   graphdist.py  -  Don Cross  -  2020-01-01
#
#   Graphs the distance between Neptune and the Sun
#   as a function of time. Shows that the movement
#   of the Solar System Barycenter causes two
#   local minima and one local maximum of the distance
#   function within one small arc of Neptune's orbit.
#

import sys
import matplotlib.pyplot as plt
import astronomy

def GraphPlanetHelioDistance(body, t1, t2):
    npoints = 1200
    xlist = []
    ylist = []
    for i in range(npoints+1):
        ut = t1.ut + ((i / npoints) * (t2.ut - t1.ut))
        time = astronomy.Time(ut)
        vec = astronomy.HelioVector(body, time)
        dist = vec.Length()
        xlist.append(ut)
        ylist.append(dist)
    plt.plot(xlist, ylist, 'b.')
    plt.show()
    plt.close('all')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('USAGE:  graphdist.py planet yyyy-mm-ddThh:mm:ssZ yyyy-mm-ddThh:mm:ssZ')
        sys.exit(1)
    body = astronomy.Body[sys.argv[1]]
    time1 = astronomy.Time.Parse(sys.argv[2])
    time2 = astronomy.Time.Parse(sys.argv[3])
    GraphPlanetHelioDistance(body, time1, time2)
    sys.exit(0)
