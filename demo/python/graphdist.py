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

import astronomy
import matplotlib.pyplot as plt

def GraphPlanetHelioDistance(body, year1, year2):
    t1 = astronomy.Time.Make(year1, 1, 1, 0, 0, 0)
    t2 = astronomy.Time.Make(year2, 1, 1, 0, 0, 0)
    npoints = 1200
    xlist = []
    ylist = []
    for i in range(npoints+1):
        year = year1 + ((i / npoints) * (year2 - year1))
        ut = t1.ut + ((i / npoints) * (t2.ut - t1.ut))
        time = astronomy.Time(ut)
        vec = astronomy.HelioVector(body, time)
        dist = vec.Length()
        xlist.append(year)
        ylist.append(dist)
    plt.plot(xlist, ylist, 'b.')
    plt.show()
    plt.close('all')

if __name__ == '__main__':
    GraphPlanetHelioDistance(astronomy.Body.Neptune, 2030, 2060)
