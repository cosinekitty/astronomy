#!/usr/bin/env python3
#
#   constellation.py  -  Don Cross  -  2021-04-30 Fri
#
#   Searches for the moon moving into different constellations.
#   This is an example of doing a search for the time of transition
#   of a body from one constellation to another.
#

import sys
from astronomy import Time, Body, Constellation, GeoVector, EquatorFromVector

#------------------------------------------------------------------------------

def BodyConstellation(body, time):
    vec = GeoVector(body, time, False)
    equ = EquatorFromVector(vec)
    return Constellation(equ.ra, equ.dec)

#------------------------------------------------------------------------------

def FindConstellationChange(body, c1, t1, t2):
    # Do a binary search to find what time between t1 and t2
    # is the boundary between being in constellation c1 and some other
    # constellation. Return the tuple (tx, cx), where tx is the
    # transition time, and cx is the constellation we are moving into.

    tolerance = 0.1 / (24.0 * 3600.0)   # one tenth of a second
    while True:
        dt = t2.ut - t1.ut
        tx = t1.AddDays(dt/2)
        cx = BodyConstellation(body, tx)
        if cx.symbol == c1.symbol:
            t1 = tx
        else:
            if dt < tolerance:
                # Always end the search inside the new constellation.
                return (tx, cx)
            t2 = tx

#------------------------------------------------------------------------------

def FindConstellationChanges(body, startTime, stopTime, dayIncrement):
    t1 = startTime
    c1 = BodyConstellation(body, t1)
    while t1 < stopTime:
        t2 = t1.AddDays(dayIncrement)
        c2 = BodyConstellation(body, t2)
        if c1.symbol != c2.symbol:
            # The body moved from one constellation to another during this time step.
            # Narrow in on the exact moment by doing a binary search.
            tx, cx = FindConstellationChange(body, c1, t1, t2)
            print('{} : {} leaves {} and enters {}.'.format(tx, body.name, c1.name, cx.name))
            c1 = cx
            t1 = tx
        else:
            # No constellation change in this time step. Try again on the next time step.
            c1 = c2
            t1 = t2
    return 0

#------------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) == 1:
        startTime = Time.Now()
    elif len(sys.argv) == 2:
        startTime = Time.Parse(sys.argv[1])
    else:
        print('USAGE: {} [yyyy-mm-ddThh:mm:ssZ]'.format(sys.argv[0]))
        sys.exit(1)
    stopTime = startTime.AddDays(30.0)
    # There are 12 zodiac constellations, and the moon takes
    # about 27.3 days in its sidereal period. Therefore, there
    # are roughly 2.2 days per constellation. We will sample
    # the Moon's constellation once every 0.1 days to reduce
    # the chance of missing a brief transition through a small
    # part of a constellation.
    dayIncrement = 0.1
    rc = FindConstellationChanges(Body.Moon, startTime, stopTime, dayIncrement)
    sys.exit(rc)
