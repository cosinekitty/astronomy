#!/usr/bin/env python3
#
#   lunar_angles.py  -  Don Cross  -  2021-04-27
#
#   Searches for the next few times the Moon reaches a relative
#   ecliptic longitude with respect to another body
#   (as seen from the Earth) that is a multiple of 30 degrees.
#   This is an example of creating a custom search algorithm.
#
import sys
from astronomy import Time, Body, PairLongitude

#------------------------------------------------------------------------------

class Event:
    def __init__(self, body, time, angle):
        self.body = body
        self.time = time
        self.angle = angle

    def __lt__(self, other):
        # This makes lists of Event objects sortable chronologically.
        return self.time < other.time

    def __str__(self):
        return '{} {:<8s} {:3.0f}'.format(self.time, self.body.name, self.angle)

#------------------------------------------------------------------------------

def AdjustAngle(angle):
    # Force the angle into the half-open range (-180.0, +180.0]
    while angle <= -180.0:
        angle += 360.0
    while angle > +180.0:
        angle -= 360.0
    return angle

def Straddle(lon1, lon2, angle):
    # A pair of longitudes "straddles" an angle if
    # the angle lies between the two longitudes modulo 360 degrees.
    a1 = AdjustAngle(lon1 - angle)
    a2 = AdjustAngle(lon2 - angle)
    return a1 <= 0.0 <= a2

def AppendEvents(event_list, body, startTime, stopTime, dayIncrement):
    angle_list = [30.0*i for i in range(12)]
    t1 = startTime
    while t1 < stopTime:
        t2 = t1.AddDays(dayIncrement)
        # Calculate the relative longitude at t1 and t2.
        lon1 = PairLongitude(Body.Moon, body, t1)
        lon2 = PairLongitude(Body.Moon, body, t2)
        # Does it straddle a 30-degree boundary?
        for angle in angle_list:
            if Straddle(lon1, lon2, angle):
                # !!! Search for when the pair reaches that exact longitude.
                # !!! Hack for now: take the mean time.
                t = t1.AddDays(dayIncrement / 2.0)
                event_list.append(Event(body, t, angle))
        t1 = t2

#------------------------------------------------------------------------------

def PrintMoonChart(startTime, stopTime, dayIncrement):
    # Make a list of all other bodies with which to compare the Moon.
    otherBodyList = [Body.Sun, Body.Mercury, Body.Venus, Body.Mars, Body.Jupiter, Body.Saturn]

    # Make a list of events for each body in that list.
    event_list = []
    for body in otherBodyList:
        AppendEvents(event_list, body, startTime, stopTime, dayIncrement)

    # Sort the list chronologically.
    event_list.sort()

    # Print the list.
    for evt in event_list:
        print(evt)

    # Success!
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
    rc = PrintMoonChart(startTime, stopTime, 1.0)
    sys.exit(rc)
