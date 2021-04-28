#!/usr/bin/env python3
#
#   lunar_angles.py  -  Don Cross  -  2021-04-28
#
#   Searches for the next few times the Moon reaches a relative
#   ecliptic longitude with respect to another body
#   (as seen from the Earth) that is a multiple of 30 degrees.
#
#   This is an example of how to implement your own custom search to
#   find the time when a function ascends from a negative value through
#   zero to a positive value. The search returns the time of the zero-crossing
#   within the specified tolerance in seconds. You must already have narrowed
#   in on a bounding time range [t1, t2] that contains exactly one
#   zero-crossing, or the search will fail.
#
#   You provide a custom context (see example PairSearchContext below)
#   and a function that returns a floating point number (see LongitudeFunc below).
#   Your custom context is passed to every single call to your custom function,
#   and holds whatever constant state your function needs to calculate its return value.
#   (Some functions may not need a context; can pass None for the context in that case.)
#   Your function must accept your custom context type and a Time
#   value, and return a floating point value. Search() assumes the value will
#   increase from negative value to a positive value over the supplied time range [t1, t2].
#   If there is no zero-crossing, or the zero-crossing is not unique,
#   or the zero-crossing descends through zero instead of ascending,
#   Search() will return None, indicating a search failure.
#
#   See the following online documentation for more information about Search():
#
#   https://github.com/cosinekitty/astronomy/tree/master/source/python#searchfunc-context-t1-t2-dt_tolerance_seconds
#
import sys
from astronomy import Time, Body, PairLongitude, Search

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

#------------------------------------------------------------------------------

class PairSearchContext:
    def __init__(self, body1, body2, targetRelLon):
        self.body1 = body1
        self.body2 = body2
        self.targetRelLon = targetRelLon


def LongitudeFunc(context, time):
    lon = PairLongitude(context.body1, context.body2, time)
    return AdjustAngle(lon - context.targetRelLon)


def LongitudeSearch(body1, body2, angle, t1, t2):
    context = PairSearchContext(body1, body2, angle)
    tolerance_seconds = 0.1
    t = Search(LongitudeFunc, context, t1, t2, tolerance_seconds)
    if not t:
        raise Exception('Search failure for body={}, t1={}, t2={}'.format(body, t1, t2))
    return t

#------------------------------------------------------------------------------

def Straddle(lon1, lon2, angle):
    # A pair of longitudes "straddles" an angle if
    # the angle lies between the two longitudes modulo 360 degrees.
    a1 = AdjustAngle(lon1 - angle)
    a2 = AdjustAngle(lon2 - angle)
    return a1 <= 0.0 <= a2

#------------------------------------------------------------------------------

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
                # Search for when the pair reaches that exact longitude.
                t = LongitudeSearch(Body.Moon, body, angle, t1, t2)
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
