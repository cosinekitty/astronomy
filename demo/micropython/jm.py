#!/usr/bin/env python3
#
#   jupiter_moons.py  -  by Don Cross  -  2021-04-15
#
#   Example Python program for Astronomy Engine:
#   https://github.com/cosinekitty/astronomy
#
#   Calculates the coordinates of Jupiter and its four major moons
#   (Io, Europa, Ganymede, and Callisto) as seen from the Earth
#   at a given date and time. This program illustrates how to correct
#   for the delay caused by the time it takes for light to reach
#   the Earth from the Jupiter system.
#
import sys
from astronomy import Time, JupiterMoons, GeoVector, EquatorFromVector, Body, C_AUDAY, Vector

def PrintBody(name: str, geovec: Vector) -> None:
    # Convert the geocentric vector into equatorial coordinates.
    equ = EquatorFromVector(geovec)
    print('{:<8s}   RA {:10.6f}   DEC {:10.6f}  {:10.6f} AU'.format(name, equ.ra, equ.dec, equ.dist))

def run():
    # If date/time is provided on the command line, use it.
    # Otherwise, use the current date and time.
    if len(sys.argv) == 2:
        time = Time.Parse(sys.argv[1])
    else:
        time = Time.Now()
    time = Time.Parse("2021-04-16T00:26:18.000Z")
    print('Calculations for:', time)

    # Call GeoVector to calculate the geocentric position of Jupiter.
    # GeoVector corrects for light travel time.
    # That means it returns a vector to where Jupiter appears to be
    # in the sky, when the light left Jupiter to travel toward the
    # Earth to arrive here at the specified time. This is different from
    # where Jupiter is at that time.

    jv = GeoVector(Body.Jupiter, time, True)

    # Calculate the amount of time it took light to reach the Earth from Jupiter.
    # The distance to Jupiter (AU) divided by the speed of light (AU/day) = time in days.
    lt_days = jv.Length() / C_AUDAY
    print()
    print('It took light {:0.2f} minutes to reach the Earth from Jupiter.'.format(lt_days * 24.0 * 60.0))
    print()

    # The JupiterMoons function calculates positions of Jupiter's moons without
    # correcting for light travel time. Correct for light travel by backdating
    # by the given amount of light travel time.
    backdate = time.AddDays(-lt_days)

    jm = JupiterMoons(backdate)

    # Tricky: I'm "cheating" a little bit below by adding Vector `jv`
    # to StateVector `jm.<moon>`, to result in a Vector position for each moon.
    # This works because StateVector has all the fields that Vector has,
    # plus the velocity components (vx, vy, vz).
    # This is alarming to type purists, but just another normal day of
    # "duck typing" for Pythonistas.

    PrintBody('Jupiter',    jv)
    PrintBody('Io',         jv + jm.io.Position())
    PrintBody('Europa',     jv + jm.europa.Position())
    PrintBody('Ganymede',   jv + jm.ganymede.Position())
    PrintBody('Callisto',   jv + jm.callisto.Position())
    print()

run()
