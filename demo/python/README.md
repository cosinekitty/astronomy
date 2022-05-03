# Astronomy Engine examples in Python

---

### [Camera](camera.py)
Suppose you want to photograph the Moon,
and you want to know what it will look like in the photo.
Given a location on the Earth, and a date/time,
this program calculates the orientation of the sunlit
side of the Moon with respect to the top of your
photo image. It assumes the camera faces directly
toward the Moon's azimuth and tilts upward to its
altitude angle above the horizon.

### [Constellation](constellation.py)
This demo finds what constellation the Moon
is in at a given time. It also shows how to do a binary
search to find the moment in time when the Moon moves
across the border between constellations.

### [Culmination](culminate.py)
Finds when the Sun, Moon, and planets reach their highest position in the sky on a given date,
as seen by an observer at a specified location on the Earth.
Culmination is also the moment a body crosses the *meridian*, the imaginary semicircle
in the sky that passes from due north on the horizon, through the zenith (straight up),
and then toward due south on the horizon.

### [Galactic to Horizontal Converter](galactic.py)
A demonstration of how to convert galactic coordinates to horizontal coordinates.
This could be useful for backyard radio astronomers who know the galactic
coordinates of a distant radio source and want to aim a radio dish at it.
Given the galactic coordinates, the geographic coordinates of the observer,
and the date and time of the observation, this program shows how to
obtain the altitude and azimuth to aim the dish at the radio source.

### [Horizon Intersection](horizon.py)
This is a more advanced example. It shows how to use coordinate
transforms to find where the ecliptic intersects with an observer's
horizon at a given date and time.

### [Jupiter's Moons](jupiter_moons.py)
Calculates the coordinates of Jupiter and its four major moons
(Io, Europa, Ganymede, and Callisto) as seen from the Earth
at a given date and time. This program illustrates how to correct
for the delay caused by the time it takes for light to reach
the Earth from the Jupiter system.

### [Lunar Angles](lunar_angles.py)
This is an example of how to implement your own custom search function
using Astronomy Engine. This program searches for the next few times
the Moon reaches a relative ecliptic longitude with respect to another body
(as seen from the Earth) that is a multiple of 30 degrees.

### [Lunar Eclipse](lunar_eclipse.py)
Calculates details about the first 10 partial/total lunar eclipses
after the given date and time.

### [Moon Phase Calculator](moonphase.py)
This example shows how to determine the Moon's current phase,
and how to predict when the next few quarter phases will occur.

### [Positions](positions.py)
Calculates equatorial and horizontal coordinates of the Sun, Moon, and planets.

### [Rise/Set](riseset.py)
Shows how to calculate sunrise, sunset, moonrise, and moonset times.

### [Seasons](seasons.py)
Calculates the equinoxes and solstices for a given calendar year.

### [Triangulate](triangulate.py)
Given the geographic coordinates of two observers, and angular
directions they are looking in, determines geographic coordinates
of the point they are both looking at. This example demonstrates
use of the geoid functions `VectorObserver` and `ObserverVector`
that convert between geographic coordinates and vectors.

---

# [API Reference](../../source/python/)
Complete documentation for all the functions and types available
in the Python version of Astronomy Engine.
