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

### [Culmination](culminate.py)
Finds when the Sun, Moon, and planets reach their highest position in the sky on a given date,
as seen by an observer at a specified location on the Earth.
Culmination is also the moment a body crosses the *meridian*, the imaginary semicircle
in the sky that passes from due north on the horizon, through the zenith (straight up),
and then toward due south on the horizon.

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

---

# [API Reference](../../source/python/)
Complete documentation for all the functions and types available
in the Python version of Astronomy Engine.
