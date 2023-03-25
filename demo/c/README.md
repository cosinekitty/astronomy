# Astronomy Engine examples in C

---

### [Camera](camera.c)
Suppose you want to photograph the Moon, and you want to know what it will look like in the photo.
Given a location on the Earth, and a date/time, this program calculates the orientation of the sunlit
side of the Moon with respect to the top of your photo image. It assumes the camera faces directly
toward the Moon's azimuth and tilts upward to its altitude angle above the horizon.
The angles are defined counterclockwise from the zenith, as shown here:

![](https://user-images.githubusercontent.com/11699954/227584171-1135ad6b-2584-4f71-a1e7-5b6ed9ccb87a.png)


### [Culmination](culminate.c)
Finds when the Sun, Moon, and planets reach their highest position in the sky on a given date,
as seen by an observer at a specified location on the Earth.
Culmination is also the moment a body crosses the *meridian*, the imaginary semicircle
in the sky that passes from due north on the horizon, through the zenith (straight up),
and then toward due south on the horizon.

### [Galactic to Horizontal Converter](galactic.c)
A demonstration of how to convert galactic coordinates to horizontal coordinates.
This could be useful for backyard radio astronomers who know the galactic
coordinates of a distant radio source and want to aim a radio dish at it.
Given the galactic coordinates, the geographic coordinates of the observer,
and the date and time of the observation, this program shows how to
obtain the altitude and azimuth to aim the dish at the radio source.

### [Horizon Intersection](horizon.c)
This is a more advanced example. It shows how to use coordinate
transforms to find where the ecliptic intersects with an observer's
horizon at a given date and time.

### [Lunar Eclipse](lunar_eclipse.c)
Calculates details about the first 10 partial/total lunar eclipses
after the given date and time.

### [Moon Phase Calculator](moonphase.c)
This example shows how to determine the Moon's current phase,
and how to predict when the next few quarter phases will occur.

### [Positions](positions.c)
Calculates equatorial and horizontal coordinates of the Sun, Moon, and planets.

### [Rise/Set](riseset.c)
Shows how to calculate sunrise, sunset, moonrise, and moonset times.

### [Seasons](seasons.c)
Calculates the equinoxes and solstices for a given calendar year.

### [Solar Time](solar_time.c)
An example of how to use the Sun's hour angle to calculate true solar time.

### [Triangulate](triangulate.c)
Given the geographic coordinates of two observers, and angular
directions they are looking in, determines geographic coordinates
of the point they are both looking at. This example demonstrates
use of the geoid functions `VectorObserver` and `ObserverVector`
that convert between geographic coordinates and vectors.

### [World Map](worldmap.cpp)
This C++ program generates a PNG image showing intensity of
sunlight and moonlight on a Mercator projection of the Earth.
This example is helpful for showing how to minimize per-pixel
calculations across the globe for a given time of observation.

---

# [API Reference](../../source/c/)
Complete documentation for all the functions and types available
in the C version of Astronomy Engine.
