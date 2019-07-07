# astronomy

## Time
```python
Time(self, ut)
```
Represents a date and time used for performing astronomy calculations.

### Make
```python
Time.Make(year, month, day, hour, minute, second)
```
Creates a `Time` object from a UTC calendar date and time.

:param year: The UTC 4-digit year value, e.g. 2019.
:param month: The UTC month in the range 1..12.
:param day: The UTC day of the month, in the range 1..31.
:param hour: The UTC hour, in the range 0..23.
:param minute: The UTC minute, in the range 0..59.
:param second: The real-valued UTC second, in the range [0, 60).
:rtype: `Time`

## Observer
```python
Observer(self, latitude, longitude, height=0)
```
Represents the geographic location of an observer on the surface of the Earth.

:param latitude: Geographic latitude in degrees north of the equator.
:param longitude: Geographic longitude in degrees east of the prime meridian at Greenwich, England.
:param height: Elevation above sea level in meters.

## GeoMoon
```python
GeoMoon(time)
```
Calculates the geocentric position of the Moon at a given time.

Given a time of observation, calculates the Moon's position as a vector.
The vector gives the location of the Moon's center relative to the Earth's center
with x-, y-, and z-components measured in astronomical units.

This algorithm is based on Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
It is adapted from Turbo Pascal code from the book
[Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
by Montenbruck and Pfleger.

:param time:  The date and time for which to calculate the Moon's position.
:return The Moon's position as a vector in J2000 Cartesian equatorial coordinates.
:rtype Vector

