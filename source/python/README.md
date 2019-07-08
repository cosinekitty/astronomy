# astronomy
Astronomy Engine by Don Cross

See the GitHub project page for full documentation, examples,
and other information:

https://github.com/cosinekitty/astronomy


## BodyCode
```python
BodyCode(name)
```
Finds the integer body code given the name of a body.

Parameters
----------
name: str
    The common English name of a supported celestial body.

Returns
-------
int
    If `name` is a valid body name, returns the integer value
    of the body code associated with that body.
    Otherwise, returns `BODY_INVALID`.

Example
-------

>>> astronomy.BodyCode('Mars')
3


## Time
```python
Time(self, ut)
```
Represents a date and time used for performing astronomy calculations.

All calculations performed by Astronomy Engine are based on
dates and times represented by `Time` objects.

Parameters
----------
ut : float
    UT1/UTC number of days since noon on January 1, 2000.
    See the `ut` attribute of this class for more details.

Attributes
----------
ut : float
    The floating point number of days of Universal Time since noon UTC January 1, 2000.
    Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
    not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds.
    This approximation is sufficient for the accuracy requirements of Astronomy Engine.
    Universal Time Coordinate (UTC) is the international standard for legal and civil
    timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
    UTC is kept in sync with unpredictable observed changes in the Earth's rotation
    by occasionally adding leap seconds as needed.
    UT1 is an idealized time scale based on observed rotation of the Earth, which
    gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
    large scale weather events like hurricanes, and internal seismic and convection effects.
    Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
    is adjusted by a scheduled whole number of leap seconds as needed.
    The value in `ut` is appropriate for any calculation involving the Earth's rotation,
    such as calculating rise/set times, culumination, and anything involving apparent
    sidereal time.
    Before the era of atomic timekeeping, days based on the Earth's rotation
    were often known as <i>mean solar days</i>.
tt : float
    Terrestrial Time days since noon on January 1, 2000.
    Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
    In this system, days are not based on Earth rotations, but instead by
    the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
    divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
    for changes in the Earth's rotation.
    The value in `tt` is used for calculations of movements not involving the Earth's rotation,
    such as the orbits of planets around the Sun, or the Moon around the Earth.
    Historically, Terrestrial Time has also been known by the term <i>Ephemeris Time</i> (ET).

### Make
```python
Time.Make(year, month, day, hour, minute, second)
```
Creates a `Time` object from a UTC calendar date and time.

Parameters
----------
year : int
    The UTC 4-digit year value, e.g. 2019.
month : int
    The UTC month in the range 1..12.
day : int
    The UTC day of the month, in the range 1..31.
hour : int
    The UTC hour, in the range 0..23.
minute : int
    The UTC minute, in the range 0..59.
second : float
    The real-valued UTC second, in the range [0, 60).

Returns
-------
Time

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

Parameters
----------
time : Time
    The date and time for which to calculate the Moon's position.

Returns
-------
Vector
    The Moon's position as a vector in J2000 Cartesian equatorial coordinates.


