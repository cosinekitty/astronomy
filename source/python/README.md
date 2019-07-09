---

<a name="functions"></a>
## Functions


---

<a name="BodyCode"></a>
### BodyCode(name) -> astronomy.Body
Finds the Body enumeration value, given the name of a body.

Parameters
----------
name: str
    The common English name of a supported celestial body.

Returns
-------
Body
    If `name` is a valid body name, returns the enumeration
    value associated with that body.
    Otherwise, returns `Body.Invalid`.

Example
-------

>>> astronomy.BodyCode('Mars')
<Body.Mars: 3>


---

<a name="GeoMoon"></a>
### GeoMoon(time: astronomy.Time) -> astronomy.Vector
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


---

<a name="unique"></a>
### unique(enumeration)
Class decorator for enumerations ensuring unique member values.

