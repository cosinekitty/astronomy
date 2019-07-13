---

<a name="classes"></a>
## Classes


---

<a name="Observer"></a>
### class Observer

**Represents the geographic location of an observer on the surface of the Earth.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `latitude` | Geographic latitude in degrees north of the equator. |
| `float` | `longitude` | Geographic longitude in degrees east of the prime meridian at Greenwich, England. |
| `float` | `height` | Elevation above sea level in meters. |




---

<a name="Time"></a>
### class Time

**Represents a date and time used for performing astronomy calculations.**

All calculations performed by Astronomy Engine are based on
dates and times represented by `Time` objects.


| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | UT1/UTC number of days since noon on January 1, 2000. See the `ut` attribute of this class for more details. |

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ut` | The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine. Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed. UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed. The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time. Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*. |
| `float` | `tt` | Terrestrial Time days since noon on January 1, 2000. Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation. The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth. Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET). |



#### member functions


<a name="Time.Make"></a>
### Time.Make(year, month, day, hour, minute, second)

**Creates a #Time object from a UTC calendar date and time.**

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The UTC 4-digit year value, e.g. 2019. |
| `int` | `month` | The UTC month in the range 1..12. |
| `int` | `day` | The UTC day of the month, in the range 1..31. |
| `int` | `hour` | The UTC hour, in the range 0..23. |
| `int` | `minute` | The UTC minute, in the range 0..59. |
| `float` | `second` | The real-valued UTC second, in the range [0, 60). |




---

<a name="Vector"></a>
### class Vector

**A Cartesian vector with 3 space coordinates and 1 time coordinate.**

The vector's space coordinates are measured in astronomical units (AU).
The coordinate system varies and depends on context.
The vector also includes a time stamp.
Properties
x : float
    The x-coordinate of the vector, measured in AU.
y : float
    The y-coordinate of the vector, measured in AU.
z : float
    The z-coordinate of the vector, measured in AU.
t : Time
    The date and time at which the coordinate is valid.




#### member functions


<a name="Vector.Length"></a>
### Vector.Length(self)

Returns the length of the vector in AU.




---

<a name="enumerations"></a>
## Enumerated Types


---

<a name="ApsisKind"></a>
### enum ApsisKind

**Represents whether a satellite is at a closest or farthest point in its orbit.**

An apsis is a point in a satellite's orbit that is closest to,
or farthest from, the body it orbits (its primary).
`ApsisKind` is an enumerated type that indicates which of these
two cases applies to a particular apsis event.


| Value | Description |
| --- | --- |
| `Pericenter` | The satellite is at its closest point to its primary. |
| `Apocenter` | The satellite is at its farthest point from its primary. |
| `Invalid` | A placeholder for an undefined, unknown, or invalid apsis. |



---

<a name="Body"></a>
### enum Body

**The celestial bodies supported by Astronomy Engine calculations.**

| Value | Description |
| --- | --- |
| `Invalid` | An unknown, invalid, or undefined celestial body. |
| `Mercury` | The planet Mercury. |
| `Venus` | The planet Venus. |
| `Earth` | The planet Earth. |
| `Mars` | The planet Mars. |
| `Jupiter` | The planet Jupiter. |
| `Saturn` | The planet Saturn. |
| `Uranus` | The planet Uranus. |
| `Neptune` | The planet Neptune. |
| `Pluto` | The planet Pluto. |
| `Sun` | The Sun. |
| `Moon` | The Earth's moon. |



---

<a name="Direction"></a>
### enum Direction

**Indicates whether a body is rising above or setting below the horizon.**

Specifies the direction of a rising or setting event for a body.
For example, `Direction.Rise` is used to find sunrise times,
and `Direction.Set` is used to find sunset times.


| Value | Description |
| --- | --- |
| `Rise` | First appearance of a body as it rises above the horizon. |
| `Set` | Last appearance of a body as it sinks below the horizon. |



---

<a name="Refraction"></a>
### enum Refraction

**Selects if/how to correct for atmospheric refraction.**

Some functions allow enabling or disabling atmospheric refraction
for the calculated apparent position of a celestial body
as seen by an observer on the surface of the Earth.


| Value | Description |
| --- | --- |
| `Airless` | No atmospheric refraction correction. |
| `Normal` | Recommended correction for standard atmospheric refraction. |
| `JplHorizons` | Used only for compatibility testing with JPL Horizons online tool. |


---

<a name="errors"></a>
## Error Types


---

<a name="BadVectorError"></a>
### BadVectorError

A vector magnitude is too small to have a direction in space.





---

<a name="EarthNotAllowedError"></a>
### EarthNotAllowedError

The Earth is not allowed as the celestial body in this calculation.





---

<a name="Error"></a>
### Error

Indicates an error in an astronomical calculation.





---

<a name="InternalError"></a>
### InternalError

**An internal error occured that should be reported as a bug.**

Indicates an unexpected and unrecoverable condition occurred.
If you encounter this error using Astronomy Engine, it would be very
helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
page on GitHub. Please include a copy of the stack trace, along with a description
of how to reproduce the error. This will help improve the quality of
Astronomy Engine for everyone! (Thank you in advance from the author.)





---

<a name="InvalidBodyError"></a>
### InvalidBodyError

The celestial body is not allowed for this calculation.





---

<a name="NoConvergeError"></a>
### NoConvergeError

**A numeric solver did not converge.**

Indicates that there was a failure of a numeric solver to converge.
If you encounter this error using Astronomy Engine, it would be very
helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
page on GitHub. Please include a copy of the stack trace, along with a description
of how to reproduce the error. This will help improve the quality of
Astronomy Engine for everyone! (Thank you in advance from the author.)




---

<a name="functions"></a>
## Functions


---

<a name="BodyCode"></a>
### BodyCode(name)

**Finds the Body enumeration value, given the name of a body.**

```
>>> astronomy.BodyCode('Mars')
<Body.Mars: 3>
```


| Type | Parameter | Description |
| --- | --- | --- |
| `str` | `name` | The common English name of a supported celestial body. |




---

<a name="GeoMoon"></a>
### GeoMoon(time)

**Calculates the geocentric position of the Moon at a given time.**

Given a time of observation, calculates the Moon's position as a vector.
The vector gives the location of the Moon's center relative to the Earth's center
with x-, y-, and z-components measured in astronomical units.
This algorithm is based on Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
It is adapted from Turbo Pascal code from the book
[Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
by Montenbruck and Pfleger.


| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's position. |



