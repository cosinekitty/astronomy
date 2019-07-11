---

<a name="classes"></a>
## Classes


---

<a name="Observer"></a>
### class Observer

**Represents the geographic location of an observer on the surface of the Earth.**




---

<a name="Time"></a>
### class Time

**Represents a date and time used for performing astronomy calculations.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | UT1/UTC number of days since noon on January 1, 2000. See the `ut` attribute of this class for more details. |

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ut` | The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine. Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed. UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed. The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time. Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*. |
| `float` | `tt` | Terrestrial Time days since noon on January 1, 2000. Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation. The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth. Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET). |



---

<a name="functions"></a>
## Functions


---

<a name="BodyCode"></a>
### BodyCode(name)

**Finds the Body enumeration value, given the name of a body.**

| Type | Parameter | Description |
| --- | --- | --- |
| `str` | `name` | The common English name of a supported celestial body. |




---

<a name="GeoMoon"></a>
### GeoMoon(time)

**Calculates the geocentric position of the Moon at a given time.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Moon's position. |



