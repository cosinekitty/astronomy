---

<a name="functions"></a>
## Functions


---

<a name="ApsisKind"></a>
### class ApsisKind




---

<a name="BadVectorError"></a>
### class BadVectorError




---

<a name="Body"></a>
### class Body




---

<a name="Direction"></a>
### class Direction




---

<a name="EarthNotAllowedError"></a>
### class EarthNotAllowedError




---

<a name="Error"></a>
### class Error




---

<a name="IntEnum"></a>
### class IntEnum




---

<a name="InternalError"></a>
### class InternalError




---

<a name="InvalidBodyError"></a>
### class InvalidBodyError




---

<a name="NoConvergeError"></a>
### class NoConvergeError




---

<a name="Observer"></a>
### class Observer

**Represents the geographic location of an observer on the surface of the Earth.**




---

<a name="Refraction"></a>
### class Refraction




---

<a name="Time"></a>
### class Time

**Represents a date and time used for performing astronomy calculations.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | UT1/UTC number of days since noon on January 1, 2000. See the `ut` attribute of this class for more details. |




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




---

<a name="unique"></a>
### unique(enumeration)



