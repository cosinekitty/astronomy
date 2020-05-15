# Astronomy Engine (Python)

This is the complete programming reference for the Python version of
Astronomy Engine. Supports Python 3. Does NOT support Python 2.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start
To get started quickly, here are some [examples](../../demo/python/).

---

## Contents

- [Topic Index](#topics)
- [Classes](#classes)
- [Enumerated Types](#enumerations)
- [Error Types](#errors)
- [Functions](#functions)

---

<a name="topics"></a>
## Topic Index

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#HelioVector) | Calculates vector with respect to the center of the Sun. |
| [GeoVector](#GeoVector)     | Calculates vector with respect to the center of the Earth. |
| [Equator](#Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Ecliptic)       | Converts J2000 equatorial coordinates to J2000 ecliptic coordinates. |
| [EclipticLongitude](#EclipticLongitude) | Calculates ecliptic longitude of a body in the J2000 system. |
| [Horizon](#Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [LongitudeFromSun](#LongitudeFromSun) | Calculates a body's apparent ecliptic longitude difference from the Sun, as seen by an observer on the Earth. |

### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchHourAngle](#SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Planet perihelion and aphelion

| Function | Description |
| -------- | ----------- |
| [SearchPlanetApsis](#SearchPlanetApsis) | Finds the next perihelion or aphelion of a planet after a specified date. |
| [NextPlanetApsis](#NextPlanetApsis) | Given an already-found apsis, finds the next perihelion or aphelion of a planet. |

### Visual magnitude and elongation

| Function | Description |
| -------- | ----------- |
| [Illumination](#Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| Function | Description |
| -------- | ----------- |
| [SearchRelativeLongitude](#SearchRelativeLongitude) | Finds oppositions and conjunctions of planets. |

### Equinoxes, solstices, and apparent solar motion

| Function | Description |
| -------- | ----------- |
| [SearchSunLongitude](#SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the *true equator of date* system. |
| [Seasons](#Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#SunPosition) | Calculates the Sun's apparent ecliptic coordinates as seen from the Earth. |

### Coordinate transforms

The following four orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = Equatorial J2000**: Uses the Earth's equator on January 1, 2000, at noon UTC.
- **EQD = Equator of-date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECL = Ecliptic**: Uses the mean plane of the Earth's orbit around the Sun. The x-axis is referenced against the J2000 equinox.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.

| Function | Description |
| -------- | ----------- |
| [RotateVector](#RotateVector) | Applies a rotation matrix to a vector, yielding a vector in another orientation system. |
| [InverseRotation](#InverseRotation) | Given a rotation matrix, finds the inverse rotation matrix that does the opposite transformation. |
| [CombineRotation](#CombineRotation) | Given two rotation matrices, returns a rotation matrix that combines them into a net transformation. |
| [VectorFromSphere](#VectorFromSphere) | Converts spherical coordinates to Cartesian coordinates. |
| [SphereFromVector](#SphereFromVector) | Converts Cartesian coordinates to spherical coordinates. |
| [VectorFromEquator](#VectorFromEquator) | Given angular equatorial coordinates, calculates equatorial vector. |
| [EquatorFromVector](#EquatorFromVector) | Given an equatorial vector, calculates equatorial angular coordinates. |
| [VectorFromHorizon](#VectorFromHorizon) | Given apparent angular horizontal coordinates, calculates horizontal vector. |
| [HorizonFromVector](#HorizonFromVector) | Given a vector in horizontal orientation, calculates horizontal angular coordinates. |
| [Rotation_EQD_EQJ](#Rotation_EQD_EQJ) | Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ). |
| [Rotation_EQD_ECL](#Rotation_EQD_ECL) | Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL). |
| [Rotation_EQD_HOR](#Rotation_EQD_HOR) | Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR). |
| [Rotation_EQJ_EQD](#Rotation_EQJ_EQD) | Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD). |
| [Rotation_EQJ_ECL](#Rotation_EQJ_ECL) | Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL). |
| [Rotation_EQJ_HOR](#Rotation_EQJ_HOR) | Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR). |
| [Rotation_ECL_EQD](#Rotation_ECL_EQD) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD). |
| [Rotation_ECL_EQJ](#Rotation_ECL_EQJ) | Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ). |
| [Rotation_ECL_HOR](#Rotation_ECL_HOR) | Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR). |
| [Rotation_HOR_EQD](#Rotation_HOR_EQD) | Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD). |
| [Rotation_HOR_EQJ](#Rotation_HOR_EQJ) | Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ). |
| [Rotation_HOR_ECL](#Rotation_HOR_ECL) | Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL). |

---

---

<a name="classes"></a>
## Classes

---

<a name="Apsis"></a>
### class Apsis

**An event where a satellite is closest to or farthest from the body it orbits.**

For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
event where the orbiting body reaches its closest or farthest point from the primary body.
The closest approach is called *pericenter* and the farthest point is *apocenter*.
More specific terminology is common for particular orbiting bodies.
The Moon's closest approach to the Earth is called *perigee* and its furthest
point is called *apogee*. The closest approach of a planet to the Sun is called
*perihelion* and the furthest point is called *aphelion*.
This data structure is returned by [`SearchLunarApsis`](#SearchLunarApsis) and [`NextLunarApsis`](#NextLunarApsis)
to iterate through consecutive alternating perigees and apogees.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the apsis. |
| [`ApsisKind`](#ApsisKind) | `kind` | Whether this is a pericenter or apocenter event. |
| `float` | `dist_au` | The distance between the centers of the bodies in astronomical units. |
| `float` | `dist_km` | The distance between the centers of the bodies in kilometers. |

---

<a name="ConstellationInfo"></a>
### class ConstellationInfo

**Reports the constellation that a given celestial point lies within.**

The [`Constellation`](#Constellation) function returns this struct
to report which constellation corresponds with a given point in the sky.
Constellations are defined with respect to the B1875 equatorial system
per IAU standard. Although `Constellation` requires J2000 equatorial
coordinates, the struct contains converted B1875 coordinates for reference.

| Type | Attribute | Description |
| --- | --- | --- |
| `string` | `symbol` | 3-character mnemonic symbol for the constellation, e.g. "Ori". |
| `string` | `name` | Full name of constellation, e.g. "Orion". |
| `float` | `ra1875` | Right ascension expressed in B1875 coordinates. |
| `float` | `dec1875` | Declination expressed in B1875 coordinates. |

---

<a name="EclipticCoordinates"></a>
### class EclipticCoordinates

**Ecliptic angular and Cartesian coordinates.**

Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ex` | Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane. |
| `float` | `ey` | Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox. |
| `float` | `ez` | Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north. |
| `float` | `elat` | Latitude in degrees north (positive) or south (negative) of the ecliptic plane. |
| `float` | `elon` | Longitude in degrees around the ecliptic plane prograde from the equinox. |

---

<a name="ElongationEvent"></a>
### class ElongationEvent

**Contains information about the visibility of a celestial body at a given date and time.**

See the [`Elongation`](#Elongation) function for more detailed information about the members of this class.
See also [`SearchMaxElongation`](#SearchMaxElongation) for how to search for maximum elongation events.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| [`Visibility`](#Visibility) | `visibility` | Whether the body is best seen in the morning or the evening. |
| `float` | `elongation` | The angle in degrees between the body and the Sun, as seen from the Earth. |
| `float` | `ecliptic_separation` | The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth. |

---

<a name="Equatorial"></a>
### class Equatorial

**Equatorial angular coordinates**

Coordinates of a celestial body as seen from the Earth.
Can be geocentric or topocentric, depending on context.
The coordinates are oriented with respect to the Earth's
equator projected onto the sky.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `ra` | Right ascension in sidereal hours. |
| `float` | `dec` | Declination in degrees. |
| `float` | `dist` | Distance to the celestial body in AU. |

---

<a name="HorizontalCoordinates"></a>
### class HorizontalCoordinates

**Coordinates of a celestial body as seen by a topocentric observer.**

Contains horizontal and equatorial coordinates as seen by an observer
on or near the surface of the Earth (a topocentric observer).
All coordinates are optionally corrected for atmospheric refraction.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `azimuth` | The compass direction laterally around the observer's horizon, measured in degrees. North is 0 degrees, east is 90 degrees, south is 180 degrees, etc. |
| `float` | `altitude` | The angle in degrees above (positive) or below (negative) the observer's horizon. |
| `float` | `ra` | The right ascension in sidereal hours. |
| `float` | `dec` | The declination in degrees. |

---

<a name="HourAngleEvent"></a>
### class HourAngleEvent

**Information about a celestial body crossing a specific hour angle.**

Returned by the function [`SearchHourAngle`](#SearchHourAngle) to report information about
a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time when the body crosses the specified hour angle. |
| [`HorizontalCoordinates`](#HorizontalCoordinates) | `hor` | Apparent coordinates of the body at the time it crosses the specified hour angle. |

---

<a name="IlluminationInfo"></a>
### class IlluminationInfo

**Information about the brightness and illuminated shape of a celestial body.**

Returned by functions [`Illumination`](#Illumination) and [`SearchPeakMagnitude`](#SearchPeakMagnitude)
to report the visual magnitude and illuminated fraction of a celestial
body at a given date and time.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| `float` | `mag` | The visual magnitude of the body. Smaller values are brighter. |
| `float` | `phase_angle` | The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. |
| `float` | `phase_fraction` | A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth. |
| `float` | `helio_dist` | The distance between the Sun and the body at the observation time, in AU. |
| `float` | `ring_tilt` | For Saturn, the tilt angle in degrees of its rings as seen from Earth. When the `ring_tilt` is very close to 0, it means the rings are edge-on as seen from observers on the Earth, and are thus very difficult to see. For bodies other than Saturn, `ring_tilt` is `None`. |

---

<a name="MoonQuarter"></a>
### class MoonQuarter

**A lunar quarter event along with its date and time.**

An object of this type represents one of the four major
lunar phases that appear on calendars:
new moon, first quarter, full moon, or third quarter.
Along with the `quarter` attribute that specifies the
type of quarter, it contains a `time` field that indicates
when the lunar quarter event happens.

| Type | Attribute | Description |
| --- | --- | --- |
| `int` | `quarter` | 0=new moon, 1=first quarter, 2=full moon, 3=third quarter. |
| [`Time`](#Time) | `time` | The date and time of the lunar quarter. |

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

<a name="RotationMatrix"></a>
### class RotationMatrix

Contains a rotation matrix that can be used to transform one
coordinate system into another.

| Type | Parameter | Description |
| --- | --- | --- |
| `float[3][3]` | `rot` | A normalized 3x3 rotation matrix. |

---

<a name="SeasonInfo"></a>
### class SeasonInfo

**The dates and times of changes of season for a given calendar year.**

Call [`Seasons`](#Seasons) to calculate this data structure for a given year.

| Type | Attribute | Description |
| --- | --- | --- |
| [`Time`](#Time) | `mar_equinox` | The date and time of the March equinox for the specified year. |
| [`Time`](#Time) | `jun_solstice` | The date and time of the June solstice for the specified year. |
| [`Time`](#Time) | `sep_equinox` | The date and time of the September equinox for the specified year. |
| [`Time`](#Time) | `dec_solstice` | The date and time of the December solstice for the specified year. |

---

<a name="Spherical"></a>
### class Spherical

**Holds spherical coordinates: latitude, longitude, distance.**

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `lat` | The latitude angle: -90..+90 degrees. |
| `float` | `lon` | The longitude angle: 0..360 degrees. |
| `float` | `dist` | Distance in AU. |

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

<a name="Time.AddDays"></a>
### Time.AddDays(self, days)

**Calculates the sum or difference of a [`Time`](#Time) with a specified real-valued number of days.**

Sometimes we need to adjust a given [`Time`](#Time) value by a certain amount of time.
This function adds the given real number of days in `days` to the date and time
in the calling object.
More precisely, the result's Universal Time field `ut` is exactly adjusted by `days`
and the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time,
according to the historical and predictive Delta-T model provided by the
[United States Naval Observatory](http://maia.usno.navy.mil/ser7/).
The value of the calling object is not modified. This function creates a brand new
[`Time`](#Time) object and returns it.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `days` | A floating point number of days by which to adjust `time`. May be negative, 0, or positive. |

### Returns: [`Time`](#Time)

<a name="Time.Make"></a>
### Time.Make(year, month, day, hour, minute, second)

**Creates a [`Time`](#Time) object from a UTC calendar date and time.**

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The UTC 4-digit year value, e.g. 2019. |
| `int` | `month` | The UTC month in the range 1..12. |
| `int` | `day` | The UTC day of the month, in the range 1..31. |
| `int` | `hour` | The UTC hour, in the range 0..23. |
| `int` | `minute` | The UTC minute, in the range 0..59. |
| `float` | `second` | The real-valued UTC second, in the range [0, 60). |

### Returns: [`Time`](#Time)

<a name="Time.Now"></a>
### Time.Now()

**Returns the computer's current date and time in the form of a [`Time`](#Time) object.**

Uses the computer's system clock to find the current UTC date and time.
Converts that date and time to a [`Time`](#Time) value and returns the result.
Callers can pass this value to other Astronomy Engine functions to
calculate current observational conditions.

### Returns: [`Time`](#Time)

<a name="Time.Parse"></a>
### Time.Parse(text)

**Creates a [`Time`](#Time) object from a string of the form 'yyyy-mm-ddThh:mm:ss.sssZ'**

Parses a UTC date and time from a string and returns a [`Time`](#Time) object.
Permits a subset of ISO 8601 format.
The year, month, and day are required.
Hours, minutes, seconds, and fractions of a second are optional.
If time is specified, there must be a 'T' between the date and the time
and a 'Z' at the end of the time.

| Type | Parameter | Description |
| --- | --- | --- |
| `string` | `text` | A string of the following formats: `yyyy-mm-dd` `yyyy-mm-ddThh:mmZ` `yyyy-mm-ddThh:mm:ssZ` `yyyy-mm-ddThh:mm:ss.sssZ` |

### Returns: [`Time`](#Time)

<a name="Time.Utc"></a>
### Time.Utc(self)

**Returns the UTC date and time as a `datetime` object.**

Uses the standard [`datetime`](https://docs.python.org/3/library/datetime.html) class
to represent the date and time in this Time object.

### Returns: `datetime`

---

<a name="Vector"></a>
### class Vector

**A Cartesian vector with 3 space coordinates and 1 time coordinate.**

The vector's space coordinates are measured in astronomical units (AU).
The coordinate system varies and depends on context.
The vector also includes a time stamp.

| Type | Attribute | Description |
| --- | --- | --- |
| `float` | `x` | The x-coordinate of the vector, measured in AU. |
| `float` | `y` | The y-coordinate of the vector, measured in AU. |
| `float` | `z` | The z-coordinate of the vector, measured in AU. |
| [`Time`](#Time) | `t` | The date and time at which the coordinate is valid. |

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
| `EMB` | The Earth/Moon Barycenter. |
| `SSB` | The Solar System Barycenter. |

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

<a name="Visibility"></a>
### enum Visibility

**Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.**

| Value | Description |
| --- | --- |
| `Morning` | The body is best visible in the morning, before sunrise. |
| `Evening` | The body is best visible in the evening, after sunset. |

---

<a name="errors"></a>
## Error Types

---

<a name="BadTimeError"></a>
### BadTimeError

Cannot calculate Pluto position for this date/time.

---

<a name="BadVectorError"></a>
### BadVectorError

A vector magnitude is too small to have a direction in space.

---

<a name="DateTimeFormatError"></a>
### DateTimeFormatError

The syntax of a UTC date/time string was not valid, or it contains invalid values.

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

<a name="AngleFromSun"></a>
### AngleFromSun(body, time)

**Returns the angle between the given body and the Sun, as seen from the Earth.**

This function calculates the angular separation between the given body and the Sun,
as seen from the center of the Earth. This angle is helpful for determining how
easy it is to see the body away from the glare of the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose angle from the Sun is to be measured. Not allowed to be `Body.Earth`. |
| [`Time`](#Time) | `time` | The time at which the observation is made. |

### Returns: `float`
A numeric value indicating the angle in degrees between the Sun
and the specified body as seen from the center of the Earth.

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

### Returns: [`Body`](#Body)
If `name` is a valid body name, returns the enumeration
value associated with that body.
Otherwise, returns `Body.Invalid`.

---

<a name="CombineRotation"></a>
### CombineRotation(a, b)

**Creates a rotation based on applying one rotation followed by another.**

Given two rotation matrices, returns a combined rotation matrix that is
equivalent to rotating based on the first matrix, followed by the second.
b : RotationMatrix
    The second rotation to apply.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `a` | The first rotation to apply. |

### Returns: [`RotationMatrix`](#RotationMatrix)
The combined rotation matrix.

---

<a name="Constellation"></a>
### Constellation(ra, dec)

**Determines the constellation that contains the given point in the sky.**

Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
constellation that contains that point.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ra` | The right ascension (RA) of a point in the sky, using the J2000 equatorial system. |
| `float` | `dec` | The declination (DEC) of a point in the sky, using the J2000 equatorial system. |

### Returns: [`ConstellationInfo`](#ConstellationInfo)
A structure that contains the 3-letter abbreviation and full name
of the constellation that contains the given (ra,dec), along with
the converted B1875 (ra,dec) for that point.

---

<a name="DeltaT_EspenakMeeus"></a>
### DeltaT_EspenakMeeus(ut)

**The default Delta T function used by Astronomy Engine.**

Espenak and Meeus use a series of piecewise polynomials to
approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
This is the default Delta T function used by Astronomy Engine.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `ut` | The floating point number of days since noon UTC on January 1, 2000. |

### Returns: `float`
The estimated difference TT-UT on the given date, expressed in seconds.

---

<a name="Ecliptic"></a>
### Ecliptic(equ)

**Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.**

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
which are relative to the plane of the Earth's orbit around the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Equatorial`](#Equatorial) | `equ` | Equatorial coordinates in the J2000 frame of reference. |

### Returns: [`EclipticCoordinates`](#EclipticCoordinates)
Ecliptic coordinates in the J2000 frame of reference.

---

<a name="EclipticLongitude"></a>
### EclipticLongitude(body, time)

**Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.**

This function calculates the angle around the plane of the Earth's orbit
of a celestial body, as seen from the center of the Sun.
The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).
time : Time
    The date and time at which the body's ecliptic longitude is to be calculated.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body other than the Sun. |

### Returns: `float`
An angular value in degrees indicating the ecliptic longitude of the body.

---

<a name="Elongation"></a>
### Elongation(body, time)

**Determines visibility of a celestial body relative to the Sun, as seen from the Earth.**

This function returns an [`ElongationEvent`](#ElongationEvent) object, which provides the following
information about the given celestial body at the given time:
- `visibility` is an enumerated type that specifies whether the body is more
  easily seen in the morning before sunrise, or in the evening after sunset.
- `elongation` is the angle in degrees between two vectors: one from the center
  of the Earth to the center of the Sun, the other from the center of the Earth
  to the center of the specified body. This angle indicates how far away the body
  is from the glare of the Sun. The elongation angle is always in the range [0, 180].
- `ecliptic_separation` is the absolute value of the difference between the body's
  ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center
  of the Earth. This angle measures around the plane of the Earth's orbit, and ignores
  how far above or below that plane the body is.
  The ecliptic separation is measured in degrees and is always in the range [0, 180].
time : Time
    The date and time of the observation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose visibility is to be calculated. |

### Returns: [`ElongationEvent`](#ElongationEvent)

---

<a name="Equator"></a>
### Equator(body, time, observer, ofdate, aberration)

**Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.**

Calculates topocentric equatorial coordinates in one of two different systems:
J2000 or true-equator-of-date, depending on the value of the `ofdate` parameter.
Equatorial coordinates include right ascension, declination, and distance in astronomical units.
This function corrects for light travel time: it adjusts the apparent location
of the observed body based on how long it takes for light to travel from the body to the Earth.
This function corrects for *topocentric parallax*, meaning that it adjusts for the
angular shift depending on where the observer is located on the Earth. This is most
significant for the Moon, because it is so close to the Earth. However, parallax corection
has a small effect on the apparent positions of other bodies.
Correction for aberration is optional, using the `aberration` parameter.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body to be observed. Not allowed to be `Body.Earth`. |
| [`Time`](#Time) | `time` | The date and time at which the observation takes place. |
| [`Observer`](#Observer) | `observer` | A location on or near the surface of the Earth. |
| `bool` | `ofdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. If `True`, returns coordinates using the equator and equinox of date. If `False`, returns coordinates converted to the J2000 system. |
| `bool` | `aberration` | If `True`, corrects for aberration of light based on the motion of the Earth with respect to the heliocentric origin. If `False`, does not correct for aberration. |

### Returns: [`Equatorial`](#Equatorial)
Equatorial coordinates in the specified frame of reference.

---

<a name="EquatorFromVector"></a>
### EquatorFromVector(vec)

**Given an equatorial vector, calculates equatorial angular coordinates.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vec` | A vector in an equatorial coordinate system. |

### Returns: [`Equatorial`](#Equatorial)
Angular coordinates expressed in the same equatorial system as `vec`.

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

### Returns: [`Vector`](#Vector)
The Moon's position as a vector in J2000 Cartesian equatorial coordinates.

---

<a name="GeoVector"></a>
### GeoVector(body, time, aberration)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Earth as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.
If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
the year range 1700..2200, this function will raise an exception.
Unlike [`HelioVector`](#HelioVector), this function always corrects for light travel time.
This means the position of the body is "back-dated" by the amount of time it takes
light to travel from that body to an observer on the Earth.
Also, the position can optionally be corrected for
[aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
causing the apparent direction of the body to be shifted due to transverse
movement of the Earth with respect to the rays of light coming from that body.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets. |
| [`Time`](#Time) | `time` | The date and time for which to calculate the position. |
| `bool` | `aberration` | A boolean value indicating whether to correct for aberration. |

### Returns: [`Vector`](#Vector)
A geocentric position vector of the center of the given body.

---

<a name="HelioDistance"></a>
### HelioDistance(body, time)

**Calculates the distance between a body and the Sun at a given time.**

Given a date and time, this function calculates the distance between
the center of `body` and the center of the Sun.
For the planets Mercury through Neptune, this function is significantly
more efficient than calling [`HelioVector`](#HelioVector) followed by taking the length
of the resulting vector.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric distance: the Sun, Moon, or any of the planets. |
| [`Time`](#Time) | `time` | The date and time for which to calculate the heliocentric distance. |

### Returns: `float`
The heliocentric distance in AU.

---

<a name="HelioVector"></a>
### HelioVector(body, time)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Sun as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.
The position is not corrected for light travel time or aberration.
This is different from the behavior of [`GeoVector`](#GeoVector).
If given an invalid value for `body`, or the body is `Body.Pluto` and `time` is outside
the year range 1700..2200, this function raise an exception.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose heliocentric position is to be calculated: The Sun, Moon, EMB, SSB, or any of the planets. |
| [`Time`](#Time) | `time` | The time at which to calculate the heliocentric position. |

### Returns: [`Vector`](#Vector)
A heliocentric position vector of the center of the given body
at the given time.

---

<a name="Horizon"></a>
### Horizon(time, observer, ra, dec, refraction)

**Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.**

Given a date and time, the geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
this function returns horizontal coordinates (azimuth and altitude angles) for the body
relative to the horizon at the geographic location.
The right ascension `ra` and declination `dec` passed in must be *equator of date*
coordinates, based on the Earth's true equator at the date and time of the observation.
Otherwise the resulting horizontal coordinates will be inaccurate.
Equator of date coordinates can be obtained by calling [`Equator`](#Equator), passing in
`True` as its `ofdate` parameter. It is also recommended to enable
aberration correction by passing in `True` for the `aberration` parameter.
This function optionally corrects for atmospheric refraction.
For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
correct for optical lensing of the Earth's atmosphere that causes objects
to appear somewhat higher above the horizon than they actually are.
However, callers may choose to avoid this correction by passing in `Refraction.Airless`.
If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
in the [`HorizontalCoordinates`](#HorizontalCoordinates) object returned by this function will all be corrected for refraction.
If refraction is disabled, none of these four coordinates will be corrected; in that case,
the right ascension and declination in the returned object will be numerically identical
to the respective `ra` and `dec` values passed in.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to find horizontal coordinates. |
| [`Observer`](#Observer) | `observer` | The location of the observer for which to find horizontal coordinates. |
| `float` | `ra` | Right ascension in sidereal hours of the celestial object, referred to the mean equinox of date for the J2000 epoch. |
| `float` | `dec` | Declination in degrees of the celestial object, referred to the mean equator of date for the J2000 epoch. Positive values are north of the celestial equator and negative values are south of it. |
| [`Refraction`](#Refraction) | `refraction` | The option for selecting whether to correct for atmospheric lensing. If `Refraction.Normal`, a well-behaved refraction model is used. If `Refraction.None`, no refraction correct is performed. `Refraction.JplHorizons` is used only for compatibility testing with the JPL Horizons online tool. |

### Returns: [`HorizontalCoordinates`](#HorizontalCoordinates)
The horizontal coordinates (altitude and azimuth), along with
equatorial coordinates (right ascension and declination), all
optionally corrected for atmospheric refraction. See remarks above
for more details.

---

<a name="HorizonFromVector"></a>
### HorizonFromVector(vector, refraction)

**Converts Cartesian coordinates to horizontal coordinates.**

Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
*IMPORTANT:* This function differs from `SphereFromVector` in two ways:
- `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
  from north (e.g., west = +90), but this function represents a clockwise rotation
  (e.g., east = +90). The difference is because `SphereFromVector` is intended
  to preserve the vector "right-hand rule", while this function defines azimuth in a more
  traditional way as used in navigation and cartography.
- This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.
The returned object contains the azimuth in `lon`.
It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
The altitude is stored in `lat`.
The distance to the observed object is stored in `dist`,
and is expressed in astronomical units (AU).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vector` | Cartesian vector to be converted to horizontal angular coordinates. |
| [`Refraction`](#Refraction) | `refraction` | See comments in the [`RefractionAngle`](#RefractionAngle) function. |

---

<a name="Illumination"></a>
### Illumination(body, time)

**Finds visual magnitude, phase angle, and other illumination information about a celestial body.**

This function calculates information about how bright a celestial body appears from the Earth,
reported as visual magnitude, which is a smaller (or even negative) number for brighter objects,
and a larger number for dimmer objects.
For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
of the body appears illuminated as seen from the Earth. For example, when the phase angle is
near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
of 90 degrees means the body appears "half full".
For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
so it doesn't have a phase angle.
When the body is Saturn, the returned object contains a field `ring_tilt` that holds
the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
0 for all bodies other than Saturn.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`Time`](#Time) | `time` | The date and time of the observation. |

### Returns: [`IlluminationInfo`](#IlluminationInfo)

---

<a name="InverseRefractionAngle"></a>
### InverseRefractionAngle(refraction, bent_altitude)

**Calculates the inverse of an atmospheric refraction angle.**

Given an observed altitude angle that includes atmospheric refraction,
calculate the negative angular correction to obtain the unrefracted
altitude. This is useful for cases where observed horizontal
coordinates are to be converted to another orientation system,
but refraction first must be removed from the observed position.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | `Refraction.Normal` - corrects for atmospheric refraction (recommended). `Refraction.Airless` - no correction is performed. `Refraction.JplHorizons` - For JPL Horizons compatibility testing only. |
| `float` | `bent_altitude` | The apparent altitude that includes atmospheric refraction. |

### Returns: `float`
The angular adjustment in degrees, to be added to the
altitude angle to correct for atmospheric lensing.
This will be less than or equal to zero.

---

<a name="InverseRotation"></a>
### InverseRotation(rotation)

**Calculates the inverse of a rotation matrix.**

Given a rotation matrix that performs some coordinate transform,
this function returns the matrix that reverses that trasnform.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | The rotation matrix to be inverted. |

### Returns: [`RotationMatrix`](#RotationMatrix)
The inverse rotation matrix.

---

<a name="LongitudeFromSun"></a>
### LongitudeFromSun(body, time)

**Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.**

This function can be used to determine where a planet appears around the ecliptic plane
(the plane of the Earth's orbit around the Sun) as seen from the Earth,
relative to the Sun's apparent position.
The angle starts at 0 when the body and the Sun are at the same ecliptic longitude
as seen from the Earth. The angle increases in the prograde direction
(the direction that the planets orbit the Sun and the Moon orbits the Earth).
When the angle is 180 degrees, it means the Sun and the body appear on opposite sides
of the sky for an Earthly observer. When `body` is a planet whose orbit around the
Sun is farther than the Earth's, 180 degrees indicates opposition. For the Moon,
it indicates a full moon.
The angle keeps increasing up to 360 degrees as the body's apparent prograde
motion continues relative to the Sun. When the angle reaches 360 degrees, it starts
over at 0 degrees.
Values between 0 and 180 degrees indicate that the body is visible in the evening sky
after sunset.  Values between 180 degrees and 360 degrees indicate that the body
is visible in the morning sky before sunrise.
time : Time
    The date and time of the observation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body for which to find longitude from the Sun. |

### Returns: `float`
An angle in degrees in the range [0, 360).

---

<a name="MoonPhase"></a>
### MoonPhase(time)

**Returns the Moon's phase as an angle from 0 to 360 degrees.**

This function determines the phase of the Moon using its apparent
ecliptic longitude relative to the Sun, as seen from the center of the Earth.
Certain values of the angle have conventional definitions:
- 0 = new moon
- 90 = first quarter
- 180 = full moon
- 270 = third quarter

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |

### Returns: `float`

---

<a name="NextLunarApsis"></a>
### NextLunarApsis(apsis)

**Finds the next lunar perigee or apogee in a series.**

This function requires an [`Apsis`](#Apsis) value obtained from a call to
[`SearchLunarApsis`](#SearchLunarApsis) or `NextLunarApsis`.
Given an apogee event, this function finds the next perigee event,
and vice versa.
See [`SearchLunarApsis`](#SearchLunarApsis) for more details.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Apsis`](#Apsis) | `apsis` |  |

### Returns: [`Apsis`](#Apsis)

---

<a name="NextMoonQuarter"></a>
### NextMoonQuarter(mq)

**Continues searching for lunar quarters from a previous search.**

After calling [`SearchMoonQuarter`](#SearchMoonQuarter), this function can be called
one or more times to continue finding consecutive lunar quarters.
This function finds the next consecutive moon quarter event after
the one passed in as the parameter `mq`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`MoonQuarter`](#MoonQuarter) | `mq` | A value returned by a prior call to [`SearchMoonQuarter`](#SearchMoonQuarter) or [`NextMoonQuarter`](#NextMoonQuarter). |

### Returns: [`MoonQuarter`](#MoonQuarter)

---

<a name="NextPlanetApsis"></a>
### NextPlanetApsis(body, apsis)

**Finds the next planetary perihelion or aphelion event in a series.**

This function requires an [`Apsis`](#Apsis) value obtained from a call
to [`SearchPlanetApsis`](#SearchPlanetApsis) or `NextPlanetApsis`.
Given an aphelion event, this function finds the next perihelion event, and vice versa.
See [`SearchPlanetApsis`](#SearchPlanetApsis) for more details.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. Must match the body passed into the call that produced the `apsis` parameter. |
| [`Apsis`](#Apsis) | `apsis` | An apsis event obtained from a call to [`SearchPlanetApsis`](#SearchPlanetApsis) or `NextPlanetApsis`. |

### Returns: [`Apsis`](#Apsis)

---

<a name="RefractionAngle"></a>
### RefractionAngle(refraction, altitude)

**Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.**

Given an altitude angle and a refraction option, calculates
the amount of "lift" caused by atmospheric refraction.
This is the number of degrees higher in the sky an object appears
due to lensing of the Earth's atmosphere.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | The option for selecting whether to correct for atmospheric lensing. If `Refraction.Normal`, a well-behaved refraction model is used. If `Refraction.Airless`, no refraction correct is performed. `Refraction.JplHorizons` is used only for compatibility testing with the JPL Horizons online tool. |
| `float` | `altitude` | The number of degrees above (positive) or below (negative) the horizon an object is, before being corrected for refraction. |

### Returns: `float`
The number of additional degrees of altitude an object appears
to have, due to atmospheric refraction, depending on the
option selected by the `refraction` parameter.

---

<a name="RotateVector"></a>
### RotateVector(rotation, vector)

**Applies a rotation to a vector, yielding a rotated vector.**

This function transforms a vector in one orientation to a vector
in another orientation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | A rotation matrix that specifies how the orientation of the vector is to be changed. |
| [`Vector`](#Vector) | `vector` | The vector whose orientation is to be changed. |

### Returns: [`Vector`](#Vector)
A vector in the orientation specified by `rotation`.

---

<a name="Rotation_ECL_EQD"></a>
### Rotation_ECL_EQD(time)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of date.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired equator. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to EQD.

---

<a name="Rotation_ECL_EQJ"></a>
### Rotation_ECL_EQJ()

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQJ = equatorial system, using equator at J2000 epoch.

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to EQJ.

---

<a name="Rotation_ECL_HOR"></a>
### Rotation_ECL_HOR(time, observer)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts ECL to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_EQD_ECL"></a>
### Rotation_EQD_ECL(time)

**Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of date.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the source equator. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD to ECL.

---

<a name="Rotation_EQD_EQJ"></a>
### Rotation_EQD_EQJ(time)

**Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator defines the source orientation. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD at `time` to EQJ.

---

<a name="Rotation_EQD_HOR"></a>
### Rotation_EQD_HOR(time, observer)

**Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value
to a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's location. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQD to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_EQJ_ECL"></a>
### Rotation_EQJ_ECL()

**Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECL = ecliptic system, using equator at J2000 epoch.

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to ECL.

---

<a name="Rotation_EQJ_EQD"></a>
### Rotation_EQJ_EQD(time)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator defines the target orientation. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to EQD at `time`.

---

<a name="Rotation_EQJ_HOR"></a>
### Rotation_EQJ_HOR(time, observer)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: HOR = horizontal system.
Use [`HorizonFromVector`](#HorizonFromVector) to convert the return value to
a traditional altitude/azimuth pair.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
The components of the horizontal vector are:
x = north, y = west, z = zenith (straight up from the observer).
These components are chosen so that the "right-hand rule" works for the vector
and so that north represents the direction where azimuth = 0.

---

<a name="Rotation_HOR_ECL"></a>
### Rotation_HOR_ECL(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the horizontal observation. |
| [`Observer`](#Observer) | `observer` | The location of the horizontal observer. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to ECL.

---

<a name="Rotation_HOR_EQD"></a>
### Rotation_HOR_EQD(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to EQD at `time` and for `observer`.

---

<a name="Rotation_HOR_EQJ"></a>
### Rotation_HOR_EQJ(time, observer)

**Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQJ = equatorial system, using equator at the J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time of the observation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that define's the observer's horizon. |

### Returns: [`RotationMatrix`](#RotationMatrix)
A rotation matrix that converts HOR to EQD at `time` and for `observer`.

---

<a name="Search"></a>
### Search(func, context, t1, t2, dt_tolerance_seconds)

**Searches for a time at which a function's value increases through zero.**

Certain astronomy calculations involve finding a time when an event occurs.
Often such events can be defined as the root of a function:
the time at which the function's value becomes zero.
`Search` finds the *ascending root* of a function: the time at which
the function's value becomes zero while having a positive slope. That is, as time increases,
the function transitions from a negative value, through zero at a specific moment,
to a positive value later. The goal of the search is to find that specific moment.
The search function is specified by two parameters: `func` and `context`.
The `func` parameter is a function itself that accepts a time
and a context containing any other arguments needed to evaluate the function.
The `context` parameter supplies that context for the given search.
As an example, a caller may wish to find the moment a celestial body reaches a certain
ecliptic longitude. In that case, the caller might create a type (class, tuple, whatever)
that contains a [`Body`](#Body) member to specify the body and a numeric value to hold the target longitude.
A different function might use a completely different context type.
Every time it is called, `func` returns a `float` value or it raises an exception.
If `func` raises an exception, the search immediately fails and the exception is
propagated back to the caller.
Otherwise, the search proceeds until it either finds the ascending root or fails for some reason.
The search calls `func` repeatedly to rapidly narrow in on any ascending
root within the time window specified by `t1` and `t2`. The search never
reports a solution outside this time window.
`Search` uses a combination of bisection and quadratic interpolation
to minimize the number of function calls. However, it is critical that the
supplied time window be small enough that there cannot be more than one root
(ascedning or descending) within it; otherwise the search can fail.
Beyond that, it helps to make the time window as small as possible, ideally
such that the function itself resembles a smooth parabolic curve within that window.
If an ascending root is not found, or more than one root
(ascending and/or descending) exists within the window `t1`..`t2`,
`Search` will return `None` to indicate a normal search failure.
If the search does not converge within 20 iterations, it will raise
an [`Error`](#Error) exception.
context : object
    An arbitrary data structure needed to be passed to the function `func`
    every time it is called.
t1 : float
    The lower time bound of the search window.
    See remarks above for more details.
t2 : float
    The upper time bound of the search window.
    See remarks above for more details.
dt_tolerance_seconds : float
    Specifies an amount of time in seconds within which a bounded ascending root
    is considered accurate enough to stop. A typical value is 1 second.

| Type | Parameter | Description |
| --- | --- | --- |
| `function(context, Time)` | `func` | A function that takes an arbitrary context parameter and a [`Time`](#Time) parameter. Returns a float value.  See remarks above for more details. |

### Returns: [`Time`](#Time) or `None`
If the search is successful, returns a #Time object that is within
`dt_tolerance_seconds` of an ascending root.
In this case, the returned time value will always be within the
inclusive range [`t1`, `t2`].
If there is no ascending root, or there is more than one ascending root,
the function returns `None`.

---

<a name="SearchHourAngle"></a>
### SearchHourAngle(body, observer, hourAngle, startTime)

**Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.**

The *hour angle* of a celestial body indicates its position in the sky with respect
to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
the number of hours that have passed since the most recent time that the body has culminated,
or reached its highest point.
This function searches for the next time a celestial body reaches the given hour angle
after the date and time specified by `startTime`.
To find when a body culminates, pass 0 for `hourAngle`.
To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.
Note that, especially close to the Earth's poles, a body as seen on a given day
may always be above the horizon or always below the horizon, so the caller cannot
assume that a culminating object is visible nor that an object is below the horizon
at its minimum altitude.
On success, the function reports the date and time, along with the horizontal coordinates
of the body at that time, as seen by the given observer.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body, which can the Sun, the Moon, or any planet other than the Earth. |
| [`Observer`](#Observer) | `observer` | Indicates a location on or near the surface of the Earth where the observer is located. |
| `float` | `hourAngle` | An hour angle value in the range [0.0, 24.0) indicating the number of sidereal hours after the body's most recent culmination. |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |

### Returns: [`HourAngleEvent`](#HourAngleEvent)

---

<a name="SearchLunarApsis"></a>
### SearchLunarApsis(startTime)

**Finds the time of the first lunar apogee or perigee after the given time.**

Given a date and time to start the search in `startTime`, this function finds
the next date and time that the center of the Moon reaches the closest or
farthest point in its orbit with respect to the center of the Earth, whichever
comes first after `startTime`.  The return value (of type [`Apsis`](#Apsis)) also
contains an indicator of whether the event is apogee or perigee.
The closest point is called *perigee* and the farthest point is called *apogee*.
The word *apsis* refers to either event.
To iterate through consecutive alternating perigee and apogee events,
call [`SearchLunarApsis`](#SearchLunarApsis) once, then use the return value to call [`NextLunarApsis`](#NextLunarApsis).
After that, keep feeding the previous return value from `NextLunarApsis` into
another call of `NextLunarApsis` as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time at which to start searching for the next perigee or apogee. |

### Returns: [`Apsis`](#Apsis)

---

<a name="SearchMaxElongation"></a>
### SearchMaxElongation(body, startTime)

**Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.**

Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
a telescope without atmospheric interference, are when these planets reach maximum elongation.
These are events where the planets reach the maximum angle from the Sun as seen from the Earth.
This function solves for those times, reporting the next maximum elongation event's date and time,
the elongation value itself, the relative longitude with the Sun, and whether the planet is best
observed in the morning or evening. See [`ElongationEvent`](#ElongationEvent) for more details about the returned object.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception. To find the best viewing opportunities for planets farther from the Sun than the Earth is (Mars through Pluto), use [`SearchRelativeLongitude`](#SearchRelativeLongitude) to find the next opposition event. |
| [`Time`](#Time) | `startTime` | The date and time at which to begin the search. The maximum elongation event found will always be the first one that occurs after this date and time. |

### Returns: [`ElongationEvent`](#ElongationEvent)

---

<a name="SearchMoonPhase"></a>
### SearchMoonPhase(targetLon, startTime, limitDays)

**Searches for the time that the Moon reaches a specified phase.**

Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
longitude with respect to the Sun's geocentric ecliptic longitude.
When the Moon and the Sun have the same longitude, that is defined as a new moon.
When their longitudes are 180 degrees apart, that is defined as a full moon.
This function searches for any value of the lunar phase expressed as an
angle in degrees in the range [0, 360).
If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
it is much easier to call the functions [`SearchMoonQuarter`](#SearchMoonQuarter) and [`NextMoonQuarter`](#NextMoonQuarter).
This function is useful for finding general phase angles outside those four quarters.

| Type | Parameter | Description |
| --- | --- | --- |
| `float` | `targetLon` | The difference in geocentric longitude between the Sun and Moon that specifies the lunar phase being sought. This can be any value in the range [0, 360).  Certain values have conventional names: 0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter. |
| [`Time`](#Time) | `startTime` | The beginning of the time window in which to search for the Moon reaching the specified phase. |
| `float` | `limitDays` | The number of days after `startTime` that limits the time window for the search. |

### Returns: [`Time`](#Time) or `None`

---

<a name="SearchMoonQuarter"></a>
### SearchMoonQuarter(startTime)

**Finds the first lunar quarter after the specified date and time.**

A lunar quarter is one of the following four lunar phase events:
new moon, first quarter, full moon, third quarter.
This function finds the lunar quarter that happens soonest
after the specified date and time.
To continue iterating through consecutive lunar quarters, call this function once,
followed by calls to [`NextMoonQuarter`](#NextMoonQuarter) as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |

### Returns: [`MoonQuarter`](#MoonQuarter)

---

<a name="SearchPeakMagnitude"></a>
### SearchPeakMagnitude(body, startTime)

**Searches for the date and time Venus will next appear brightest as seen from the Earth.**

This function searches for the date and time Venus appears brightest as seen from the Earth.
Currently only Venus is supported for the `body` parameter, though this could change in the future.
Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see
from the Earth, so peak magnitude events have little practical value for that planet.
Planets other than Venus and Mercury reach peak magnitude at opposition, which can
be found using [`SearchRelativeLongitude`](#SearchRelativeLongitude).
The Moon reaches peak magnitude at full moon, which can be found using
[`SearchMoonQuarter`](#SearchMoonQuarter) or [`SearchMoonPhase`](#SearchMoonPhase).
The Sun reaches peak magnitude at perihelion, which occurs each year in January.
However, the difference is minor and has little practical value.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Currently only `Body.Venus` is allowed. Any other value results in an exception. See remarks above for more details. |
| [`Time`](#Time) | `startTime` | The date and time to start searching for the next peak magnitude event. |

### Returns: [`IlluminationInfo`](#IlluminationInfo)

---

<a name="SearchPlanetApsis"></a>
### SearchPlanetApsis(body, startTime)

**Finds the next planet perihelion or aphelion, after a given time.**

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the specified planet reaches the closest or farthest point
in its orbit with respect to the center of the Sun, whichever comes first after `startTime`.
The closest point is called *perihelion* and the farthest point is called *aphelion*.
The word *apsis* refers to either event.
To iterate through consecutive alternating perihelion and aphelion events,
call `SearchPlanetApsis` once, then use the return value to call [`NextPlanetApsis`](#NextPlanetApsis).
After that, keep feeding the previous return value from `NextPlanetApsis`
into another call of `NextPlanetApsis` as many times as desired.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. |
| [`Time`](#Time) | `startTime` | The date and time at which to start searching for the next perihelion or aphelion. |

### Returns: [`Apsis`](#Apsis)

---

<a name="SearchRelativeLongitude"></a>
### SearchRelativeLongitude(body, targetRelLon, startTime)

**Searches for when the Earth and another planet are separated by a certain ecliptic longitude.**

Searches for the time when the Earth and another planet are separated by a specified angle
in ecliptic longitude, as seen from the Sun.
A relative longitude is the angle between two bodies measured in the plane of the
Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic
plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane,
and the angle measured around that plane from one body to the other in the direction
the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees.
This is the relative longitude.
Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
this function searches for the next time that the relative longitude measured from the
planet to the Earth is `targetRelLon`.
Certain astronomical events are defined in terms of relative longitude between
the Earth and another planet:
- When the relative longitude is 0 degrees, it means both planets are in the same
  direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus),
  this is known as *inferior conjunction*, a time when the other planet becomes very
  difficult to see because of being lost in the Sun's glare.
  (The only exception is in the rare event of a transit, when we see the silhouette
  of the planet passing between the Earth and the Sun.)
- When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
  this is known as *opposition*. Opposition is when the planet is closest to the Earth,
  and also when it is visible for most of the night, so it is considered the best time
  to observe the planet.
- When the relative longitude is 180 degrees, it means the other planet is on the opposite
  side of the Sun from the Earth.  This is called *superior conjunction*.  Like inferior
  conjunction, the planet is very difficult to see from the Earth.
  Superior conjunction is possible for any planet other than the Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A planet other than the Earth. If `body` is not a planet, or if it is `Body.Earth`, an error occurs. |
| `float` | `targetRelLon` | The desired relative longitude, expressed in degrees. Must be in the range [0, 360). |
| [`Time`](#Time) | `startTime` | The date and time at which to begin the search. |

### Returns: [`Time`](#Time)
The date and time of the relative longitude event.

---

<a name="SearchRiseSet"></a>
### SearchRiseSet(body, observer, direction, startTime, limitDays)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.**

This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
Rise time is when the body first starts to be visible above the horizon.
For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
Set time is the moment when the body appears to vanish below the horizon.
This function corrects for typical atmospheric refraction, which causes celestial
bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).
Note that rise or set may not occur in every 24 hour period.
For example, near the Earth's poles, there are long periods of time where
the Sun stays below the horizon, never rising.
Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
significant amount during each rotation of the Earth.
Therefore callers must not assume that the function will always succeed.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`Observer`](#Observer) | `observer` | The location where observation takes place. |
| [`Direction`](#Direction) | `direction` | Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time. |
| [`Time`](#Time) | `startTime` | The date and time at which to start the search. |
| `float` | `limitDays` | Limits how many days to search for a rise or set time. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |

### Returns: [`Time`](#Time) or `None`
If the rise or set time is found within the specified time window,
this function returns that time. Otherwise, it returns `None`.

---

<a name="SearchSunLongitude"></a>
### SearchSunLongitude(targetLon, startTime, limitDays)

**Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.**

This function finds the moment in time, if any exists in the given time window,
that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.
This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call [`Seasons`](#Seasons)
to calculate all equinoxes and solstices for a given calendar year.
The function searches the window of time specified by `startTime` and `startTime+limitDays`.
The search will return `None` if the Sun never reaches the longitude `targetLon` or
if the window is so large that the longitude ranges more than 180 degrees within it.
It is recommended to keep the window smaller than 10 days when possible.
targetLon : float
     The desired ecliptic longitude in degrees, relative to the true equinox of date.
     This may be any value in the range [0, 360), although certain values have
     conventional meanings:
     0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice.
startTime : Time
     The date and time for starting the search for the desired longitude event.
limitDays : float
     The real-valued number of days, which when added to `startTime`, limits the
     range of time over which the search looks.
     It is recommended to keep this value between 1 and 10 days.
     See remarks above for more details.

### Returns: [`Time`](#Time) or `None`

---

<a name="Seasons"></a>
### Seasons(year)

**Finds both equinoxes and both solstices for a given calendar year.**

The changes of seasons are defined by solstices and equinoxes.
Given a calendar year number, this function calculates the
March and September equinoxes and the June and December solstices.
The equinoxes are the moments twice each year when the plane of the
Earth's equator passes through the center of the Sun. In other words,
the Sun's declination is zero at both equinoxes.
The March equinox defines the beginning of spring in the northern hemisphere
and the beginning of autumn in the southern hemisphere.
The September equinox defines the beginning of autumn in the northern hemisphere
and the beginning of spring in the southern hemisphere.
The solstices are the moments twice each year when one of the Earth's poles
is most tilted toward the Sun. More precisely, the Sun's declination reaches
its minimum value at the December solstice, which defines the beginning of
winter in the northern hemisphere and the beginning of summer in the southern
hemisphere. The Sun's declination reaches its maximum value at the June solstice,
which defines the beginning of summer in the northern hemisphere and the beginning
of winter in the southern hemisphere.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The calendar year number for which to calculate equinoxes and solstices. The value may be any integer, but only the years 1800 through 2100 have been validated for accuracy: unit testing against data from the United States Naval Observatory confirms that all equinoxes and solstices for that range of years are within 2 minutes of the correct time. |

### Returns: [`SeasonInfo`](#SeasonInfo)

---

<a name="SphereFromVector"></a>
### SphereFromVector(vector)

**Converts Cartesian coordinates to spherical coordinates.**

Given a Cartesian vector, returns latitude, longitude, and distance.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Vector`](#Vector) | `vector` | Cartesian vector to be converted to spherical coordinates. |

### Returns: [`Spherical`](#Spherical)
Spherical coordinates that are equivalent to the given vector.

---

<a name="SunPosition"></a>
### SunPosition(time)

**Calculates geocentric ecliptic coordinates for the Sun.**

This function calculates the position of the Sun as seen from the Earth.
The returned value includes both Cartesian and spherical coordinates.
The x-coordinate and longitude values in the returned object are based
on the *true equinox of date*: one of two points in the sky where the instantaneous
plane of the Earth's equator at the given date and time (the *equatorial plane*)
intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
By convention, the apparent location of the Sun at the March equinox is chosen
as the longitude origin and x-axis direction, instead of the one for September.
`SunPosition` corrects for precession and nutation of the Earth's axis
in order to obtain the exact equatorial plane at the given time.
This function can be used for calculating changes of seasons: equinoxes and solstices.
In fact, the function [`Seasons`](#Seasons) does use this function for that purpose.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Time`](#Time) | `time` | The date and time for which to calculate the Sun's position. |

### Returns: [`EclipticCoordinates`](#EclipticCoordinates)
The ecliptic coordinates of the Sun using the Earth's true equator of date.

---

<a name="VectorFromEquator"></a>
### VectorFromEquator(equ, time)

**Given angular equatorial coordinates in `equ`, calculates equatorial vector.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Equatorial`](#Equatorial) | `equ` | Angular equatorial coordinates to be converted to a vector. |
| [`Time`](#Time) | `time` | The date and time of the observation. This is needed because the returned vector object requires a valid time value when passed to certain other functions. |

### Returns: [`Vector`](#Vector)
A vector in the equatorial system.

---

<a name="VectorFromHorizon"></a>
### VectorFromHorizon(sphere, time, refraction)

**Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | A structure that contains apparent horizontal coordinates: `lat` holds the refracted azimuth angle, `lon` holds the azimuth in degrees clockwise from north, and `dist` holds the distance from the observer to the object in AU. |
| [`Time`](#Time) | `time` | The date and time of the observation. This is needed because the returned vector object requires a valid time value when passed to certain other functions. |
| [`Refraction`](#Refraction) | `refraction` | See remarks in function [`RefractionAngle`](#RefractionAngle). |

### Returns: [`Vector`](#Vector)
A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).

---

<a name="VectorFromSphere"></a>
### VectorFromSphere(sphere, time)

**Converts spherical coordinates to Cartesian coordinates.**

Given spherical coordinates and a time at which they are valid,
returns a vector of Cartesian coordinates. The returned value
includes the time, as required by all `Time` objects.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | Spherical coordinates to be converted. |
| [`Time`](#Time) | `time` | The time that should be included in the returned vector. |

### Returns: [`Vector`](#Vector)
The vector form of the supplied spherical coordinates.

