# Astronomy Engine (C/C++)

This is the complete programming reference for the C version of
Astronomy Engine. It can be used directly from C++ programs also.
Other programming languages are supported.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start
To get started quickly, here are some [examples](../../demo/c/).

---

## Contents

- [Topic Index](#topics)
- [Functions](#functions)
- [Enumerated Types](#enums)
- [Structures](#structs)
- [Type Definitions](#typedefs)

---

<a name="topics"></a>
## Topic Index

### Dates and times

| Function | Description |
| -------- | ----------- |
| [CurrentTime](#Astronomy_CurrentTime) | Obtains the current date and time of the computer's clock in the form of an [`astro_time_t`](#astro_time_t) that can be used for astronomy calculations. |
| [MakeTime](#Astronomy_MakeTime) | Converts a UTC calendar date and time given as separate numeric parameters into an [`astro_time_t`](#astro_time_t) that can be used for astronomy calculations. |
| [AddDays](#Astronomy_AddDays) | Adds or subtracts an amount of time to an [`astro_time_t`](#astro_time_t) to get another [`astro_time_t`](#astro_time_t). |
| [TimeFromUtc](#Astronomy_TimeFromUtc) | Converts UTC calendar date and time from an [`astro_utc_t`](#astro_utc_t) structure to an [`astro_time_t`](#astro_time_t) structure that can be used for astronomy calculations. |
| [UtcFromTime](#Astronomy_UtcFromTime) | Converts an astronomical [`astro_time_t`](#astro_time_t) time value to an [`astro_utc_t`](#astro_utc_t) structure that can be used for displaying a UTC calendar date and time. |

### Celestial bodies

| Function | Description |
| -------- | ----------- |
| [BodyCode](#Astronomy_BodyCode) | Converts the English name of a celestial body to its equivalent [`astro_body_t`](#astro_body_t) enumeration value. |
| [BodyName](#Astronomy_BodyName) | Converts an [`astro_body_t`](#astro_body_t) enumeration value to its equivalent English name as a string. |

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#Astronomy_HelioVector) | Calculates vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy_GeoVector)     | Calculates vector with respect to the center of the Earth. |
| [Equator](#Astronomy_Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Astronomy_Ecliptic)       | Converts J2000 equatorial coordinates to J2000 ecliptic coordinates. |
| [EclipticLongitude](#Astronomy_EclipticLongitude) | Calculates ecliptic longitude of a body in the J2000 system. |
| [Horizon](#Astronomy_Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [LongitudeFromSun](#Astronomy_LongitudeFromSun) | Calculates a body's apparent ecliptic longitude difference from the Sun, as seen by an observer on the Earth. |

### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#Astronomy_SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchHourAngle](#Astronomy_SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#Astronomy_MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#Astronomy_SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#Astronomy_SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy_NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#Astronomy_SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy_NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Visual magnitude and elongation

| Function | Description |
| -------- | ----------- |
| [Illumination](#Astronomy_Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#Astronomy_SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#Astronomy_AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Astronomy_Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#Astronomy_SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| Function | Description |
| -------- | ----------- |
| [SearchRelativeLongitude](#Astronomy_SearchRelativeLongitude) | Finds oppositions and conjunctions of planets. |

### Equinoxes, solstices, and apparent solar motion

| Function | Description |
| -------- | ----------- |
| [SearchSunLongitude](#Astronomy_SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the *true equator of date* system. |
| [Seasons](#Astronomy_Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#Astronomy_SunPosition) | Calculates the Sun's apparent ecliptic coordinates as seen from the Earth. |

---


<a name="functions"></a>
## Functions



---

<a name="Astronomy_AddDays"></a>
### Astronomy_AddDays(time, days) &#8658; [`astro_time_t`](#astro_time_t)

**Calculates the sum or difference of an [`astro_time_t`](#astro_time_t) with a specified floating point number of days.** 



Sometimes we need to adjust a given [`astro_time_t`](#astro_time_t) value by a certain amount of time. This function adds the given real number of days in `days` to the date and time in `time`.

More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time, according to the historical and predictive Delta-T model provided by the [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).

The value stored in `time` will not be modified; it is passed by value.



**Returns:**  A date and time that is conceptually equal to `time + days`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  A date and time for which to calculate an adjusted date and time.  | 
| `double` | `days` |  A floating point number of days by which to adjust `time`. May be negative, 0, or positive.  | 




---

<a name="Astronomy_AngleFromSun"></a>
### Astronomy_AngleFromSun(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Returns the angle between the given body and the Sun, as seen from the Earth.** 



This function calculates the angular separation between the given body and the Sun, as seen from the center of the Earth. This angle is helpful for determining how easy it is to see the body away from the glare of the Sun.



**Returns:**  If successful, the returned structure contains `ASTRO_SUCCESS` in the `status` field and `angle` holds the angle in degrees between the Sun and the specified body as seen from the center of the Earth. If an error occurs, the `status` field contains a value other than `ASTRO_SUCCESS` that indicates the error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body whose angle from the Sun is to be measured. Not allowed to be `BODY_EARTH`. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The time at which the observation is made. | 




---

<a name="Astronomy_BodyCode"></a>
### Astronomy_BodyCode(name) &#8658; [`astro_body_t`](#astro_body_t)

**Returns the [`astro_body_t`](#astro_body_t) value corresponding to the given English name.** 





**Returns:**  If `name` is one of the strings (case-sensitive) listed above, the returned value is the corresponding [`astro_body_t`](#astro_body_t) value, otherwise it is [`BODY_INVALID`](#BODY_INVALID). 



| Type | Parameter | Description |
| --- | --- | --- |
| `const char *` | `name` |  One of the following strings: Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto.  | 




---

<a name="Astronomy_BodyName"></a>
### Astronomy_BodyName(body) &#8658; `const char *`

**Finds the name of a celestial body.** 





**Returns:**  The English-language name of the celestial body, or "" if the body is not valid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body whose name is to be found.  | 




---

<a name="Astronomy_CombineRotation"></a>
### Astronomy_CombineRotation(a, b) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Creates a rotation based on applying one rotation followed by another.** 



Given two rotation matrices, returns a combined rotation matrix that is equivalent to rotating based on the first matrix, followed by the second.



**Returns:**  The combined rotation matrix. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `a` |  The first rotation to apply. | 
| [`astro_rotation_t`](#astro_rotation_t) | `b` |  The second rotation to apply. | 




---

<a name="Astronomy_CurrentTime"></a>
### Astronomy_CurrentTime() &#8658; [`astro_time_t`](#astro_time_t)

**Returns the computer's current date and time in the form of an [`astro_time_t`](#astro_time_t).** 



Uses the computer's system clock to find the current UTC date and time with 1-second granularity. Converts that date and time to an [`astro_time_t`](#astro_time_t) value and returns the result. Callers can pass this value to other Astronomy Engine functions to calculate current observational conditions. 

---

<a name="Astronomy_Ecliptic"></a>
### Astronomy_Ecliptic(equ) &#8658; [`astro_ecliptic_t`](#astro_ecliptic_t)

**Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.** 



Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates, which are relative to the plane of the Earth's orbit around the Sun.



**Returns:**  Ecliptic coordinates in the J2000 frame of reference. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `equ` |  Equatorial coordinates in the J2000 frame of reference. You can call [`Astronomy_GeoVector`](#Astronomy_GeoVector) to obtain suitable equatorial coordinates. | 




---

<a name="Astronomy_EclipticLongitude"></a>
### Astronomy_EclipticLongitude(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.** 



This function calculates the angle around the plane of the Earth's orbit of a celestial body, as seen from the center of the Sun. The angle is measured prograde (in the direction of the Earth's orbit around the Sun) in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).



**Returns:**  On success, returns a structure whose `status` is `ASTRO_SUCCESS` and whose `angle` holds the ecliptic longitude in degrees. On failure, `status` holds a value other than `ASTRO_SUCCESS`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body other than the Sun. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time at which the body's ecliptic longitude is to be calculated. | 




---

<a name="Astronomy_Elongation"></a>
### Astronomy_Elongation(body, time) &#8658; [`astro_elongation_t`](#astro_elongation_t)

**Determines visibility of a celestial body relative to the Sun, as seen from the Earth.** 



This function returns an [`astro_elongation_t`](#astro_elongation_t) structure, which provides the following information about the given celestial body at the given time:



- `visibility` is an enumerated type that specifies whether the body is more easily seen in the morning before sunrise, or in the evening after sunset.
- `elongation` is the angle in degrees between two vectors: one from the center of the Earth to the center of the Sun, the other from the center of the Earth to the center of the specified body. This angle indicates how far away the body is from the glare of the Sun. The elongation angle is always in the range [0, 180].
- `ecliptic_separation` is the absolute value of the difference between the body's ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures around the plane of the Earth's orbit, and ignores how far above or below that plane the body is. The ecliptic separation is measured in degrees and is always in the range [0, 180].




**Returns:**  If successful, the `status` field in the returned structure contains `ASTRO_SUCCESS` and all the other fields in the structure are valid. On failure, `status` contains some other value as an error code and the other fields contain invalid values. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body whose visibility is to be calculated. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_Equator"></a>
### Astronomy_Equator(body, time, observer, equdate, aberration) &#8658; [`astro_equatorial_t`](#astro_equatorial_t)

**Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.** 



Calculates topocentric equatorial coordinates in one of two different systems: J2000 or true-equator-of-date, depending on the value of the `equdate` parameter. Equatorial coordinates include right ascension, declination, and distance in astronomical units.

This function corrects for light travel time: it adjusts the apparent location of the observed body based on how long it takes for light to travel from the body to the Earth.

This function corrects for *topocentric parallax*, meaning that it adjusts for the angular shift depending on where the observer is located on the Earth. This is most significant for the Moon, because it is so close to the Earth. However, parallax corection has a small effect on the apparent positions of other bodies.

Correction for aberration is optional, using the `aberration` parameter.



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body to be observed. Not allowed to be `BODY_EARTH`.  | 
| [`astro_time_t *`](#astro_time_t *) | `time` |  The date and time at which the observation takes place.  | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location on or near the surface of the Earth.  | 
| [`astro_equator_date_t`](#astro_equator_date_t) | `equdate` |  Selects the date of the Earth's equator in which to express the equatorial coordinates.  | 
| [`astro_aberration_t`](#astro_aberration_t) | `aberration` |  Selects whether or not to correct for aberration.  | 




---

<a name="Astronomy_GeoMoon"></a>
### Astronomy_GeoMoon(time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates the geocentric position of the Moon at a given time.** 



Given a time of observation, calculates the Moon's position as a vector. The vector gives the location of the Moon's center relative to the Earth's center with x-, y-, and z-components measured in astronomical units.

This algorithm is based on Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954, which in turn derives from E. W. Brown's lunar theories from the early twentieth century. It is adapted from Turbo Pascal code from the book [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210) by Montenbruck and Pfleger.



**Returns:**  The Moon's position as a vector in J2000 Cartesian equatorial coordinates. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the Moon's position.  | 




---

<a name="Astronomy_GeoVector"></a>
### Astronomy_GeoVector(body, time, aberration) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.** 



This function calculates the position of the given celestial body as a vector, using the center of the Earth as the origin. The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

If given an invalid value for `body`, or the body is `BODY_PLUTO` and the `time` is outside the year range 1700..2200, this function will fail. The caller should always check the `status` field inside the returned [`astro_vector_t`](#astro_vector_t) for `ASTRO_SUCCESS` (success) or any other value (failure) before trusting the resulting vector.

Unlike [`Astronomy_HelioVector`](#Astronomy_HelioVector), this function always corrects for light travel time. This means the position of the body is "back-dated" by the amount of time it takes light to travel from that body to an observer on the Earth.

Also, the position can optionally be corrected for [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect causing the apparent direction of the body to be shifted due to transverse movement of the Earth with respect to the rays of light coming from that body.



**Returns:**  A geocentric position vector of the center of the given body. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the position.  | 
| [`astro_aberration_t`](#astro_aberration_t) | `aberration` |  `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected.  | 




---

<a name="Astronomy_HelioVector"></a>
### Astronomy_HelioVector(body, time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.** 



This function calculates the position of the given celestial body as a vector, using the center of the Sun as the origin. The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

The position is not corrected for light travel time or aberration. This is different from the behavior of [`Astronomy_GeoVector`](#Astronomy_GeoVector).

If given an invalid value for `body`, or the body is `BODY_PLUTO` and the `time` is outside the year range 1700..2200, this function will fail. The caller should always check the `status` field inside the returned [`astro_vector_t`](#astro_vector_t) for `ASTRO_SUCCESS` (success) or any other value (failure) before trusting the resulting vector.



**Returns:**  A heliocentric position vector of the center of the given body. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the position.  | 




---

<a name="Astronomy_Horizon"></a>
### Astronomy_Horizon(time, observer, ra, dec, refraction) &#8658; [`astro_horizon_t`](#astro_horizon_t)

**Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.** 



Given a date and time, the geographic location of an observer on the Earth, and equatorial coordinates (right ascension and declination) of a celestial body, this function returns horizontal coordinates (azimuth and altitude angles) for the body relative to the horizon at the geographic location.

The right ascension `ra` and declination `dec` passed in must be *equator of date* coordinates, based on the Earth's true equator at the date and time of the observation. Otherwise the resulting horizontal coordinates will be inaccurate. Equator of date coordinates can be obtained by calling [`Astronomy_Equator`](#Astronomy_Equator), passing in `EQUATOR_OF_DATE` as its `equdate` parameter. It is also recommended to enable aberration correction by passing in `ABERRATION` as the `aberration` parameter.

This function optionally corrects for atmospheric refraction. For most uses, it is recommended to pass `REFRACTION_NORMAL` in the `refraction` parameter to correct for optical lensing of the Earth's atmosphere that causes objects to appear somewhat higher above the horizon than they actually are. However, callers may choose to avoid this correction by passing in `REFRACTION_NONE`. If refraction correction is enabled, the azimuth, altitude, right ascension, and declination in the [`astro_horizon_t`](#astro_horizon_t) structure returned by this function will all be corrected for refraction. If refraction is disabled, none of these four coordinates will be corrected; in that case, the right ascension and declination in the returned structure will be numerically identical to the respective `ra` and `dec` values passed in.



**Returns:**  The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t *`](#astro_time_t *) | `time` |  The date and time of the observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of the observer. | 
| `double` | `ra` |  The right ascension of the body in sidereal hours. See remarks above for more details. | 
| `double` | `dec` |  The declination of the body in degrees. See remarks above for more details. | 
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  Selects whether to correct for atmospheric refraction, and if so, which model to use. The recommended value for most uses is `REFRACTION_NORMAL`. See remarks above for more details. | 




---

<a name="Astronomy_HorizonFromVector"></a>
### Astronomy_HorizonFromVector(vector, refraction) &#8658; [`astro_spherical_t`](#astro_spherical_t)

**Converts Cartesian coordinates to horizontal coordinates.** 



Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

*IMPORTANT:* This function differs from [`Astronomy_SphereFromVector`](#Astronomy_SphereFromVector) in two ways:

- `Astronomy_SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise from north (e.g., west = +90), but this function represents a clockwise rotation (e.g., east = +90). The difference is because `Astronomy_SphereFromVector` is intended to preserve the vector "right-hand rule", while this function defines azimuth in a more traditional way as used in navigation and cartography.
- This function optionally corrects for atmospheric refraction, while `Astronomy_SphereFromVector` does not.


The returned structure will contain the azimuth in `lon`. It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

The altitude will be stored in `lat`.

The distance to the observed object is stored in `dist`, and is expressed in astronomical units (AU).



**Returns:**  If successful, `status` hold `ASTRO_SUCCESS` and the other fields are valid as described above. Otherwise `status` holds an error code and the other fields are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `vector` |  Cartesian vector to be converted to horizontal coordinates. | 
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  `REFRACTION_NORMAL`: correct altitude for atmospheric refraction (recommended). `REFRACTION_NONE`: no atmospheric refraction correction is performed. `REFRACTION_JPLHOR`: for JPL Horizons compatibility testing only; not recommended for normal use. | 




---

<a name="Astronomy_Illumination"></a>
### Astronomy_Illumination(body, time) &#8658; [`astro_illum_t`](#astro_illum_t)

**Finds visual magnitude, phase angle, and other illumination information about a celestial body.** 



This function calculates information about how bright a celestial body appears from the Earth, reported as visual magnitude, which is a smaller (or even negative) number for brighter objects and a larger number for dimmer objects.

For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction of the body appears illuminated as seen from the Earth. For example, when the phase angle is near zero, it means the body appears "full" as seen from the Earth. A phase angle approaching 180 degrees means the body appears as a thin crescent as seen from the Earth. A phase angle of 90 degrees means the body appears "half full". For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it, so it doesn't have a phase angle.

When the body is Saturn, the returned structure contains a field `ring_tilt` that holds the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds 0 for all bodies other than Saturn.



**Returns:**  On success, the `status` field of the return structure holds `ASTRO_SUCCESS` and the other structure fields are valid. Any other value indicates an error, in which case the remaining structure fields are not valid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Moon, or any planet other than the Earth. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_InverseRotation"></a>
### Astronomy_InverseRotation(rotation) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates the inverse of a rotation matrix.** 



Given a rotation matrix that performs some coordinate transform, this function returns the matrix that reverses that trasnform.



**Returns:**  A rotation matrix that performs the opposite transformation. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `rotation` |  The rotation matrix to be inverted. | 




---

<a name="Astronomy_LongitudeFromSun"></a>
### Astronomy_LongitudeFromSun(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.** 



This function can be used to determine where a planet appears around the ecliptic plane (the plane of the Earth's orbit around the Sun) as seen from the Earth, relative to the Sun's apparent position.

The angle starts at 0 when the body and the Sun are at the same ecliptic longitude as seen from the Earth. The angle increases in the prograde direction (the direction that the planets orbit the Sun and the Moon orbits the Earth).

When the angle is 180 degrees, it means the Sun and the body appear on opposite sides of the sky for an Earthly observer. When `body` is a planet whose orbit around the Sun is farther than the Earth's, 180 degrees indicates opposition. For the Moon, it indicates a full moon.

The angle keeps increasing up to 360 degrees as the body's apparent prograde motion continues relative to the Sun. When the angle reaches 360 degrees, it starts over at 0 degrees.

Values between 0 and 180 degrees indicate that the body is visible in the evening sky after sunset. Values between 180 degrees and 360 degrees indicate that the body is visible in the morning sky before sunrise.



**Returns:**  On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the `angle` field holds a value in the range [0, 360). On failure, the `status` field contains some other value indicating an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body for which to find longitude from the Sun. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_MakeObserver"></a>
### Astronomy_MakeObserver(latitude, longitude, height) &#8658; [`astro_observer_t`](#astro_observer_t)

**Creates an observer object that represents a location on or near the surface of the Earth.** 



Some Astronomy Engine functions calculate values pertaining to an observer on the Earth. These functions require a value of type [`astro_observer_t`](#astro_observer_t) that represents the location of such an observer.



**Returns:**  An observer object that can be passed to astronomy functions that require a geographic location. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `latitude` |  The geographic latitude of the observer in degrees north (positive) or south (negative) of the equator.  | 
| `double` | `longitude` |  The geographic longitude of the observer in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.  | 
| `double` | `height` |  The height of the observer in meters above mean sea level.  | 




---

<a name="Astronomy_MakeTime"></a>
### Astronomy_MakeTime(year, month, day, hour, minute, second) &#8658; [`astro_time_t`](#astro_time_t)

**Creates an [`astro_time_t`](#astro_time_t) value from a given calendar date and time.** 



Given a UTC calendar date and time, calculates an [`astro_time_t`](#astro_time_t) value that can be passed to other Astronomy Engine functions for performing various calculations relating to that date and time.

It is the caller's responsibility to ensure that the parameter values are correct. The parameters are not checked for validity, and this function never returns any indication of an error. Invalid values, for example passing in February 31, may cause unexpected return values.



**Returns:**  An [`astro_time_t`](#astro_time_t) value that represents the given calendar date and time. 



| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` |  The UTC calendar year, e.g. 2019.  | 
| `int` | `month` |  The UTC calendar month in the range 1..12.  | 
| `int` | `day` |  The UTC calendar day in the range 1..31.  | 
| `int` | `hour` |  The UTC hour of the day in the range 0..23.  | 
| `int` | `minute` |  The UTC minute in the range 0..59.  | 
| `double` | `second` |  The UTC floating-point second in the range [0, 60). | 




---

<a name="Astronomy_MoonPhase"></a>
### Astronomy_MoonPhase(time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Returns the Moon's phase as an angle from 0 to 360 degrees.** 



This function determines the phase of the Moon using its apparent ecliptic longitude relative to the Sun, as seen from the center of the Earth. Certain values of the angle have conventional definitions:



- 0 = new moon
- 90 = first quarter
- 180 = full moon
- 270 = third quarter




**Returns:**  On success, the function returns the angle as described above in the `angle` field and `ASTRO_SUCCESS` in the `status` field. The function should always succeed, but it is a good idea for callers to check the `status` field in the returned structure. Any other value in `status` indicates a failure that should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_NextLunarApsis"></a>
### Astronomy_NextLunarApsis(apsis) &#8658; [`astro_apsis_t`](#astro_apsis_t)

**Finds the next lunar perigee or apogee event in a series.** 



This function requires an [`astro_apsis_t`](#astro_apsis_t) value obtained from a call to [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) or `Astronomy_NextLunarApsis`. Given an apogee event, this function finds the next perigee event, and vice versa.

See [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) for more details.



**Returns:**  Same as the return value for [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_apsis_t`](#astro_apsis_t) | `apsis` |  An apsis event obtained from a call to [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) or `Astronomy_NextLunarApsis`. See [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) for more details. | 




---

<a name="Astronomy_NextMoonQuarter"></a>
### Astronomy_NextMoonQuarter(mq) &#8658; [`astro_moon_quarter_t`](#astro_moon_quarter_t)

**Continues searching for lunar quarters from a previous search.** 



After calling [`Astronomy_SearchMoonQuarter`](#Astronomy_SearchMoonQuarter), this function can be called one or more times to continue finding consecutive lunar quarters. This function finds the next consecutive moon quarter event after the one passed in as the parameter `mq`.



**Returns:**  If `mq` is valid, this function should always succeed, indicated by the `status` field in the returned structure holding `ASTRO_SUCCESS`. Any other value indicates an internal error, which (after confirming that `mq` is valid) should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues). To be safe, calling code should always check the `status` field for errors. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_moon_quarter_t`](#astro_moon_quarter_t) | `mq` |  A value returned by a prior call to [`Astronomy_SearchMoonQuarter`](#Astronomy_SearchMoonQuarter) or [`Astronomy_NextMoonQuarter`](#Astronomy_NextMoonQuarter). | 




---

<a name="Astronomy_Refraction"></a>
### Astronomy_Refraction(refraction, altitude) &#8658; `double`

**Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.** 





**Returns:**  The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  The option selecting which refraction correction to use. If `REFRACTION_NORMAL`, uses a well-behaved refraction model that works well for all valid values (-90 to +90) of `altitude`. If `REFRACTION_JPLHOR`, this function returns a compatible value with the JPL Horizons tool. If any other value (including `REFRACTION_NONE`), this function returns 0. | 
| `double` | `altitude` |  An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90. | 




---

<a name="Astronomy_RotateVector"></a>
### Astronomy_RotateVector(rotation, vector) &#8658; [`astro_vector_t`](#astro_vector_t)

**Applies a rotation to a vector, yielding a rotated vector.** 



This function transforms a vector in one orientation to a vector in another orientation.



**Returns:**  A vector in the orientation specified by `rotation`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `rotation` |  A rotation matrix that specifies how the orientation of the vector is to be changed. | 
| [`astro_vector_t`](#astro_vector_t) | `vector` |  The vector whose orientation is to be changed. | 




---

<a name="Astronomy_Rotation_ECL_EQD"></a>
### Astronomy_Rotation_ECL_EQD(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of date.



**Returns:**  A rotation matrix that converts ECL to EQD. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the desired equator. | 




---

<a name="Astronomy_Rotation_ECL_EQJ"></a>
### Astronomy_Rotation_ECL_EQJ() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: EQJ = equatorial system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts ECL to EQJ. 



---

<a name="Astronomy_Rotation_ECL_HOR"></a>
### Astronomy_Rotation_ECL_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: HOR = horizontal system.

Use [`Astronomy_HorizonFromVector`](#Astronomy_HorizonFromVector) to convert the return value to a traditional altitude/azimuth pair.



**Returns:**  A rotation matrix that converts ECL to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the desired horizontal orientation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_EQD_ECL"></a>
### Astronomy_Rotation_EQD_ECL(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of date. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQD to ECL. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the source equator. | 




---

<a name="Astronomy_Rotation_EQD_EQJ"></a>
### Astronomy_Rotation_EQD_EQJ(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of the specified date/time. Target: EQJ = equatorial system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQD at `time` to EQJ. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time at which the Earth's equator defines the source orientation. | 




---

<a name="Astronomy_Rotation_EQD_HOR"></a>
### Astronomy_Rotation_EQD_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of the specified date/time. Target: HOR = horizontal system.

Use [`Astronomy_HorizonFromVector`](#Astronomy_HorizonFromVector) to convert the return value to a traditional altitude/azimuth pair.



**Returns:**  A rotation matrix that converts EQD to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time at which the Earth's equator applies. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_EQJ_ECL"></a>
### Astronomy_Rotation_EQJ_ECL() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQJ to ECL. 



---

<a name="Astronomy_Rotation_EQJ_EQD"></a>
### Astronomy_Rotation_EQJ_EQD(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of the specified date/time.



**Returns:**  A rotation matrix that converts EQJ to EQD at `time`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time at which the Earth's equator defines the target orientation. | 




---

<a name="Astronomy_Rotation_EQJ_HOR"></a>
### Astronomy_Rotation_EQJ_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using the equator at the J2000 epoch. Target: HOR = horizontal system.

Use [`Astronomy_HorizonFromVector`](#Astronomy_HorizonFromVector) to convert the return value to a traditional altitude/azimuth pair.



**Returns:**  A rotation matrix that converts EQJ to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the desired horizontal orientation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_HOR_ECL"></a>
### Astronomy_Rotation_HOR_ECL(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts HOR to ECL. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the horizontal observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The location of the horizontal observer. | 




---

<a name="Astronomy_Rotation_HOR_EQD"></a>
### Astronomy_Rotation_HOR_EQD(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Source: EQD = equatorial system, using equator of the specified date/time.



**Returns:**  A rotation matrix that converts HOR to EQD at `time` and for `observer`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time at which the Earth's equator applies. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_HOR_EQJ"></a>
### Astronomy_Rotation_HOR_EQJ(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Source: EQJ = equatorial system, using equator at the J2000 epoch.



**Returns:**  A rotation matrix that converts HOR to EQD at `time` and for `observer`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Search"></a>
### Astronomy_Search(func, context, t1, t2, dt_tolerance_seconds) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for a time at which a function's value increases through zero.** 



Certain astronomy calculations involve finding a time when an event occurs. Often such events can be defined as the root of a function: the time at which the function's value becomes zero.

`Astronomy_Search` finds the *ascending root* of a function: the time at which the function's value becomes zero while having a positive slope. That is, as time increases, the function transitions from a negative value, through zero at a specific moment, to a positive value later. The goal of the search is to find that specific moment.

The search function is specified by two parameters: `func` and `context`. The `func` parameter is a pointer to the function itself, which accepts a time and a context containing any other arguments needed to evaluate the function. The `context` parameter supplies that context for the given search. As an example, a caller may wish to find the moment a celestial body reaches a certain ecliptic longitude. In that case, the caller might create a structure that contains an [`astro_body_t`](#astro_body_t) member to specify the body and a `double` to hold the target longitude. The function would cast the pointer `context` passed in as a pointer to that structure type. It could subtract the target longitude from the actual longitude at a given time; thus the difference would equal zero at the moment in time the planet reaches the desired longitude.

The `func` returns an [`astro_func_result_t`](#astro_func_result_t) structure every time it is called. If the returned structure has a value of `status` other than `ASTRO_SUCCESS`, the search immediately fails and reports that same error code in the `status` returned by `Astronomy_Search`. Otherwise, `status` is `ASTRO_SUCCESS` and `value` is the value of the function, and the search proceeds until it either finds the ascending root or fails for some reason.

The search calls `func` repeatedly to rapidly narrow in on any ascending root within the time window specified by `t1` and `t2`. The search never reports a solution outside this time window.

`Astronomy_Search` uses a combination of bisection and quadratic interpolation to minimize the number of function calls. However, it is critical that the supplied time window be small enough that there cannot be more than one root (ascedning or descending) within it; otherwise the search can fail. Beyond that, it helps to make the time window as small as possible, ideally such that the function itself resembles a smooth parabolic curve within that window.

If an ascending root is not found, or more than one root (ascending and/or descending) exists within the window `t1`..`t2`, the search will fail with status code `ASTRO_SEARCH_FAILURE`.

If the search does not converge within 20 iterations, it will fail with status code `ASTRO_NO_CONVERGE`.



**Returns:**  If successful, the returned structure has `status` equal to `ASTRO_SUCCESS` and `time` set to a value within `dt_tolerance_seconds` of an ascending root. On success, the `time` value will always be in the inclusive range [`t1`, `t2`]. If the search fails, `status` will be set to a value other than `ASTRO_SUCCESS`. See the remarks above for more details. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_search_func_t`](#astro_search_func_t) | `func` |  The function for which to find the time of an ascending root. See remarks above for more details. | 
| `void *` | `context` |  Any ancillary data needed by the function `func` to calculate a value. The data type varies depending on the function passed in. For example, the function may involve a specific celestial body that must be specified somehow. | 
| [`astro_time_t`](#astro_time_t) | `t1` |  The lower time bound of the search window. See remarks above for more details. | 
| [`astro_time_t`](#astro_time_t) | `t2` |  The upper time bound of the search window. See remarks above for more details. | 
| `double` | `dt_tolerance_seconds` |  Specifies an amount of time in seconds within which a bounded ascending root is considered accurate enough to stop. A typical value is 1 second. | 




---

<a name="Astronomy_SearchHourAngle"></a>
### Astronomy_SearchHourAngle(body, observer, hourAngle, startTime) &#8658; [`astro_hour_angle_t`](#astro_hour_angle_t)

**Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.** 



The *hour angle* of a celestial body indicates its position in the sky with respect to the Earth's rotation. The hour angle depends on the location of the observer on the Earth. The hour angle is 0 when the body reaches its highest angle above the horizon in a given day. The hour angle increases by 1 unit for every sidereal hour that passes after that point, up to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates the number of hours that have passed since the most recent time that the body has culminated, or reached its highest point.

This function searches for the next time a celestial body reaches the given hour angle after the date and time specified by `startTime`. To find when a body culminates, pass 0 for `hourAngle`. To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.

Note that, especially close to the Earth's poles, a body as seen on a given day may always be above the horizon or always below the horizon, so the caller cannot assume that a culminating object is visible nor that an object is below the horizon at its minimum altitude.

On success, the function reports the date and time, along with the horizontal coordinates of the body at that time, as seen by the given observer.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the other structure fields are valid. Otherwise, `status` holds some other value that indicates an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body, which can the Sun, the Moon, or any planet other than the Earth. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  Indicates a location on or near the surface of the Earth where the observer is located. Call [`Astronomy_MakeObserver`](#Astronomy_MakeObserver) to create an observer structure. | 
| `double` | `hourAngle` |  An hour angle value in the range [0, 24) indicating the number of sidereal hours after the body's most recent culmination. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 




---

<a name="Astronomy_SearchLunarApsis"></a>
### Astronomy_SearchLunarApsis(startTime) &#8658; [`astro_apsis_t`](#astro_apsis_t)

**Finds the date and time of the Moon's closest distance (perigee) or farthest distance (apogee) with respect to the Earth.** 



Given a date and time to start the search in `startTime`, this function finds the next date and time that the center of the Moon reaches the closest or farthest point in its orbit with respect to the center of the Earth, whichever comes first after `startTime`.

The closest point is called *perigee* and the farthest point is called *apogee*. The word *apsis* refers to either event.

To iterate through consecutive alternating perigee and apogee events, call `Astronomy_SearchLunarApsis` once, then use the return value to call [`Astronomy_NextLunarApsis`](#Astronomy_NextLunarApsis). After that, keep feeding the previous return value from `Astronomy_NextLunarApsis` into another call of `Astronomy_NextLunarApsis` as many times as desired.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`, `time` holds the date and time of the next lunar apsis, `kind` holds either `APSIS_PERICENTER` for perigee or `APSIS_APOCENTER` for apogee, and the distance values `dist_au` (astronomical units) and `dist_km` (kilometers) are valid. If the function fails, `status` holds some value other than `ASTRO_SUCCESS` that indicates what went wrong, and the other structure fields are invalid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start searching for the next perigee or apogee. | 




---

<a name="Astronomy_SearchMaxElongation"></a>
### Astronomy_SearchMaxElongation(body, startTime) &#8658; [`astro_elongation_t`](#astro_elongation_t)

**Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.** 



Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is. Mercury especially is almost always impossible to see because it gets lost in the Sun's glare. The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through a telescope without atmospheric interference, are when these planets reach maximum elongation. These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

This function solves for those times, reporting the next maximum elongation event's date and time, the elongation value itself, the relative longitude with the Sun, and whether the planet is best observed in the morning or evening. See [`Astronomy_Elongation`](#Astronomy_Elongation) for more details about the returned structure.



**Returns:**  If successful, the `status` field of the returned structure will be `ASTRO_SUCCESS` and the other structure fields will be valid. Otherwise, `status` will contain some other value indicating an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  Either `BODY_MERCURY` or `BODY_VENUS`. Any other value will fail with the error `ASTRO_INVALID_BODY`. To find the best viewing opportunites for planets farther from the Sun than the Earth is (Mars through Pluto) use [`Astronomy_SearchRelativeLongitude`](#Astronomy_SearchRelativeLongitude) to find the next opposition event. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to begin the search. The maximum elongation event found will always be the first one that occurs after this date and time. | 




---

<a name="Astronomy_SearchMoonPhase"></a>
### Astronomy_SearchMoonPhase(targetLon, startTime, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for the time that the Moon reaches a specified phase.** 



Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic longitude with respect to the Sun's geocentric ecliptic longitude. When the Moon and the Sun have the same longitude, that is defined as a new moon. When their longitudes are 180 degrees apart, that is defined as a full moon.

This function searches for any value of the lunar phase expressed as an angle in degrees in the range [0, 360).

If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter) it is much easier to call the functions [`Astronomy_SearchMoonQuarter`](#Astronomy_SearchMoonQuarter) and [`Astronomy_NextMoonQuarter`](#Astronomy_NextMoonQuarter). This function is useful for finding general phase angles outside those four quarters.



**Returns:**  On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the `time` field holds the date and time when the Moon reaches the target longitude. On failure, `status` holds some other value as an error code. One possible error code is `ASTRO_NO_MOON_QUARTER` if `startTime` and `limitDays` do not enclose the desired event. See remarks in [`Astronomy_Search`](#Astronomy_Search) for other possible error codes. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `targetLon` |  The difference in geocentric longitude between the Sun and Moon that specifies the lunar phase being sought. This can be any value in the range [0, 360). Certain values have conventional names: 0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The beginning of the time window in which to search for the Moon reaching the specified phase. | 
| `double` | `limitDays` |  The number of days after `startTime` that limits the time window for the search. | 




---

<a name="Astronomy_SearchMoonQuarter"></a>
### Astronomy_SearchMoonQuarter(startTime) &#8658; [`astro_moon_quarter_t`](#astro_moon_quarter_t)

**Finds the first lunar quarter after the specified date and time.** 



A lunar quarter is one of the following four lunar phase events: new moon, first quarter, full moon, third quarter. This function finds the lunar quarter that happens soonest after the specified date and time.

To continue iterating through consecutive lunar quarters, call this function once, followed by calls to [`Astronomy_NextMoonQuarter`](#Astronomy_NextMoonQuarter) as many times as desired.



**Returns:**  This function should always succeed, indicated by the `status` field in the returned structure holding `ASTRO_SUCCESS`. Any other value indicates an internal error, which should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues). To be safe, calling code should always check the `status` field for errors. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 




---

<a name="Astronomy_SearchPeakMagnitude"></a>
### Astronomy_SearchPeakMagnitude(body, startTime) &#8658; [`astro_illum_t`](#astro_illum_t)

**Searches for the date and time Venus will next appear brightest as seen from the Earth.** 



This function searches for the date and time Venus appears brightest as seen from the Earth. Currently only Venus is supported for the `body` parameter, though this could change in the future. Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth, so peak magnitude events have little practical value for that planet. Planets other than Venus and Mercury reach peak magnitude at opposition, which can be found using [`Astronomy_SearchRelativeLongitude`](#Astronomy_SearchRelativeLongitude). The Moon reaches peak magnitude at full moon, which can be found using [`Astronomy_SearchMoonQuarter`](#Astronomy_SearchMoonQuarter) or [`Astronomy_SearchMoonPhase`](#Astronomy_SearchMoonPhase). The Sun reaches peak magnitude at perihelion, which occurs each year in January. However, the difference is minor and has little practical value.



**Returns:**  See documentation about the return value from [`Astronomy_Illumination`](#Astronomy_Illumination). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  Currently only `BODY_VENUS` is allowed. Any other value results in the error `ASTRO_INVALID_BODY`. See remarks above for more details. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time to start searching for the next peak magnitude event. | 




---

<a name="Astronomy_SearchRelativeLongitude"></a>
### Astronomy_SearchRelativeLongitude(body, targetRelLon, startTime) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for the time when the Earth and another planet are separated by a specified angle in ecliptic longitude, as seen from the Sun.** 



A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around that plane from one body to the other in the direction the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees. This is the relative longitude.

Given a planet other than the Earth in `body` and a time to start the search in `startTime`, this function searches for the next time that the relative longitude measured from the planet to the Earth is `targetRelLon`.

Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:



- When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*, a time when the other planet becomes very difficult to see because of being lost in the Sun's glare. (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing between the Earth and the Sun.)
- When the relative longitude is 0 degrees and the other planet orbits farther from the Sun, this is known as *opposition*. Opposition is when the planet is closest to the Earth, and also when it is visible for most of the night, so it is considered the best time to observe the planet.
- When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.




**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and `time` will hold the date and time of the relative longitude event. Otherwise `status` will hold some other value that indicates an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A planet other than the Earth. If `body` is not a planet other than the Earth, an error occurs. | 
| `double` | `targetRelLon` |  The desired relative longitude, expressed in degrees. Must be in the range [0, 360). | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to begin the search. | 




---

<a name="Astronomy_SearchRiseSet"></a>
### Astronomy_SearchRiseSet(body, observer, direction, startTime, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.** 



This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth. Rise time is when the body first starts to be visible above the horizon. For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon. Set time is the moment when the body appears to vanish below the horizon.

This function corrects for typical atmospheric refraction, which causes celestial bodies to appear higher above the horizon than they would if the Earth had no atmosphere. It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).

Note that rise or set may not occur in every 24 hour period. For example, near the Earth's poles, there are long periods of time where the Sun stays below the horizon, never rising. Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day. This is because the Moon sets nearly an hour later each day due to orbiting the Earth a significant amount during each rotation of the Earth. Therefore callers must not assume that the function will always succeed.



**Returns:**  On success, the `status` field in the returned structure contains `ASTRO_SUCCESS` and the `time` field contains the date and time of the rise or set time as requested. If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the rise or set event does not occur within `limitDays` days of `startTime`. This is a normal condition, not an error. Any other value of `status` indicates an error of some kind. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Moon, or any planet other than the Earth. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The location where observation takes place. You can create an observer structure by calling [`Astronomy_MakeObserver`](#Astronomy_MakeObserver). | 
| [`astro_direction_t`](#astro_direction_t) | `direction` |  Either `DIRECTION_RISE` to find a rise time or `DIRECTION_SET` to find a set time. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 
| `double` | `limitDays` |  Limits how many days to search for a rise or set time. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. | 




---

<a name="Astronomy_SearchSunLongitude"></a>
### Astronomy_SearchSunLongitude(targetLon, startTime, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.** 



This function finds the moment in time, if any exists in the given time window, that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

This function can be used to determine equinoxes and solstices. However, it is usually more convenient and efficient to call [`Astronomy_Seasons`](#Astronomy_Seasons) to calculate all equinoxes and solstices for a given calendar year.

The function searches the window of time specified by `startTime` and `startTime+limitDays`. The search will return an error if the Sun never reaches the longitude `targetLon` or if the window is so large that the longitude ranges more than 180 degrees within it. It is recommended to keep the window smaller than 10 days when possible.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the `time` field will contain the date and time the Sun reaches the target longitude. Any other value indicates an error. See remarks in [`Astronomy_Search`](#Astronomy_Search) (which this function calls) for more information about possible error codes. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `targetLon` |  The desired ecliptic longitude in degrees, relative to the true equinox of date. This may be any value in the range [0, 360), although certain values have conventional meanings: 0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for the desired longitude event. | 
| `double` | `limitDays` |  The real-valued number of days, which when added to `startTime`, limits the range of time over which the search looks. It is recommended to keep this value between 1 and 10 days. See remarks above for more details. | 




---

<a name="Astronomy_Seasons"></a>
### Astronomy_Seasons(year) &#8658; [`astro_seasons_t`](#astro_seasons_t)

**Finds both equinoxes and both solstices for a given calendar year.** 



The changes of seasons are defined by solstices and equinoxes. Given a calendar year number, this function calculates the March and September equinoxes and the June and December solstices.

The equinoxes are the moments twice each year when the plane of the Earth's equator passes through the center of the Sun. In other words, the Sun's declination is zero at both equinoxes. The March equinox defines the beginning of spring in the northern hemisphere and the beginning of autumn in the southern hemisphere. The September equinox defines the beginning of autumn in the northern hemisphere and the beginning of spring in the southern hemisphere.

The solstices are the moments twice each year when one of the Earth's poles is most tilted toward the Sun. More precisely, the Sun's declination reaches its minimum value at the December solstice, which defines the beginning of winter in the northern hemisphere and the beginning of summer in the southern hemisphere. The Sun's declination reaches its maximum value at the June solstice, which defines the beginning of summer in the northern hemisphere and the beginning of winter in the southern hemisphere.



**Returns:**  The times of the four seasonal changes in the given calendar year. This function should always succeed. However, to be safe, callers should check the `status` field of the returned structure to make sure it contains `ASTRO_SUCCESS`. Any failures indicate a bug in the algorithm and should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues). 



| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` |  The calendar year number for which to calculate equinoxes and solstices. The value may be any integer, but only the years 1800 through 2100 have been validated for accuracy: unit testing against data from the United States Naval Observatory confirms that all equinoxes and solstices for that range of years are within 2 minutes of the correct time. | 




---

<a name="Astronomy_SphereFromVector"></a>
### Astronomy_SphereFromVector(vector) &#8658; [`astro_spherical_t`](#astro_spherical_t)

**Converts Cartesian coordinates to spherical coordinates.** 



Given a Cartesian vector, returns latitude, longitude, and distance.



**Returns:**  Spherical coordinates that are equivalent to the given vector. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `vector` |  Cartesian vector to be converted to spherical coordinates. | 




---

<a name="Astronomy_SunPosition"></a>
### Astronomy_SunPosition(time) &#8658; [`astro_ecliptic_t`](#astro_ecliptic_t)

**Calculates geocentric ecliptic coordinates for the Sun.** 



This function calculates the position of the Sun as seen from the Earth. The returned value includes both Cartesian and spherical coordinates. The x-coordinate and longitude values in the returned structure are based on the *true equinox of date*: one of two points in the sky where the instantaneous plane of the Earth's equator at the given date and time (the *equatorial plane*) intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*). By convention, the apparent location of the Sun at the March equinox is chosen as the longitude origin and x-axis direction, instead of the one for September.

`Astronomy_SunPosition` corrects for precession and nutation of the Earth's axis in order to obtain the exact equatorial plane at the given time.

This function can be used for calculating changes of seasons: equinoxes and solstices. In fact, the function [`Astronomy_Seasons`](#Astronomy_Seasons) does use this function for that purpose.



**Returns:**  The ecliptic coordinates of the Sun using the Earth's true equator of date. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the Sun's position. | 




---

<a name="Astronomy_TimeFromDays"></a>
### Astronomy_TimeFromDays(ut) &#8658; [`astro_time_t`](#astro_time_t)

**Converts a J2000 day value to an [`astro_time_t`](#astro_time_t) value.** 



This function can be useful for reproducing an [`astro_time_t`](#astro_time_t) structure from its `ut` field only.



**Returns:**  An [`astro_time_t`](#astro_time_t) value for the given `ut` value. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` |  The floating point number of days since noon UTC on January 1, 2000. | 




---

<a name="Astronomy_TimeFromUtc"></a>
### Astronomy_TimeFromUtc(utc) &#8658; [`astro_time_t`](#astro_time_t)

**Creates an [`astro_time_t`](#astro_time_t) value from a given calendar date and time.** 



This function is similar to [`Astronomy_MakeTime`](#Astronomy_MakeTime), only it receives a UTC calendar date and time in the form of an [`astro_utc_t`](#astro_utc_t) structure instead of as separate numeric parameters. Astronomy_TimeFromUtc is the inverse of [`Astronomy_UtcFromTime`](#Astronomy_UtcFromTime).



**Returns:**  A value that can be used for astronomical calculations for the given date and time. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_utc_t`](#astro_utc_t) | `utc` |  The UTC calendar date and time to be converted to [`astro_time_t`](#astro_time_t).  | 




---

<a name="Astronomy_UtcFromTime"></a>
### Astronomy_UtcFromTime(time) &#8658; [`astro_utc_t`](#astro_utc_t)

**Determines the calendar year, month, day, and time from an [`astro_time_t`](#astro_time_t) value.** 



After calculating the date and time of an astronomical event in the form of an [`astro_time_t`](#astro_time_t) value, it is often useful to display the result in a human-readable form. This function converts the linear time scales in the `ut` field of [`astro_time_t`](#astro_time_t) into a calendar date and time: year, month, day, hours, minutes, and seconds, expressed in UTC.



**Returns:**  A date and time broken out into conventional year, month, day, hour, minute, and second. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The astronomical time value to be converted to calendar date and time.  | 




---

<a name="Astronomy_VectorFromSphere"></a>
### Astronomy_VectorFromSphere(sphere, time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Converts spherical coordinates to Cartesian coordinates.** 



Given spherical coordinates and a time at which they are valid, returns a vector of Cartesian coordinates. The returned value includes the time, as required by the type [`astro_vector_t`](#astro_vector_t).



**Returns:**  The vector form of the supplied spherical coordinates. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_spherical_t`](#astro_spherical_t) | `sphere` |  Spherical coordinates to be converted. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The time that should be included in the return value. | 




---

<a name="Astronomy_VectorLength"></a>
### Astronomy_VectorLength(vector) &#8658; `double`

**Calculates the length of the given vector.** 



Calculates the non-negative length of the given vector. The length is expressed in the same units as the vector's components, usually astronomical units (AU).



**Returns:**  The length of the vector. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `vector` |  The vector whose length is to be calculated.  | 



<a name="enums"></a>
## Enumerated Types



---

<a name="astro_aberration_t"></a>
### `astro_aberration_t`

**Aberration calculation options.** 



[Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect causing the apparent direction of an observed body to be shifted due to transverse movement of the Earth with respect to the rays of light coming from that body. This angular correction can be anywhere from 0 to about 20 arcseconds, depending on the position of the observed body relative to the instantaneous velocity vector of the Earth.

Some Astronomy Engine functions allow optional correction for aberration by passing in a value of this enumerated type.

Aberration correction is useful to improve accuracy of coordinates of apparent locations of bodies seen from the Earth. However, because aberration affects not only the observed body (such as a planet) but the surrounding stars, aberration may be unhelpful (for example) for determining exactly when a planet crosses from one constellation to another. 

| Enum Value | Description |
| --- | --- |
| `ABERRATION` |  Request correction for aberration.  |
| `NO_ABERRATION` |  Do not correct for aberration.  |



---

<a name="astro_apsis_kind_t"></a>
### `astro_apsis_kind_t`

**The type of apsis: pericenter (closest approach) or apocenter (farthest distance).** 



| Enum Value | Description |
| --- | --- |
| `APSIS_PERICENTER` |  The body is at its closest approach to the object it orbits.  |
| `APSIS_APOCENTER` |  The body is at its farthest distance from the object it orbits.  |
| `APSIS_INVALID` |  Undefined or invalid apsis.  |



---

<a name="astro_body_t"></a>
### `astro_body_t`

**A celestial body.** 



| Enum Value | Description |
| --- | --- |
| `BODY_INVALID` |  An invalid or undefined celestial body.  |
| `BODY_MERCURY` |  Mercury  |
| `BODY_VENUS` |  Venus  |
| `BODY_EARTH` |  Earth  |
| `BODY_MARS` |  Mars  |
| `BODY_JUPITER` |  Jupiter  |
| `BODY_SATURN` |  Saturn  |
| `BODY_URANUS` |  Uranus  |
| `BODY_NEPTUNE` |  Neptune  |
| `BODY_PLUTO` |  Pluto  |
| `BODY_SUN` |  Sun  |
| `BODY_MOON` |  Moon  |



---

<a name="astro_direction_t"></a>
### `astro_direction_t`

**Selects whether to search for a rise time or a set time.** 



The [`Astronomy_SearchRiseSet`](#Astronomy_SearchRiseSet) function finds the rise or set time of a body depending on the value of its `direction` parameter. 

| Enum Value | Description |
| --- | --- |
| `DIRECTION_RISE` |  Search for the time a body begins to rise above the horizon.  |
| `DIRECTION_SET` |  Search for the time a body finishes sinking below the horizon.  |



---

<a name="astro_equator_date_t"></a>
### `astro_equator_date_t`

**Selects the date on which the Earth's equator to be used for representing equatorial coordinates.** 



The Earth's equator is not always in the same plane due to precession and nutation.

Sometimes it is useful to have a fixed plane of reference for equatorial coordinates across different calendar dates. In these cases, a fixed *epoch*, or reference time, is helpful. Astronomy Engine provides the J2000 epoch for such cases. This refers to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.

For some other purposes, it is more helpful to represent coordinates using the Earth's equator exactly as it is on that date. For example, when calculating rise/set times or horizontal coordinates, it is most accurate to use the orientation of the Earth's equator at that same date and time. For these uses, Astronomy Engine allows *of-date* calculations. 

| Enum Value | Description |
| --- | --- |
| `EQUATOR_J2000` |  Represent equatorial coordinates in the J2000 epoch.  |
| `EQUATOR_OF_DATE` |  Represent equatorial coordinates using the Earth's equator at the given date and time.  |



---

<a name="astro_refraction_t"></a>
### `astro_refraction_t`

**Selects whether to correct for atmospheric refraction, and if so, how.** 



| Enum Value | Description |
| --- | --- |
| `REFRACTION_NONE` |  No atmospheric refraction correction (airless).  |
| `REFRACTION_NORMAL` |  Recommended correction for standard atmospheric refraction.  |
| `REFRACTION_JPLHOR` |  Used only for compatibility testing with JPL Horizons online tool.  |



---

<a name="astro_status_t"></a>
### `astro_status_t`

**Indicates success/failure of an Astronomy Engine function call.** 



| Enum Value | Description |
| --- | --- |
| `ASTRO_SUCCESS` |  The operation was successful.  |
| `ASTRO_NOT_INITIALIZED` |  A placeholder that can be used for data that is not yet initialized.  |
| `ASTRO_INVALID_BODY` |  The celestial body was not valid. Different sets of bodies are supported depending on the function.  |
| `ASTRO_NO_CONVERGE` |  A numeric solver failed to converge. This should not happen unless there is a bug in Astronomy Engine.  |
| `ASTRO_BAD_TIME` |  Cannot calculate Pluto's position outside the year range 1700..2200.  |
| `ASTRO_BAD_VECTOR` |  Vector magnitude is too small to be normalized into a unit vector.  |
| `ASTRO_SEARCH_FAILURE` |  Search was not able to find an ascending root crossing of the function in the specified time interval.  |
| `ASTRO_EARTH_NOT_ALLOWED` |  The Earth cannot be treated as a celestial body seen from an observer on the Earth itself.  |
| `ASTRO_NO_MOON_QUARTER` |  No lunar quarter occurs inside the specified time range.  |
| `ASTRO_WRONG_MOON_QUARTER` |  Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter.  |
| `ASTRO_INTERNAL_ERROR` |  A self-check failed inside the code somewhere, indicating a bug needs to be fixed.  |
| `ASTRO_INVALID_PARAMETER` |  A parameter value passed to a function was not valid.  |



---

<a name="astro_visibility_t"></a>
### `astro_visibility_t`

**Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.** 



| Enum Value | Description |
| --- | --- |
| `VISIBLE_MORNING` |  The body is best visible in the morning, before sunrise.  |
| `VISIBLE_EVENING` |  The body is best visible in the evening, after sunset.  |


<a name="structs"></a>
## Structures



---

<a name="astro_angle_result_t"></a>
### `astro_angle_result_t`

**An angular value expressed in degrees.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `angle` |  An angle expressed in degrees.  |


---

<a name="astro_apsis_t"></a>
### `astro_apsis_t`

**An apsis event: pericenter (closest approach) or apocenter (farthest distance).** 



For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an event where the orbiting body reaches its closest or farthest point from the primary body. The closest approach is called *pericenter* and the farthest point is *apocenter*.

More specific terminology is common for particular orbiting bodies. The Moon's closest approach to the Earth is called *perigee* and its furthest point is called *apogee*. The closest approach of a planet to the Sun is called *perihelion* and the furthest point is called *aphelion*.

This data structure is returned by [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) and [`Astronomy_NextLunarApsis`](#Astronomy_NextLunarApsis) to iterate through consecutive alternating perigees and apogees. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the apsis.  |
| [`astro_apsis_kind_t`](#astro_apsis_kind_t) | `kind` |  Whether this is a pericenter or apocenter event.  |
| `double` | `dist_au` |  The distance between the centers of the bodies in astronomical units.  |
| `double` | `dist_km` |  The distance between the centers of the bodies in kilometers.  |


---

<a name="astro_ecliptic_t"></a>
### `astro_ecliptic_t`

**Ecliptic angular and Cartesian coordinates.** 



Coordinates of a celestial body as seen from the center of the Sun (heliocentric), oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic). 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `ex` |  Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane.  |
| `double` | `ey` |  Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox.  |
| `double` | `ez` |  Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north.  |
| `double` | `elat` |  Latitude in degrees north (positive) or south (negative) of the ecliptic plane.  |
| `double` | `elon` |  Longitude in degrees around the ecliptic plane prograde from the equinox.  |


---

<a name="astro_elongation_t"></a>
### `astro_elongation_t`

**Contains information about the visibility of a celestial body at a given date and time. See [`Astronomy_Elongation`](#Astronomy_Elongation) for more detailed information about the members of this structure. See also [`Astronomy_SearchMaxElongation`](#Astronomy_SearchMaxElongation) for how to search for maximum elongation events.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation.  |
| [`astro_visibility_t`](#astro_visibility_t) | `visibility` |  Whether the body is best seen in the morning or the evening.  |
| `double` | `elongation` |  The angle in degrees between the body and the Sun, as seen from the Earth.  |
| `double` | `ecliptic_separation` |  The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.  |


---

<a name="astro_equatorial_t"></a>
### `astro_equatorial_t`

**Equatorial angular coordinates.** 



Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context), oriented with respect to the projection of the Earth's equator onto the sky. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `ra` |  right ascension in sidereal hours.  |
| `double` | `dec` |  declination in degrees  |
| `double` | `dist` |  distance to the celestial body in AU.  |


---

<a name="astro_func_result_t"></a>
### `astro_func_result_t`

**A real value returned by a function whose ascending root is to be found.** 



When calling [`Astronomy_Search`](#Astronomy_Search), the caller must pass in a callback function compatible with the function-pointer type [`astro_search_func_t`](#astro_search_func_t) whose ascending root is to be found. That callback function must return [`astro_func_result_t`](#astro_func_result_t). If the function call is successful, it will set `status` to [`ASTRO_SUCCESS`](#ASTRO_SUCCESS) and `value` to the numeric value appropriate for the given date and time. If the call fails for some reason, it should set `status` to an appropriate error value other than `ASTRO_SUCCESS`; in the error case, to guard against any possible misuse of `value`, it is recommended to set `value` to `NAN`, though this is not strictly necessary. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `value` |  The value returned by a function whose ascending root is to be found.  |


---

<a name="astro_horizon_t"></a>
### `astro_horizon_t`

**Coordinates of a celestial body as seen by a topocentric observer.** 



Contains horizontal and equatorial coordinates seen by an observer on or near the surface of the Earth (a topocentric observer). Optionally corrected for atmospheric refraction. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `azimuth` |  Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West.  |
| `double` | `altitude` |  Angle in degrees above (positive) or below (negative) the observer's horizon.  |
| `double` | `ra` |  Right ascension in sidereal hours.  |
| `double` | `dec` |  Declination in degrees.  |


---

<a name="astro_hour_angle_t"></a>
### `astro_hour_angle_t`

**Information about a celestial body crossing a specific hour angle.** 



Returned by the function [`Astronomy_SearchHourAngle`](#Astronomy_SearchHourAngle) to report information about a celestial body crossing a certain hour angle as seen by a specified topocentric observer. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time when the body crosses the specified hour angle.  |
| [`astro_horizon_t`](#astro_horizon_t) | `hor` |  Apparent coordinates of the body at the time it crosses the specified hour angle.  |


---

<a name="astro_illum_t"></a>
### `astro_illum_t`

**Information about the brightness and illuminated shape of a celestial body.** 



Returned by the functions [`Astronomy_Illumination`](#Astronomy_Illumination) and [`Astronomy_SearchPeakMagnitude`](#Astronomy_SearchPeakMagnitude) to report the visual magnitude and illuminated fraction of a celestial body at a given date and time. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation.  |
| `double` | `mag` |  The visual magnitude of the body. Smaller values are brighter.  |
| `double` | `phase_angle` |  The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth.  |
| `double` | `helio_dist` |  The distance between the Sun and the body at the observation time.  |
| `double` | `ring_tilt` |  For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.  |


---

<a name="astro_moon_quarter_t"></a>
### `astro_moon_quarter_t`

**A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `int` | `quarter` |  0=new moon, 1=first quarter, 2=full moon, 3=third quarter.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the lunar quarter.  |


---

<a name="astro_observer_t"></a>
### `astro_observer_t`

**The location of an observer on (or near) the surface of the Earth.** 



This structure is passed to functions that calculate phenomena as observed from a particular place on the Earth.

You can create this structure directly, or you can call the convenience function [`Astronomy_MakeObserver`](#Astronomy_MakeObserver)# to create one for you. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `latitude` |  Geographic latitude in degrees north (positive) or south (negative) of the equator.  |
| `double` | `longitude` |  Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.  |
| `double` | `height` |  The height above (positive) or below (negative) sea level, expressed in meters.  |


---

<a name="astro_rotation_t"></a>
### `astro_rotation_t`

**Contains a rotation matrix that can be used to transform one coordinate system to another.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `rot` |  A normalized 3x3 rotation matrix.  |


---

<a name="astro_search_result_t"></a>
### `astro_search_result_t`

**The result of a search for an astronomical event.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The time at which a searched-for event occurs.  |


---

<a name="astro_seasons_t"></a>
### `astro_seasons_t`

**The dates and times of changes of season for a given calendar year. Call [`Astronomy_Seasons`](#Astronomy_Seasons) to calculate this data structure for a given year.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `mar_equinox` |  The date and time of the March equinox for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `jun_solstice` |  The date and time of the June soltice for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `sep_equinox` |  The date and time of the September equinox for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `dec_solstice` |  The date and time of the December solstice for the specified year.  |


---

<a name="astro_spherical_t"></a>
### `astro_spherical_t`

**Spherical coordinates: latitude, longitude, distance.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `lat` |  The latitude angle: -90..+90 degrees.  |
| `double` | `lon` |  The longitude angle: 0..360 degrees.  |
| `double` | `dist` |  Distance in AU.  |


---

<a name="astro_time_t"></a>
### `astro_time_t`

**A date and time used for astronomical calculations.** 



This type is of fundamental importance to Astronomy Engine. It is used to represent dates and times for all astronomical calculations. It is also included in the values returned by many Astronomy Engine functions.

To create a valid [`astro_time_t`](#astro_time_t) value from scratch, call [`Astronomy_MakeTime`](#Astronomy_MakeTime) (for a given calendar date and time) or [`Astronomy_CurrentTime`](#Astronomy_CurrentTime) (for the system's current date and time).

To adjust an existing [`astro_time_t`](#astro_time_t) by a certain real number of days, call [`Astronomy_AddDays`](#Astronomy_AddDays).

The [`astro_time_t`](#astro_time_t) type contains `ut` to represent Universal Time (UT1/UTC) and `tt` to represent Terrestrial Time (TT, also known as *ephemeris time*). The difference `tt-ut` is known as *&Delta;T*, and is obtained from a model provided by the [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).

Both `tt` and `ut` are necessary for performing different astronomical calculations. Indeed, certain calculations (such as rise/set times) require both time scales. See the documentation for the `ut` and `tt` fields for more detailed information.

In cases where [`astro_time_t`](#astro_time_t) is included in a structure returned by a function that can fail, the astro_status_t field `status` will contain a value other than [`ASTRO_SUCCESS`](#ASTRO_SUCCESS); in that case the `ut` and `tt` will hold `NAN` (not a number). In general, when there is an error code stored in a struct field `status`, the caller should ignore all other values in that structure, including the `ut` and `tt` inside [`astro_time_t`](#astro_time_t). 

| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `ut` | **UT1/UTC number of days since noon on January 1, 2000.**  The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine. Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed. UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed. The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time. Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*.  |
| `double` | `tt` | **Terrestrial Time days since noon on January 1, 2000.**  Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation. The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth. Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).  |
| `double` | `psi` | **For internal use only. Used to optimize Earth tilt calculations.**  |
| `double` | `eps` | **For internal use only. Used to optimize Earth tilt calculations.**  |


---

<a name="astro_utc_t"></a>
### `astro_utc_t`

**A calendar date and time expressed in UTC.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| `int` | `year` |  The year value, e.g. 2019.  |
| `int` | `month` |  The month value: 1=January, 2=February, ..., 12=December.  |
| `int` | `day` |  The day of the month in the range 1..31.  |
| `int` | `hour` |  The hour of the day in the range 0..23.  |
| `int` | `minute` |  The minute of the hour in the range 0..59.  |
| `double` | `second` |  The floating point number of seconds in the range [0,60).  |


---

<a name="astro_vector_t"></a>
### `astro_vector_t`

**A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `x` |  The Cartesian x-coordinate of the vector in AU.  |
| `double` | `y` |  The Cartesian y-coordinate of the vector in AU.  |
| `double` | `z` |  The Cartesian z-coordinate of the vector in AU.  |
| [`astro_time_t`](#astro_time_t) | `t` |  The date and time at which this vector is valid.  |

<a name="typedefs"></a>
## Type Definitions



---

<a name="astro_search_func_t"></a>
### `astro_search_func_t`

`typedef astro_func_result_t(*  astro_search_func_t) (void *context, astro_time_t time);`

**A pointer to a function that is to be passed as a callback to [`Astronomy_Search`](#Astronomy_Search).** 



The function [`Astronomy_Search`](#Astronomy_Search) numerically solves for the time that a given event occurs. An event is defined as the time when an arbitrary function transitions between having a negative value and a non-negative value. This transition is called an *ascending root*.

The type astro_search_func_t represents such a callback function that accepts a custom `context` pointer and an [`astro_time_t`](#astro_time_t) representing the time to probe. The function returns an [`astro_func_result_t`](#astro_func_result_t) that contains either a real number in `value` or an error code in `status` that aborts the search.

The `context` points to some data whose type varies depending on the callback function. It can contain any auxiliary parameters (other than time) needed to evaluate the function. For example, a function may pertain to a specific celestial body, in which case `context` may point to a value of type astro_body_t. The `context` parameter is supplied by the caller of [`Astronomy_Search`](#Astronomy_Search), which passes it along to every call to the callback function. If the caller of `Astronomy_Search` knows that the callback function does not need a context, it is safe to pass `NULL` as the context pointer. 