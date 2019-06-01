# Astronomy Engine (C/C++)

This is the complete programming reference for the C version of 
[Astronomy Engine](../../). It can be used directly from C++ programs also.
Other programming languages are supported. See the [home page](../../) for more info.

---

## Quick Start
To get started quickly, here are some [examples](../../demo/c/).

---

## Topic Index

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#Astronomy_HelioVector) | Calculates vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy_GeoVector)     | Calculates vector with respect to the center of the Earth. |
| [Equator](#Astronomy_Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Astronomy_Ecliptic)       | Calculates ecliptic latitude, longitude, and Cartesian coordinates. |
| [Horizon](#Astronomy_Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |

### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#Astronomy_SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchHourAngle](#Astronomy_SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#Astronomy_MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonQuarter](#Astronomy_SearchMoonQuarter) | Find the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy_NextMoonQuarter) | Find the next quarter moon phase after a previous one that has been found. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#Astronomy_SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy_NextLunarApsis) | Given an already-found apsis, find the next perigee or apogee of the Moon. |

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
| [SearchRelativeLongitude](#Astronomy_SearchRelativeLongitude) | Find oppositions and conjunctions of planets. |

### Equinoxes and solstices

| Function | Description |
| -------- | ----------- |
| [Seasons](#Astronomy_Seasons) | Finds the equinoxes and solstices for a given calendar year. |

---


## Functions



---

<a name="Astronomy_AddDays"></a>
### Astronomy_AddDays(time, days) &#8658; [`astro_time_t`](#astro_time_t)

**Calculates the sum or difference of an [`astro_time_t`](#astro_time_t) with a specified floating point number of days.** 



Sometimes we need to adjust a given [`astro_time_t`](#astro_time_t) value by a certain amount of time. This function adds the given real number of days in `days` to the date and time in `time`.

More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time, according to the historical and predictive Delta-T model provided by the [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).

The value stored in `time` will not be modified; it is passed by value.



**Returns:**  A date and time that is conceptually equal to `time + days`. 



---

<a name="Astronomy_AngleFromSun"></a>
### Astronomy_AngleFromSun(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)



---

<a name="Astronomy_BodyCode"></a>
### Astronomy_BodyCode(name) &#8658; [`astro_body_t`](#astro_body_t)

**Returns the [`astro_body_t`](#astro_body_t) value corresponding to the given English name.** 





**Returns:**  If `name` is one of the strings (case-sensitive) listed above, the returned value is the corresponding [`astro_body_t`](#astro_body_t) value, otherwise it is [`BODY_INVALID`](#BODY_INVALID). 



---

<a name="Astronomy_BodyName"></a>
### Astronomy_BodyName(body) &#8658; `const char *`

**Finds the name of a celestial body.** 





**Returns:**  The English-language name of the celestial body, or "" if the body is not valid. 



---

<a name="Astronomy_CurrentTime"></a>
### Astronomy_CurrentTime() &#8658; [`astro_time_t`](#astro_time_t)

**Returns the computer's current date and time in the form of an [`astro_time_t`](#astro_time_t).** 



Uses the computer's system clock to find the current UTC date and time with 1-second granularity. Converts that date and time to an [`astro_time_t`](#astro_time_t) value and returns the result. Callers can pass this value to other Astronomy Engine functions to calculate current observational conditions. 

---

<a name="Astronomy_Ecliptic"></a>
### Astronomy_Ecliptic(equ) &#8658; [`astro_ecliptic_t`](#astro_ecliptic_t)



---

<a name="Astronomy_EclipticLongitude"></a>
### Astronomy_EclipticLongitude(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)



---

<a name="Astronomy_Elongation"></a>
### Astronomy_Elongation(body, time) &#8658; [`astro_elongation_t`](#astro_elongation_t)



---

<a name="Astronomy_Equator"></a>
### Astronomy_Equator(body, time, observer, ofdate, aberration) &#8658; [`astro_equatorial_t`](#astro_equatorial_t)



---

<a name="Astronomy_GeoMoon"></a>
### Astronomy_GeoMoon(time) &#8658; [`astro_vector_t`](#astro_vector_t)



---

<a name="Astronomy_GeoVector"></a>
### Astronomy_GeoVector(body, time, correct_aberration) &#8658; [`astro_vector_t`](#astro_vector_t)



---

<a name="Astronomy_HelioVector"></a>
### Astronomy_HelioVector(body, time) &#8658; [`astro_vector_t`](#astro_vector_t)



---

<a name="Astronomy_Horizon"></a>
### Astronomy_Horizon(time, observer, ra, dec, refraction) &#8658; [`astro_horizon_t`](#astro_horizon_t)



---

<a name="Astronomy_Illumination"></a>
### Astronomy_Illumination(body, time) &#8658; [`astro_illum_t`](#astro_illum_t)



---

<a name="Astronomy_LongitudeFromSun"></a>
### Astronomy_LongitudeFromSun(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)



---

<a name="Astronomy_MakeObserver"></a>
### Astronomy_MakeObserver(latitude, longitude, height) &#8658; [`astro_observer_t`](#astro_observer_t)



---

<a name="Astronomy_MakeTime"></a>
### Astronomy_MakeTime(year, month, day, hour, minute, second) &#8658; [`astro_time_t`](#astro_time_t)

**Creates an [`astro_time_t`](#astro_time_t) value from a given calendar date and time.** 



Given a UTC calendar date and time, calculates an [`astro_time_t`](#astro_time_t) value that can be passed to other Astronomy Engine functions for performing various calculations relating to that date and time.

It is the caller's responsibility to ensure that the parameter values are correct. The parameters are not checked for validity, and this function never returns any indication of an error. Invalid values, for example passing in February 31, may cause unexpected return values.



**Returns:**  An [`astro_time_t`](#astro_time_t) value that represents the given calendar date and time. 



---

<a name="Astronomy_MoonPhase"></a>
### Astronomy_MoonPhase(time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)



---

<a name="Astronomy_NextLunarApsis"></a>
### Astronomy_NextLunarApsis(apsis) &#8658; [`astro_apsis_t`](#astro_apsis_t)



---

<a name="Astronomy_NextMoonQuarter"></a>
### Astronomy_NextMoonQuarter(mq) &#8658; [`astro_moon_quarter_t`](#astro_moon_quarter_t)



---

<a name="Astronomy_Search"></a>
### Astronomy_Search(func, context, t1, t2, dt_tolerance_seconds) &#8658; [`astro_search_result_t`](#astro_search_result_t)



---

<a name="Astronomy_SearchHourAngle"></a>
### Astronomy_SearchHourAngle(body, observer, hourAngle, dateStart) &#8658; [`astro_hour_angle_t`](#astro_hour_angle_t)



---

<a name="Astronomy_SearchLunarApsis"></a>
### Astronomy_SearchLunarApsis(startTime) &#8658; [`astro_apsis_t`](#astro_apsis_t)



---

<a name="Astronomy_SearchMaxElongation"></a>
### Astronomy_SearchMaxElongation(body, startDate) &#8658; [`astro_elongation_t`](#astro_elongation_t)



---

<a name="Astronomy_SearchMoonPhase"></a>
### Astronomy_SearchMoonPhase(targetLon, dateStart, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)



---

<a name="Astronomy_SearchMoonQuarter"></a>
### Astronomy_SearchMoonQuarter(dateStart) &#8658; [`astro_moon_quarter_t`](#astro_moon_quarter_t)



---

<a name="Astronomy_SearchPeakMagnitude"></a>
### Astronomy_SearchPeakMagnitude(body, startDate) &#8658; [`astro_illum_t`](#astro_illum_t)



---

<a name="Astronomy_SearchRelativeLongitude"></a>
### Astronomy_SearchRelativeLongitude(body, targetRelLon, startDate) &#8658; [`astro_search_result_t`](#astro_search_result_t)



---

<a name="Astronomy_SearchRiseSet"></a>
### Astronomy_SearchRiseSet(body, observer, direction, dateStart, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)



---

<a name="Astronomy_SearchSunLongitude"></a>
### Astronomy_SearchSunLongitude(targetLon, dateStart, limitDays) &#8658; [`astro_search_result_t`](#astro_search_result_t)



---

<a name="Astronomy_Seasons"></a>
### Astronomy_Seasons(calendar_year) &#8658; [`astro_seasons_t`](#astro_seasons_t)



---

<a name="Astronomy_SunPosition"></a>
### Astronomy_SunPosition(time) &#8658; [`astro_ecliptic_t`](#astro_ecliptic_t)



---

<a name="Astronomy_TimeFromUtc"></a>
### Astronomy_TimeFromUtc(utc) &#8658; [`astro_time_t`](#astro_time_t)



---

<a name="Astronomy_UtcFromTime"></a>
### Astronomy_UtcFromTime(time) &#8658; [`astro_utc_t`](#astro_utc_t)



---

<a name="Astronomy_VectorLength"></a>
### Astronomy_VectorLength(vector) &#8658; `double`

**Calculates the length of the given vector.** 




## Enumerated Types



---

<a name="astro_apsis_kind_t"></a>


---

<a name="astro_body_t"></a>


---

<a name="astro_refraction_t"></a>


---

<a name="astro_status_t"></a>


---

<a name="astro_visibility_t"></a>

## Structures



---

<a name="astro_angle_result_t"></a>
#### `astro_angle_result_t`

**An angular value expressed in degrees.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `angle` |  An angle expressed in degrees.  |


---

<a name="astro_apsis_t"></a>
#### `astro_apsis_t`

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
#### `astro_ecliptic_t`

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
#### `astro_elongation_t`

**Contains information about the visibility of a celestial body at a given date and time.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation.  |
| [`astro_visibility_t`](#astro_visibility_t) | `visibility` |  Whether the body is best seen in the morning or the evening.  |
| `double` | `elongation` |  The angle in degrees between the body and the Sun, as seen from the Earth.  |
| `double` | `relative_longitude` |  The difference between the ecliptic longitudes of the body and the Sun.  |


---

<a name="astro_equatorial_t"></a>
#### `astro_equatorial_t`

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
#### `astro_func_result_t`

**A real value returned by a function whose ascending root is to be found.** 



When calling [`Astronomy_Search`](#Astronomy_Search), the caller must pass in a callback function compatible with the function-pointer type astro_search_func_t whose ascending root is to be found. That callback function must return [`astro_func_result_t`](#astro_func_result_t). If the function call is successful, it will set `status` to [`ASTRO_SUCCESS`](#ASTRO_SUCCESS) and `value` to the numeric value appropriate for the given date and time. If the call fails for some reason, it should set `status` to an appropriate error value other than `ASTRO_SUCCESS`; in the error case, to guard against any possible misuse of `value`, it is recommended to set `value` to `NAN`, though this is not strictly necessary. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `value` |  The value returned by a function whose ascending root is to be found.  |


---

<a name="astro_horizon_t"></a>
#### `astro_horizon_t`

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
#### `astro_hour_angle_t`

**Information about a celestial body crossing a specific hour angle.** 



Returned by the function [`Astronomy_SearchHourAngle`](#Astronomy_SearchHourAngle) to report information about a celestial body crossing a certain hour angle as seen by a specified topocentric observer. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time when the body crosses the specified hour angle.  |
| [`astro_horizon_t`](#astro_horizon_t) | `hor` |  Apparent coordinates of the body at the time it crosses the specified hour angle.  |


---

<a name="astro_illum_t"></a>
#### `astro_illum_t`

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
#### `astro_moon_quarter_t`

**A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `int` | `quarter` |  0=new moon, 1=first quarter, 2=full moon, 3=third quarter.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the lunar quarter.  |


---

<a name="astro_observer_t"></a>
#### `astro_observer_t`

**The location of an observer on (or near) the surface of the Earth.** 



This structure is passed to functions that calculate phenomena as observed from a particular place on the Earth. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `latitude` |  Geographic latitude in degrees north (positive) or south (negative) of the equator.  |
| `double` | `longitude` |  Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.  |
| `double` | `height` |  The height above (positive) or below (negative) sea level, expressed in meters.  |


---

<a name="astro_search_result_t"></a>
#### `astro_search_result_t`

**The result of a search for an astronomical event.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The time at which a searched-for event occurs.  |


---

<a name="astro_seasons_t"></a>
#### `astro_seasons_t`

**The dates and times of changes of season for a given calendar year.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `mar_equinox` |  The date and time of the March equinox for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `jun_solstice` |  The date and time of the June soltice for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `sep_equinox` |  The date and time of the September equinox for the specified year.  |
| [`astro_time_t`](#astro_time_t) | `dec_solstice` |  The date and time of the December solstice for the specified year.  |


---

<a name="astro_time_t"></a>
#### `astro_time_t`

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


---

<a name="astro_utc_t"></a>
#### `astro_utc_t`

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
#### `astro_vector_t`

**A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  ASTRO_SUCCESS if this struct is valid; otherwise an error code.  |
| `double` | `x` |  The Cartesian x-coordinate of the vector in AU.  |
| `double` | `y` |  The Cartesian y-coordinate of the vector in AU.  |
| `double` | `z` |  The Cartesian z-coordinate of the vector in AU.  |
| [`astro_time_t`](#astro_time_t) | `t` |  The date and time at which this vector is valid.  |

## Type Definitions



---

<a name="astro_search_func_t"></a>
