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


---

<a name="Astronomy_AngleFromSun"></a>


---

<a name="Astronomy_BodyCode"></a>


---

<a name="Astronomy_BodyName"></a>


---

<a name="Astronomy_CurrentTime"></a>


---

<a name="Astronomy_Ecliptic"></a>


---

<a name="Astronomy_EclipticLongitude"></a>


---

<a name="Astronomy_Elongation"></a>


---

<a name="Astronomy_Equator"></a>


---

<a name="Astronomy_GeoMoon"></a>


---

<a name="Astronomy_GeoVector"></a>


---

<a name="Astronomy_HelioVector"></a>


---

<a name="Astronomy_Horizon"></a>


---

<a name="Astronomy_Illumination"></a>


---

<a name="Astronomy_LongitudeFromSun"></a>


---

<a name="Astronomy_MakeObserver"></a>


---

<a name="Astronomy_MakeTime"></a>


---

<a name="Astronomy_MoonPhase"></a>


---

<a name="Astronomy_NextLunarApsis"></a>


---

<a name="Astronomy_NextMoonQuarter"></a>


---

<a name="Astronomy_Search"></a>


---

<a name="Astronomy_SearchHourAngle"></a>


---

<a name="Astronomy_SearchLunarApsis"></a>


---

<a name="Astronomy_SearchMaxElongation"></a>


---

<a name="Astronomy_SearchMoonPhase"></a>


---

<a name="Astronomy_SearchMoonQuarter"></a>


---

<a name="Astronomy_SearchPeakMagnitude"></a>


---

<a name="Astronomy_SearchRelativeLongitude"></a>


---

<a name="Astronomy_SearchRiseSet"></a>


---

<a name="Astronomy_SearchSunLongitude"></a>


---

<a name="Astronomy_Seasons"></a>


---

<a name="Astronomy_SunPosition"></a>


---

<a name="Astronomy_TimeFromUtc"></a>


---

<a name="Astronomy_UtcFromTime"></a>


---

<a name="Astronomy_VectorLength"></a>

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
#### `astro_angle_result_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `angle` | `double` | <p>An angle expressed in degrees. </p> |


---

<a name="astro_apsis_t"></a>
#### `astro_apsis_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `time` | `astro_time_t` | <p>The date and time of the apsis. </p> |
| `kind` | `astro_apsis_kind_t` | <p>Whether this is a pericenter or apocenter event. </p> |
| `dist_au` | `double` | <p>The distance between the centers of the bodies in astronomical units. </p> |
| `dist_km` | `double` | <p>The distance between the centers of the bodies in kilometers. </p> |


---

<a name="astro_ecliptic_t"></a>
#### `astro_ecliptic_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `ex` | `double` | <p>Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane. </p> |
| `ey` | `double` | <p>Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox. </p> |
| `ez` | `double` | <p>Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north. </p> |
| `elat` | `double` | <p>Latitude in degrees north (positive) or south (negative) of the ecliptic plane. </p> |
| `elon` | `double` | <p>Longitude in degrees around the ecliptic plane prograde from the equinox. </p> |


---

<a name="astro_elongation_t"></a>
#### `astro_elongation_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `time` | `astro_time_t` | <p>The date and time of the observation. </p> |
| `visibility` | `astro_visibility_t` | <p>Whether the body is best seen in the morning or the evening. </p> |
| `elongation` | `double` | <p>The angle in degrees between the body and the Sun, as seen from the Earth. </p> |
| `relative_longitude` | `double` | <p>The difference between the ecliptic longitudes of the body and the Sun. </p> |


---

<a name="astro_equatorial_t"></a>
#### `astro_equatorial_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `ra` | `double` | <p>right ascension in sidereal hours. </p> |
| `dec` | `double` | <p>declination in degrees </p> |
| `dist` | `double` | <p>distance to the celestial body in AU. </p> |


---

<a name="astro_func_result_t"></a>
#### `astro_func_result_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `value` | `double` | <p>The value returned by a function whose ascending root is to be found. </p> |


---

<a name="astro_horizon_t"></a>
#### `astro_horizon_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `azimuth` | `double` | <p>Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. </p> |
| `altitude` | `double` | <p>Angle in degrees above (positive) or below (negative) the observer's horizon. </p> |
| `ra` | `double` | <p>Right ascension in sidereal hours. </p> |
| `dec` | `double` | <p>Declination in degrees. </p> |


---

<a name="astro_hour_angle_t"></a>
#### `astro_hour_angle_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `time` | `astro_time_t` | <p>The date and time when the body crosses the specified hour angle. </p> |
| `hor` | `astro_horizon_t` | <p>Apparent coordinates of the body at the time it crosses the specified hour angle. </p> |


---

<a name="astro_illum_t"></a>
#### `astro_illum_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `time` | `astro_time_t` | <p>The date and time of the observation. </p> |
| `mag` | `double` | <p>The visual magnitude of the body. Smaller values are brighter. </p> |
| `phase_angle` | `double` | <p>The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. </p> |
| `helio_dist` | `double` | <p>The distance between the Sun and the body at the observation time. </p> |
| `ring_tilt` | `double` | <p>For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0. </p> |


---

<a name="astro_moon_quarter_t"></a>
#### `astro_moon_quarter_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `quarter` | `int` | <p>0=new moon, 1=first quarter, 2=full moon, 3=third quarter. </p> |
| `time` | `astro_time_t` | <p>The date and time of the lunar quarter. </p> |


---

<a name="astro_observer_t"></a>
#### `astro_observer_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `latitude` | `double` | <p>Geographic latitude in degrees north (positive) or south (negative) of the equator. </p> |
| `longitude` | `double` | <p>Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. </p> |
| `height` | `double` | <p>The height above (positive) or below (negative) sea level, expressed in meters. </p> |


---

<a name="astro_search_result_t"></a>
#### `astro_search_result_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `time` | `astro_time_t` | <p>The time at which a searched-for event occurs. </p> |


---

<a name="astro_seasons_t"></a>
#### `astro_seasons_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `mar_equinox` | `astro_time_t` | <p>The date and time of the March equinox for the specified year. </p> |
| `jun_solstice` | `astro_time_t` | <p>The date and time of the June soltice for the specified year. </p> |
| `sep_equinox` | `astro_time_t` | <p>The date and time of the September equinox for the specified year. </p> |
| `dec_solstice` | `astro_time_t` | <p>The date and time of the December solstice for the specified year. </p> |


---

<a name="astro_time_t"></a>
#### `astro_time_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `ut` | `double` | <b><p>UT1/UTC number of days since noon on January 1, 2000. </p></b><p>The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine.</p><p>Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed.</p><p>UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed.</p><p>The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time.</p><p>Before the era of atomic timekeeping, days based on the Earth's rotation were often known as <i>mean solar days</i>. </p> |
| `tt` | `double` | <b><p>Terrestrial Time days since noon on January 1, 2000. </p></b><p>Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed <a href="https://physics.nist.gov/cuu/Units/second.html>SI seconds</a> divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation.</p><p>The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth.</p><p>Historically, Terrestrial Time has also been known by the term <i>Ephemeris Time</i> (ET). </p> |


---

<a name="astro_utc_t"></a>
#### `astro_utc_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `year` | `int` | <p>The year value, e.g. 2019. </p> |
| `month` | `int` | <p>The month value: 1=January, 2=February, ..., 12=December. </p> |
| `day` | `int` | <p>The day of the month in the range 1..31. </p> |
| `hour` | `int` | <p>The hour of the day in the range 0..23. </p> |
| `minute` | `int` | <p>The minute of the hour in the range 0..59. </p> |
| `second` | `double` | <p>The floating point number of seconds in the range [0,60). </p> |


---

<a name="astro_vector_t"></a>
#### `astro_vector_t` (structure type)

| Member | Type | Description |
| ------ | ---- | ----------- |
| `status` | `astro_status_t` | <p>ASTRO_SUCCESS if this struct is valid; otherwise an error code. </p> |
| `x` | `double` | <p>The Cartesian x-coordinate of the vector in AU. </p> |
| `y` | `double` | <p>The Cartesian y-coordinate of the vector in AU. </p> |
| `z` | `double` | <p>The Cartesian z-coordinate of the vector in AU. </p> |
| `t` | `astro_time_t` | <p>The date and time at which this vector is valid. </p> |

## Type Definitions



---

<a name="astro_search_func_t"></a>
