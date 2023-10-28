# Astronomy Engine (C/C++)

This is the complete programming reference for the C version of
Astronomy Engine. It can be used directly from C++ programs also.
Other programming languages are supported.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start
To include Astronomy Engine in your own C or C++ program, all you need are the
files `astronomy.h` and `astronomy.c` from this directory.

To get started quickly, here are some [examples](../../demo/c/).

---

## Contents

- [Topic Index](#topics)
- [Functions](#functions)
- [Constants](#constants)
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
| [FormatTime](#Astronomy_FormatTime) | Formats an [`astro_time_t`](#astro_time_t) value as an ISO 8601 string. |

### Celestial bodies

| Function | Description |
| -------- | ----------- |
| [BodyCode](#Astronomy_BodyCode) | Converts the English name of a celestial body to its equivalent [`astro_body_t`](#astro_body_t) enumeration value. |
| [BodyName](#Astronomy_BodyName) | Converts an [`astro_body_t`](#astro_body_t) enumeration value to its equivalent English name as a string. |

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#Astronomy_HelioVector) | Calculates body position vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy_GeoVector)     | Calculates body position vector with respect to the center of the Earth. |
| [Equator](#Astronomy_Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Astronomy_Ecliptic)       | Converts J2000 mean equator (EQJ) coordinates to true ecliptic of date (ECT) coordinates. |
| [EclipticLongitude](#Astronomy_EclipticLongitude) | Calculates true ecliptic of date (ECT) longitude for a body. |
| [Horizon](#Astronomy_Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [PairLongitude](#Astronomy_PairLongitude) | Calculates the difference in apparent ecliptic longitude between two bodies, as seen from the Earth. |
| [BaryState](#Astronomy_BaryState) | Calculates the barycentric position and velocity vectors of the Sun or a planet. |

### Geographic helper functions

| Function | Description |
| -------- | ----------- |
| [ObserverVector](#Astronomy_ObserverVector) | Calculates a vector from the center of the Earth to an observer on the Earth's surface. |
| [VectorObserver](#Astronomy_VectorObserver) | Calculates the geographic coordinates for a geocentric equatorial vector. |


### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSetEx](#Astronomy_SearchRiseSetEx) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchAltitude](#Astronomy_SearchAltitude) | Finds time when a body reaches a given altitude above or below the horizon. Useful for finding civil, nautical, or astronomical twilight. |
| [SearchHourAngleEx](#Astronomy_SearchHourAngleEx) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#Astronomy_MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#Astronomy_SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#Astronomy_SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy_NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Eclipses and Transits

| Function | Description |
| -------- | ----------- |
| [SearchLunarEclipse](#Astronomy_SearchLunarEclipse) | Search for the first lunar eclipse after a given date. |
| [NextLunarEclipse](#Astronomy_NextLunarEclipse) | Continue searching for more lunar eclipses. |
| [SearchGlobalSolarEclipse](#Astronomy_SearchGlobalSolarEclipse) | Search for the first solar eclipse that is visible anywhere in the world after a given date. |
| [NextGlobalSolarEclipse](#Astronomy_NextGlobalSolarEclipse) | Continue searching for more global solar eclipses. |
| [SearchLocalSolarEclipse](#Astronomy_SearchLocalSolarEclipse) | Search for the first solar eclipse as seen at a particular location after a given date. |
| [NextLocalSolarEclipse](#Astronomy_NextLocalSolarEclipse) | Continue searching for more local solar eclipses. |
| [SearchTransit](#Astronomy_SearchTransit) | Search for the next transit of Mercury or Venus. |
| [NextTransit](#Astronomy_NextTransit) | Continue searching for transits of Mercury or Venus. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#Astronomy_SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy_NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Planet perihelion and aphelion

| Function | Description |
| -------- | ----------- |
| [SearchPlanetApsis](#Astronomy_SearchPlanetApsis) | Finds the next perihelion or aphelion of a planet after a specified date. |
| [NextPlanetApsis](#Astronomy_NextPlanetApsis) | Given an already-found apsis, finds the next perihelion or aphelion of a planet. |

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
| [SearchSunLongitude](#Astronomy_SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the true ecliptic of date (ECT) system. |
| [Seasons](#Astronomy_Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#Astronomy_SunPosition) | Calculates the Sun's apparent true ecliptic of date (ECT) coordinates as seen from the Earth. |

### Coordinate transforms

The following orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = J2000 Mean Equator**: Uses the Earth's mean equator (corrected for precession but ignoring nutation) on January 1, 2000, at noon UTC. This moment in time is called J2000.
- **EQD = True Equator of Date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECL = J2000 Mean Ecliptic**: Uses the plane of the Earth's orbit around the Sun at J2000. The x-axis is referenced against the J2000 mean equinox.
- **ECT = True Ecliptic of Date**: Uses the true (corrected for precession and nutation) orbital plane of the Earth on the given date. The x-axis is referenced against the true equinox for that date.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.
- **GAL = Galactic**: Based on the IAU 1958 definition of galactic coordinates.

| Function | Description |
| -------- | ----------- |
| [RotateVector](#Astronomy_RotateVector) | Applies a rotation matrix to a vector, yielding a vector in another orientation system. |
| [InverseRotation](#Astronomy_InverseRotation) | Given a rotation matrix, finds the inverse rotation matrix that does the opposite transformation. |
| [CombineRotation](#Astronomy_CombineRotation) | Given two rotation matrices, returns a rotation matrix that combines them into a net transformation. |
| [IdentityMatrix](#Astronomy_IdentityMatrix) | Returns a 3x3 identity matrix, which can be used to form other rotation matrices. |
| [Pivot](#Astronomy_Pivot) | Transforms a rotation matrix by pivoting it around a given axis by a given angle. |
| [VectorFromSphere](#Astronomy_VectorFromSphere) | Converts spherical coordinates to Cartesian coordinates. |
| [SphereFromVector](#Astronomy_SphereFromVector) | Converts Cartesian coordinates to spherical coordinates. |
| [EquatorFromVector](#Astronomy_EquatorFromVector) | Given an equatorial vector, calculates equatorial angular coordinates. |
| [VectorFromHorizon](#Astronomy_VectorFromHorizon) | Given apparent angular horizontal coordinates, calculates horizontal vector. |
| [HorizonFromVector](#Astronomy_HorizonFromVector) | Given a vector in horizontal orientation, calculates horizontal angular coordinates. |
| [Rotation_EQD_EQJ](#Astronomy_Rotation_EQD_EQJ) | Calculates a rotation matrix from true equator of date (EQD) to J2000 mean equator (EQJ). |
| [Rotation_EQD_ECT](#Astronomy_Rotation_EQD_ECT) | Calculates a rotation matrix from true equator of date (EQD) to true ecliptic of date (ECT). |
| [Rotation_EQD_ECL](#Astronomy_Rotation_EQD_ECL) | Calculates a rotation matrix from true equator of date (EQD) to J2000 mean ecliptic (ECL). |
| [Rotation_EQD_HOR](#Astronomy_Rotation_EQD_HOR) | Calculates a rotation matrix from true equator of date (EQD) to horizontal (HOR). |
| [Rotation_EQJ_EQD](#Astronomy_Rotation_EQJ_EQD) | Calculates a rotation matrix from J2000 mean equator (EQJ) to true equator of date (EQD). |
| [Rotation_EQJ_ECT](#Astronomy_Rotation_EQJ_ECT) | Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL). |
| [Rotation_EQJ_ECL](#Astronomy_Rotation_EQJ_ECL) | Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT). |
| [Rotation_EQJ_HOR](#Astronomy_Rotation_EQJ_HOR) | Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR). |
| [Rotation_ECT_EQD](#Astronomy_Rotation_ECT_EQD) | Calculates a rotation matrix from true ecliptic of date (ECT) to true equator of date (EQD). |
| [Rotation_ECT_EQJ](#Astronomy_Rotation_ECT_EQJ) | Calculates a rotation matrix from true ecliptic of date (ECT) J2000 mean equator (EQJ). |
| [Rotation_ECL_EQD](#Astronomy_Rotation_ECL_EQD) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to true true equator of date (EQD). |
| [Rotation_ECL_EQJ](#Astronomy_Rotation_ECL_EQJ) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ). |
| [Rotation_ECL_HOR](#Astronomy_Rotation_ECL_HOR) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR). |
| [Rotation_HOR_EQD](#Astronomy_Rotation_HOR_EQD) | Calculates a rotation matrix from horizontal (HOR) to true equator of date (EQD). |
| [Rotation_HOR_EQJ](#Astronomy_Rotation_HOR_EQJ) | Calculates a rotation matrix from horizontal (HOR) to J2000 mean equator (EQJ). |
| [Rotation_HOR_ECL](#Astronomy_Rotation_HOR_ECL) | Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL). |
| [Rotation_EQJ_GAL](#Astronomy_Rotation_EQJ_GAL) | Calculates a rotation matrix from J2000 mean equator (EQJ) to galactic (GAL). |
| [Rotation_GAL_EQJ](#Astronomy_Rotation_EQJ_GAL) | Calculates a rotation matrix from galactic (GAL) to J2000 mean equator (EQJ). |

### Gravitational simulation of small bodies

Astronomy Engine provides a generic gravity simulator that allows you to
model the trajectories of one or more small bodies like asteroids,
comets, or coasting spacecraft. If you know an initial position vector
and velocity vector for a small body, the gravity simulator can incrementally
simulate the pull of gravity on it from the Sun and planets, to calculate its
movement through the Solar System.

| Function | Description |
| -------- | ----------- |
| [GravSimInit](#Astronomy_GravSimInit) | Creates a gravity simulator object. |
| [GravSimFree](#Astronomy_GravSimFree) | Releases memory allocated to a gravity simulator object. |
| [GravSimUpdate](#Astronomy_GravSimUpdate) | Advances the gravity simulation by a small time step. |
| [GravSimSwap](#Astronomy_GravSimSwap) | Exchanges the current time step with the previous time step. |
| [GravSimTime](#Astronomy_GravSimTime) | Returns the time of the current simulation step. |
| [GravSimBodyState](#Astronomy_GravSimBodyState) | Get the position and velocity of a Solar System body included in the simulation. |
| [GravSimNumBodies](#Astronomy_GravSimNumBodies) | Returns the number of small bodies represented in this simulation. |
| [GravSimOrigin](#Astronomy_GravSimOrigin) | Returns the body whose center is the coordinate origin that small bodies are referenced to. |

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

<a name="Astronomy_AngleBetween"></a>
### Astronomy_AngleBetween(a, b) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Calculates the angle between two vectors.** 



Given a pair of vectors, this function returns the angle in degrees between the two vectors in 3D space. The angle is measured in the plane that contains both vectors.



**Returns:**  On success, the `status` field holds `ASTRO_SUCCESS` and `angle` holds a number of degrees in the range [0, 180]. If either vector has a zero magnitude or contains NAN (not a number) components, the `status` will hold the error code `ASTRO_BAD_VECTOR`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `a` |  The first vector. | 
| [`astro_vector_t`](#astro_vector_t) | `b` |  The second vector. | 




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

<a name="Astronomy_Atmosphere"></a>
### Astronomy_Atmosphere(elevationMeters) &#8658; [`astro_atmosphere_t`](#astro_atmosphere_t)

**Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.** 



This function calculates idealized values of pressure, temperature, and density using the U.S. Standard Atmosphere (1976) model.

See: [https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml](https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml)[https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf](https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf)[https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf](https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf)



**Returns:**  astro_atmosphere_tp0 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `elevationMeters` |  The elevation above sea level at which to calculate atmospheric variables. The value must be at least -500 to +100000, or the function will fail with status `ASTRO_INVALID_PARAMETER`. | 




---

<a name="Astronomy_BackdatePosition"></a>
### Astronomy_BackdatePosition(time, observerBody, targetBody, aberration) &#8658; [`astro_vector_t`](#astro_vector_t)

**Solve for light travel time correction of apparent position.** 



When observing a distant object, for example Jupiter as seen from Earth, the amount of time it takes for light to travel from the object to the observer can significantly affect the object's apparent position.

This function solves the light travel time correction for the apparent relative position vector of a target body as seen by an observer body at a given observation time.

For geocentric calculations, [`Astronomy_GeoVector`](#Astronomy_GeoVector) also includes light travel time correction, but the time `t` embedded in its returned vector refers to the observation time, not the backdated time that light left the observed body. Thus `Astronomy_BackdatePosition` provides direct access to the light departure time for callers that need it.

For a more generalized light travel correction solver, see [`Astronomy_CorrectLightTravel`](#Astronomy_CorrectLightTravel).



**Returns:**  On success, the position vector at the solved backdated time. The returned vector will hold `ASTRO_SUCCESS` in its `status` field, the backdated time in its `t` field, along with the apparent relative position. If an error occurs, `status` will hold an error code and the remaining fields should be ignored. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The time of observation.  | 
| [`astro_body_t`](#astro_body_t) | `observerBody` |  The body to be used as the observation location.  | 
| [`astro_body_t`](#astro_body_t) | `targetBody` |  The body to be observed.  | 
| [`astro_aberration_t`](#astro_aberration_t) | `aberration` |  `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected. | 




---

<a name="Astronomy_BaryState"></a>
### Astronomy_BaryState(body, time) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates barycentric position and velocity vectors for the given body.** 



Given a body and a time, calculates the barycentric position and velocity vectors for the center of that body at that time. The vectors are expressed in J2000 mean equator coordinates (EQJ).



**Returns:**  A structure that contains barycentric position and velocity vectors. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body whose barycentric state vector is to be calculated. Supported values are `BODY_SUN`, `BODY_MOON`, `BODY_EMB`, `BODY_SSB`, and all planets: `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`, `BODY_JUPITER`, `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`.  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate position and velocity.  | 




---

<a name="Astronomy_BodyCode"></a>
### Astronomy_BodyCode(name) &#8658; [`astro_body_t`](#astro_body_t)

**Returns the [`astro_body_t`](#astro_body_t) value corresponding to the given English name.** 





**Returns:**  If `name` is one of the listed strings (case-sensitive), the returned value is the corresponding [`astro_body_t`](#astro_body_t) value, otherwise it is `BODY_INVALID`. 



| Type | Parameter | Description |
| --- | --- | --- |
| `const char *` | `name` |  One of the following strings: Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, EMB, SSB.  | 




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

<a name="Astronomy_Constellation"></a>
### Astronomy_Constellation(ra, dec) &#8658; [`astro_constellation_t`](#astro_constellation_t)

**Determines the constellation that contains the given point in the sky.** 



Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the constellation that contains that point.



**Returns:**  If successful, `status` holds `ASTRO_SUCCESS`, `symbol` holds a pointer to a 3-character string like "Ori", and `name` holds a pointer to the full constellation name like "Orion". 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ra` |  The right ascension (RA) of a point in the sky, using the J2000 equatorial system. | 
| `double` | `dec` |  The declination (DEC) of a point in the sky, using the J2000 equatorial system. | 




---

<a name="Astronomy_CorrectLightTravel"></a>
### Astronomy_CorrectLightTravel(context, func, time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Solve for light travel time of a vector function.** 



When observing a distant object, for example Jupiter as seen from Earth, the amount of time it takes for light to travel from the object to the observer can significantly affect the object's apparent position. This function is a generic solver that figures out how long in the past light must have left the observed object to reach the observer at the specified observation time. It uses a context/function pair as a generic interface that expresses an arbitrary position vector as a function of time.

This function repeatedly calls `func`, passing `context` and a series of time estimates in the past. Then `func` must return a relative position vector between the observer and the target. `Astronomy_CorrectLightTravel` keeps calling `func` with more and more refined estimates of the time light must have left the target to arrive at the observer.

For common use cases, it is simpler to use [`Astronomy_BackdatePosition`](#Astronomy_BackdatePosition) for calculating the light travel time correction of one body observing another body.

For geocentric calculations, [`Astronomy_GeoVector`](#Astronomy_GeoVector) also backdates the returned position vector for light travel time, only it returns the observation time in the returned vector's `t` field rather than the backdated time.



**Returns:**  The position vector returned by `func` at the solved backdated time. On success, the vector will hold `ASTRO_SUCCESS` in its `status` field, the backdated time in its `t` field, along with the apparent relative position. If an error occurs, `status` will hold an error code and the remaining fields should be ignored. 



| Type | Parameter | Description |
| --- | --- | --- |
| `void *` | `context` |  Holds any parameters needed by `func`.  | 
| [`astro_position_func_t`](#astro_position_func_t) | `func` |  Pointer to a function that returns a relative position vector as a function of time.  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The observation time for which to solve for light travel delay.  | 




---

<a name="Astronomy_CurrentTime"></a>
### Astronomy_CurrentTime() &#8658; [`astro_time_t`](#astro_time_t)

**Returns the computer's current date and time in the form of an [`astro_time_t`](#astro_time_t).** 



Uses the computer's system clock to find the current UTC date and time. Converts that date and time to an [`astro_time_t`](#astro_time_t) value and returns the result. Callers can pass this value to other Astronomy Engine functions to calculate current observational conditions.

On supported platforms (Linux/Unix, Mac, Windows), the time is measured with microsecond resolution.

On unsupported platforms, a compiler error will occur due to lack of microsecond resolution support. However, if whole second resolution is good enough for your application, you can define the preprocessor symbol `ASTRONOMY_ENGINE_WHOLE_SECOND` to use the portable function `time(NULL)`. Alternatively, if you do not need to use `Astronomy_CurrentTime`, you can define the preprocessor symbol `ASTRONOMY_ENGINE_NO_CURRENT_TIME` to exclude this function from your code. 

---

<a name="Astronomy_DefineStar"></a>
### Astronomy_DefineStar(body, ra, dec, distanceLightYears) &#8658; [`astro_status_t`](#astro_status_t)

**Assign equatorial coordinates to a user-defined star.** 



Some Astronomy Engine functions allow their `body` parameter to be a user-defined fixed point in the sky, loosely called a "star". This function assigns a right ascension, declination, and distance to one of the eight user-defined stars `BODY_STAR1` .. `BODY_STAR8`.

Stars are not valid until defined. Once defined, they retain their definition until re-defined by another call to `Astronomy_DefineStar`.



**Returns:**  `ASTRO_SUCCESS` indicates the star has been defined. Any other value indicates an error, in which case no change has taken place to any of the star definitions. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  One of the eight user-defined star identifiers: `BODY_STAR1` .. `BODY_STAR8`. | 
| `double` | `ra` |  The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of sidereal hours, and must be within the half-open range [0, 24). | 
| `double` | `dec` |  The declination to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of degrees north (positive) or south (negative) of the J2000 equator, and must be within the closed range [-90, +90]. | 
| `double` | `distanceLightYears` |  The distance between the star and the Sun, expressed in light-years. This value is used to calculate the tiny parallax shift as seen by an observer on Earth. If you don't know the distance to the star, using a large value like 1000 will generally work well. The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations. | 




---

<a name="Astronomy_DeltaT_EspenakMeeus"></a>
### Astronomy_DeltaT_EspenakMeeus(ut) &#8658; `double`

**The default Delta T function used by Astronomy Engine.** 



Espenak and Meeus use a series of piecewise polynomials to approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses". See: [https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html) This is the default Delta T function used by Astronomy Engine.



**Returns:**  The estimated difference TT-UT on the given date, expressed in seconds. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` |  The floating point number of days since noon UTC on January 1, 2000. | 




---

<a name="Astronomy_DeltaT_JplHorizons"></a>
### Astronomy_DeltaT_JplHorizons(ut) &#8658; `double`

**A Delta T function that approximates the one used by the JPL Horizons tool.** 



In order to support unit tests based on data generated by the JPL Horizons online tool, I had to reverse engineer their Delta T function by generating a table that contained it. The main difference between their tool and the Espenak/Meeus function is that they stop extrapolating the Earth's deceleration after the year 2017.



**Returns:**  The estimated difference TT-UT on the given date, expressed in seconds. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` |  The floating point number of days since noon UTC on January 1, 2000. | 




---

<a name="Astronomy_Ecliptic"></a>
### Astronomy_Ecliptic(eqj) &#8658; [`astro_ecliptic_t`](#astro_ecliptic_t)

**Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.** 



Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC on 1 January 2000), this function converts those coordinates to true ecliptic coordinates that are relative to the plane of the Earth's orbit around the Sun on that date.



**Returns:**  Spherical and vector coordinates expressed in true ecliptic coordinates of date (ECT). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `eqj` |  Equatorial coordinates in the EQJ frame of reference. You can call [`Astronomy_GeoVector`](#Astronomy_GeoVector) to obtain suitable equatorial coordinates. | 




---

<a name="Astronomy_EclipticGeoMoon"></a>
### Astronomy_EclipticGeoMoon(time) &#8658; [`astro_spherical_t`](#astro_spherical_t)

**Calculates spherical ecliptic geocentric position of the Moon.** 



Given a time of observation, calculates the Moon's geocentric position in ecliptic spherical coordinates. Provides the ecliptic latitude and longitude in degrees, and the geocentric distance in astronomical units (AU).

The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and equatorial plane at the specified time. This means the Earth's equator is corrected for precession and nutation, and the plane of the Earth's orbit is corrected for gradual obliquity drift.

This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954, which in turn derives from E. W. Brown's lunar theories from the early twentieth century. It is adapted from Turbo Pascal code from the book [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210) by Montenbruck and Pfleger.

To calculate a J2000 mean equator vector instead, use [`Astronomy_GeoMoon`](#Astronomy_GeoMoon).



**Returns:**  The Moon's position expressed in ecliptic coordinates using the true equinox of date (ECT). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the Moon's position.  | 




---

<a name="Astronomy_EclipticLongitude"></a>
### Astronomy_EclipticLongitude(body, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Calculates heliocentric ecliptic longitude of a body.** 



This function calculates the angle around the plane of the Earth's orbit of a celestial body, as seen from the center of the Sun. The angle is measured prograde (in the direction of the Earth's orbit around the Sun) in degrees from the true equinox of date. The ecliptic longitude is always in the range [0, 360).



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



**Returns:**  Topocentric equatorial coordinates of the celestial body. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body to be observed. Not allowed to be `BODY_EARTH`.  | 
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the observation takes place.  | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location on or near the surface of the Earth.  | 
| [`astro_equator_date_t`](#astro_equator_date_t) | `equdate` |  Selects the date of the Earth's equator in which to express the equatorial coordinates.  | 
| [`astro_aberration_t`](#astro_aberration_t) | `aberration` |  Selects whether or not to correct for aberration.  | 




---

<a name="Astronomy_EquatorFromVector"></a>
### Astronomy_EquatorFromVector(vector) &#8658; [`astro_equatorial_t`](#astro_equatorial_t)

**Given an equatorial vector, calculates equatorial angular coordinates.** 





**Returns:**  Angular coordinates expressed in the same equatorial system as `vector`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `vector` |  A vector in an equatorial coordinate system. | 




---

<a name="Astronomy_FormatTime"></a>
### Astronomy_FormatTime(time, format, text, size) &#8658; [`astro_status_t`](#astro_status_t)

**Formats an [`astro_time_t`](#astro_time_t) value as an ISO 8601 string.** 



Given an [`astro_time_t`](#astro_time_t) value `time`, formats it as an ISO 8601 string to the resolution specified by the `format` parameter. The result is stored in the `text` buffer whose capacity in bytes is specified by `size`.



**Returns:**  `ASTRO_SUCCESS` on success; otherwise an error as described in the parameter notes. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time whose civil time `time.ut` is to be formatted as an ISO 8601 string. If the civil time is outside the year range -999999 to +999999, the function fails and returns `ASTRO_BAD_TIME`. Years prior to 1583 are treated as if they are using the modern Gregorian calendar, even when the Julian calendar was actually in effect. The year before 1 AD, commonly known as 1 BC, is represented by the value 0. The year 2 BC is represented by -1, etc. | 
| [`astro_time_format_t`](#astro_time_format_t) | `format` |  Specifies the resolution to which the date and time should be formatted, as explained at [`astro_time_format_t`](#astro_time_format_t). If the value of `format` is not recognized, the function fails and returns `ASTRO_INVALID_PARAMETER`. | 
| `char *` | `text` |  A pointer to a text buffer to receive the output. If `text` is `NULL`, this function returns `ASTRO_INVALID_PARAMETER`. If the function fails for any reason, and `text` is not `NULL`, and `size` is greater than 0, the `text` buffer is set to an empty string. | 
| `size_t` | `size` |  The size in bytes of the buffer pointed to by `text`. The buffer must be large enough to accomodate the output format selected by the `format` parameter, as specified at [`astro_time_format_t`](#astro_time_format_t). If `size` is too small to hold the string as specified by `format`, the `text` buffer is set to `""` (if possible) and the function returns `ASTRO_BUFFER_TOO_SMALL`. A buffer that is `TIME_TEXT_BYTES` (28) bytes or larger is always large enough for this function. | 




---

<a name="Astronomy_GeoEmbState"></a>
### Astronomy_GeoEmbState(time) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates the geocentric position and velocity of the Earth/Moon barycenter.** 



Given a time of observation, calculates the geocentric position and velocity vectors of the Earth/Moon barycenter (EMB). The position (x, y, z) components are expressed in AU (astronomical units). The velocity (vx, vy, vz) components are expressed in AU/day. The coordinates are oriented with respect to the Earth's equator at the J2000 epoch. In Astronomy Engine, this orientation is called EQJ.



**Returns:**  The EMB's position and velocity vectors in geocentric J2000 equatorial coordinates. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the EMB vectors.  | 




---

<a name="Astronomy_GeoMoon"></a>
### Astronomy_GeoMoon(time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates equatorial geocentric position of the Moon at a given time.** 



Given a time of observation, calculates the Moon's position as a vector. The vector gives the location of the Moon's center relative to the Earth's center with x-, y-, and z-components measured in astronomical units. The coordinates are oriented with respect to the Earth's equator at the J2000 epoch. In Astronomy Engine, this orientation is called EQJ.

This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954, which in turn derives from E. W. Brown's lunar theories from the early twentieth century. It is adapted from Turbo Pascal code from the book [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210) by Montenbruck and Pfleger.

To calculate ecliptic spherical coordinates instead, see [`Astronomy_EclipticGeoMoon`](#Astronomy_EclipticGeoMoon).



**Returns:**  The Moon's position as a vector in J2000 Cartesian equatorial (EQJ) coordinates. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the Moon's position.  | 




---

<a name="Astronomy_GeoMoonState"></a>
### Astronomy_GeoMoonState(time) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates equatorial geocentric position and velocity of the Moon at a given time.** 



Given a time of observation, calculates the Moon's position and velocity vectors. The position and velocity are of the Moon's center relative to the Earth's center. The position (x, y, z) components are expressed in AU (astronomical units). The velocity (vx, vy, vz) components are expressed in AU/day. The coordinates are oriented with respect to the Earth's equator at the J2000 epoch. In Astronomy Engine, this orientation is called EQJ.

If you need the Moon's position only, and not its velocity, it is much more efficient to use [`Astronomy_GeoMoon`](#Astronomy_GeoMoon) instead.



**Returns:**  The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the Moon's position and velocity.  | 




---

<a name="Astronomy_GeoVector"></a>
### Astronomy_GeoVector(body, time, aberration) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.** 



This function calculates the position of the given celestial body as a vector, using the center of the Earth as the origin. The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

If given an invalid value for `body`, this function will fail. The caller should always check the `status` field inside the returned [`astro_vector_t`](#astro_vector_t) for `ASTRO_SUCCESS` (success) or any other value (failure) before trusting the resulting vector.

Unlike [`Astronomy_HelioVector`](#Astronomy_HelioVector), this function corrects for light travel time. This means the position of the body is "back-dated" by the amount of time it takes light to travel from that body to an observer on the Earth.

Also, the position can optionally be corrected for [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect causing the apparent direction of the body to be shifted due to transverse movement of the Earth with respect to the rays of light coming from that body.



**Returns:**  A geocentric position vector of the center of the given body. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets. Can also be a star defined by [`Astronomy_DefineStar`](#Astronomy_DefineStar).  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the position.  | 
| [`astro_aberration_t`](#astro_aberration_t) | `aberration` |  `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected.  | 




---

<a name="Astronomy_GravSimBodyState"></a>
### Astronomy_GravSimBodyState(sim, body) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Get the position and velocity of a Solar System body included in the simulation.** 



In order to simulate the movement of small bodies through the Solar System, the simulator needs to calculate the state vectors for the Sun and planets.

If an application wants to know the positions of one or more of the planets in addition to the small bodies, this function provides a way to obtain their state vectors. This is provided for the sake of efficiency, to avoid redundant calculations.



**Returns:**  If the given body is part of the set of calculated bodies (Sun and planets), returns the current time step's state vector for that body, expressed in the coordinate system that was specified by the `originBody` parameter to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit). Success is indicated by the returned structure's `status` field holding `ASTRO_SUCCESS`. Any other `status` value indicates an error, meaning the returned state vector is invalid. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_grav_sim_t">astro_grav_sim_t</a> *</code> | `sim` |  A gravity simulator object created by a successful call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit). | 
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune. | 




---

<a name="Astronomy_GravSimFree"></a>
### Astronomy_GravSimFree(sim) &#8658; `void`

**Releases memory allocated to a gravity simulator object.** 



To avoid memory leaks, any successful call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit) must be paired with a matching call to `Astronomy_GravSimFree`.



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_grav_sim_t">astro_grav_sim_t</a> *</code> | `sim` |  A gravity simulator object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit).  | 




---

<a name="Astronomy_GravSimInit"></a>
### Astronomy_GravSimInit(simOut, originBody, time, numBodies, bodyStateArray) &#8658; [`astro_status_t`](#astro_status_t)

**Allocate and initialize a gravity step simulator.** 



Prepares to simulate a series of incremental time steps, simulating the movement of zero or more small bodies through the Solar System acting under gravitational attraction from the Sun and planets.

After calling this function, you can call [`Astronomy_GravSimUpdate`](#Astronomy_GravSimUpdate) as many times as desired to advance the simulation by small time steps.

If this function succeeds (returns `ASTRO_SUCCESS`), `sim` will be set to a dynamically allocated object. The caller is then responsible for eventually calling [`Astronomy_GravSimFree`](#Astronomy_GravSimFree) to release the memory.



**Returns:**  `ASTRO_SUCCESS` on success, with `*sim` set to a non-NULL value. Otherwise an error code with `*sim` set to NULL. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_grav_sim_t">astro_grav_sim_t</a> **</code> | `simOut` |  The address of a pointer to store the newly allocated simulation object. The type [`astro_grav_sim_t`](#astro_grav_sim_t) is an opaque type, so its internal structure is not documented. | 
| [`astro_body_t`](#astro_body_t) | `originBody` |  Specifies the origin of the reference frame. All position vectors and velocity vectors will use `originBody` as the origin of the coordinate system. This origin applies to all the input vectors provided in the `bodyStateArray` parameter of this function, along with all output vectors returned by [`Astronomy_GravSimUpdate`](#Astronomy_GravSimUpdate). Most callers will want to provide one of the following: `BODY_SUN` for heliocentric coordinates, `BODY_SSB` for solar system barycentric coordinates, or `BODY_EARTH` for geocentric coordinates. Note that the gravity simulator does not correct for light travel time; all state vectors are tied to a Newtonian "instantaneous" time. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The initial time at which to start the simulation. | 
| `int` | `numBodies` |  The number of small bodies to be simulated. This may be any non-negative integer. | 
| `const astro_state_vector_t *` | `bodyStateArray` |  An array of initial state vectors (positions and velocities) of the small bodies to be simulated. The caller must know the positions and velocities of the small bodies at an initial moment in time. Their positions and velocities are expressed with respect to `originBody`, using J2000 mean equator orientation (EQJ). Positions are expressed in astronomical units (AU). Velocities are expressed in AU/day. All the times embedded within the state vectors must be exactly equal to `time`, or this function will fail with the error `ASTRO_INCONSISTENT_TIMES`. | 




---

<a name="Astronomy_GravSimNumBodies"></a>
### Astronomy_GravSimNumBodies(sim) &#8658; `int`

**Returns the number of small bodies represented in this simulation.** 



When a simulation is created by a call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit), the caller specifies the number of small bodies. This function returns that same number, which may be convenient for a caller, so that it does not need to track the body count separately.



| Type | Parameter | Description |
| --- | --- | --- |
| `const astro_grav_sim_t *` | `sim` |  A gravity simulator object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit).  | 




---

<a name="Astronomy_GravSimOrigin"></a>
### Astronomy_GravSimOrigin(sim) &#8658; [`astro_body_t`](#astro_body_t)

**Returns the body whose center is the coordinate origin that small bodies are referenced to.** 



When a simulation is created by a call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit), the caller specifies an `originBody` to indicate the coordinate origin used to represent the small bodies being simulated. This function returns that same [`astro_body_t`](#astro_body_t) value.



| Type | Parameter | Description |
| --- | --- | --- |
| `const astro_grav_sim_t *` | `sim` |  A gravity simulator object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit).  | 




---

<a name="Astronomy_GravSimSwap"></a>
### Astronomy_GravSimSwap(sim) &#8658; `void`

**Exchange the current time step with the previous time step.** 



Sometimes it is helpful to "explore" various times near a given simulation time step, while repeatedly returning to the original time step. For example, when backdating a position for light travel time, the caller may wish to repeatedly try different amounts of backdating. When the backdating solver has converged, the caller wants to leave the simulation in its original state.

This function allows a single "undo" of a simulation, and does so very efficiently.

Usually this function will be called immediately after a matching call to [`Astronomy_GravSimUpdate`](#Astronomy_GravSimUpdate). It has the effect of rolling back the most recent update. If called twice in a row, it reverts the swap and thus has no net effect.

[`Astronomy_GravSimInit`](#Astronomy_GravSimInit) initializes the current state and previous state to be identical. Both states represent the `time` parameter that was passed into the initializer. Therefore, `Astronomy_GravSimSwap` will have no effect from the caller's point of view when passed a simulator that has not yet been updated by a call to [`Astronomy_GravSimUpdate`](#Astronomy_GravSimUpdate).



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_grav_sim_t">astro_grav_sim_t</a> *</code> | `sim` |  A gravity simulator object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit).  | 




---

<a name="Astronomy_GravSimTime"></a>
### Astronomy_GravSimTime(sim) &#8658; [`astro_time_t`](#astro_time_t)

**Returns the time of the current simulation step.** 





| Type | Parameter | Description |
| --- | --- | --- |
| `const astro_grav_sim_t *` | `sim` |  A gravity simulator object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit).  | 




---

<a name="Astronomy_GravSimUpdate"></a>
### Astronomy_GravSimUpdate(sim, time, numBodies, bodyStateArray) &#8658; [`astro_status_t`](#astro_status_t)

**Advances a gravity simulation by a small time step.** 





**Returns:**  `ASTRO_SUCCESS` if the calculation was successful. Otherwise, an error code if something went wrong, in which case the simulation should be considered "broken". This means there is no reliable output in `bodyStateArray` and that no more calculations can be performed with `sim`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_grav_sim_t">astro_grav_sim_t</a> *</code> | `sim` |  A simulation object that was created by a prior call to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit). | 
| [`astro_time_t`](#astro_time_t) | `time` |  A time that is a small increment away from the current simulation time. It is up to the developer to figure out an appropriate time increment. Depending on the trajectories, a smaller or larger increment may be needed for the desired accuracy. Some experimentation may be needed. Generally, bodies that stay in the outer Solar System and move slowly can use larger time steps. Bodies that pass into the inner Solar System and move faster will need a smaller time step to maintain accuracy. The `time` value may be after or before the current simulation time to move forward or backward in time. | 
| `int` | `numBodies` |  The number of bodies whose state vectors are to be updated. This is the number of elements in the `bodyStateArray`. This parameter is passed as a sanity check, and must be equal to the value passed to [`Astronomy_GravSimInit`](#Astronomy_GravSimInit) when `sim` was created. | 
| <code><a href="#astro_state_vector_t">astro_state_vector_t</a> *</code> | `bodyStateArray` |  An array big enough to hold `numBodies` state vectors, to receive the updated positions and velocities of the simulated small bodies. Alternatively, `bodyStateArray` may be NULL if the output of this simulation step is not needed. This makes the call slightly faster. | 




---

<a name="Astronomy_HelioDistance"></a>
### Astronomy_HelioDistance(body, time) &#8658; [`astro_func_result_t`](#astro_func_result_t)

**Calculates the distance from a body to the Sun at a given time.** 



Given a date and time, this function calculates the distance between the center of `body` and the center of the Sun. For the planets Mercury through Neptune, this function is significantly more efficient than calling [`Astronomy_HelioVector`](#Astronomy_HelioVector) followed by [`Astronomy_VectorLength`](#Astronomy_VectorLength).



**Returns:**  If successful, an [`astro_func_result_t`](#astro_func_result_t) structure whose `status` is `ASTRO_SUCCESS` and whose `value` holds the heliocentric distance in AU. Otherwise, `status` reports an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body for which to calculate a heliocentric distance: the Sun, Moon, any of the planets, or a user-defined star. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the heliocentric distance. | 




---

<a name="Astronomy_HelioState"></a>
### Astronomy_HelioState(body, time) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates heliocentric position and velocity vectors for the given body.** 



Given a body and a time, calculates the position and velocity vectors for the center of that body at that time, relative to the center of the Sun. The vectors are expressed in J2000 mean equator coordinates (EQJ). If you need the position vector only, it is more efficient to call [`Astronomy_HelioVector`](#Astronomy_HelioVector). The Sun's center is a non-inertial frame of reference. In other words, the Sun experiences acceleration due to gravitational forces, mostly from the larger planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum, kinetic energy, or other quantities that require a non-accelerating frame of reference, consider using [`Astronomy_BaryState`](#Astronomy_BaryState) instead.



**Returns:**  A structure that contains heliocentric position and velocity vectors. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The celestial body whose heliocentric state vector is to be calculated. Supported values are `BODY_SUN`, `BODY_MOON`, `BODY_EMB`, `BODY_SSB`, and all planets: `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`, `BODY_JUPITER`, `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`. Also allowed to be a user-defined star created by [`Astronomy_DefineStar`](#Astronomy_DefineStar).  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate position and velocity.  | 




---

<a name="Astronomy_HelioVector"></a>
### Astronomy_HelioVector(body, time) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.** 



This function calculates the position of the given celestial body as a vector, using the center of the Sun as the origin. The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

The position is not corrected for light travel time or aberration. This is different from the behavior of [`Astronomy_GeoVector`](#Astronomy_GeoVector).

If given an invalid value for `body`, this function will fail. The caller should always check the `status` field inside the returned [`astro_vector_t`](#astro_vector_t) for `ASTRO_SUCCESS` (success) or any other value (failure) before trusting the resulting vector.



**Returns:**  A heliocentric position vector of the center of the given body. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  A body for which to calculate a heliocentric position: the Sun, Moon, any of the planets, the Solar System Barycenter (SSB), or the Earth Moon Barycenter (EMB). Can also be a star defined by [`Astronomy_DefineStar`](#Astronomy_DefineStar).  | 
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
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of the observer. | 
| `double` | `ra` |  The right ascension of the body in sidereal hours. See function remarks for more details. | 
| `double` | `dec` |  The declination of the body in degrees. See function remarks for more details. | 
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  Selects whether to correct for atmospheric refraction, and if so, which model to use. The recommended value for most uses is `REFRACTION_NORMAL`. See function remarks for more details. | 




---

<a name="Astronomy_HorizonFromVector"></a>
### Astronomy_HorizonFromVector(vector, refraction) &#8658; [`astro_spherical_t`](#astro_spherical_t)

**Converts Cartesian coordinates to horizontal coordinates.** 



Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

*IMPORTANT:* This function differs from [`Astronomy_SphereFromVector`](#Astronomy_SphereFromVector) in two ways:

- `Astronomy_SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise from north (e.g., west = +90), but this function represents a clockwise rotation (e.g., east = +90). The difference is because `Astronomy_SphereFromVector` is intended to preserve the vector "right-hand rule", while this function defines azimuth in a more traditional way as used in navigation and cartography.
- This function optionally corrects for atmospheric refraction, while `Astronomy_SphereFromVector` does not.


The returned structure contains the azimuth in `lon`. It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

The altitude is stored in `lat`.

The distance to the observed object is stored in `dist`, and is expressed in astronomical units (AU).



**Returns:**  If successful, `status` holds `ASTRO_SUCCESS` and the other fields are valid as described in the function remarks. Otherwise `status` holds an error code and the other fields are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_vector_t`](#astro_vector_t) | `vector` |  Cartesian vector to be converted to horizontal coordinates. | 
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  `REFRACTION_NORMAL`: correct altitude for atmospheric refraction (recommended). `REFRACTION_NONE`: no atmospheric refraction correction is performed. `REFRACTION_JPLHOR`: for JPL Horizons compatibility testing only; not recommended for normal use. | 




---

<a name="Astronomy_HourAngle"></a>
### Astronomy_HourAngle(body, time, observer) &#8658; [`astro_func_result_t`](#astro_func_result_t)

**Finds the hour angle of a body for a given observer and time.** 



The *hour angle* of a celestial body indicates its position in the sky with respect to the Earth's rotation. The hour angle depends on the location of the observer on the Earth. The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day. The hour angle increases by 1 unit for every sidereal hour that passes after that point, up to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates the number of hours that have passed since the most recent time that the body has culminated, or reached its highest point.



**Returns:**  [`astro_func_result_t`](#astro_func_result_t) If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and `value` holds the hour angle in the half-open range [0, 24). Otherwise, `status` is an error code that indicates failure. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The body whose observed hour angle is to be found. | 
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The time of the observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location where the observation takes place. | 




---

<a name="Astronomy_IdentityMatrix"></a>
### Astronomy_IdentityMatrix() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Creates an identity rotation matrix.** 



Returns a rotation matrix that has no effect on orientation. This matrix can be the starting point for other operations, such as using a series of calls to [`Astronomy_Pivot`](#Astronomy_Pivot) to create a custom rotation matrix.



**Returns:**  The identity matrix. 



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

<a name="Astronomy_InverseRefraction"></a>
### Astronomy_InverseRefraction(refraction, bent_altitude) &#8658; `double`

**Calculates the inverse of an atmospheric refraction angle.** 



Given an observed altitude angle that includes atmospheric refraction, calculates the negative angular correction to obtain the unrefracted altitude. This is useful for cases where observed horizontal coordinates are to be converted to another orientation system, but refraction first must be removed from the observed position.



**Returns:**  The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing. This will be less than or equal to zero. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  The option selecting which refraction correction to use. See [`Astronomy_Refraction`](#Astronomy_Refraction). | 
| `double` | `bent_altitude` |  The apparent altitude that includes atmospheric refraction. | 




---

<a name="Astronomy_InverseRotation"></a>
### Astronomy_InverseRotation(rotation) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates the inverse of a rotation matrix.** 



Given a rotation matrix that performs some coordinate transform, this function returns the matrix that reverses that transform.



**Returns:**  A rotation matrix that performs the opposite transformation. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `rotation` |  The rotation matrix to be inverted. | 




---

<a name="Astronomy_JupiterMoons"></a>
### Astronomy_JupiterMoons(time) &#8658; [`astro_jupiter_moons_t`](#astro_jupiter_moons_t)

**Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.** 



Calculates position and velocity vectors for Jupiter's moons Io, Europa, Ganymede, and Callisto, at the given date and time. The vectors are jovicentric (relative to the center of Jupiter). Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ). The position components are expressed in astronomical units (AU), and the velocity components are in AU/day.

To convert to heliocentric position vectors, call [`Astronomy_HelioVector`](#Astronomy_HelioVector) with `BODY_JUPITER` to get Jupiter's heliocentric position, then add the jovicentric positions.

Likewise, you can call [`Astronomy_GeoVector`](#Astronomy_GeoVector) with `BODY_JUPITER` to convert to geocentric positions. However, you will have to manually correct for light travel time from the Jupiter system to Earth to figure out what time to pass to `Astronomy_JupiterMoons` to get an accurate picture of how Jupiter and its moons look from Earth.



**Returns:**  Position vectors of Jupiter's largest 4 moons, as described above. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate the position vectors.  | 




---

<a name="Astronomy_LagrangePoint"></a>
### Astronomy_LagrangePoint(point, time, major_body, minor_body) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.** 



Given a more massive "major" body and a much less massive "minor" body, calculates one of the five Lagrange points in relation to the minor body's orbit around the major body. The parameter `point` is an integer that selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body. 2 = the Lagrange point on the far side of the minor body. 3 = the Lagrange point on the far side of the major body. 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position. 5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The function returns the state vector for the selected Lagrange point in J2000 mean equator coordinates (EQJ), with respect to the center of the major body.

To calculate Sun/Earth Lagrange points, pass in `BODY_SUN` for `major_body` and `BODY_EMB` (Earth/Moon barycenter) for `minor_body`. For Lagrange points of the Sun and any other planet, pass in that planet (e.g. `BODY_JUPITER`) for `minor_body`. To calculate Earth/Moon Lagrange points, pass in `BODY_EARTH` and `BODY_MOON` for the major and minor bodies respectively.

In some cases, it may be more efficient to call [`Astronomy_LagrangePointFast`](#Astronomy_LagrangePointFast), especially when the state vectors have already been calculated, or are needed for some other purpose.



**Returns:**  The position and velocity of the selected Lagrange point with respect to the major body's center. 



| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` |  A value 1..5 that selects which of the Lagrange points to calculate.  | 
| [`astro_time_t`](#astro_time_t) | `time` |  The time at which the Lagrange point is to be calculated.  | 
| [`astro_body_t`](#astro_body_t) | `major_body` |  The more massive of the co-orbiting bodies: `BODY_SUN` or `BODY_EARTH`.  | 
| [`astro_body_t`](#astro_body_t) | `minor_body` |  The less massive of the co-orbiting bodies. See main remarks.  | 




---

<a name="Astronomy_LagrangePointFast"></a>
### Astronomy_LagrangePointFast(point, major_state, major_mass, minor_state, minor_mass) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates one of the 5 Lagrange points from body masses and state vectors.** 



Given a more massive "major" body and a much less massive "minor" body, calculates one of the five Lagrange points in relation to the minor body's orbit around the major body. The parameter `point` is an integer that selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body. 2 = the Lagrange point on the far side of the minor body. 3 = the Lagrange point on the far side of the major body. 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position. 5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The caller passes in the state vector and mass for both bodies. The state vectors can be in any orientation and frame of reference. The body masses are expressed as GM products, where G = the universal gravitation constant and M = the body's mass. Thus the units for `major_mass` and `minor_mass` must be au^3/day^2. Use [`Astronomy_MassProduct`](#Astronomy_MassProduct) to obtain GM values for various solar system bodies.

The function returns the state vector for the selected Lagrange point using the same orientation as the state vector parameters `major_state` and `minor_state`, and the position and velocity components are with respect to the major body's center.

Consider calling [`Astronomy_LagrangePoint`](#Astronomy_LagrangePoint), instead of this function, for simpler usage in most cases.



**Returns:**  The position and velocity of the selected Lagrange point with respect to the major body's center. 



| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` |  A value 1..5 that selects which of the Lagrange points to calculate.  | 
| [`astro_state_vector_t`](#astro_state_vector_t) | `major_state` |  The state vector of the major (more massive) of the pair of bodies.  | 
| `double` | `major_mass` |  The mass product GM of the major body.  | 
| [`astro_state_vector_t`](#astro_state_vector_t) | `minor_state` |  The state vector of the minor (less massive) of the pair of bodies.  | 
| `double` | `minor_mass` |  The mass product GM of the minor body.  | 




---

<a name="Astronomy_Libration"></a>
### Astronomy_Libration(time) &#8658; [`astro_libration_t`](#astro_libration_t)

**Calculates the Moon's libration angles at a given moment in time.** 



Libration is an observed back-and-forth wobble of the portion of the Moon visible from the Earth. It is caused by the imperfect tidal locking of the Moon's fixed rotation rate, compared to its variable angular speed of orbit around the Earth.

This function calculates a pair of perpendicular libration angles, one representing rotation of the Moon in ecliptic longitude `elon`, the other in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.

This function also returns the geocentric position of the Moon expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the distance `dist_km` between the centers of the Earth and Moon expressed in kilometers, and the apparent angular diameter of the Moon `diam_deg`.



**Returns:**  The Moon's ecliptic position and libration angles as seen from the Earth. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time for which to calculate libration angles.  | 




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

<a name="Astronomy_MassProduct"></a>
### Astronomy_MassProduct(body) &#8658; `double`

**Returns the product of mass and universal gravitational constant of a Solar System body.** 



For problems involving the gravitational interactions of Solar System bodies, it is helpful to know the product GM, where G = the universal gravitational constant and M = the mass of the body. In practice, GM is known to a higher precision than either G or M alone, and thus using the product results in the most accurate results. This function returns the product GM in the units au^3/day^2, or 0 for invalid bodies. The values come from page 10 of a [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).



**Returns:**  The mass product of the given body in au^3/day^2. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The body for which to find the GM product.  | 




---

<a name="Astronomy_MoonPhase"></a>
### Astronomy_MoonPhase(time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Returns the Moon's phase as an angle from 0 to 360 degrees.** 



This function determines the phase of the Moon using its apparent ecliptic longitude relative to the Sun, as seen from the center of the Earth. Certain values of the angle have conventional definitions:



- 0 = new moon
- 90 = first quarter
- 180 = full moon
- 270 = third quarter




**Returns:**  On success, the function returns the angle as described in the function remarks in the `angle` field and `ASTRO_SUCCESS` in the `status` field. The function should always succeed, but it is a good idea for callers to check the `status` field in the returned structure. Any other value in `status` indicates a failure that should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_NextGlobalSolarEclipse"></a>
### Astronomy_NextGlobalSolarEclipse(prevEclipseTime) &#8658; [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t)

**Searches for the next global solar eclipse in a series.** 



After using [`Astronomy_SearchGlobalSolarEclipse`](#Astronomy_SearchGlobalSolarEclipse) to find the first solar eclipse in a series, you can call this function to find the next consecutive solar eclipse. Pass in the `peak` value from the [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t) returned by the previous call to `Astronomy_SearchGlobalSolarEclipse` or `Astronomy_NextGlobalSolarEclipse` to find the next solar eclipse.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields are as described in [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t). Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `prevEclipseTime` |  A date and time near a new moon. Solar eclipse search will start at the next new moon. | 




---

<a name="Astronomy_NextLocalSolarEclipse"></a>
### Astronomy_NextLocalSolarEclipse(prevEclipseTime, observer) &#8658; [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t)

**Searches for the next local solar eclipse in a series.** 



After using [`Astronomy_SearchLocalSolarEclipse`](#Astronomy_SearchLocalSolarEclipse) to find the first solar eclipse in a series, you can call this function to find the next consecutive solar eclipse. Pass in the `peak` value from the [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t) returned by the previous call to `Astronomy_SearchLocalSolarEclipse` or `Astronomy_NextLocalSolarEclipse` to find the next solar eclipse.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields are as described in [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t). Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `prevEclipseTime` |  A date and time near a new moon. Solar eclipse search will start at the next new moon. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of the observer. | 




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

<a name="Astronomy_NextLunarEclipse"></a>
### Astronomy_NextLunarEclipse(prevEclipseTime) &#8658; [`astro_lunar_eclipse_t`](#astro_lunar_eclipse_t)

**Searches for the next lunar eclipse in a series.** 



After using [`Astronomy_SearchLunarEclipse`](#Astronomy_SearchLunarEclipse) to find the first lunar eclipse in a series, you can call this function to find the next consecutive lunar eclipse. Pass in the `peak` value from the [`astro_lunar_eclipse_t`](#astro_lunar_eclipse_t) returned by the previous call to `Astronomy_SearchLunarEclipse` or `Astronomy_NextLunarEclipse` to find the next lunar eclipse.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields will be valid. Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `prevEclipseTime` |  A date and time near a full moon. Lunar eclipse search will start at the next full moon. | 




---

<a name="Astronomy_NextMoonNode"></a>
### Astronomy_NextMoonNode(prevNode) &#8658; [`astro_node_event_t`](#astro_node_event_t)

**Searches for the next time when the Moon's center crosses through the ecliptic plane.** 



Call [`Astronomy_SearchMoonNode`](#Astronomy_SearchMoonNode) to find the first of a series of nodes. Then call `Astronomy_NextMoonNode` to find as many more consecutive nodes as desired.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the other fields are as documented in [`astro_node_event_t`](#astro_node_event_t). Otherwise, `status` holds an error code and the other structure members are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_node_event_t`](#astro_node_event_t) | `prevNode` |  The previous node found from calling [`Astronomy_SearchMoonNode`](#Astronomy_SearchMoonNode) or `Astronomy_NextMoonNode`. | 




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

<a name="Astronomy_NextPlanetApsis"></a>
### Astronomy_NextPlanetApsis(body, apsis) &#8658; [`astro_apsis_t`](#astro_apsis_t)

**Finds the next planetary perihelion or aphelion event in a series.** 



This function requires an [`astro_apsis_t`](#astro_apsis_t) value obtained from a call to [`Astronomy_SearchPlanetApsis`](#Astronomy_SearchPlanetApsis) or `Astronomy_NextPlanetApsis`. Given an aphelion event, this function finds the next perihelion event, and vice versa.

See [`Astronomy_SearchPlanetApsis`](#Astronomy_SearchPlanetApsis) for more details.



**Returns:**  Same as the return value for [`Astronomy_SearchPlanetApsis`](#Astronomy_SearchPlanetApsis). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The planet for which to find the next perihelion/aphelion event. Not allowed to be `BODY_SUN` or `BODY_MOON`. Must match the body passed into the call that produced the `apsis` parameter. | 
| [`astro_apsis_t`](#astro_apsis_t) | `apsis` |  An apsis event obtained from a call to [`Astronomy_SearchPlanetApsis`](#Astronomy_SearchPlanetApsis) or `Astronomy_NextPlanetApsis`. | 




---

<a name="Astronomy_NextTransit"></a>
### Astronomy_NextTransit(body, prevTransitTime) &#8658; [`astro_transit_t`](#astro_transit_t)

**Searches for another transit of Mercury or Venus.** 



After calling [`Astronomy_SearchTransit`](#Astronomy_SearchTransit) to find a transit of Mercury or Venus, this function finds the next transit after that. Keep calling this function as many times as you want to keep finding more transits.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the other fields are as documented in [`astro_transit_t`](#astro_transit_t). Otherwise, `status` holds an error code and the other structure members are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The planet whose transit is to be found. Must be `BODY_MERCURY` or `BODY_VENUS`. | 
| [`astro_time_t`](#astro_time_t) | `prevTransitTime` |  A date and time near the previous transit. | 




---

<a name="Astronomy_ObserverGravity"></a>
### Astronomy_ObserverGravity(latitude, height) &#8658; `double`

**Calculates the gravitational acceleration experienced by an observer on the Earth.** 



This function implements the WGS 84 Ellipsoidal Gravity Formula. The result is a combination of inward gravitational acceleration with outward centrifugal acceleration, as experienced by an observer in the Earth's rotating frame of reference. The resulting value increases toward the Earth's poles and decreases toward the equator, consistent with changes of the weight measured by a spring scale of a fixed mass moved to different latitudes and heights on the Earth.



**Returns:**  The effective gravitational acceleration expressed in meters per second squared [m/s^2]. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `latitude` |  The latitude of the observer in degrees north or south of the equator. By formula symmetry, positive latitudes give the same answer as negative latitudes, so the sign does not matter. | 
| `double` | `height` |  The height above the sea level geoid in meters. No range checking is done; however, accuracy is only valid in the range 0 to 100000 meters. | 




---

<a name="Astronomy_ObserverState"></a>
### Astronomy_ObserverState(time, observer, equdate) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.** 



This function calculates position and velocity vectors of an observer on or near the surface of the Earth, expressed in equatorial coordinates. It takes into account the rotation of the Earth at the given time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in `equdate` to select either `EQUATOR_J2000` for using J2000 coordinates, or `EQUATOR_OF_DATE` for using coordinates relative to the Earth's equator at the specified time.

The returned position vector has components expressed in astronomical units (AU). To convert to kilometers, multiply the `x`, `y`, and `z` values by the constant value [`KM_PER_AU`](#KM_PER_AU).

The returned velocity vector is measured in AU/day.

If you need the position only, and not the velocity, [`Astronomy_ObserverVector`](#Astronomy_ObserverVector) is slightly more efficient than this function.



**Returns:**  If successful, the returned state vector holds `ASTRO_SUCCESS` in its `status` field, and the position (x, y, z) and velocity (vx, vy, vz) vectors are valid. Otherwise, `status` holds an error code. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time for which to calculate the observer's geocentric state vector. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of a point on or near the surface of the Earth. | 
| [`astro_equator_date_t`](#astro_equator_date_t) | `equdate` |  Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `time` as the orientation. | 




---

<a name="Astronomy_ObserverVector"></a>
### Astronomy_ObserverVector(time, observer, equdate) &#8658; [`astro_vector_t`](#astro_vector_t)

**Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.** 



This function calculates a vector from the center of the Earth to a point on or near the surface of the Earth, expressed in equatorial coordinates. It takes into account the rotation of the Earth at the given time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in `equdate` to select either `EQUATOR_J2000` for using J2000 coordinates, or `EQUATOR_OF_DATE` for using coordinates relative to the Earth's equator at the specified time.

The returned vector has components expressed in astronomical units (AU). To convert to kilometers, multiply the `x`, `y`, and `z` values by the constant value [`KM_PER_AU`](#KM_PER_AU).

The inverse of this function is also available: [`Astronomy_VectorObserver`](#Astronomy_VectorObserver).



**Returns:**  If successful, the returned vector holds `ASTRO_SUCCESS` in its `status` field, and is an equatorial vector from the center of the Earth to the specified location on (or near) the Earth's surface. Otherwise, `status` holds an error code. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time for which to calculate the observer's position vector. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of a point on or near the surface of the Earth. | 
| [`astro_equator_date_t`](#astro_equator_date_t) | `equdate` |  Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `time` as the orientation. | 




---

<a name="Astronomy_PairLongitude"></a>
### Astronomy_PairLongitude(body1, body2, time) &#8658; [`astro_angle_result_t`](#astro_angle_result_t)

**Returns one body's ecliptic longitude with respect to another, as seen from the Earth.** 



This function determines where one body appears around the ecliptic plane (the plane of the Earth's orbit around the Sun) as seen from the Earth, relative to the another body's apparent position. The function returns an angle in the half-open range [0, 360) degrees. The value is the ecliptic longitude of `body1` relative to the ecliptic longitude of `body2`.

The angle is 0 when the two bodies are at the same ecliptic longitude as seen from the Earth. The angle increases in the prograde direction (the direction that the planets orbit the Sun and the Moon orbits the Earth).

When the angle is 180 degrees, it means the two bodies appear on opposite sides of the sky for an Earthly observer.

Neither `body1` nor `body2` is allowed to be `BODY_EARTH`. If this happens, the function fails with the error code `ASTRO_EARTH_NOT_ALLOWED`.



**Returns:**  On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the `angle` field holds a value in the range [0, 360). On failure, the `status` field contains some other value indicating an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body1` |  The first body, whose longitude is to be found relative to the second body. | 
| [`astro_body_t`](#astro_body_t) | `body2` |  The second body, relative to which the longitude of the first body is to be found. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. | 




---

<a name="Astronomy_Pivot"></a>
### Astronomy_Pivot(rotation, axis, angle) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Re-orients a rotation matrix by pivoting it by an angle around one of its axes.** 



Given a rotation matrix, a selected coordinate axis, and an angle in degrees, this function pivots the rotation matrix by that angle around that coordinate axis.

For example, if you have rotation matrix that converts ecliptic coordinates (ECL) to horizontal coordinates (HOR), but you really want to convert ECL to the orientation of a telescope camera pointed at a given body, you can use `Astronomy_Pivot` twice: (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the western axis by the body's altitude angle. The resulting rotation matrix will then reorient ECL coordinates to the orientation of your telescope camera.



**Returns:**  If successful, the return value will have `ASTRO_SUCCESS` in the `status` field, along with a pivoted rotation matrix. Otherwise, `status` holds an appropriate error code and the rotation matrix is invalid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `rotation` |  The input rotation matrix. | 
| `int` | `axis` |  An integer that selects which coordinate axis to rotate around: 0 = x, 1 = y, 2 = z. Any other value will fail with the error code `ASTRO_INVALID_PARAMETER` in the `status` field of the return value. | 
| `double` | `angle` |  An angle in degrees indicating the amount of rotation around the specified axis. Positive angles indicate rotation counterclockwise as seen from the positive direction along that axis, looking towards the origin point of the orientation system. If `angle` is NAN or infinite, the function will fail with the error code `ASTRO_INVALID_PARAMETER`. Any finite number of degrees is allowed, but best precision will result from keeping `angle` in the range [-360, +360]. | 




---

<a name="Astronomy_PlanetOrbitalPeriod"></a>
### Astronomy_PlanetOrbitalPeriod(body) &#8658; `double`

**Returns the average number of days it takes for a planet to orbit the Sun.** 





**Returns:**  The mean orbital period of the body, or 0.0 if the `body` parameter is not valid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  One of the planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, or Pluto.  | 




---

<a name="Astronomy_Refraction"></a>
### Astronomy_Refraction(refraction, altitude) &#8658; `double`

**Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.** 



Given an altitude angle and a refraction option, calculates the amount of "lift" caused by atmospheric refraction. This is the number of degrees higher in the sky an object appears due to the lensing of the Earth's atmosphere. This function works best near sea level. To correct for higher elevations, call [`Astronomy_Atmosphere`](#Astronomy_Atmosphere) for that elevation and multiply the refraction angle by the resulting relative density.



**Returns:**  The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  The option selecting which refraction correction to use. If `REFRACTION_NORMAL`, uses a well-behaved refraction model that works well for all valid values (-90 to +90) of `altitude`. If `REFRACTION_JPLHOR`, this function returns a compatible value with the JPL Horizons tool. If any other value (including `REFRACTION_NONE`), this function returns 0. | 
| `double` | `altitude` |  An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90. | 




---

<a name="Astronomy_Reset"></a>
### Astronomy_Reset() &#8658; `void`

**Frees up all dynamic memory allocated by Astronomy Engine.** 



Astronomy Engine uses dynamic memory allocation in only one place: it makes calculation of Pluto's orbit more efficient by caching 11 KB segments and recycling them. To force purging this cache and freeing all the dynamic memory, you can call this function at any time. It is always safe to call, although it will slow down the very next calculation of Pluto's position for a nearby time value. Calling this function before your program exits is optional, but it will be helpful for leak-checkers like valgrind. 

---

<a name="Astronomy_RotateState"></a>
### Astronomy_RotateState(rotation, state) &#8658; [`astro_state_vector_t`](#astro_state_vector_t)

**Applies a rotation to a state vector, yielding a rotated vector.** 



This function transforms a state vector in one orientation to a vector in another orientation.



**Returns:**  A state vector in the orientation specified by `rotation`. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_rotation_t`](#astro_rotation_t) | `rotation` |  A rotation matrix that specifies how the orientation of the state vector is to be changed. | 
| [`astro_state_vector_t`](#astro_state_vector_t) | `state` |  The state vector whose orientation is to be changed. Both the position and velocity components are transformed. | 




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

<a name="Astronomy_RotationAxis"></a>
### Astronomy_RotationAxis(body, time) &#8658; [`astro_axis_t`](#astro_axis_t)

**Calculates information about a body's rotation axis at a given time.** 



Calculates the orientation of a body's rotation axis, along with the rotation angle of its prime meridian, at a given moment in time.

This function uses formulas standardized by the IAU Working Group on Cartographics and Rotational Elements 2015 report, as described in the following document:

[https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf)

See [`astro_axis_t`](#astro_axis_t) for more detailed information.



**Returns:**  [`astro_axis_t`](#astro_axis_t)



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The body whose rotation axis is to be found. The supported bodies are: `BODY_SUN`, `BODY_MOON`, `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`, `BODY_JUPITER`, `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`. | 
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The time for which the body's rotation axis is to be found. | 




---

<a name="Astronomy_Rotation_ECL_EQD"></a>
### Astronomy_Rotation_ECL_EQD(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of date.



**Returns:**  A rotation matrix that converts ECL to EQD. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the desired equator. | 




---

<a name="Astronomy_Rotation_ECL_EQJ"></a>
### Astronomy_Rotation_ECL_EQJ() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: EQJ = equatorial system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts ECL to EQJ. 



---

<a name="Astronomy_Rotation_ECL_HOR"></a>
### Astronomy_Rotation_ECL_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: HOR = horizontal system.



**Returns:**  A rotation matrix that converts ECL to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the desired horizontal orientation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_ECT_EQD"></a>
### Astronomy_Rotation_ECT_EQD(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Returns a rotation matrix from true ecliptic of date (ECT) to equator of date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECT = true ecliptic of date. Target: EQD = equator of date.



**Returns:**  A rotation matrix that converts ECT to EQD. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the ecliptic/equator conversion. | 




---

<a name="Astronomy_Rotation_ECT_EQJ"></a>
### Astronomy_Rotation_ECT_EQJ(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECT = ecliptic system, using true equinox of the specified date/time. Target: EQJ = equatorial system, using mean equator at J2000 epoch.



**Returns:**  A rotation matrix that converts ECT to EQJ at `time`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator defines the target orientation. | 




---

<a name="Astronomy_Rotation_EQD_ECL"></a>
### Astronomy_Rotation_EQD_ECL(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of date. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQD to ECL. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the source equator. | 




---

<a name="Astronomy_Rotation_EQD_ECT"></a>
### Astronomy_Rotation_EQD_ECT(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Returns a rotation matrix from equator of date (EQD) to true ecliptic of date (ECT).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equator of date. Target: ECT = true ecliptic of date.



**Returns:**  A rotation matrix that converts EQD to ECT. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the equator/ecliptic conversion. | 




---

<a name="Astronomy_Rotation_EQD_EQJ"></a>
### Astronomy_Rotation_EQD_EQJ(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean equator (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of the specified date/time. Target: EQJ = equatorial system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQD at `time` to EQJ. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator defines the source orientation. | 




---

<a name="Astronomy_Rotation_EQD_HOR"></a>
### Astronomy_Rotation_EQD_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of the specified date/time. Target: HOR = horizontal system.



**Returns:**  A rotation matrix that converts EQD to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator applies. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_EQJ_ECL"></a>
### Astronomy_Rotation_EQJ_ECL() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts EQJ to ECL. 



---

<a name="Astronomy_Rotation_EQJ_ECT"></a>
### Astronomy_Rotation_EQJ_ECT(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using mean equator at J2000 epoch. Target: ECT = ecliptic system, using true equinox of the specified date/time.



**Returns:**  A rotation matrix that converts EQJ to ECT at `time`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator defines the target orientation. | 




---

<a name="Astronomy_Rotation_EQJ_EQD"></a>
### Astronomy_Rotation_EQJ_EQD(time) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of the specified date/time.



**Returns:**  A rotation matrix that converts EQJ to EQD at `time`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator defines the target orientation. | 




---

<a name="Astronomy_Rotation_EQJ_GAL"></a>
### Astronomy_Rotation_EQJ_GAL() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Returns a rotation matrix from J2000 mean ecliptic (EQJ) to galactic (GAL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using the equator at the J2000 epoch. Target: GAL = galactic system (IAU 1958 definition).



**Returns:**  A rotation matrix that converts EQJ to GAL. 



---

<a name="Astronomy_Rotation_EQJ_HOR"></a>
### Astronomy_Rotation_EQJ_HOR(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using the equator at the J2000 epoch. Target: HOR = horizontal system.



**Returns:**  A rotation matrix that converts EQJ to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the desired horizontal orientation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_GAL_EQJ"></a>
### Astronomy_Rotation_GAL_EQJ() &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Returns a rotation matrix from ecliptic galactic (GAL) to J2000 (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: GAL = galactic system (IAU 1958 definition). Target: EQJ = equatorial system, using the equator at the J2000 epoch.



**Returns:**  A rotation matrix that converts GAL to EQJ. 



---

<a name="Astronomy_Rotation_HOR_ECL"></a>
### Astronomy_Rotation_HOR_ECL(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system. Target: ECL = ecliptic system, using equator at J2000 epoch.



**Returns:**  A rotation matrix that converts HOR to ECL. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the horizontal observation. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The location of the horizontal observer. | 




---

<a name="Astronomy_Rotation_HOR_EQD"></a>
### Astronomy_Rotation_HOR_EQD(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Target: EQD = equatorial system, using equator of the specified date/time.



**Returns:**  A rotation matrix that converts HOR to EQD at `time` and for `observer`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time at which the Earth's equator applies. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  A location near the Earth's mean sea level that defines the observer's horizon. | 




---

<a name="Astronomy_Rotation_HOR_EQJ"></a>
### Astronomy_Rotation_HOR_EQJ(time, observer) &#8658; [`astro_rotation_t`](#astro_rotation_t)

**Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).** 



This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Target: EQJ = equatorial system, using equator at the J2000 epoch.



**Returns:**  A rotation matrix that converts HOR to EQJ at `time` and for `observer`. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time of the observation. | 
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



**Returns:**  If successful, the returned structure has `status` equal to `ASTRO_SUCCESS` and `time` set to a value within `dt_tolerance_seconds` of an ascending root. On success, the `time` value will always be in the inclusive range [`t1`, `t2`]. If the search fails, `status` will be set to a value other than `ASTRO_SUCCESS`. See function remarks for more details. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_search_func_t`](#astro_search_func_t) | `func` |  The function for which to find the time of an ascending root. See function remarks for more details. | 
| `void *` | `context` |  Any ancillary data needed by the function `func` to calculate a value. The data type varies depending on the function passed in. For example, the function may involve a specific celestial body that must be specified somehow. | 
| [`astro_time_t`](#astro_time_t) | `t1` |  The lower time bound of the search window. See function remarks for more details. | 
| [`astro_time_t`](#astro_time_t) | `t2` |  The upper time bound of the search window. See function remarks for more details. | 
| `double` | `dt_tolerance_seconds` |  Specifies an amount of time in seconds within which a bounded ascending root is considered accurate enough to stop. A typical value is 1 second. | 




---

<a name="Astronomy_SearchAltitude"></a>
### Astronomy_SearchAltitude(body, observer, direction, startTime, limitDays, altitude) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Finds the next time the center of a body passes through a given altitude.** 



Finds when the center of the given body ascends or descends through a given altitude angle, as seen by an observer at the specified location on the Earth. By using the appropriate combination of `direction` and `altitude` parameters, this function can be used to find when civil, nautical, or astronomical twilight begins (dawn) or ends (dusk).

Civil dawn begins before sunrise when the Sun ascends through 6 degrees below the horizon. To find civil dawn, pass `DIRECTION_RISE` for `direction` and -6 for `altitude`.

Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon. To find civil dusk, pass `DIRECTION_SET` for `direction` and -6 for `altitude`.

Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.

Astronomical twilight uses -18 degrees as the `altitude` value.

By convention for twilight time calculations, the altitude is not corrected for atmospheric refraction. This is because the target altitudes are below the horizon, and refraction is not directly observable.

`Astronomy_SearchAltitude` is not intended to find rise/set times of a body for two reasons: (1) Rise/set times of the Sun or Moon are defined by their topmost visible portion, not their centers. (2) Rise/set times are affected significantly by atmospheric refraction. Therefore, it is better to use [`Astronomy_SearchRiseSetEx`](#Astronomy_SearchRiseSetEx) to find rise/set times, which corrects for both of these considerations.

`Astronomy_SearchAltitude` will not work reliably for altitudes at or near the body's maximum or minimum altitudes. To find the time a body reaches minimum or maximum altitude angles, use [`Astronomy_SearchHourAngleEx`](#Astronomy_SearchHourAngleEx).



**Returns:**  On success, the `status` field in the returned structure contains `ASTRO_SUCCESS` and the `time` field contains the date and time of the requested altitude event. If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the altitude event does not occur within `limitDays` days of `startTime`. This is a normal condition, not an error. Any other value of `status` indicates an error of some kind. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy_DefineStar`](#Astronomy_DefineStar). | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The location where observation takes place. You can create an observer structure by calling [`Astronomy_MakeObserver`](#Astronomy_MakeObserver). | 
| [`astro_direction_t`](#astro_direction_t) | `direction` |  Either `DIRECTION_RISE` to find when the body ascends through the altitude, or `DIRECTION_SET` for when the body descends through the altitude. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 
| `double` | `limitDays` |  Limits how many days to search for the body reaching the altitude angle, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit the search to the same day, you can use a value of 1 day. In cases where you want to find the altitude event no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. | 
| `double` | `altitude` |  The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range [-90, +90]. | 




---

<a name="Astronomy_SearchGlobalSolarEclipse"></a>
### Astronomy_SearchGlobalSolarEclipse(startTime) &#8658; [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t)

**Searches for a solar eclipse visible anywhere on the Earth's surface.** 



This function finds the first solar eclipse that occurs after `startTime`. A solar eclipse may be partial, annular, or total. See [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t) for more information. To find a series of solar eclipses, call this function once, then keep calling [`Astronomy_NextGlobalSolarEclipse`](#Astronomy_NextGlobalSolarEclipse) as many times as desired, passing in the `peak` value returned from the previous call.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields are as described in [`astro_global_solar_eclipse_t`](#astro_global_solar_eclipse_t). Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for a solar eclipse. | 




---

<a name="Astronomy_SearchHourAngleEx"></a>
### Astronomy_SearchHourAngleEx(body, observer, hourAngle, startTime, direction) &#8658; [`astro_hour_angle_t`](#astro_hour_angle_t)

**Searches for the time when the center of a body reaches a specified hour angle as seen by an observer on the Earth.** 



The *hour angle* of a celestial body indicates its position in the sky with respect to the Earth's rotation. The hour angle depends on the location of the observer on the Earth. The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day. The hour angle increases by 1 unit for every sidereal hour that passes after that point, up to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates the number of hours that have passed since the most recent time that the body has culminated, or reached its highest point.

This function searches for the next or previous time a celestial body reaches the given hour angle relative to the date and time specified by `startTime`. To find when a body culminates, pass 0 for `hourAngle`. To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.

Note that, especially close to the Earth's poles, a body as seen on a given day may always be above the horizon or always below the horizon, so the caller cannot assume that a culminating object is visible nor that an object is below the horizon at its minimum altitude.

On success, the function reports the date and time, along with the horizontal coordinates of the body at that time, as seen by the given observer.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the other structure fields are valid. Otherwise, `status` holds some other value that indicates an error condition. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy_DefineStar`](#Astronomy_DefineStar). | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  Indicates a location on or near the surface of the Earth where the observer is located. Call [`Astronomy_MakeObserver`](#Astronomy_MakeObserver) to create an observer structure. | 
| `double` | `hourAngle` |  An hour angle value in the range [0, 24) indicating the number of sidereal hours after the body's most recent culmination. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 
| `int` | `direction` |  The direction in time to perform the search: a positive value searches forward in time, a negative value searches backward in time. The function will fail with `ASTRO_INVALID_PARAMETER` if `direction` is zero. | 




---

<a name="Astronomy_SearchLocalSolarEclipse"></a>
### Astronomy_SearchLocalSolarEclipse(startTime, observer) &#8658; [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t)

**Searches for a solar eclipse visible at a specific location on the Earth's surface.** 



This function finds the first solar eclipse that occurs after `startTime`. A solar eclipse may be partial, annular, or total. See [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t) for more information. To find a series of solar eclipses, call this function once, then keep calling [`Astronomy_NextLocalSolarEclipse`](#Astronomy_NextLocalSolarEclipse) as many times as desired, passing in the `peak` value returned from the previous call.

IMPORTANT: An eclipse reported by this function might be partly or completely invisible to the observer due to the time of day.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields are as described in [`astro_local_solar_eclipse_t`](#astro_local_solar_eclipse_t). Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for a solar eclipse. | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The geographic location of the observer. | 




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

<a name="Astronomy_SearchLunarEclipse"></a>
### Astronomy_SearchLunarEclipse(startTime) &#8658; [`astro_lunar_eclipse_t`](#astro_lunar_eclipse_t)

**Searches for a lunar eclipse.** 



This function finds the first lunar eclipse that occurs after `startTime`. A lunar eclipse may be penumbral, partial, or total. See [`astro_lunar_eclipse_t`](#astro_lunar_eclipse_t) for more information. To find a series of lunar eclipses, call this function once, then keep calling [`Astronomy_NextLunarEclipse`](#Astronomy_NextLunarEclipse) as many times as desired, passing in the `peak` value returned from the previous call.



**Returns:**  If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS` and the remaining structure fields will be valid. Any other value indicates an error. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for a lunar eclipse. | 




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

<a name="Astronomy_SearchMoonNode"></a>
### Astronomy_SearchMoonNode(startTime) &#8658; [`astro_node_event_t`](#astro_node_event_t)

**Searches for a time when the Moon's center crosses through the ecliptic plane.** 



Searches for the first ascending or descending node of the Moon after `startTime`. An ascending node is when the Moon's center passes through the ecliptic plane (the plane of the Earth's orbit around the Sun) from south to north. A descending node is when the Moon's center passes through the ecliptic plane from north to south. Nodes indicate possible times of solar or lunar eclipses, if the Moon also happens to be in the correct phase (new or full, respectively).

Call `Astronomy_SearchMoonNode` to find the first of a series of nodes. Then call [`Astronomy_NextMoonNode`](#Astronomy_NextMoonNode) to find as many more consecutive nodes as desired.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS` and the other fields are as documented in [`astro_node_event_t`](#astro_node_event_t). Otherwise, `status` holds an error code and the other structure members are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for an ascending or descending node of the Moon. | 




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
| `double` | `limitDays` |  The number of days away from `startTime` that limits the time window for the search. If the value is negative, the search is performed into the past from `startTime`. Otherwise, the search is performed into the future from `startTime`. | 




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
| [`astro_body_t`](#astro_body_t) | `body` |  Currently only `BODY_VENUS` is allowed. Any other value results in the error `ASTRO_INVALID_BODY`. See function remarks for more details. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time to start searching for the next peak magnitude event. | 




---

<a name="Astronomy_SearchPlanetApsis"></a>
### Astronomy_SearchPlanetApsis(body, startTime) &#8658; [`astro_apsis_t`](#astro_apsis_t)

**Finds the date and time of a planet's perihelion (closest approach to the Sun) or aphelion (farthest distance from the Sun) after a given time.** 



Given a date and time to start the search in `startTime`, this function finds the next date and time that the center of the specified planet reaches the closest or farthest point in its orbit with respect to the center of the Sun, whichever comes first after `startTime`.

The closest point is called *perihelion* and the farthest point is called *aphelion*. The word *apsis* refers to either event.

To iterate through consecutive alternating perihelion and aphelion events, call `Astronomy_SearchPlanetApsis` once, then use the return value to call [`Astronomy_NextPlanetApsis`](#Astronomy_NextPlanetApsis). After that, keep feeding the previous return value from `Astronomy_NextPlanetApsis` into another call of `Astronomy_NextPlanetApsis` as many times as desired.



**Returns:**  If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`, `time` holds the date and time of the next planetary apsis, `kind` holds either `APSIS_PERICENTER` for perihelion or `APSIS_APOCENTER` for aphelion, and the distance values `dist_au` (astronomical units) and `dist_km` (kilometers) are valid. If the function fails, `status` holds some value other than `ASTRO_SUCCESS` that indicates what went wrong, and the other structure fields are invalid. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The planet for which to find the next perihelion/aphelion event. Not allowed to be `BODY_SUN` or `BODY_MOON`. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start searching for the next perihelion or aphelion. | 




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

<a name="Astronomy_SearchRiseSetEx"></a>
### Astronomy_SearchRiseSetEx(body, observer, direction, startTime, limitDays, metersAboveGround) &#8658; [`astro_search_result_t`](#astro_search_result_t)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.** 



This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth. Rise time is when the body first starts to be visible above the horizon. For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon. Set time is the moment when the body appears to vanish below the horizon. Therefore, this function adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).

This function corrects for a typical value of atmospheric refraction, which causes celestial bodies to appear higher above the horizon than they would if the Earth had no atmosphere. Astronomy Engine uses a correction of 34 arcminutes. Real-world refraction varies based on air temperature, pressure, and humidity; such weather-based conditions are outside the scope of Astronomy Engine.

Note that rise or set may not occur in every 24 hour period. For example, near the Earth's poles, there are long periods of time where the Sun stays below the horizon, never rising. Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day. This is because the Moon sets nearly an hour later each day due to orbiting the Earth a significant amount during each rotation of the Earth. Therefore callers must not assume that the function will always succeed.



**Returns:**  On success, the `status` field in the returned structure contains `ASTRO_SUCCESS` and the `time` field contains the date and time of the rise or set time as requested. If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the rise or set event does not occur within `limitDays` days of `startTime`. This is a normal condition, not an error. Any other value of `status` indicates an error of some kind. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy_DefineStar`](#Astronomy_DefineStar). | 
| [`astro_observer_t`](#astro_observer_t) | `observer` |  The location where observation takes place. You can create an observer structure by calling [`Astronomy_MakeObserver`](#Astronomy_MakeObserver). | 
| [`astro_direction_t`](#astro_direction_t) | `direction` |  Either `DIRECTION_RISE` to find a rise time or `DIRECTION_SET` to find a set time. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time at which to start the search. | 
| `double` | `limitDays` |  Limits how many days to search for a rise or set time, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. | 
| `double` | `metersAboveGround` |  Usually the observer is located at ground level. Then this parameter should be zero. But if the observer is significantly higher than ground level, for example in an airplane, this parameter should be a positive number indicating how far above the ground the observer is. An error occurs if `metersAboveGround` is negative. | 




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
| `double` | `limitDays` |  The real-valued number of days, which when added to `startTime`, limits the range of time over which the search looks. It is recommended to keep this value between 1 and 10 days. See function remarks for more details. | 




---

<a name="Astronomy_SearchTransit"></a>
### Astronomy_SearchTransit(body, startTime) &#8658; [`astro_transit_t`](#astro_transit_t)

**Searches for the first transit of Mercury or Venus after a given date.** 



Finds the first transit of Mercury or Venus after a specified date. A transit is when an inferior planet passes between the Sun and the Earth so that the silhouette of the planet is visible against the Sun in the background. To continue the search, pass the `finish` time in the returned structure to [`Astronomy_NextTransit`](#Astronomy_NextTransit).



**Returns:**  If successful, the `status` field in the returned structure hold `ASTRO_SUCCESS` and the other fields are as documented in [`astro_transit_t`](#astro_transit_t). Otherwise, `status` holds an error code and the other structure members are undefined. 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_body_t`](#astro_body_t) | `body` |  The planet whose transit is to be found. Must be `BODY_MERCURY` or `BODY_VENUS`. | 
| [`astro_time_t`](#astro_time_t) | `startTime` |  The date and time for starting the search for a transit. | 




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

<a name="Astronomy_SetDeltaTFunction"></a>
### Astronomy_SetDeltaTFunction(func) &#8658; `void`

**Changes the function Astronomy Engine uses to calculate Delta T.** 



Most programs should not call this function. It is for advanced use cases only. By default, Astronomy Engine uses the function [`Astronomy_DeltaT_EspenakMeeus`](#Astronomy_DeltaT_EspenakMeeus) to estimate changes in the Earth's rotation rate over time. However, for the sake of unit tests that compare calculations against external data sources that use alternative models for Delta T, it is sometimes useful to replace the Delta T model to match. This function allows replacing the Delta T model with any other desired model.



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_deltat_func`](#astro_deltat_func) | `func` |  A pointer to a function to convert UT values to DeltaT values.  | 




---

<a name="Astronomy_SiderealTime"></a>
### Astronomy_SiderealTime(time) &#8658; `double`

**Calculates Greenwich Apparent Sidereal Time (GAST).** 



Given a date and time, this function calculates the rotation of the Earth, represented by the equatorial angle of the Greenwich prime meridian with respect to distant stars (not the Sun, which moves relative to background stars by almost one degree per day). This angle is called Greenwich Apparent Sidereal Time (GAST). GAST is measured in sidereal hours in the half-open range [0, 24). When GAST = 0, it means the prime meridian is aligned with the of-date equinox, corrected at that time for precession and nutation of the Earth's axis. In this context, the "equinox" is the direction in space where the Earth's orbital plane (the ecliptic) intersects with the plane of the Earth's equator, at the location on the Earth's orbit of the (seasonal) March equinox. As the Earth rotates, GAST increases from 0 up to 24 sidereal hours, then starts over at 0. To convert to degrees, multiply the return value by 15.



**Returns:**  {number} 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_time_t">astro_time_t</a> *</code> | `time` |  The date and time for which to find GAST. The parameter is passed by address because it can be modified by the call: As an optimization, this function caches the sidereal time value in `time`, unless it has already been cached, in which case the cached value is reused. If the `time` pointer is NULL, this function returns a NAN value. | 




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

<a name="Astronomy_TerrestrialTime"></a>
### Astronomy_TerrestrialTime(tt) &#8658; [`astro_time_t`](#astro_time_t)

**Converts a terrestrial time value into an [`astro_time_t`](#astro_time_t) value.** 



This function can be used in rare cases where a time must be based on Terrestrial Time (TT) rather than Universal Time (UT). Most developers will want to call [`Astronomy_TimeFromDays`](#Astronomy_TimeFromDays) instead of this function, because usually time is based on civil time adjusted by leap seconds to match the Earth's rotation, rather than the uniformly flowing TT used to calculate solar system dynamics. In rare cases where the caller already knows TT, this function is provided to create an [`astro_time_t`](#astro_time_t) value that can be passed to Astronomy Engine functions.



**Returns:**  An [`astro_time_t`](#astro_time_t) value for the given `tt` value. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `tt` |  The floating point number of days of uniformly flowing Terrestrial Time since the J2000 epoch. | 




---

<a name="Astronomy_TimeFromDays"></a>
### Astronomy_TimeFromDays(ut) &#8658; [`astro_time_t`](#astro_time_t)

**Converts a J2000 day value to an [`astro_time_t`](#astro_time_t) value.** 



This function can be useful for reproducing an [`astro_time_t`](#astro_time_t) structure from its `ut` field only.



**Returns:**  An [`astro_time_t`](#astro_time_t) value for the given `ut` value. 



| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` |  The floating point number of days since noon UTC on January 1, 2000. This time is based on UTC/UT1 civil time. See [`Astronomy_TerrestrialTime`](#Astronomy_TerrestrialTime) if you instead want to create a time value based on atomic Terrestrial Time (TT). | 




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

<a name="Astronomy_VectorFromHorizon"></a>
### Astronomy_VectorFromHorizon(sphere, time, refraction) &#8658; [`astro_vector_t`](#astro_vector_t)

**Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.** 





**Returns:**  A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up). 



| Type | Parameter | Description |
| --- | --- | --- |
| [`astro_spherical_t`](#astro_spherical_t) | `sphere` |  A structure that contains apparent horizontal coordinates: `lat` holds the refracted altitude angle, `lon` holds the azimuth in degrees clockwise from north, and `dist` holds the distance from the observer to the object in AU. | 
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation. This is needed because the returned [`astro_vector_t`](#astro_vector_t) structure requires a valid time value when passed to certain other functions. | 
| [`astro_refraction_t`](#astro_refraction_t) | `refraction` |  The refraction option used to model atmospheric lensing. See [`Astronomy_Refraction`](#Astronomy_Refraction). This specifies how refraction is to be removed from the altitude stored in `sphere.lat`. | 




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




---

<a name="Astronomy_VectorObserver"></a>
### Astronomy_VectorObserver(vector, equdate) &#8658; [`astro_observer_t`](#astro_observer_t)

**Calculates the geographic location corresponding to an equatorial vector.** 



This is the inverse function of [`Astronomy_ObserverVector`](#Astronomy_ObserverVector). Given a geocentric equatorial vector, it returns the geographic latitude, longitude, and elevation for that vector.



**Returns:**  The geographic latitude, longitude, and elevation above sea level that corresponds to the given equatorial vector. 



| Type | Parameter | Description |
| --- | --- | --- |
| <code><a href="#astro_vector_t">astro_vector_t</a> *</code> | `vector` |  The geocentric equatorial position vector for which to find geographic coordinates. The components are expressed in Astronomical Units (AU). You can calculate AU by dividing kilometers by the constant [`KM_PER_AU`](#KM_PER_AU). The time `vector.t` determines the Earth's rotation. The caller must set `vector.t` to a valid time. The vector is passed by reference (using a pointer) so that nutation calculations can be cached inside `vector.t` as an optimization. | 
| [`astro_equator_date_t`](#astro_equator_date_t) | `equdate` |  Selects the date of the Earth's equator in which `vector` is expressed. The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by `vector.t`. Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `vector.t` as the orientation. | 



<a name="constants"></a>
## Constants



---

<a name="AU_PER_LY"></a>
### `AU_PER_LY`

**The number of astronomical units per light-year.** 



```C
#define AU_PER_LY  63241.07708807546
```



---

<a name="CALLISTO_RADIUS_KM"></a>
### `CALLISTO_RADIUS_KM`

**The mean radius of Jupiter's moon Callisto, expressed in kilometers.** 



```C
#define CALLISTO_RADIUS_KM  2410.3
```



---

<a name="C_AUDAY"></a>
### `C_AUDAY`

**The speed of light in AU/day.** 



```C
#define C_AUDAY  173.1446326846693
```



---

<a name="DEG2RAD"></a>
### `DEG2RAD`

**The factor to convert degrees to radians = pi/180.** 



```C
#define DEG2RAD  0.017453292519943296
```



---

<a name="EARTH_EQUATORIAL_RADIUS_KM"></a>
### `EARTH_EQUATORIAL_RADIUS_KM`

**The equatorial radius of the Earth, expressed in kilometers.** 



```C
#define EARTH_EQUATORIAL_RADIUS_KM  6378.1366
```



---

<a name="EARTH_FLATTENING"></a>
### `EARTH_FLATTENING`

**The Earth's polar radius divided by its equatorial radius.** 



```C
#define EARTH_FLATTENING  0.996647180302104
```



---

<a name="EARTH_POLAR_RADIUS_KM"></a>
### `EARTH_POLAR_RADIUS_KM`

**The polar radius of the Earth, expressed in kilometers.** 



```C
#define EARTH_POLAR_RADIUS_KM  (EARTH_EQUATORIAL_RADIUS_KM * EARTH_FLATTENING)
```



---

<a name="EUROPA_RADIUS_KM"></a>
### `EUROPA_RADIUS_KM`

**The mean radius of Jupiter's moon Europa, expressed in kilometers.** 



```C
#define EUROPA_RADIUS_KM  1560.8
```



---

<a name="GANYMEDE_RADIUS_KM"></a>
### `GANYMEDE_RADIUS_KM`

**The mean radius of Jupiter's moon Ganymede, expressed in kilometers.** 



```C
#define GANYMEDE_RADIUS_KM  2631.2
```



---

<a name="HOUR2RAD"></a>
### `HOUR2RAD`

**The factor to convert sidereal hours to radians = pi/12.** 



```C
#define HOUR2RAD  0.2617993877991494365
```



---

<a name="IO_RADIUS_KM"></a>
### `IO_RADIUS_KM`

**The mean radius of Jupiter's moon Io, expressed in kilometers.** 



```C
#define IO_RADIUS_KM  1821.6
```



---

<a name="JUPITER_EQUATORIAL_RADIUS_KM"></a>
### `JUPITER_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Jupiter, expressed in kilometers.** 



```C
#define JUPITER_EQUATORIAL_RADIUS_KM  71492.0
```



---

<a name="JUPITER_MEAN_RADIUS_KM"></a>
### `JUPITER_MEAN_RADIUS_KM`

**The volumetric mean radius of Jupiter, expressed in kilometers.** 



```C
#define JUPITER_MEAN_RADIUS_KM  69911.0
```



---

<a name="JUPITER_POLAR_RADIUS_KM"></a>
### `JUPITER_POLAR_RADIUS_KM`

**The polar radius of Jupiter, expressed in kilometers.** 



```C
#define JUPITER_POLAR_RADIUS_KM  66854.0
```



---

<a name="KM_PER_AU"></a>
### `KM_PER_AU`

**The number of kilometers in one astronomical unit (AU).** 



```C
#define KM_PER_AU  1.4959787069098932e+8
```



---

<a name="MARS_EQUATORIAL_RADIUS_KM"></a>
### `MARS_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Mars, expressed in kilometers.** 



```C
#define MARS_EQUATORIAL_RADIUS_KM  3396.2
```



---

<a name="MARS_POLAR_RADIUS_KM"></a>
### `MARS_POLAR_RADIUS_KM`

**The polar radius of Mars, expressed in kilometers.** 



```C
#define MARS_POLAR_RADIUS_KM  3376.2
```



---

<a name="MERCURY_EQUATORIAL_RADIUS_KM"></a>
### `MERCURY_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Mercury, expressed in kilometers.** 



```C
#define MERCURY_EQUATORIAL_RADIUS_KM  2440.5
```



---

<a name="MERCURY_POLAR_RADIUS_KM"></a>
### `MERCURY_POLAR_RADIUS_KM`

**The polar radius of Mercury, expressed in kilometers.** 



```C
#define MERCURY_POLAR_RADIUS_KM  2438.3
```



---

<a name="MOON_EQUATORIAL_RADIUS_KM"></a>
### `MOON_EQUATORIAL_RADIUS_KM`

**The equatorial radius of the Moon, expressed in kilometers.** 



```C
#define MOON_EQUATORIAL_RADIUS_KM  1738.1
```



---

<a name="MOON_POLAR_RADIUS_KM"></a>
### `MOON_POLAR_RADIUS_KM`

**The polar radius of the Moon, expressed in kilometers.** 



```C
#define MOON_POLAR_RADIUS_KM  1736.0
```



---

<a name="NEPTUNE_EQUATORIAL_RADIUS_KM"></a>
### `NEPTUNE_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Neptune, expressed in kilometers.** 



```C
#define NEPTUNE_EQUATORIAL_RADIUS_KM  24764.0
```



---

<a name="NEPTUNE_POLAR_RADIUS_KM"></a>
### `NEPTUNE_POLAR_RADIUS_KM`

**The polar radius of Neptune, expressed in kilometers.** 



```C
#define NEPTUNE_POLAR_RADIUS_KM  24341.0
```



---

<a name="PLUTO_RADIUS_KM"></a>
### `PLUTO_RADIUS_KM`

**The mean radius of Pluto, expressed in kilometers. Pluto is nearly spherical.** 



```C
#define PLUTO_RADIUS_KM  1188.3
```



---

<a name="RAD2DEG"></a>
### `RAD2DEG`

**The factor to convert radians to degrees = 180/pi.** 



```C
#define RAD2DEG  57.295779513082321
```



---

<a name="RAD2HOUR"></a>
### `RAD2HOUR`

**The factor to convert radians to sidereal hours = 12/pi.** 



```C
#define RAD2HOUR  3.819718634205488
```



---

<a name="SATURN_EQUATORIAL_RADIUS_KM"></a>
### `SATURN_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Saturn, expressed in kilometers.** 



```C
#define SATURN_EQUATORIAL_RADIUS_KM  60268.0
```



---

<a name="SATURN_POLAR_RADIUS_KM"></a>
### `SATURN_POLAR_RADIUS_KM`

**The polar radius of Saturn, expressed in kilometers.** 



```C
#define SATURN_POLAR_RADIUS_KM  54364.0
```



---

<a name="SUN_RADIUS_KM"></a>
### `SUN_RADIUS_KM`

**The mean radius of the Sun's photosphere, expressed in kilometers. The Sun is nearly spherical.** 



```C
#define SUN_RADIUS_KM  695700.0
```



---

<a name="TIME_TEXT_BYTES"></a>
### `TIME_TEXT_BYTES`

**The smallest number of characters that is always large enough for [`Astronomy_FormatTime`](#Astronomy_FormatTime).** 



```C
#define TIME_TEXT_BYTES  28
```



---

<a name="URANUS_EQUATORIAL_RADIUS_KM"></a>
### `URANUS_EQUATORIAL_RADIUS_KM`

**The equatorial radius of Uranus, expressed in kilometers.** 



```C
#define URANUS_EQUATORIAL_RADIUS_KM  25559.0
```



---

<a name="URANUS_POLAR_RADIUS_KM"></a>
### `URANUS_POLAR_RADIUS_KM`

**The polar radius of Uranus, expressed in kilometers.** 



```C
#define URANUS_POLAR_RADIUS_KM  24973.0
```



---

<a name="VENUS_RADIUS_KM"></a>
### `VENUS_RADIUS_KM`

**The mean radius of Venus, expressed in kilometers. Venus is nearly spherical.** 



```C
#define VENUS_RADIUS_KM  6051.8
```


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
| `BODY_EMB` |  Earth/Moon Barycenter  |
| `BODY_SSB` |  Solar System Barycenter  |
| `BODY_STAR1` |  user-defined star #1  |
| `BODY_STAR2` |  user-defined star #2  |
| `BODY_STAR3` |  user-defined star #3  |
| `BODY_STAR4` |  user-defined star #4  |
| `BODY_STAR5` |  user-defined star #5  |
| `BODY_STAR6` |  user-defined star #6  |
| `BODY_STAR7` |  user-defined star #7  |
| `BODY_STAR8` |  user-defined star #8  |



---

<a name="astro_direction_t"></a>
### `astro_direction_t`

**Selects whether to search for a rise time or a set time.** 



The [`Astronomy_SearchRiseSetEx`](#Astronomy_SearchRiseSetEx) function finds the rise or set time of a body depending on the value of its `direction` parameter. 

| Enum Value | Description |
| --- | --- |
| `DIRECTION_RISE` |  Search for the time a body begins to rise above the horizon.  |
| `DIRECTION_SET` |  Search for the time a body finishes sinking below the horizon.  |



---

<a name="astro_eclipse_kind_t"></a>
### `astro_eclipse_kind_t`

**The different kinds of lunar/solar eclipses.** 



| Enum Value | Description |
| --- | --- |
| `ECLIPSE_NONE` |  No eclipse found.  |
| `ECLIPSE_PENUMBRAL` |  A penumbral lunar eclipse. (Never used for a solar eclipse.)  |
| `ECLIPSE_PARTIAL` |  A partial lunar/solar eclipse.  |
| `ECLIPSE_ANNULAR` |  An annular solar eclipse. (Never used for a lunar eclipse.)  |
| `ECLIPSE_TOTAL` |  A total lunar/solar eclipse.  |



---

<a name="astro_equator_date_t"></a>
### `astro_equator_date_t`

**Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.** 



The Earth's equator is not always in the same plane due to precession and nutation.

Sometimes it is useful to have a fixed plane of reference for equatorial coordinates across different calendar dates. In these cases, a fixed *epoch*, or reference time, is helpful. Astronomy Engine provides the J2000 epoch for such cases. This refers to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.

For some other purposes, it is more helpful to represent coordinates using the Earth's equator exactly as it is on that date. For example, when calculating rise/set times or horizontal coordinates, it is most accurate to use the orientation of the Earth's equator at that same date and time. For these uses, Astronomy Engine allows *of-date* calculations. 

| Enum Value | Description |
| --- | --- |
| `EQUATOR_J2000` |  Represent equatorial coordinates in the J2000 epoch.  |
| `EQUATOR_OF_DATE` |  Represent equatorial coordinates using the Earth's equator at the given date and time.  |



---

<a name="astro_node_kind_t"></a>
### `astro_node_kind_t`

**Indicates whether a crossing through the ecliptic plane is ascending or descending.** 



| Enum Value | Description |
| --- | --- |
| `INVALID_NODE` |  Placeholder value for a missing or invalid node.  |
| `ASCENDING_NODE` |  The body passes through the ecliptic plane from south to north.  |
| `DESCENDING_NODE` |  The body passes through the ecliptic plane from north to south.  |



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
| `ASTRO_BAD_TIME` |  The provided date/time is outside the range allowed by this function.  |
| `ASTRO_BAD_VECTOR` |  Vector magnitude is too small to be normalized into a unit vector.  |
| `ASTRO_SEARCH_FAILURE` |  Search was not able to find an ascending root crossing of the function in the specified time interval.  |
| `ASTRO_EARTH_NOT_ALLOWED` |  The Earth cannot be treated as a celestial body seen from an observer on the Earth itself.  |
| `ASTRO_NO_MOON_QUARTER` |  No lunar quarter occurs inside the specified time range.  |
| `ASTRO_WRONG_MOON_QUARTER` |  Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter.  |
| `ASTRO_INTERNAL_ERROR` |  A self-check failed inside the code somewhere, indicating a bug needs to be fixed.  |
| `ASTRO_INVALID_PARAMETER` |  A parameter value passed to a function was not valid.  |
| `ASTRO_FAIL_APSIS` |  Special-case logic for finding Neptune/Pluto apsis failed.  |
| `ASTRO_BUFFER_TOO_SMALL` |  A provided buffer's size is too small to receive the requested data.  |
| `ASTRO_OUT_OF_MEMORY` |  An attempt to allocate memory failed.  |
| `ASTRO_INCONSISTENT_TIMES` |  The provided initial state vectors did not have matching times.  |



---

<a name="astro_time_format_t"></a>
### `astro_time_format_t`

**Selects the output format of the function [`Astronomy_FormatTime`](#Astronomy_FormatTime).** 



| Enum Value | Description |
| --- | --- |
| `TIME_FORMAT_DAY` |  Truncate to UTC calendar date only, e.g. `2020-12-31`. Buffer size must be at least 11 characters.  |
| `TIME_FORMAT_MINUTE` |  Round to nearest UTC minute, e.g. `2020-12-31T15:47Z`. Buffer size must be at least 18 characters.  |
| `TIME_FORMAT_SECOND` |  Round to nearest UTC second, e.g. `2020-12-31T15:47:32Z`. Buffer size must be at least 21 characters.  |
| `TIME_FORMAT_MILLI` |  Round to nearest UTC millisecond, e.g. `2020-12-31T15:47:32.397Z`. Buffer size must be at least 25 characters.  |



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
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `angle` |  An angle expressed in degrees.  |


---

<a name="astro_apsis_t"></a>
### `astro_apsis_t`

**An apsis event: pericenter (closest approach) or apocenter (farthest distance).** 



For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an event where the orbiting body reaches its closest or farthest point from the primary body. The closest approach is called *pericenter* and the farthest point is *apocenter*.

More specific terminology is common for particular orbiting bodies. The Moon's closest approach to the Earth is called *perigee* and its farthest point is called *apogee*. The closest approach of a planet to the Sun is called *perihelion* and the furthest point is called *aphelion*.

This data structure is returned by [`Astronomy_SearchLunarApsis`](#Astronomy_SearchLunarApsis) and [`Astronomy_NextLunarApsis`](#Astronomy_NextLunarApsis) to iterate through consecutive alternating perigees and apogees. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the apsis.  |
| [`astro_apsis_kind_t`](#astro_apsis_kind_t) | `kind` |  Whether this is a pericenter or apocenter event.  |
| `double` | `dist_au` |  The distance between the centers of the bodies in astronomical units.  |
| `double` | `dist_km` |  The distance between the centers of the bodies in kilometers.  |


---

<a name="astro_atmosphere_t"></a>
### `astro_atmosphere_t`

**Information about idealized atmospheric variables at a given elevation.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `pressure` |  Atmospheric pressure in pascals  |
| `double` | `temperature` |  Atmospheric temperature in kelvins  |
| `double` | `density` |  Atmospheric density relative to sea level  |


---

<a name="astro_axis_t"></a>
### `astro_axis_t`

**Information about a body's rotation axis at a given time.** 



This structure is returned by [`Astronomy_RotationAxis`](#Astronomy_RotationAxis) to report the orientation of a body's rotation axis at a given moment in time. The axis is specified by the direction in space that the body's north pole points, using angular equatorial coordinates in the J2000 system (EQJ).

Thus `ra` is the right ascension, and `dec` is the declination, of the body's north pole vector at the given moment in time. The north pole of a body is defined as the pole that lies on the north side of the [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane), regardless of the body's direction of rotation.

The `spin` field indicates the angular position of a prime meridian arbitrarily recommended for the body by the International Astronomical Union (IAU).

The fields `ra`, `dec`, and `spin` correspond to the variables α0, δ0, and W, respectively, from [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).

The field `north` is a unit vector pointing in the direction of the body's north pole. It is expressed in the equatorial J2000 system (EQJ). 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `ra` |  The J2000 right ascension of the body's north pole direction, in sidereal hours.  |
| `double` | `dec` |  The J2000 declination of the body's north pole direction, in degrees.  |
| `double` | `spin` |  Rotation angle of the body's prime meridian, in degrees.  |
| [`astro_vector_t`](#astro_vector_t) | `north` |  A J2000 dimensionless unit vector pointing in the direction of the body's north pole.  |


---

<a name="astro_constellation_t"></a>
### `astro_constellation_t`

**Reports the constellation that a given celestial point lies within.** 



The [`Astronomy_Constellation`](#Astronomy_Constellation) function returns this struct to report which constellation corresponds with a given point in the sky. Constellations are defined with respect to the B1875 equatorial system per IAU standard. Although `Astronomy.Constellation` requires J2000 equatorial coordinates, the struct contains converted B1875 coordinates for reference. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `const char *` | `symbol` |  3-character mnemonic symbol for the constellation, e.g. "Ori".  |
| `const char *` | `name` |  Full name of constellation, e.g. "Orion".  |
| `double` | `ra_1875` |  Right ascension expressed in B1875 coordinates.  |
| `double` | `dec_1875` |  Declination expressed in B1875 coordinates.  |


---

<a name="astro_eclipse_event_t"></a>
### `astro_eclipse_event_t`

**Holds a time and the observed altitude of the Sun at that time.** 



When reporting a solar eclipse observed at a specific location on the Earth (a "local" solar eclipse), a series of events occur. In addition to the time of each event, it is important to know the altitude of the Sun, because each event may be invisible to the observer if the Sun is below the horizon.

If `altitude` is negative, the event is theoretical only; it would be visible if the Earth were transparent, but the observer cannot actually see it. If `altitude` is positive but less than a few degrees, visibility will be impaired by atmospheric interference (sunrise or sunset conditions). 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the event.  |
| `double` | `altitude` |  The angular altitude of the center of the Sun above/below the horizon, at `time`, corrected for atmospheric refraction and expressed in degrees.  |


---

<a name="astro_ecliptic_t"></a>
### `astro_ecliptic_t`

**Ecliptic angular and Cartesian coordinates.** 



Coordinates of a celestial body as seen from the center of the Sun (heliocentric), oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic). 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_vector_t`](#astro_vector_t) | `vec` |  Cartesian ecliptic vector: x=equinox, y=90 degrees prograde in ecliptic plane, z=northward perpendicular to ecliptic.  |
| `double` | `elat` |  Latitude in degrees north (positive) or south (negative) of the ecliptic plane.  |
| `double` | `elon` |  Longitude in degrees around the ecliptic plane prograde from the equinox.  |


---

<a name="astro_elongation_t"></a>
### `astro_elongation_t`

**Contains information about the visibility of a celestial body at a given date and time. See [`Astronomy_Elongation`](#Astronomy_Elongation) for more detailed information about the members of this structure. See also [`Astronomy_SearchMaxElongation`](#Astronomy_SearchMaxElongation) for how to search for maximum elongation events.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation.  |
| [`astro_visibility_t`](#astro_visibility_t) | `visibility` |  Whether the body is best seen in the morning or the evening.  |
| `double` | `elongation` |  The angle in degrees between the body and the Sun, as seen from the Earth.  |
| `double` | `ecliptic_separation` |  The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.  |


---

<a name="astro_equatorial_t"></a>
### `astro_equatorial_t`

**Equatorial angular and cartesian coordinates.** 



Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context), oriented with respect to the projection of the Earth's equator onto the sky. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `ra` |  right ascension in sidereal hours.  |
| `double` | `dec` |  declination in degrees  |
| `double` | `dist` |  distance to the celestial body in AU.  |
| [`astro_vector_t`](#astro_vector_t) | `vec` |  equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north.  |


---

<a name="astro_func_result_t"></a>
### `astro_func_result_t`

**A real value returned by a function whose ascending root is to be found.** 



When calling [`Astronomy_Search`](#Astronomy_Search), the caller must pass in a callback function compatible with the function-pointer type [`astro_search_func_t`](#astro_search_func_t) whose ascending root is to be found. That callback function must return [`astro_func_result_t`](#astro_func_result_t). If the function call is successful, it will set `status` to `ASTRO_SUCCESS` and `value` to the numeric value appropriate for the given date and time. If the call fails for some reason, it should set `status` to an appropriate error value other than `ASTRO_SUCCESS`; in the error case, to guard against any possible misuse of `value`, it is recommended to set `value` to `NAN`, though this is not strictly necessary. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `value` |  The value returned by a function whose ascending root is to be found.  |


---

<a name="astro_global_solar_eclipse_t"></a>
### `astro_global_solar_eclipse_t`

**Reports the time and geographic location of the peak of a solar eclipse.** 



Returned by [`Astronomy_SearchGlobalSolarEclipse`](#Astronomy_SearchGlobalSolarEclipse) or [`Astronomy_NextGlobalSolarEclipse`](#Astronomy_NextGlobalSolarEclipse) to report information about a solar eclipse event. If a solar eclipse is found, `status` holds `ASTRO_SUCCESS` and `kind`, `peak`, and `distance` have valid values. The `latitude` and `longitude` are set only for total and annular eclipses (see more below). If `status` holds any value other than `ASTRO_SUCCESS`, it is an error code; in that case, `kind` holds `ECLIPSE_NONE` and all the other fields are undefined.

The eclipse is classified as partial, annular, or total, depending on the maximum amount of the Sun's disc obscured, as seen at the peak location on the surface of the Earth.

The `kind` field thus holds `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`. A total eclipse is when the peak observer sees the Sun completely blocked by the Moon. An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface to completely block the Sun; instead, the Sun takes on a ring-shaped appearance. A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth observes either a total or annular eclipse.

If `kind` is `ECLIPSE_TOTAL` or `ECLIPSE_ANNULAR`, the `latitude` and `longitude` fields give the geographic coordinates of the center of the Moon's shadow projected onto the daytime side of the Earth at the instant of the eclipse's peak. If `kind` has any other value, `latitude` and `longitude` are undefined and should not be used.

For total or annular eclipses, the `obscuration` field holds the fraction (0, 1] of the Sun's apparent disc area that is blocked from view by the Moon's silhouette, as seen by an observer located at the geographic coordinates `latitude`, `longitude` at the darkest time `peak`. The value will always be 1 for total eclipses, and less than 1 for annular eclipses. For partial eclipses, `obscuration` is undefined and should not be used. This is because there is little practical use for an obscuration value of a partial eclipse without supplying a particular observation location. Developers who wish to find an obscuration value for partial solar eclipses should therefore use [`Astronomy_SearchLocalSolarEclipse`](#Astronomy_SearchLocalSolarEclipse) and provide the geographic coordinates of an observer. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_eclipse_kind_t`](#astro_eclipse_kind_t) | `kind` |  The type of solar eclipse found.  |
| `double` | `obscuration` |  The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only).  |
| [`astro_time_t`](#astro_time_t) | `peak` |  The date and time when the solar eclipse is darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.  |
| `double` | `distance` |  The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.  |
| `double` | `latitude` |  The geographic latitude at the center of the peak eclipse shadow.  |
| `double` | `longitude` |  The geographic longitude at the center of the peak eclipse shadow.  |


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



Returned by the function [`Astronomy_SearchHourAngleEx`](#Astronomy_SearchHourAngleEx) to report information about a celestial body crossing a certain hour angle as seen by a specified topocentric observer. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time when the body crosses the specified hour angle.  |
| [`astro_horizon_t`](#astro_horizon_t) | `hor` |  Apparent coordinates of the body at the time it crosses the specified hour angle.  |


---

<a name="astro_illum_t"></a>
### `astro_illum_t`

**Information about the brightness and illuminated shape of a celestial body.** 



Returned by the functions [`Astronomy_Illumination`](#Astronomy_Illumination) and [`Astronomy_SearchPeakMagnitude`](#Astronomy_SearchPeakMagnitude) to report the visual magnitude and illuminated fraction of a celestial body at a given date and time. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the observation.  |
| `double` | `mag` |  The visual magnitude of the body. Smaller values are brighter.  |
| `double` | `phase_angle` |  The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth.  |
| `double` | `phase_fraction` |  A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth.  |
| `double` | `helio_dist` |  The distance between the Sun and the body at the observation time.  |
| `double` | `ring_tilt` |  For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.  |


---

<a name="astro_jupiter_moons_t"></a>
### `astro_jupiter_moons_t`

**Holds the positions and velocities of Jupiter's major 4 moons.** 



The [`Astronomy_JupiterMoons`](#Astronomy_JupiterMoons) function returns this struct to report position and velocity vectors for Jupiter's largest 4 moons Io, Europa, Ganymede, and Callisto. Each position vector is relative to the center of Jupiter. Both position and velocity are oriented in the EQJ system (that is, using Earth's equator at the J2000 epoch.) The positions are expressed in astronomical units (AU), and the velocities in AU/day. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_state_vector_t`](#astro_state_vector_t) | `io` |  Jovicentric position and velocity of Io.  |
| [`astro_state_vector_t`](#astro_state_vector_t) | `europa` |  Jovicentric position and velocity of Europa.  |
| [`astro_state_vector_t`](#astro_state_vector_t) | `ganymede` |  Jovicentric position and velocity of Ganymede.  |
| [`astro_state_vector_t`](#astro_state_vector_t) | `callisto` |  Jovicentric position and velocity of Callisto.  |


---

<a name="astro_libration_t"></a>
### `astro_libration_t`

**Lunar libration angles, returned by [`Astronomy_Libration`](#Astronomy_Libration).** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `elat` |  Sub-Earth libration ecliptic latitude angle, in degrees.  |
| `double` | `elon` |  Sub-Earth libration ecliptic longitude angle, in degrees.  |
| `double` | `mlat` |  Moon's geocentric ecliptic latitude, in degrees.  |
| `double` | `mlon` |  Moon's geocentric ecliptic longitude, in degrees.  |
| `double` | `dist_km` |  Distance between the centers of the Earth and Moon in kilometers.  |
| `double` | `diam_deg` |  The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.  |


---

<a name="astro_local_solar_eclipse_t"></a>
### `astro_local_solar_eclipse_t`

**Information about a solar eclipse as seen by an observer at a given time and geographic location.** 



Returned by [`Astronomy_SearchLocalSolarEclipse`](#Astronomy_SearchLocalSolarEclipse) or [`Astronomy_NextLocalSolarEclipse`](#Astronomy_NextLocalSolarEclipse) to report information about a solar eclipse as seen at a given geographic location. If a solar eclipse is found, `status` holds `ASTRO_SUCCESS` and the other fields are set. If `status` holds any other value, it is an error code and the other fields are undefined.

When a solar eclipse is found, it is classified as partial, annular, or total. The `kind` field thus holds `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`. A partial solar eclipse is when the Moon does not line up directly enough with the Sun to completely block the Sun's light from reaching the observer. An annular eclipse occurs when the Moon's disc is completely visible against the Sun but the Moon is too far away to completely block the Sun's light; this leaves the Sun with a ring-like appearance. A total eclipse occurs when the Moon is close enough to the Earth and aligned with the Sun just right to completely block all sunlight from reaching the observer.

The `obscuration` field reports what fraction of the Sun's disc appears blocked by the Moon when viewed by the observer at the peak eclipse time. This is a value that ranges from 0 (no blockage) to 1 (total eclipse). The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses. The value will be exactly 1 for total eclipses. Obscuration gives an indication of how dark the eclipse appears.

There are 5 "event" fields, each of which contains a time and a solar altitude. Field `peak` holds the date and time of the center of the eclipse, when it is at its peak. The fields `partial_begin` and `partial_end` are always set, and indicate when the eclipse begins/ends. If the eclipse reaches totality or becomes annular, `total_begin` and `total_end` indicate when the total/annular phase begins/ends. When an event field is valid, the caller must also check its `altitude` field to see whether the Sun is above the horizon at that time. See [`astro_eclipse_kind_t`](#astro_eclipse_kind_t) for more information. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_eclipse_kind_t`](#astro_eclipse_kind_t) | `kind` |  The type of solar eclipse found: `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`.  |
| `double` | `obscuration` |  The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak.  |
| [`astro_eclipse_event_t`](#astro_eclipse_event_t) | `partial_begin` |  The time and Sun altitude at the beginning of the eclipse.  |
| [`astro_eclipse_event_t`](#astro_eclipse_event_t) | `total_begin` |  If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid.  |
| [`astro_eclipse_event_t`](#astro_eclipse_event_t) | `peak` |  The time and Sun altitude when the eclipse reaches its peak.  |
| [`astro_eclipse_event_t`](#astro_eclipse_event_t) | `total_end` |  If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid.  |
| [`astro_eclipse_event_t`](#astro_eclipse_event_t) | `partial_end` |  The time and Sun altitude at the end of the eclipse.  |


---

<a name="astro_lunar_eclipse_t"></a>
### `astro_lunar_eclipse_t`

**Information about a lunar eclipse.** 



Returned by [`Astronomy_SearchLunarEclipse`](#Astronomy_SearchLunarEclipse) or [`Astronomy_NextLunarEclipse`](#Astronomy_NextLunarEclipse) to report information about a lunar eclipse event. If a lunar eclipse is found, `status` holds `ASTRO_SUCCESS` and the other fields are set. If `status` holds any other value, it is an error code and the other fields are undefined.

When a lunar eclipse is found, it is classified as penumbral, partial, or total. Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed by the Earth's penumbra; no part of the Moon touches the Earth's umbra. Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra. Total eclipses occur when the entire Moon passes into the Earth's umbra.

The `kind` field thus holds `ECLIPSE_PENUMBRAL`, `ECLIPSE_PARTIAL`, or `ECLIPSE_TOTAL`, depending on the kind of lunar eclipse found.

The `obscuration` field holds a value in the range [0, 1] that indicates what fraction of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak. This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses, the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.

Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.

Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase of the eclipse, which is half of the amount of time the eclipse spends in each phase (expressed in minutes), or 0 if the eclipse never reaches that phase. By converting from minutes to days, and subtracting/adding with `center`, the caller may determine the date and time of the beginning/end of each eclipse phase. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_eclipse_kind_t`](#astro_eclipse_kind_t) | `kind` |  The type of lunar eclipse found.  |
| `double` | `obscuration` |  The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra.  |
| [`astro_time_t`](#astro_time_t) | `peak` |  The time of the eclipse at its peak.  |
| `double` | `sd_penum` |  The semi-duration of the penumbral phase in minutes.  |
| `double` | `sd_partial` |  The semi-duration of the partial phase in minutes, or 0.0 if none.  |
| `double` | `sd_total` |  The semi-duration of the total phase in minutes, or 0.0 if none.  |


---

<a name="astro_moon_quarter_t"></a>
### `astro_moon_quarter_t`

**A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `int` | `quarter` |  0=new moon, 1=first quarter, 2=full moon, 3=third quarter.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The date and time of the lunar quarter.  |


---

<a name="astro_node_event_t"></a>
### `astro_node_event_t`

**Information about an ascending or descending node of a body.** 



This structure is returned by [`Astronomy_SearchMoonNode`](#Astronomy_SearchMoonNode) and [`Astronomy_NextMoonNode`](#Astronomy_NextMoonNode) to report information about the center of the Moon passing through the ecliptic plane. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The time when the body passes through the ecliptic plane.  |
| [`astro_node_kind_t`](#astro_node_kind_t) | `kind` |  Either `ASCENDING_NODE` or `DESCENDING_NODE`, depending on the direction of the ecliptic plane crossing.  |


---

<a name="astro_observer_t"></a>
### `astro_observer_t`

**The location of an observer on (or near) the surface of the Earth.** 



This structure is passed to functions that calculate phenomena as observed from a particular place on the Earth.

You can create this structure directly, or you can call the convenience function [`Astronomy_MakeObserver`](#Astronomy_MakeObserver) to create one for you. 

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
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `rot` |  A normalized 3x3 rotation matrix.  |


---

<a name="astro_search_result_t"></a>
### `astro_search_result_t`

**The result of a search for an astronomical event.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `time` |  The time at which a searched-for event occurs.  |


---

<a name="astro_seasons_t"></a>
### `astro_seasons_t`

**The dates and times of changes of season for a given calendar year. Call [`Astronomy_Seasons`](#Astronomy_Seasons) to calculate this data structure for a given year.** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
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
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `lat` |  The latitude angle: -90..+90 degrees.  |
| `double` | `lon` |  The longitude angle: 0..360 degrees.  |
| `double` | `dist` |  Distance in AU.  |


---

<a name="astro_state_vector_t"></a>
### `astro_state_vector_t`

**A state vector that contains a position (AU) and velocity (AU/day).** 



| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `x` |  The Cartesian position x-coordinate of the vector in AU.  |
| `double` | `y` |  The Cartesian position y-coordinate of the vector in AU.  |
| `double` | `z` |  The Cartesian position z-coordinate of the vector in AU.  |
| `double` | `vx` |  The Cartesian velocity x-coordinate of the vector in AU/day.  |
| `double` | `vy` |  The Cartesian velocity y-coordinate of the vector in AU/day.  |
| `double` | `vz` |  The Cartesian velocity z-coordinate of the vector in AU/day.  |
| [`astro_time_t`](#astro_time_t) | `t` |  The date and time at which this state vector is valid.  |


---

<a name="astro_time_t"></a>
### `astro_time_t`

**A date and time used for astronomical calculations.** 



This type is of fundamental importance to Astronomy Engine. It is used to represent dates and times for all astronomical calculations. It is also included in the values returned by many Astronomy Engine functions.

To create a valid [`astro_time_t`](#astro_time_t) value from scratch, call [`Astronomy_MakeTime`](#Astronomy_MakeTime) (for a given calendar date and time) or [`Astronomy_CurrentTime`](#Astronomy_CurrentTime) (for the system's current date and time).

To adjust an existing [`astro_time_t`](#astro_time_t) by a certain real number of days, call [`Astronomy_AddDays`](#Astronomy_AddDays).

The [`astro_time_t`](#astro_time_t) type contains `ut` to represent Universal Time (UT1/UTC) and `tt` to represent Terrestrial Time (TT, also known as *ephemeris time*). The difference `tt-ut` is known as *&Delta;T*, using a best-fit piecewise model devised by [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).

Both `tt` and `ut` are necessary for performing different astronomical calculations. Indeed, certain calculations (such as rise/set times) require both time scales. See the documentation for the `ut` and `tt` fields for more detailed information.

In cases where `[`astro_time_t`](#astro_time_t)` is included in a structure returned by a function that can fail, the `astro_status_t` field `status` will contain a value other than `ASTRO_SUCCESS`; in that case the `ut` and `tt` will hold `NAN` (not a number). In general, when there is an error code stored in a struct field `status`, the caller should ignore all other values in that structure, including the `ut` and `tt` inside `[`astro_time_t`](#astro_time_t)`. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| `double` | `ut` | **UT1/UTC number of days since noon on January 1, 2000.**  The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine. Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed. UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed. The value in `ut` is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time. Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*.  |
| `double` | `tt` | **Terrestrial Time days since noon on January 1, 2000.**  Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments for changes in the Earth's rotation. The value in `tt` is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth. Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).  |
| `double` | `psi` | **For internal use only. Used to optimize Earth tilt calculations.**  |
| `double` | `eps` | **For internal use only. Used to optimize Earth tilt calculations.**  |
| `double` | `st` | **For internal use only. Lazy-caches sidereal time (Earth rotation).**  |


---

<a name="astro_transit_t"></a>
### `astro_transit_t`

**Information about a transit of Mercury or Venus, as seen from the Earth.** 



Returned by [`Astronomy_SearchTransit`](#Astronomy_SearchTransit) or [`Astronomy_NextTransit`](#Astronomy_NextTransit) to report information about a transit of Mercury or Venus. A transit is when Mercury or Venus passes between the Sun and Earth so that the other planet is seen in silhouette against the Sun.

The `start` field reports the moment in time when the planet first becomes visible against the Sun in its background. The `peak` field reports when the planet is most aligned with the Sun, as seen from the Earth. The `finish` field reports the last moment when the planet is visible against the Sun in its background.

The calculations are performed from the point of view of a geocentric observer. 

| Type | Member | Description |
| ---- | ------ | ----------- |
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| [`astro_time_t`](#astro_time_t) | `start` |  Date and time at the beginning of the transit.  |
| [`astro_time_t`](#astro_time_t) | `peak` |  Date and time of the peak of the transit.  |
| [`astro_time_t`](#astro_time_t) | `finish` |  Date and time at the end of the transit.  |
| `double` | `separation` |  Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`.  |


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
| [`astro_status_t`](#astro_status_t) | `status` |  `ASTRO_SUCCESS` if this struct is valid; otherwise an error code.  |
| `double` | `x` |  The Cartesian x-coordinate of the vector in AU.  |
| `double` | `y` |  The Cartesian y-coordinate of the vector in AU.  |
| `double` | `z` |  The Cartesian z-coordinate of the vector in AU.  |
| [`astro_time_t`](#astro_time_t) | `t` |  The date and time at which this vector is valid.  |

<a name="typedefs"></a>
## Type Definitions



---

<a name="astro_deltat_func"></a>
### `astro_deltat_func`

`typedef double(* astro_deltat_func) (double ut);`

**A pointer to a function that calculates Delta T.** 



Delta T is the discrepancy between times measured using an atomic clock and times based on observations of the Earth's rotation, which is gradually slowing down over time. Delta T = TT - UT, where TT = Terrestrial Time, based on atomic time, and UT = Universal Time, civil time based on the Earth's rotation. Astronomy Engine defaults to using a Delta T function defined by Espenak and Meeus in their "Five Millennium Canon of Solar Eclipses". See: [https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html)

---

<a name="astro_grav_sim_t"></a>
### `astro_grav_sim_t`

`typedef struct astro_grav_sim_s astro_grav_sim_t;`

**A data type used for managing simulation of the gravitational forces on a small body.** 



This is an opaque data type used to hold the internal state of a numeric integrator used to calculate the trajectory of a small body moving through the Solar System. 

---

<a name="astro_position_func_t"></a>
### `astro_position_func_t`

`typedef astro_vector_t(* astro_position_func_t) (void *context, astro_time_t time);`

**A function for which to solve a light-travel time problem.** 



The function [`Astronomy_CorrectLightTravel`](#Astronomy_CorrectLightTravel) solves a generalized problem of deducing how far in the past light must have left a target object to be seen by an observer at a specified time. This function pointer type expresses an arbitrary position vector as function of time. Such a function must be passed to `Astronomy_CorrectLightTravel`. 

---

<a name="astro_search_func_t"></a>
### `astro_search_func_t`

`typedef astro_func_result_t(* astro_search_func_t) (void *context, astro_time_t time);`

**A pointer to a function that is to be passed as a callback to [`Astronomy_Search`](#Astronomy_Search).** 



The function [`Astronomy_Search`](#Astronomy_Search) numerically solves for the time that a given event occurs. An event is defined as the time when an arbitrary function transitions between having a negative value and a non-negative value. This transition is called an *ascending root*.

The type astro_search_func_t represents such a callback function that accepts a custom `context` pointer and an [`astro_time_t`](#astro_time_t) representing the time to probe. The function returns an [`astro_func_result_t`](#astro_func_result_t) that contains either a real number in `value` or an error code in `status` that aborts the search.

The `context` points to some data whose type varies depending on the callback function. It can contain any auxiliary parameters (other than time) needed to evaluate the function. For example, a function may pertain to a specific celestial body, in which case `context` may point to a value of type astro_body_t. The `context` parameter is supplied by the caller of [`Astronomy_Search`](#Astronomy_Search), which passes it along to every call to the callback function. If the caller of `Astronomy_Search` knows that the callback function does not need a context, it is safe to pass `NULL` as the context pointer. 