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
| [Rotation_EQJ_ECT](#Astronomy_Rotation_EQJ_ECT) | Calculates a rotation matrix from J2000 mean equator (EQJ) to  true ecliptic of date (ECT). |
| [Rotation_EQJ_ECL](#Astronomy_Rotation_EQJ_ECL) | Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL). |
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

