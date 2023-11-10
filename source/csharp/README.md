# Astronomy Engine (C#)

This is the complete programming reference for the C# version of Astronomy Engine.
Other programming languages are supported.
See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start

[![NuGet](https://img.shields.io/nuget/v/CosineKitty.AstronomyEngine)](https://www.nuget.org/packages/CosineKitty.AstronomyEngine/)

[Astronomy Engine is available as a NuGet package](https://www.nuget.org/packages/CosineKitty.AstronomyEngine/)
that targets .NET Standard 2.0. This means you can use the same NuGet package in your
.NET Framework 4+ or .NET Core 5+ projects.

Alternatively, you can include Astronomy Engine in your own C# program
by copying the source file `astronomy.cs` from this directory into your own project.

To get started quickly, here are some [examples](../../demo/csharp/).

---

## Contents

- [Topic Index](#topics)
- [Constants](#constants)
- [Functions](#functions)
- [Types](#types)

---

<a name="topics"></a>
## Topic Index

### Position of Sun, Moon, and planets

| Function | Description |
| -------- | ----------- |
| [HelioVector](#Astronomy.HelioVector) | Calculates body position vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy.GeoVector)     | Calculates body position vector with respect to the center of the Earth. |
| [Equator](#Astronomy.Equator)         | Calculates right ascension and declination. |
| [EquatorialToEcliptic](#Astronomy.EquatorialToEcliptic) | Converts J2000 mean equator (EQJ) coordinates to true ecliptic of date (ECT) coordinates. |
| [EclipticLongitude](#Astronomy.EclipticLongitude) | Calculates true ecliptic of date (ECT) longitude for a body. |
| [Horizon](#Astronomy.Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |
| [PairLongitude](#Astronomy.PairLongitude) | Calculates the difference in apparent ecliptic longitude between two bodies, as seen from the Earth. |
| [BaryState](#Astronomy.BaryState) | Calculates the barycentric position and velocity vectors of the Sun or a planet. |

### Geographic helper functions

| Function | Description |
| -------- | ----------- |
| [ObserverVector](#Astronomy.ObserverVector) | Calculates a vector from the center of the Earth to an observer on the Earth's surface. |
| [VectorObserver](#Astronomy.VectorObserver) | Calculates the geographic coordinates for a geocentric equatorial vector. |
### Rise, set, and culmination times

| Function | Description |
| -------- | ----------- |
| [SearchRiseSet](#Astronomy.SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchAltitude](#Astronomy.SearchAltitude) | Finds time when a body reaches a given altitude above or below the horizon. Useful for finding civil, nautical, or astronomical twilight. |
| [SearchHourAngle](#Astronomy.SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| Function | Description |
| -------- | ----------- |
| [MoonPhase](#Astronomy.MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonPhase](#Astronomy.SearchMoonPhase) | Finds the next instance of the Moon reaching a specific ecliptic longitude separation from the Sun. |
| [SearchMoonQuarter](#Astronomy.SearchMoonQuarter) | Finds the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy.NextMoonQuarter) | Finds the next quarter moon phase after a previous one that has been found. |

### Eclipses and Transits

| Function | Description |
| -------- | ----------- |
| [SearchLunarEclipse](#Astronomy.SearchLunarEclipse) | Search for the first lunar eclipse after a given date. |
| [NextLunarEclipse](#Astronomy.NextLunarEclipse) | Continue searching for more lunar eclipses. |
| [SearchGlobalSolarEclipse](#Astronomy.SearchGlobalSolarEclipse) | Search for the first solar eclipse after a given date that is visible anywhere on the Earth. |
| [NextGlobalSolarEclipse](#Astronomy.NextGlobalSolarEclipse) | Continue searching for solar eclipses visible anywhere on the Earth. |
| [SearchLocalSolarEclipse](#Astronomy.SearchLocalSolarEclipse) | Search for the first solar eclipse after a given date that is visible at a particular location on the Earth. |
| [NextLocalSolarEclipse](#Astronomy.NextLocalSolarEclipse) | Continue searching for solar eclipses visible at a particular location on the Earth. |
| [SearchTransit](#Astronomy.SearchTransit) | Search for the next transit of Mercury or Venus. |
| [NextTransit](#Astronomy.NextTransit) | Continue searching for transits of Mercury or Venus. |

### Lunar perigee and apogee

| Function | Description |
| -------- | ----------- |
| [SearchLunarApsis](#Astronomy.SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy.NextLunarApsis) | Given an already-found apsis, finds the next perigee or apogee of the Moon. |

### Planet perihelion and aphelion

| Function | Description |
| -------- | ----------- |
| [SearchPlanetApsis](#Astronomy.SearchPlanetApsis) | Finds the next perihelion or aphelion of a planet after a specified date. |
| [NextPlanetApsis](#Astronomy.NextPlanetApsis) | Given an already-found apsis, finds the next perihelion or aphelion of a planet. |

### Visual magnitude and elongation

| Function | Description |
| -------- | ----------- |
| [Illumination](#Astronomy.Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#Astronomy.SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#Astronomy.AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Astronomy.Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#Astronomy.SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| Function | Description |
| -------- | ----------- |
| [SearchRelativeLongitude](#Astronomy.SearchRelativeLongitude) | Finds oppositions and conjunctions of planets. |

### Equinoxes, solstices, and apparent solar motion

| Function | Description |
| -------- | ----------- |
| [SearchSunLongitude](#Astronomy.SearchSunLongitude) | Finds the next time the Sun reaches a specified apparent ecliptic longitude in the true ecliptic of date (ECT) system. |
| [Seasons](#Astronomy.Seasons) | Finds the equinoxes and solstices for a given calendar year. |
| [SunPosition](#Astronomy.SunPosition) | Calculates the Sun's apparent true ecliptic of date (ECT) coordinates as seen from the Earth. |

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
| [RotateVector](#Astronomy.RotateVector) | Applies a rotation matrix to a vector, yielding a vector in another orientation system. |
| [InverseRotation](#Astronomy.InverseRotation) | Given a rotation matrix, finds the inverse rotation matrix that does the opposite transformation. |
| [CombineRotation](#Astronomy.CombineRotation) | Given two rotation matrices, returns a rotation matrix that combines them into a net transformation. |
| [IdentityMatrix](#Astronomy.IdentityMatrix) | Returns a 3x3 identity matrix, which can be used to form other rotation matrices. |
| [Pivot](#Astronomy.Pivot) | Transforms a rotation matrix by pivoting it around a given axis by a given angle. |
| [VectorFromSphere](#Astronomy.VectorFromSphere) | Converts spherical coordinates to Cartesian coordinates. |
| [SphereFromVector](#Astronomy.SphereFromVector) | Converts Cartesian coordinates to spherical coordinates. |
| [EquatorFromVector](#Astronomy.EquatorFromVector) | Given an equatorial vector, calculates equatorial angular coordinates. |
| [VectorFromHorizon](#Astronomy.VectorFromHorizon) | Given apparent angular horizontal coordinates, calculates horizontal vector. |
| [HorizonFromVector](#Astronomy.HorizonFromVector) | Given a vector in horizontal orientation, calculates horizontal angular coordinates. |
| [Rotation_EQD_EQJ](#Astronomy.Rotation_EQD_EQJ) | Calculates a rotation matrix from true equator of date (EQD) to J2000 mean equator (EQJ). |
| [Rotation_EQD_ECT](#Astronomy.Rotation_EQD_ECT) | Calculates a rotation matrix from true equator of date (EQD) to true ecliptic of date (ECT). |
| [Rotation_EQD_ECL](#Astronomy.Rotation_EQD_ECL) | Calculates a rotation matrix from true equator of date (EQD) to J2000 mean ecliptic (ECL). |
| [Rotation_EQD_HOR](#Astronomy.Rotation_EQD_HOR) | Calculates a rotation matrix from true equator of date (EQD) to horizontal (HOR). |
| [Rotation_EQJ_EQD](#Astronomy.Rotation_EQJ_EQD) | Calculates a rotation matrix from J2000 mean equator (EQJ) to true equator of date (EQD). |
| [Rotation_EQJ_ECT](#Astronomy.Rotation_EQJ_ECT) | Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT). |
| [Rotation_EQJ_ECL](#Astronomy.Rotation_EQJ_ECL) | Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL). |
| [Rotation_EQJ_HOR](#Astronomy.Rotation_EQJ_HOR) | Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR). |
| [Rotation_ECT_EQD](#Astronomy.Rotation_ECT_EQD) | Calculates a rotation matrix from true ecliptic of date (ECT) to true equator of date (EQD). |
| [Rotation_ECT_EQJ](#Astronomy.Rotation_ECT_EQJ) | Calculates a rotation matrix from true ecliptic of date (ECT) J2000 mean equator (EQJ). |
| [Rotation_ECL_EQD](#Astronomy.Rotation_ECL_EQD) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to true equator of date (EQD). |
| [Rotation_ECL_EQJ](#Astronomy.Rotation_ECL_EQJ) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ). |
| [Rotation_ECL_HOR](#Astronomy.Rotation_ECL_HOR) | Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR). |
| [Rotation_HOR_EQD](#Astronomy.Rotation_HOR_EQD) | Calculates a rotation matrix from horizontal (HOR) to true equator of date (EQD). |
| [Rotation_HOR_EQJ](#Astronomy.Rotation_HOR_EQJ) | Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ). |
| [Rotation_HOR_ECL](#Astronomy.Rotation_HOR_ECL) | Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL). |
| [Rotation_EQJ_GAL](#Astronomy.Rotation_EQJ_GAL) | Calculates a rotation matrix from J2000 mean equator (EQJ) to galactic (GAL). |
| [Rotation_GAL_EQJ](#Astronomy.Rotation_GAL_EQJ) | Calculates a rotation matrix from galactic (GAL) to J2000 mean equator (EQJ). |

### Gravitational simulation of small bodies

Astronomy Engine provides a [GravitySimulator](#GravitySimulator) class
that allows you to model the trajectories of one or more small bodies like asteroids,
comets, or coasting spacecraft. If you know an initial position vector
and velocity vector for a small body, the gravity simulator can incrementally
simulate the pull of gravity on it from the Sun and planets, to calculate its
movement through the Solar System.

---

<a name="constants"></a>
## Constants

---

<a name="Astronomy.AU_PER_LY"></a>
### `const double Astronomy.AU_PER_LY = 63241.07708807546;`

**The number of astronomical units in one light-year.**


---

<a name="Astronomy.C_AUDAY"></a>
### `const double Astronomy.C_AUDAY = 173.1446326846693;`

**The speed of light in AU/day.**


---

<a name="Astronomy.CALLISTO_RADIUS_KM"></a>
### `const double Astronomy.CALLISTO_RADIUS_KM = 2410.3;`

**The mean radius of Jupiter's moon Callisto, expressed in kilometers.**


---

<a name="Astronomy.DEG2RAD"></a>
### `const double Astronomy.DEG2RAD = 0.017453292519943295;`

**The factor to convert degrees to radians = pi/180.**


---

<a name="Astronomy.EUROPA_RADIUS_KM"></a>
### `const double Astronomy.EUROPA_RADIUS_KM = 1560.8;`

**The mean radius of Jupiter's moon Europa, expressed in kilometers.**


---

<a name="Astronomy.GANYMEDE_RADIUS_KM"></a>
### `const double Astronomy.GANYMEDE_RADIUS_KM = 2631.2;`

**The mean radius of Jupiter's moon Ganymede, expressed in kilometers.**


---

<a name="Astronomy.HOUR2RAD"></a>
### `const double Astronomy.HOUR2RAD = 0.26179938779914946;`

**The factor to convert sidereal hours to radians = pi/12.**


---

<a name="Astronomy.IO_RADIUS_KM"></a>
### `const double Astronomy.IO_RADIUS_KM = 1821.6;`

**The mean radius of Jupiter's moon Io, expressed in kilometers.**


---

<a name="Astronomy.JUPITER_EQUATORIAL_RADIUS_KM"></a>
### `const double Astronomy.JUPITER_EQUATORIAL_RADIUS_KM = 71492;`

**The equatorial radius of Jupiter, expressed in kilometers.**


---

<a name="Astronomy.JUPITER_MEAN_RADIUS_KM"></a>
### `const double Astronomy.JUPITER_MEAN_RADIUS_KM = 69911;`

**The volumetric mean radius of Jupiter, expressed in kilometers.**


---

<a name="Astronomy.JUPITER_POLAR_RADIUS_KM"></a>
### `const double Astronomy.JUPITER_POLAR_RADIUS_KM = 66854;`

**The polar radius of Jupiter, expressed in kilometers.**


---

<a name="Astronomy.KM_PER_AU"></a>
### `const double Astronomy.KM_PER_AU = 149597870.69098932;`

**The number of kilometers in one astronomical unit (AU).**


---

<a name="Astronomy.RAD2DEG"></a>
### `const double Astronomy.RAD2DEG = 57.29577951308232;`

**The factor to convert radians to degrees = 180/pi.**


---

<a name="Astronomy.RAD2HOUR"></a>
### `const double Astronomy.RAD2HOUR = 3.819718634205488;`

**The factor to convert radians to sidereal hours = 12/pi.**


---

<a name="functions"></a>
## Functions

---

<a name="Astronomy.AngleBetween"></a>
### Astronomy.AngleBetween(a, b) &#8658; `double`

**Calculates the angle in degrees between two vectors.**

Given a pair of vectors, this function returns the angle in degrees
between the two vectors in 3D space.
The angle is measured in the plane that contains both vectors.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `a` | The first of a pair of vectors between which to measure an angle. |
| [`AstroVector`](#AstroVector) | `b` | The second of a pair of vectors between which to measure an angle. |

**Returns:** The angle between the two vectors expressed in degrees. The value is in the range [0, 180].

<a name="Astronomy.AngleFromSun"></a>
### Astronomy.AngleFromSun(body, time) &#8658; `double`

**Returns the angle between the given body and the Sun, as seen from the Earth.**

This function calculates the angular separation between the given body and the Sun,
as seen from the center of the Earth. This angle is helpful for determining how
easy it is to see the body away from the glare of the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose angle from the Sun is to be measured. Not allowed to be `Body.Earth`. |
| [`AstroTime`](#AstroTime) | `time` | The time at which the observation is made. |

**Returns:** Returns the angle in degrees between the Sun and the specified body as seen from the center of the Earth.

<a name="Astronomy.Atmosphere"></a>
### Astronomy.Atmosphere(elevationMeters) &#8658; [`AtmosphereInfo`](#AtmosphereInfo)

**Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.**

This function calculates idealized values of pressure, temperature, and density
using the U.S. Standard Atmosphere (1976) model.
1. COESA, U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, DC, 1976.
2. Jursa, A. S., Ed., Handbook of Geophysics and the Space Environment, Air Force Geophysics Laboratory, 1985.
See:
https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml
https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `elevationMeters` | The elevation above sea level at which to calculate atmospheric variables. Must be in the range -500 to +100000, or an exception will occur. |

<a name="Astronomy.BackdatePosition"></a>
### Astronomy.BackdatePosition(time, observerBody, targetBody, aberration) &#8658; [`AstroVector`](#AstroVector)

**Solve for light travel time correction of apparent position.**

When observing a distant object, for example Jupiter as seen from Earth,
the amount of time it takes for light to travel from the object to the
observer can significantly affect the object's apparent position.

This function solves the light travel time correction for the apparent
relative position vector of a target body as seen by an observer body
at a given observation time.

For geocentric calculations, [`Astronomy.GeoVector`](#Astronomy.GeoVector) also includes light
travel time correction, but the time `t` embedded in its returned vector
refers to the observation time, not the backdated time that light left
the observed body. Thus `BackdatePosition` provides direct
access to the light departure time for callers that need it.

For a more generalized light travel correction solver, see [`Astronomy.CorrectLightTravel`](#Astronomy.CorrectLightTravel).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The time of observation. |
| [`Body`](#Body) | `observerBody` | The body to be used as the observation location. |
| [`Body`](#Body) | `targetBody` | The body to be observed. |
| [`Aberration`](#Aberration) | `aberration` | `Aberration.Corrected` to correct for aberration, or `Aberration.None` to leave uncorrected. |

**Returns:** The position vector at the solved backdated time. Its `t` field holds the time that light left the observed body to arrive at the observer at the observation time.

<a name="Astronomy.BaryState"></a>
### Astronomy.BaryState(body, time) &#8658; [`StateVector`](#StateVector)

**Calculates barycentric position and velocity vectors for the given body.**

Given a body and a time, calculates the barycentric position and velocity
vectors for the center of that body at that time.
The vectors are expressed in J2000 mean equator coordinates (EQJ).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose barycentric state vector is to be calculated. Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets: `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate position and velocity. |

**Returns:** A structure that contains barycentric position and velocity vectors.

<a name="Astronomy.CombineRotation"></a>
### Astronomy.CombineRotation(a, b) &#8658; [`RotationMatrix`](#RotationMatrix)

**Creates a rotation based on applying one rotation followed by another.**

Given two rotation matrices, returns a combined rotation matrix that is
equivalent to rotating based on the first matrix, followed by the second.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `a` | The first rotation to apply. |
| [`RotationMatrix`](#RotationMatrix) | `b` | The second rotation to apply. |

**Returns:** The combined rotation matrix.

<a name="Astronomy.Constellation"></a>
### Astronomy.Constellation(ra, dec) &#8658; [`ConstellationInfo`](#ConstellationInfo)

**Determines the constellation that contains the given point in the sky.**

Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
constellation that contains that point.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ra` | The right ascension (RA) of a point in the sky, using the J2000 equatorial system (EQJ). |
| `double` | `dec` | The declination (DEC) of a point in the sky, using the J2000 equatorial system (EQJ). |

**Returns:** A structure that contains the 3-letter abbreviation and full name of the constellation that contains the given (ra,dec), along with the converted B1875 (ra,dec) for that point.

<a name="Astronomy.CorrectLightTravel"></a>
### Astronomy.CorrectLightTravel(func, time) &#8658; [`AstroVector`](#AstroVector)

**Solve for light travel time of a vector function.**

When observing a distant object, for example Jupiter as seen from Earth,
the amount of time it takes for light to travel from the object to the
observer can significantly affect the object's apparent position.
This function is a generic solver that figures out how long in the
past light must have left the observed object to reach the observer
at the specified observation time. It uses [`IPositionFunction`](#IPositionFunction)
to express an arbitrary position vector as a function of time.

This function repeatedly calls `func.Position`, passing a series of time
estimates in the past. Then `func.Position` must return a relative state vector between
the observer and the target. `CorrectLightTravel` keeps calling
`func.Position` with more and more refined estimates of the time light must have
left the target to arrive at the observer.

For common use cases, it is simpler to use [`Astronomy.BackdatePosition`](#Astronomy.BackdatePosition)
for calculating the light travel time correction of one body observing another body.

For geocentric calculations, [`Astronomy.GeoVector`](#Astronomy.GeoVector) also backdates the returned
position vector for light travel time, only it returns the observation time in
the returned vector's `t` field rather than the backdated time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`IPositionFunction`](#IPositionFunction) | `func` | An arbitrary position vector as a function of time. |
| [`AstroTime`](#AstroTime) | `time` | The observation time for which to solve for light travel delay. |

**Returns:** The position vector at the solved backdated time. The `t` field holds the time that light left the observed body to arrive at the observer at the observation time.

<a name="Astronomy.DefineStar"></a>
### Astronomy.DefineStar(body, ra, dec, distanceLightYears) &#8658; `void`

**Assign equatorial coordinates to a user-defined star.**

Some Astronomy Engine functions allow their `body` parameter to
be a user-defined fixed point in the sky, loosely called a "star".
This function assigns a right ascension, declination, and distance
to one of the eight user-defined stars `Body.Star1`..`Body.Star8`.

Stars are not valid until defined. Once defined, they retain their
definition until re-defined by another call to `DefineStar`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the eight user-defined star identifiers: `Body.Star1`, `Body.Star2`, ..., `Body.Star8`. |
| `double` | `ra` | The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of sidereal hours, and must be within the half-open range [0, 24). |
| `double` | `dec` | The declination to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of degrees north (positive) or south (negative) of the J2000 equator, and must be within the closed range [-90, +90]. |
| `double` | `distanceLightYears` | The distance between the star and the Sun, expressed in light-years. This value is used to calculate the tiny parallax shift as seen by an observer on Earth. If you don't know the distance to the star, using a large value like 1000 will generally work well. The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations. |

<a name="Astronomy.DeltaT_EspenakMeeus"></a>
### Astronomy.DeltaT_EspenakMeeus(ut) &#8658; `double`

**The default Delta T function used by Astronomy Engine.**

Espenak and Meeus use a series of piecewise polynomials to
approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
This is the default Delta T function used by Astronomy Engine.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` | The floating point number of days since noon UTC on January 1, 2000. |

**Returns:** The estimated difference TT-UT on the given date, expressed in seconds.

<a name="Astronomy.EclipticGeoMoon"></a>
### Astronomy.EclipticGeoMoon(time) &#8658; [`Spherical`](#Spherical)

**Calculates spherical ecliptic geocentric position of the Moon.**

Given a time of observation, calculates the Moon's geocentric position
in ecliptic spherical coordinates. Provides the ecliptic latitude and
longitude in degrees, and the geocentric distance in astronomical units (AU).

The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and
equatorial plane at the specified time. This means the Earth's equator
is corrected for precession and nutation, and the plane of the Earth's
orbit is corrected for gradual obliquity drift.

This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
It is adapted from Turbo Pascal code from the book
[Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
by Montenbruck and Pfleger.

To calculate a J2000 mean equator vector instead, use [`Astronomy.GeoMoon`](#Astronomy.GeoMoon).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the Moon's position. |

<a name="Astronomy.EclipticLongitude"></a>
### Astronomy.EclipticLongitude(body, time) &#8658; `double`

**Calculates heliocentric ecliptic longitude of a body.**

This function calculates the angle around the plane of the Earth's orbit
of a celestial body, as seen from the center of the Sun.
The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
in degrees from the true equinox of date. The ecliptic longitude is always in the range [0, 360).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body other than the Sun. |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the body's ecliptic longitude is to be calculated. |

**Returns:** Returns the ecliptic longitude in degrees of the given body at the given time.

<a name="Astronomy.Elongation"></a>
### Astronomy.Elongation(body, time) &#8658; [`ElongationInfo`](#ElongationInfo)

**Determines visibility of a celestial body relative to the Sun, as seen from the Earth.**

This function returns an [`ElongationInfo`](#ElongationInfo) structure, which provides the following
information about the given celestial body at the given time:

- `visibility` is an enumerated type that specifies whether the body is more easily seen
   in the morning before sunrise, or in the evening after sunset.

- `elongation` is the angle in degrees between two vectors: one from the center of the Earth to the
   center of the Sun, the other from the center of the Earth to the center of the specified body.
   This angle indicates how far away the body is from the glare of the Sun.
   The elongation angle is always in the range [0, 180].

- `ecliptic_separation` is the absolute value of the difference between the body's ecliptic longitude
  and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures
  around the plane of the Earth's orbit, and ignores how far above or below that plane the body is.
  The ecliptic separation is measured in degrees and is always in the range [0, 180].

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose visibility is to be calculated. |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |

**Returns:** Returns a valid [`ElongationInfo`](#ElongationInfo) structure, or throws an exception if there is an error.

<a name="Astronomy.Equator"></a>
### Astronomy.Equator(body, time, observer, equdate, aberration) &#8658; [`Equatorial`](#Equatorial)

**Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.**

Calculates topocentric equatorial coordinates in one of two different systems:
J2000 or true-equator-of-date, depending on the value of the `equdate` parameter.
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
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the observation takes place. |
| [`Observer`](#Observer) | `observer` | A location on or near the surface of the Earth. |
| [`EquatorEpoch`](#EquatorEpoch) | `equdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. |
| [`Aberration`](#Aberration) | `aberration` | Selects whether or not to correct for aberration. |

**Returns:** Topocentric equatorial coordinates of the celestial body.

<a name="Astronomy.EquatorFromVector"></a>
### Astronomy.EquatorFromVector(vector) &#8658; [`Equatorial`](#Equatorial)

**Given an equatorial vector, calculates equatorial angular coordinates.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `vector` | A vector in an equatorial coordinate system. |

**Returns:** Angular coordinates expressed in the same equatorial system as `vector`.

<a name="Astronomy.EquatorialToEcliptic"></a>
### Astronomy.EquatorialToEcliptic(eqj) &#8658; [`Ecliptic`](#Ecliptic)

**Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.**

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
on 1 January 2000), this function converts those coordinates to true ecliptic coordinates of date,
which are relative to the plane of the Earth's orbit around the Sun.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `eqj` | Equatorial coordinates in the J2000 frame of reference. You can call [`Astronomy.GeoVector`](#Astronomy.GeoVector) to obtain suitable equatorial coordinates. |

**Returns:** Spherical and vector coordinates expressed in true ecliptic coordinates of date (ECT).

<a name="Astronomy.GeoEmbState"></a>
### Astronomy.GeoEmbState(time) &#8658; [`StateVector`](#StateVector)

**Calculates the geocentric position and velocity of the Earth/Moon barycenter.**

Given a time of observation, calculates the geocentric position and velocity vectors
of the Earth/Moon barycenter (EMB).
The position (x, y, z) components are expressed in AU (astronomical units).
The velocity (vx, vy, vz) components are expressed in AU/day.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the EMB vectors. |

**Returns:** The EMB's position and velocity vectors in geocentric J2000 equatorial coordinates.

<a name="Astronomy.GeoMoon"></a>
### Astronomy.GeoMoon(time) &#8658; [`AstroVector`](#AstroVector)

**Calculates equatorial geocentric position of the Moon at a given time.**

Given a time of observation, calculates the Moon's position vector.
The vector indicates the Moon's center relative to the Earth's center.
The vector components are expressed in AU (astronomical units).
The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
In Astronomy Engine, this orientation is called EQJ.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the Moon's position. |

**Returns:** The Moon's position vector in J2000 equatorial coordinates (EQJ).

<a name="Astronomy.GeoMoonState"></a>
### Astronomy.GeoMoonState(time) &#8658; [`StateVector`](#StateVector)

**Calculates equatorial geocentric position and velocity of the Moon at a given time.**

Given a time of observation, calculates the Moon's position and velocity vectors.
The position and velocity are of the Moon's center relative to the Earth's center.
The position (x, y, z) components are expressed in AU (astronomical units).
The velocity (vx, vy, vz) components are expressed in AU/day.
The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
In Astronomy Engine, this orientation is called EQJ.
If you need the Moon's position only, and not its velocity,
it is much more efficient to use [`Astronomy.GeoMoon`](#Astronomy.GeoMoon) instead.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the Moon's position and velocity. |

**Returns:** The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).

<a name="Astronomy.GeoVector"></a>
### Astronomy.GeoVector(body, time, aberration) &#8658; [`AstroVector`](#AstroVector)

**Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Earth as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.

If given an invalid value for `body`, this function will throw an exception.

Unlike [`Astronomy.HelioVector`](#Astronomy.HelioVector), this function corrects for light travel time.
This means the position of the body is "back-dated" by the amount of time it takes
light to travel from that body to an observer on the Earth.

Also, the position can optionally be corrected for
[aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
causing the apparent direction of the body to be shifted due to transverse
movement of the Earth with respect to the rays of light coming from that body.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets. |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the position. |
| [`Aberration`](#Aberration) | `aberration` | `Aberration.Corrected` to correct for aberration, or `Aberration.None` to leave uncorrected. |

**Returns:** A geocentric position vector of the center of the given body.

<a name="Astronomy.GlobalSolarEclipsesAfter"></a>
### Astronomy.GlobalSolarEclipsesAfter(startTime) &#8658; `IEnumerable<`[`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo)`>`

**Enumerates a series of global solar eclipses that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchGlobalSolarEclipse`](#Astronomy.SearchGlobalSolarEclipse) and [`Astronomy.NextGlobalSolarEclipse`](#Astronomy.NextGlobalSolarEclipse).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive solar eclipses. |

<a name="Astronomy.HelioDistance"></a>
### Astronomy.HelioDistance(body, time) &#8658; `double`

**Calculates the distance between a body and the Sun at a given time.**

Given a date and time, this function calculates the distance between
the center of `body` and the center of the Sun, expressed in AU.
For the planets Mercury through Neptune, this function is significantly
more efficient than calling [`Astronomy.HelioVector`](#Astronomy.HelioVector) followed by taking the length
of the resulting vector.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric distance: the Sun, Moon, EMB, SSB, any of the planets, or a user-defined star. |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the heliocentric distance. |

**Returns:** The heliocentric distance in AU.

<a name="Astronomy.HelioState"></a>
### Astronomy.HelioState(body, time) &#8658; [`StateVector`](#StateVector)

**Calculates heliocentric position and velocity vectors for the given body.**

Given a body and a time, calculates the position and velocity
vectors for the center of that body at that time, relative to the center of the Sun.
The vectors are expressed in J2000 mean equator coordinates (EQJ).
If you need the position vector only, it is more efficient to call [`Astronomy.HelioVector`](#Astronomy.HelioVector).
The Sun's center is a non-inertial frame of reference. In other words, the Sun
experiences acceleration due to gravitational forces, mostly from the larger
planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
kinetic energy, or other quantities that require a non-accelerating frame
of reference, consider using [`Astronomy.BaryState`](#Astronomy.BaryState) instead.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The celestial body whose heliocentric state vector is to be calculated. Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets: `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. Also allowed to be a user-defined star created by [`Astronomy.DefineStar`](#Astronomy.DefineStar). |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate position and velocity. |

**Returns:** A structure that contains heliocentric position and velocity vectors.

<a name="Astronomy.HelioVector"></a>
### Astronomy.HelioVector(body, time) &#8658; [`AstroVector`](#AstroVector)

**Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.**

This function calculates the position of the given celestial body as a vector,
using the center of the Sun as the origin.  The result is expressed as a Cartesian
vector in the J2000 equatorial system: the coordinates are based on the mean equator
of the Earth at noon UTC on 1 January 2000.

The position is not corrected for light travel time or aberration.
This is different from the behavior of [`Astronomy.GeoVector`](#Astronomy.GeoVector).

If given an invalid value for `body`, this function will throw an [`InvalidBodyException`](#InvalidBodyException).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A body for which to calculate a heliocentric position: the Sun, Moon, EMB, SSB, or any of the planets. Also allowed to be a user-defined star created by [`Astronomy.DefineStar`](#Astronomy.DefineStar). |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the position. |

**Returns:** A heliocentric position vector of the center of the given body.

<a name="Astronomy.Horizon"></a>
### Astronomy.Horizon(time, observer, ra, dec, refraction) &#8658; [`Topocentric`](#Topocentric)

**Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.**

Given a date and time, the geographic location of an observer on the Earth, and
equatorial coordinates (right ascension and declination) of a celestial body,
this function returns horizontal coordinates (azimuth and altitude angles) for the body
relative to the horizon at the geographic location.

The right ascension `ra` and declination `dec` passed in must be *equator of date*
coordinates, based on the Earth's true equator at the date and time of the observation.
Otherwise the resulting horizontal coordinates will be inaccurate.
Equator of date coordinates can be obtained by calling [`Astronomy.Equator`](#Astronomy.Equator), passing in
`EquatorEpoch.OfDate` as its `equdate` parameter. It is also recommended to enable
aberration correction by passing in `Aberration.Corrected` as the `aberration` parameter.

This function optionally corrects for atmospheric refraction.
For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
correct for optical lensing of the Earth's atmosphere that causes objects
to appear somewhat higher above the horizon than they actually are.
However, callers may choose to avoid this correction by passing in `Refraction.None`.
If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
in the [`Topocentric`](#Topocentric) structure returned by this function will all be corrected for refraction.
If refraction is disabled, none of these four coordinates will be corrected; in that case,
the right ascension and declination in the returned structure will be numerically identical
to the respective `ra` and `dec` values passed in.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |
| `double` | `ra` | The right ascension of the body in sidereal hours. See remarks above for more details. |
| `double` | `dec` | The declination of the body in degrees. See remarks above for more details. |
| [`Refraction`](#Refraction) | `refraction` | Selects whether to correct for atmospheric refraction, and if so, which model to use. The recommended value for most uses is `Refraction.Normal`. See remarks above for more details. |

**Returns:** The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.

<a name="Astronomy.HorizonFromVector"></a>
### Astronomy.HorizonFromVector(vector, refraction) &#8658; [`Spherical`](#Spherical)

**Converts Cartesian coordinates to horizontal coordinates.**

Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

*IMPORTANT:* This function differs from [`Astronomy.SphereFromVector`](#Astronomy.SphereFromVector) in two ways:
- `Astronomy.SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
  from north (e.g., west = +90), but this function represents a clockwise rotation
  (e.g., east = +90). The difference is because `Astronomy.SphereFromVector` is intended
  to preserve the vector "right-hand rule", while this function defines azimuth in a more
  traditional way as used in navigation and cartography.
- This function optionally corrects for atmospheric refraction, while `Astronomy.SphereFromVector`
  does not.

The returned structure contains the azimuth in `lon`.
It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

The altitude is stored in `lat`.

The distance to the observed object is stored in `dist`,
and is expressed in astronomical units (AU).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `vector` | Cartesian vector to be converted to horizontal coordinates. |
| [`Refraction`](#Refraction) | `refraction` | `Refraction.Normal`: correct altitude for atmospheric refraction (recommended). `Refraction.None`: no atmospheric refraction correction is performed. `Refraction.JplHor`: for JPL Horizons compatibility testing only; not recommended for normal use. |

**Returns:** Horizontal spherical coordinates as described above.

<a name="Astronomy.HourAngle"></a>
### Astronomy.HourAngle(body, time, observer) &#8658; `double`

**Finds the hour angle of a body for a given observer and time.**

The *hour angle* of a celestial body indicates its position in the sky with respect
to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
the number of hours that have passed since the most recent time that the body has culminated,
or reached its highest point.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The body whose observed hour angle is to be found. |
| [`AstroTime`](#AstroTime) | `time` | The time of the observation. |
| [`Observer`](#Observer) | `observer` | The geographic location where the observation takes place. |

**Returns:** The real-valued hour angle of the body in the half-open range [0, 24).

<a name="Astronomy.IdentityMatrix"></a>
### Astronomy.IdentityMatrix() &#8658; [`RotationMatrix`](#RotationMatrix)

**Creates an identity rotation matrix.**

Returns a rotation matrix that has no effect on orientation.
This matrix can be the starting point for other operations,
such as using a series of calls to [`Astronomy.Pivot`](#Astronomy.Pivot) to
create a custom rotation matrix.

**Returns:** The identity matrix.

<a name="Astronomy.Illumination"></a>
### Astronomy.Illumination(body, time) &#8658; [`IllumInfo`](#IllumInfo)

**Finds visual magnitude, phase angle, and other illumination information about a celestial body.**

This function calculates information about how bright a celestial body appears from the Earth,
reported as visual magnitude, which is a smaller (or even negative) number for brighter objects
and a larger number for dimmer objects.

For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
of the body appears illuminated as seen from the Earth. For example, when the phase angle is
near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
of 90 degrees means the body appears "half full".
For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
so it doesn't have a phase angle.

When the body is Saturn, the returned structure contains a field `ring_tilt` that holds
the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
0 for all bodies other than Saturn.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, or any planet other than the Earth. |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |

**Returns:** An [`IllumInfo`](#IllumInfo) structure with fields as documented above.

<a name="Astronomy.InverseRefractionAngle"></a>
### Astronomy.InverseRefractionAngle(refraction, bent_altitude) &#8658; `double`

**Calculates the inverse of an atmospheric refraction angle.**

Given an observed altitude angle that includes atmospheric refraction,
calculates the negative angular correction to obtain the unrefracted
altitude. This is useful for cases where observed horizontal
coordinates are to be converted to another orientation system,
but refraction first must be removed from the observed position.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | The option selecting which refraction correction to use. See [`Astronomy.RefractionAngle`](#Astronomy.RefractionAngle). |
| `double` | `bent_altitude` | The apparent altitude that includes atmospheric refraction. |

**Returns:** The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing. This will be less than or equal to zero.

<a name="Astronomy.InverseRotation"></a>
### Astronomy.InverseRotation(rotation) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates the inverse of a rotation matrix.**

Given a rotation matrix that performs some coordinate transform,
this function returns the matrix that reverses that transform.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | The rotation matrix to be inverted. |

**Returns:** A rotation matrix that performs the opposite transformation.

<a name="Astronomy.JupiterMoons"></a>
### Astronomy.JupiterMoons(time) &#8658; [`JupiterMoonsInfo`](#JupiterMoonsInfo)

**Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.**

Calculates position and velocity vectors for Jupiter's moons
Io, Europa, Ganymede, and Callisto, at the given date and time.
The vectors are jovicentric (relative to the center of Jupiter).
Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
The position components are expressed in astronomical units (AU), and the
velocity components are in AU/day.

To convert to heliocentric position vectors, call [`Astronomy.HelioVector`](#Astronomy.HelioVector)
with `Body.Jupiter` to get Jupiter's heliocentric position, then
add the jovicentric positions. Likewise, you can call [`Astronomy.GeoVector`](#Astronomy.GeoVector)
to convert to geocentric positions; however, you will have to manually
correct for light travel time from the Jupiter system to Earth to
figure out what time to pass to `jupiterMoons` to get an accurate picture
of how Jupiter and its moons look from Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the position vectors. |

**Returns:** Position and velocity vectors of Jupiter's largest 4 moons.

<a name="Astronomy.LagrangePoint"></a>
### Astronomy.LagrangePoint(point, time, major_body, minor_body) &#8658; [`StateVector`](#StateVector)

**Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.**

Given a more massive "major" body and a much less massive "minor" body,
calculates one of the five Lagrange points in relation to the minor body's
orbit around the major body. The parameter `point` is an integer that
selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body.
2 = the Lagrange point on the far side of the minor body.
3 = the Lagrange point on the far side of the major body.
4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The function returns the state vector for the selected Lagrange point
in J2000 mean equator coordinates (EQJ), with respect to the center of the
major body.

To calculate Sun/Earth Lagrange points, pass in `Body.Sun` for `major_body`
and `Body.EMB` (Earth/Moon barycenter) for `minor_body`.
For Lagrange points of the Sun and any other planet, pass in just that planet
(e.g. `Body.Jupiter`) for `minor_body`.
To calculate Earth/Moon Lagrange points, pass in `Body.Earth` and `Body.Moon`
for the major and minor bodies respectively.

In some cases, it may be more efficient to call [`Astronomy.LagrangePointFast`](#Astronomy.LagrangePointFast),
especially when the state vectors have already been calculated, or are needed
for some other purpose.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` | An integer 1..5 that selects which of the Lagrange points to calculate. |
| [`AstroTime`](#AstroTime) | `time` | The time for which the Lagrange point is to be calculated. |
| [`Body`](#Body) | `major_body` | The more massive of the co-orbiting bodies: `Body.Sun` or `Body.Earth`. |
| [`Body`](#Body) | `minor_body` | The less massive of the co-orbiting bodies. See main remarks. |

**Returns:** The position and velocity of the selected Lagrange point with respect to the major body's center.

<a name="Astronomy.LagrangePointFast"></a>
### Astronomy.LagrangePointFast(point, major_state, major_mass, minor_state, minor_mass) &#8658; [`StateVector`](#StateVector)

**Calculates one of the 5 Lagrange points from body masses and state vectors.**

Given a more massive "major" body and a much less massive "minor" body,
calculates one of the five Lagrange points in relation to the minor body's
orbit around the major body. The parameter `point` is an integer that
selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body.
2 = the Lagrange point on the far side of the minor body.
3 = the Lagrange point on the far side of the major body.
4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The caller passes in the state vector and mass for both bodies.
The state vectors can be in any orientation and frame of reference.
The body masses are expressed as GM products, where G = the universal
gravitation constant and M = the body's mass. Thus the units for
`major_mass` and `minor_mass` must be au^3/day^2.
Use [`Astronomy.MassProduct`](#Astronomy.MassProduct) to obtain GM values for various solar system bodies.

The function returns the state vector for the selected Lagrange point
using the same orientation as the state vector parameters `major_state` and `minor_state`,
and the position and velocity components are with respect to the major body's center.

Consider calling [`Astronomy.LagrangePoint`](#Astronomy.LagrangePoint), instead of this function, for simpler usage in most cases.

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `point` | An integer 1..5 that selects which of the Lagrange points to calculate. |
| [`StateVector`](#StateVector) | `major_state` | The state vector of the major (more massive) of the pair of bodies. |
| `double` | `major_mass` | The mass product GM of the major body. |
| [`StateVector`](#StateVector) | `minor_state` | The state vector of the minor (less massive) of the pair of bodies. |
| `double` | `minor_mass` | The mass product GM of the minor body. |

**Returns:** The position and velocity of the selected Lagrange point with respect to the major body's center.

<a name="Astronomy.Libration"></a>
### Astronomy.Libration(time) &#8658; [`LibrationInfo`](#LibrationInfo)

**Calculates the Moon's libration angles at a given moment in time.**

Libration is an observed back-and-forth wobble of the portion of the
Moon visible from the Earth. It is caused by the imperfect tidal locking
of the Moon's fixed rotation rate, compared to its variable angular speed
of orbit around the Earth.

This function calculates a pair of perpendicular libration angles,
one representing rotation of the Moon in ecliptic longitude `elon`, the other
in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.

This function also returns the geocentric position of the Moon
expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
and the apparent angular diameter of the Moon `diam_deg`.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate lunar libration. |

**Returns:** The Moon's ecliptic position and libration angles as seen from the Earth.

<a name="Astronomy.LocalSolarEclipsesAfter"></a>
### Astronomy.LocalSolarEclipsesAfter(startTime, observer) &#8658; `IEnumerable<`[`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo)`>`

**Enumerates a series of local solar eclipses that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchLocalSolarEclipse`](#Astronomy.SearchLocalSolarEclipse) and [`Astronomy.NextLocalSolarEclipse`](#Astronomy.NextLocalSolarEclipse).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive solar eclipses. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |

<a name="Astronomy.LunarApsidesAfter"></a>
### Astronomy.LunarApsidesAfter(startTime) &#8658; `IEnumerable<`[`ApsisInfo`](#ApsisInfo)`>`

**Enumerates a series of apogees/perigees that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) and [`Astronomy.NextLunarApsis`](#Astronomy.NextLunarApsis).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive lunar apsides. |

<a name="Astronomy.LunarEclipsesAfter"></a>
### Astronomy.LunarEclipsesAfter(startTime) &#8658; `IEnumerable<`[`LunarEclipseInfo`](#LunarEclipseInfo)`>`

**Enumerates a series of lunar eclipses that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchLunarEclipse`](#Astronomy.SearchLunarEclipse) and [`Astronomy.NextLunarEclipse`](#Astronomy.NextLunarEclipse).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive lunar eclipses. |

<a name="Astronomy.MassProduct"></a>
### Astronomy.MassProduct(body) &#8658; `double`

**Returns the product of mass and universal gravitational constant of a Solar System body.**

For problems involving the gravitational interactions of Solar System bodies,
it is helpful to know the product GM, where G = the universal gravitational constant
and M = the mass of the body. In practice, GM is known to a higher precision than
either G or M alone, and thus using the product results in the most accurate results.
This function returns the product GM in the units au^3/day^2.
The values come from page 10 of a
[JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The body for which to find the GM product. Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet. Any other value will cause an exception to be thrown. |

**Returns:** The mass product of the given body in au^3/day^2.

<a name="Astronomy.MoonNodesAfter"></a>
### Astronomy.MoonNodesAfter(startTime) &#8658; `IEnumerable<`[`NodeEventInfo`](#NodeEventInfo)`>`

**Enumerates a series of ascending/descending nodes of the Moon that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchMoonNode`](#Astronomy.SearchMoonNode) and [`Astronomy.NextMoonNode`](#Astronomy.NextMoonNode).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive lunar apsides. |

<a name="Astronomy.MoonPhase"></a>
### Astronomy.MoonPhase(time) &#8658; `double`

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
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |

**Returns:** The angle as described above, a value in the range 0..360 degrees.

<a name="Astronomy.MoonQuartersAfter"></a>
### Astronomy.MoonQuartersAfter(startTime) &#8658; `IEnumerable<`[`MoonQuarterInfo`](#MoonQuarterInfo)`>`

**Enumerates a series of lunar quarter phases that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchMoonQuarter`](#Astronomy.SearchMoonQuarter) and [`Astronomy.NextMoonQuarter`](#Astronomy.NextMoonQuarter).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive lunar quarter phases. |

<a name="Astronomy.NextGlobalSolarEclipse"></a>
### Astronomy.NextGlobalSolarEclipse(prevEclipseTime) &#8658; [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo)

**Searches for the next global solar eclipse in a series.**

After using [`Astronomy.SearchGlobalSolarEclipse`](#Astronomy.SearchGlobalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo) returned by the
previous call to `Astronomy.SearchGlobalSolarEclipse` or `Astronomy.NextGlobalSolarEclipse`
to find the next solar eclipse.

See [`Astronomy.GlobalSolarEclipsesAfter`](#Astronomy.GlobalSolarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `prevEclipseTime` | A date and time near a new moon. Solar eclipse search will start at the next new moon. |

<a name="Astronomy.NextLocalSolarEclipse"></a>
### Astronomy.NextLocalSolarEclipse(prevEclipseTime, observer) &#8658; [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo)

**Searches for the next local solar eclipse in a series.**

After using [`Astronomy.SearchLocalSolarEclipse`](#Astronomy.SearchLocalSolarEclipse) to find the first solar eclipse
in a series, you can call this function to find the next consecutive solar eclipse.
Pass in the `peak` value from the [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) returned by the
previous call to `Astronomy.SearchLocalSolarEclipse` or `Astronomy.NextLocalSolarEclipse`
to find the next solar eclipse.

See [`Astronomy.LocalSolarEclipsesAfter`](#Astronomy.LocalSolarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `prevEclipseTime` | A date and time near a new moon. Solar eclipse search will start at the next new moon. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |

<a name="Astronomy.NextLunarApsis"></a>
### Astronomy.NextLunarApsis(apsis) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the next lunar perigee or apogee event in a series.**

This function requires an [`ApsisInfo`](#ApsisInfo) value obtained from a call
to [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) or `Astronomy.NextLunarApsis`. Given
an apogee event, this function finds the next perigee event, and vice versa.

See [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) for more details.
See [`Astronomy.LunarApsidesAfter`](#Astronomy.LunarApsidesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`ApsisInfo`](#ApsisInfo) | `apsis` | An apsis event obtained from a call to [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) or `Astronomy.NextLunarApsis`. See [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) for more details. |

**Returns:** Same as the return value for [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis).

<a name="Astronomy.NextLunarEclipse"></a>
### Astronomy.NextLunarEclipse(prevEclipseTime) &#8658; [`LunarEclipseInfo`](#LunarEclipseInfo)

**Searches for the next lunar eclipse in a series.**

After using [`Astronomy.SearchLunarEclipse`](#Astronomy.SearchLunarEclipse) to find the first lunar eclipse
in a series, you can call this function to find the next consecutive lunar eclipse.
Pass in the `center` value from the [`LunarEclipseInfo`](#LunarEclipseInfo) returned by the
previous call to `Astronomy.SearchLunarEclipse` or `Astronomy.NextLunarEclipse`
to find the next lunar eclipse.

See [`Astronomy.LunarEclipsesAfter`](#Astronomy.LunarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `prevEclipseTime` | A date and time near a full moon. Lunar eclipse search will start at the next full moon. |

**Returns:** A [`LunarEclipseInfo`](#LunarEclipseInfo) structure containing information about the lunar eclipse.

<a name="Astronomy.NextMoonNode"></a>
### Astronomy.NextMoonNode(prevNode) &#8658; [`NodeEventInfo`](#NodeEventInfo)

**Searches for the next time when the Moon's center crosses through the ecliptic plane.**

Call [`Astronomy.SearchMoonNode`](#Astronomy.SearchMoonNode) to find the first of a series of nodes.
Then call `Astronomy.NextMoonNode` to find as many more consecutive nodes as desired.

See [`Astronomy.MoonNodesAfter`](#Astronomy.MoonNodesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`NodeEventInfo`](#NodeEventInfo) | `prevNode` | The previous node found from calling [`Astronomy.SearchMoonNode`](#Astronomy.SearchMoonNode) or `Astronomy.NextMoonNode`. |

<a name="Astronomy.NextMoonQuarter"></a>
### Astronomy.NextMoonQuarter(mq) &#8658; [`MoonQuarterInfo`](#MoonQuarterInfo)

**Continues searching for lunar quarters from a previous search.**

After calling [`Astronomy.SearchMoonQuarter`](#Astronomy.SearchMoonQuarter), this function can be called
one or more times to continue finding consecutive lunar quarters.
This function finds the next consecutive moon quarter event after
the one passed in as the parameter `mq`.

See [`Astronomy.MoonQuartersAfter`](#Astronomy.MoonQuartersAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`MoonQuarterInfo`](#MoonQuarterInfo) | `mq` | The previous moon quarter found by a call to [`Astronomy.SearchMoonQuarter`](#Astronomy.SearchMoonQuarter) or `Astronomy.NextMoonQuarter`. |

**Returns:** The moon quarter that occurs next in time after the one passed in `mq`.

<a name="Astronomy.NextPlanetApsis"></a>
### Astronomy.NextPlanetApsis(body, apsis) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the next planetary perihelion or aphelion event in a series.**

This function requires an [`ApsisInfo`](#ApsisInfo) value obtained from a call
to [`Astronomy.SearchPlanetApsis`](#Astronomy.SearchPlanetApsis) or `Astronomy.NextPlanetApsis`.
Given an aphelion event, this function finds the next perihelion event, and vice versa.

See [`Astronomy.SearchPlanetApsis`](#Astronomy.SearchPlanetApsis) for more details.
See [`Astronomy.PlanetApsidesAfter`](#Astronomy.PlanetApsidesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. Must match the body passed into the call that produced the `apsis` parameter. |
| [`ApsisInfo`](#ApsisInfo) | `apsis` | An apsis event obtained from a call to [`Astronomy.SearchPlanetApsis`](#Astronomy.SearchPlanetApsis) or `Astronomy.NextPlanetApsis`. |

**Returns:** Same as the return value for [`Astronomy.SearchPlanetApsis`](#Astronomy.SearchPlanetApsis).

<a name="Astronomy.NextTransit"></a>
### Astronomy.NextTransit(body, prevTransitTime) &#8658; [`TransitInfo`](#TransitInfo)

**Searches for another transit of Mercury or Venus.**

After calling [`Astronomy.SearchTransit`](#Astronomy.SearchTransit) to find a transit of Mercury or Venus,
this function finds the next transit after that.
Keep calling this function as many times as you want to keep finding more transits.

See [`Astronomy.TransitsAfter`](#Astronomy.TransitsAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`. |
| [`AstroTime`](#AstroTime) | `prevTransitTime` | A date and time near the previous transit. |

<a name="Astronomy.ObserverGravity"></a>
### Astronomy.ObserverGravity(latitude, height) &#8658; `double`

**Calculates the gravitational acceleration experienced by an observer on the Earth.**

This function implements the WGS 84 Ellipsoidal Gravity Formula.
The result is a combination of inward gravitational acceleration
with outward centrifugal acceleration, as experienced by an observer
in the Earth's rotating frame of reference.
The resulting value increases toward the Earth's poles and decreases
toward the equator, consistent with changes of the weight measured
by a spring scale of a fixed mass moved to different latitudes and heights
on the Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `latitude` | The latitude of the observer in degrees north or south of the equator. By formula symmetry, positive latitudes give the same answer as negative latitudes, so the sign does not matter. |
| `double` | `height` | The height above the sea level geoid in meters. No range checking is done; however, accuracy is only valid in the range 0 to 100000 meters. |

**Returns:** The effective gravitational acceleration expressed in meters per second squared [m/s^2].

<a name="Astronomy.ObserverState"></a>
### Astronomy.ObserverState(time, observer, equdate) &#8658; [`StateVector`](#StateVector)

**Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.**

This function calculates position and velocity vectors of an observer
on or near the surface of the Earth, expressed in equatorial
coordinates. It takes into account the rotation of the Earth at the given
time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in `equdate` to select either `EquatorEpoch.J2000`
for using J2000 coordinates, or `EquatorEpoch.OfDate` for using coordinates relative
to the Earth's equator at the specified time.

The returned position vector has components expressed in astronomical units (AU).
To convert to kilometers, multiply the `x`, `y`, and `z` values by
the constant value [`Astronomy.KM_PER_AU`](#Astronomy.KM_PER_AU).

The returned velocity vector is measured in AU/day.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the observer's geocentric state vector. |
| [`Observer`](#Observer) | `observer` | The geographic location of a point on or near the surface of the Earth. |
| [`EquatorEpoch`](#EquatorEpoch) | `equdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `time` as the orientation. |

**Returns:** The position and velocity of the given geographic location, relative to the center of the Earth.

<a name="Astronomy.ObserverVector"></a>
### Astronomy.ObserverVector(time, observer, equdate) &#8658; [`AstroVector`](#AstroVector)

**Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.**

This function calculates a vector from the center of the Earth to
a point on or near the surface of the Earth, expressed in equatorial
coordinates. It takes into account the rotation of the Earth at the given
time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in `equdate` to select either `EquatorEpoch.J2000`
for using J2000 coordinates, or `EquatorEpoch.OfDate` for using coordinates relative
to the Earth's equator at the specified time.

The returned vector has components expressed in astronomical units (AU).
To convert to kilometers, multiply the `x`, `y`, and `z` values by
the constant value [`Astronomy.KM_PER_AU`](#Astronomy.KM_PER_AU).

The inverse of this function is also available: [`Astronomy.VectorObserver`](#Astronomy.VectorObserver).

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the observer's position vector. |
| [`Observer`](#Observer) | `observer` | The geographic location of a point on or near the surface of the Earth. |
| [`EquatorEpoch`](#EquatorEpoch) | `equdate` | Selects the date of the Earth's equator in which to express the equatorial coordinates. The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by the `time` parameter. Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `time` as the orientation. |

**Returns:** An equatorial vector from the center of the Earth to the specified location on (or near) the Earth's surface.

<a name="Astronomy.PairLongitude"></a>
### Astronomy.PairLongitude(body1, body2, time) &#8658; `double`

**Returns one body's ecliptic longitude with respect to another, as seen from the Earth.**

This function determines where one body appears around the ecliptic plane
(the plane of the Earth's orbit around the Sun) as seen from the Earth,
relative to the another body's apparent position.
The function returns an angle in the half-open range [0, 360) degrees.
The value is the ecliptic longitude of `body1` relative to the ecliptic
longitude of `body2`.

The angle is 0 when the two bodies are at the same ecliptic longitude
as seen from the Earth. The angle increases in the prograde direction
(the direction that the planets orbit the Sun and the Moon orbits the Earth).

When the angle is 180 degrees, it means the two bodies appear on opposite sides
of the sky for an Earthly observer.

Neither `body1` nor `body2` is allowed to be `Body.Earth`.
If this happens, the function throws an exception.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body1` | The first body, whose longitude is to be found relative to the second body. |
| [`Body`](#Body) | `body2` | The second body, relative to which the longitude of the first body is to be found. |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |

**Returns:** An angle in the range [0, 360), expressed in degrees.

<a name="Astronomy.Pivot"></a>
### Astronomy.Pivot(rotation, axis, angle) &#8658; [`RotationMatrix`](#RotationMatrix)

**Re-orients a rotation matrix by pivoting it by an angle around one of its axes.**

Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
this function pivots the rotation matrix by that angle around that coordinate axis.

For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
of a telescope camera pointed at a given body, you can use `Astronomy.Pivot` twice:
(1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
western axis by the body's altitude angle. The resulting rotation matrix will then
reorient ECL coordinates to the orientation of your telescope camera.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | The input rotation matrix. |
| `int` | `axis` | An integer that selects which coordinate axis to rotate around: 0 = x, 1 = y, 2 = z. Any other value will cause an ArgumentException to be thrown. |
| `double` | `angle` | An angle in degrees indicating the amount of rotation around the specified axis. Positive angles indicate rotation counterclockwise as seen from the positive direction along that axis, looking towards the origin point of the orientation system. Any finite number of degrees is allowed, but best precision will result from keeping `angle` in the range [-360, +360]. |

**Returns:** A pivoted matrix object.

<a name="Astronomy.PlanetApsidesAfter"></a>
### Astronomy.PlanetApsidesAfter(body, startTime) &#8658; `IEnumerable<`[`ApsisInfo`](#ApsisInfo)`>`

**Enumerates a series of planet aphelia/perihelia that occur after a specified time.**

This is a convenience wrapper around the functions
[`Astronomy.SearchPlanetApsis`](#Astronomy.SearchPlanetApsis) and [`Astronomy.NextPlanetApsis`](#Astronomy.NextPlanetApsis).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find a series of consecutive aphelia/perihelia. Not allowed to be `Body.Sun` or `Body.Moon`. |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive planetary apsides. |

<a name="Astronomy.PlanetOrbitalPeriod"></a>
### Astronomy.PlanetOrbitalPeriod(body) &#8658; `double`

**Returns the average number of days it takes for a planet to orbit the Sun.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, or Pluto. |

**Returns:** The mean orbital period of the body in days.

<a name="Astronomy.RefractionAngle"></a>
### Astronomy.RefractionAngle(refraction, altitude) &#8658; `double`

**Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.**

Given an altitude angle and a refraction option, calculates
the amount of "lift" caused by atmospheric refraction.
This is the number of degrees higher in the sky an object appears
due to the lensing of the Earth's atmosphere.
This function works best near sea level.
To correct for higher elevations, call [`Astronomy.Atmosphere`](#Astronomy.Atmosphere) for that
elevation and multiply the refraction angle by the resulting relative density.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Refraction`](#Refraction) | `refraction` | The option selecting which refraction correction to use. If `Refraction.Normal`, uses a well-behaved refraction model that works well for all valid values (-90 to +90) of `altitude`. If `Refraction.JplHor`, this function returns a compatible value with the JPL Horizons tool. If any other value (including `Refraction.None`), this function returns 0. |
| `double` | `altitude` | An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90. |

**Returns:** The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.

<a name="Astronomy.RotateState"></a>
### Astronomy.RotateState(rotation, state) &#8658; [`StateVector`](#StateVector)

**Applies a rotation to a state vector, yielding a rotated state vector.**

This function transforms a state vector in one orientation to a state vector in another orientation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | A rotation matrix that specifies how the orientation of the state vector is to be changed. |
| [`StateVector`](#StateVector) | `state` | The state vector whose orientation is to be changed. |

**Returns:** A state vector in the orientation specified by `rotation`.

<a name="Astronomy.RotateVector"></a>
### Astronomy.RotateVector(rotation, vector) &#8658; [`AstroVector`](#AstroVector)

**Applies a rotation to a vector, yielding a rotated vector.**

This function transforms a vector in one orientation to a vector
in another orientation.

| Type | Parameter | Description |
| --- | --- | --- |
| [`RotationMatrix`](#RotationMatrix) | `rotation` | A rotation matrix that specifies how the orientation of the vector is to be changed. |
| [`AstroVector`](#AstroVector) | `vector` | The vector whose orientation is to be changed. |

**Returns:** A vector in the orientation specified by `rotation`.

<a name="Astronomy.Rotation_ECL_EQD"></a>
### Astronomy.Rotation_ECL_EQD(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of date.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the desired equator. |

**Returns:** A rotation matrix that converts ECL to EQD.

<a name="Astronomy.Rotation_ECL_EQJ"></a>
### Astronomy.Rotation_ECL_EQJ() &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: EQJ = equatorial system, using equator at J2000 epoch.

**Returns:** A rotation matrix that converts ECL to EQJ.

<a name="Astronomy.Rotation_ECL_HOR"></a>
### Astronomy.Rotation_ECL_HOR(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECL = ecliptic system, using equator at J2000 epoch.
Target: HOR = horizontal system.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns:** A rotation matrix that converts ECL to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0.

<a name="Astronomy.Rotation_ECT_EQD"></a>
### Astronomy.Rotation_ECT_EQD(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from true ecliptic of date (ECT) to equator of date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECT = true ecliptic of date.
Target: EQD = equator of date.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the ecliptic/equator conversion. |

**Returns:** A rotation matrix that converts ECT to EQD.

<a name="Astronomy.Rotation_ECT_EQJ"></a>
### Astronomy.Rotation_ECT_EQJ(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: ECT = ecliptic system, using true equinox of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator defines the target orientation. |

**Returns:** A rotation matrix that converts ECT to EQJ at `time`.

<a name="Astronomy.Rotation_EQD_ECL"></a>
### Astronomy.Rotation_EQD_ECL(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of date.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the source equator. |

**Returns:** A rotation matrix that converts EQD to ECL.

<a name="Astronomy.Rotation_EQD_ECT"></a>
### Astronomy.Rotation_EQD_ECT(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from equator of date (EQD) to true ecliptic of date (ECT) .**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equator of date.
Target: ECT = true ecliptic of date.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the equator/ecliptic conversion. |

**Returns:** A rotation matrix that converts EQD to ECT.

<a name="Astronomy.Rotation_EQD_EQJ"></a>
### Astronomy.Rotation_EQD_EQJ(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean equator (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: EQJ = equatorial system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator defines the source orientation. |

**Returns:** A rotation matrix that converts EQD at `time` to EQJ.

<a name="Astronomy.Rotation_EQD_HOR"></a>
### Astronomy.Rotation_EQD_HOR(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQD = equatorial system, using equator of the specified date/time.
Target: HOR = horizontal system.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns:** A rotation matrix that converts EQD to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0.

<a name="Astronomy.Rotation_EQJ_ECL"></a>
### Astronomy.Rotation_EQJ_ECL() &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECL = ecliptic system, using equator at J2000 epoch.

**Returns:** A rotation matrix that converts EQJ to ECL.

<a name="Astronomy.Rotation_EQJ_ECT"></a>
### Astronomy.Rotation_EQJ_ECT(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: ECT = ecliptic system, using true equinox of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator defines the target orientation. |

**Returns:** A rotation matrix that converts EQJ to ECT at `time`.

<a name="Astronomy.Rotation_EQJ_EQD"></a>
### Astronomy.Rotation_EQJ_EQD(time) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using equator at J2000 epoch.
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator defines the target orientation. |

**Returns:** A rotation matrix that converts EQJ to EQD at `time`.

<a name="Astronomy.Rotation_EQJ_GAL"></a>
### Astronomy.Rotation_EQJ_GAL() &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to galactic (GAL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: GAL = galactic system (IAU 1958 definition).

**Returns:** A rotation matrix that converts EQJ to GAL.

<a name="Astronomy.Rotation_EQJ_HOR"></a>
### Astronomy.Rotation_EQJ_HOR(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: EQJ = equatorial system, using the equator at the J2000 epoch.
Target: HOR = horizontal system.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the desired horizontal orientation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns:** A rotation matrix that converts EQJ to HOR at `time` and for `observer`. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0.

<a name="Astronomy.Rotation_GAL_EQJ"></a>
### Astronomy.Rotation_GAL_EQJ() &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from galactic (GAL) to J2000 mean equator (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: GAL = galactic system (IAU 1958 definition).
Target: EQJ = equatorial system, using the equator at the J2000 epoch.

**Returns:** A rotation matrix that converts GAL to EQJ.

<a name="Astronomy.Rotation_HOR_ECL"></a>
### Astronomy.Rotation_HOR_ECL(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system.
Target: ECL = ecliptic system, using equator at J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the horizontal observation. |
| [`Observer`](#Observer) | `observer` | The location of the horizontal observer. |

**Returns:** A rotation matrix that converts HOR to ECL.

<a name="Astronomy.Rotation_HOR_EQD"></a>
### Astronomy.Rotation_HOR_EQD(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQD = equatorial system, using equator of the specified date/time.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time at which the Earth's equator applies. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns:** A rotation matrix that converts HOR to EQD at `time` and for `observer`.

<a name="Astronomy.Rotation_HOR_EQJ"></a>
### Astronomy.Rotation_HOR_EQJ(time, observer) &#8658; [`RotationMatrix`](#RotationMatrix)

**Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).**

This is one of the family of functions that returns a rotation matrix
for converting from one orientation to another.
Source: HOR = horizontal system (x=North, y=West, z=Zenith).
Target: EQJ = equatorial system, using equator at the J2000 epoch.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |
| [`Observer`](#Observer) | `observer` | A location near the Earth's mean sea level that defines the observer's horizon. |

**Returns:** A rotation matrix that converts HOR to EQJ at `time` and for `observer`.

<a name="Astronomy.RotationAxis"></a>
### Astronomy.RotationAxis(body, time) &#8658; [`AxisInfo`](#AxisInfo)

**Calculates information about a body's rotation axis at a given time.**

Calculates the orientation of a body's rotation axis, along with
the rotation angle of its prime meridian, at a given moment in time.

This function uses formulas standardized by the IAU Working Group
on Cartographics and Rotational Elements 2015 report, as described
in the following document:

https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf

See [`AxisInfo`](#AxisInfo) for more detailed information.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | One of the following values: `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| [`AstroTime`](#AstroTime) | `time` | The time at which to calculate the body's rotation axis. |

**Returns:** North pole orientation and body spin angle.

<a name="Astronomy.Search"></a>
### Astronomy.Search(func, t1, t2, dt_tolerance_seconds) &#8658; [`AstroTime`](#AstroTime)

**Searches for a time at which a function's value increases through zero.**

Certain astronomy calculations involve finding a time when an event occurs.
Often such events can be defined as the root of a function:
the time at which the function's value becomes zero.

`Search` finds the *ascending root* of a function: the time at which
the function's value becomes zero while having a positive slope. That is, as time increases,
the function transitions from a negative value, through zero at a specific moment,
to a positive value later. The goal of the search is to find that specific moment.

The `func` parameter is an instance of the abstract class [`SearchContext`](#SearchContext).
As an example, a caller may wish to find the moment a celestial body reaches a certain
ecliptic longitude. In that case, the caller might derive a class that contains
a [`Body`](#Body) member to specify the body and a `double` to hold the target longitude.
It could subtract the target longitude from the actual longitude at a given time;
thus the difference would equal zero at the moment in time the planet reaches the
desired longitude.

The search calls `func.Eval` repeatedly to rapidly narrow in on any ascending
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
the search will return `null`.

If the search does not converge within 20 iterations, it will throw an exception.

| Type | Parameter | Description |
| --- | --- | --- |
| [`SearchContext`](#SearchContext) | `func` | The function for which to find the time of an ascending root. See remarks above for more details. |
| [`AstroTime`](#AstroTime) | `t1` | The lower time bound of the search window. See remarks above for more details. |
| [`AstroTime`](#AstroTime) | `t2` | The upper time bound of the search window. See remarks above for more details. |
| `double` | `dt_tolerance_seconds` | Specifies an amount of time in seconds within which a bounded ascending root is considered accurate enough to stop. A typical value is 1 second. |

**Returns:** If successful, returns an [`AstroTime`](#AstroTime) value indicating a date and time that is within `dt_tolerance_seconds` of an ascending root. If no ascending root is found, or more than one root exists in the time window `t1`..`t2`, the function returns `null`. If the search does not converge within 20 iterations, an exception is thrown.

<a name="Astronomy.SearchAltitude"></a>
### Astronomy.SearchAltitude(body, observer, direction, startTime, limitDays, altitude) &#8658; [`AstroTime`](#AstroTime)

**Finds the next time the center of a body passes through a given altitude.**

Finds when the center of the given body ascends or descends through a given
altitude angle, as seen by an observer at the specified location on the Earth.
By using the appropriate combination of `direction` and `altitude` parameters,
this function can be used to find when civil, nautical, or astronomical twilight
begins (dawn) or ends (dusk).

Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.

Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.

Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.

Astronomical twilight uses -18 degrees as the `altitude` value.

By convention for twilight time calculations, the altitude is not corrected for
atmospheric refraction. This is because the target altitudes are below the horizon,
and refraction is not directly observable.

`SearchAltitude` is not intended to find rise/set times of a body for two reasons:
(1) Rise/set times of the Sun or Moon are defined by their topmost visible portion, not their centers.
(2) Rise/set times are affected significantly by atmospheric refraction.
Therefore, it is better to use [`Astronomy.SearchRiseSet`](#Astronomy.SearchRiseSet) to find rise/set times, which
corrects for both of these considerations.

`SearchAltitude` will not work reliably for altitudes at or near the body's
maximum or minimum altitudes. To find the time a body reaches minimum or maximum altitude
angles, use [`Astronomy.SearchHourAngle`](#Astronomy.SearchHourAngle).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy.DefineStar`](#Astronomy.DefineStar). |
| [`Observer`](#Observer) | `observer` | The location where observation takes place. |
| [`Direction`](#Direction) | `direction` | Either `Direction.Rise` to find an ascending altitude event or `Direction.Set` to find a descending altitude event. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start the search. |
| `double` | `limitDays` | Limits how many days to search for the body reaching the altitude angle, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit the search to the same day, you can use a value of 1 day. In cases where you want to find the altitude event no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |
| `double` | `altitude` | The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range [-90, +90]. |

**Returns:** The date and time of the altitude event, or `null` if no such event occurs within the specified time window.

<a name="Astronomy.SearchGlobalSolarEclipse"></a>
### Astronomy.SearchGlobalSolarEclipse(startTime) &#8658; [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo)

**Searches for a solar eclipse visible anywhere on the Earth's surface.**

This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [`GlobalSolarEclipseInfo`](#GlobalSolarEclipseInfo) for more information.
To find a series of solar eclipses, call this function once,
then keep calling [`Astronomy.NextGlobalSolarEclipse`](#Astronomy.NextGlobalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.

See [`Astronomy.GlobalSolarEclipsesAfter`](#Astronomy.GlobalSolarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for a solar eclipse. |

<a name="Astronomy.SearchHourAngle"></a>
### Astronomy.SearchHourAngle(body, observer, hourAngle, startTime, direction) &#8658; [`HourAngleInfo`](#HourAngleInfo)

**Searches for the time when the center of a body reaches a specified hour angle as seen by an observer on the Earth.**

The *hour angle* of a celestial body indicates its position in the sky with respect
to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
the number of hours that have passed since the most recent time that the body has culminated,
or reached its highest point.

This function searches for the next or previous time a celestial body reaches the given hour angle
relative to the date and time specified by `startTime`.
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
| [`Body`](#Body) | `body` | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy.DefineStar`](#Astronomy.DefineStar). |
| [`Observer`](#Observer) | `observer` | Indicates a location on or near the surface of the Earth where the observer is located. |
| `double` | `hourAngle` | An hour angle value in the range [0, 24) indicating the number of sidereal hours after the body's most recent culmination. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start the search. |
| `int` | `direction` | The direction in time to perform the search: a positive value searches forward in time, a negative value searches backward in time. The function throws an exception if `direction` is zero. |

**Returns:** This function returns a valid [`HourAngleInfo`](#HourAngleInfo) object on success. If any error occurs, it throws an exception. It never returns a null value.

<a name="Astronomy.SearchLocalSolarEclipse"></a>
### Astronomy.SearchLocalSolarEclipse(startTime, observer) &#8658; [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo)

**Searches for a solar eclipse visible at a specific location on the Earth's surface.**

This function finds the first solar eclipse that occurs after `startTime`.
A solar eclipse may be partial, annular, or total.
See [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) for more information.

To find a series of solar eclipses, call this function once,
then keep calling [`Astronomy.NextLocalSolarEclipse`](#Astronomy.NextLocalSolarEclipse) as many times as desired,
passing in the `peak` value returned from the previous call.

IMPORTANT: An eclipse reported by this function might be partly or
completely invisible to the observer due to the time of day.

See [`LocalSolarEclipseInfo`](#LocalSolarEclipseInfo) for more information about this topic.
See [`Astronomy.LocalSolarEclipsesAfter`](#Astronomy.LocalSolarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for a solar eclipse. |
| [`Observer`](#Observer) | `observer` | The geographic location of the observer. |

<a name="Astronomy.SearchLunarApsis"></a>
### Astronomy.SearchLunarApsis(startTime) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the date and time of the Moon's closest distance (perigee)
or farthest distance (apogee) with respect to the Earth.**

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the Moon reaches the closest or farthest point
in its orbit with respect to the center of the Earth, whichever comes first
after `startTime`.

The closest point is called *perigee* and the farthest point is called *apogee*.
The word *apsis* refers to either event.

To iterate through consecutive alternating perigee and apogee events, call `Astronomy.SearchLunarApsis`
once, then use the return value to call [`Astronomy.NextLunarApsis`](#Astronomy.NextLunarApsis). After that,
keep feeding the previous return value from `Astronomy.NextLunarApsis` into another
call of `Astronomy.NextLunarApsis` as many times as desired.

See [`Astronomy.LunarApsidesAfter`](#Astronomy.LunarApsidesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start searching for the next perigee or apogee. |

**Returns:** Returns an [`ApsisInfo`](#ApsisInfo) structure containing information about the next lunar apsis.

<a name="Astronomy.SearchLunarEclipse"></a>
### Astronomy.SearchLunarEclipse(startTime) &#8658; [`LunarEclipseInfo`](#LunarEclipseInfo)

**Searches for a lunar eclipse.**

This function finds the first lunar eclipse that occurs after `startTime`.
A lunar eclipse may be penumbral, partial, or total.
See [`LunarEclipseInfo`](#LunarEclipseInfo) for more information.
To find a series of lunar eclipses, call this function once,
then keep calling [`Astronomy.NextLunarEclipse`](#Astronomy.NextLunarEclipse) as many times as desired,
passing in the `center` value returned from the previous call.

See [`Astronomy.LunarEclipsesAfter`](#Astronomy.LunarEclipsesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for a lunar eclipse. |

**Returns:** A [`LunarEclipseInfo`](#LunarEclipseInfo) structure containing information about the lunar eclipse.

<a name="Astronomy.SearchMaxElongation"></a>
### Astronomy.SearchMaxElongation(body, startTime) &#8658; [`ElongationInfo`](#ElongationInfo)

**Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.**

Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
a telescope without atmospheric interference, are when these planets reach maximum elongation.
These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

This function solves for those times, reporting the next maximum elongation event's date and time,
the elongation value itself, the relative longitude with the Sun, and whether the planet is best
observed in the morning or evening. See [`Astronomy.Elongation`](#Astronomy.Elongation) for more details about the returned structure.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception. To find the best viewing opportunites for planets farther from the Sun than the Earth is (Mars through Pluto) use [`Astronomy.SearchRelativeLongitude`](#Astronomy.SearchRelativeLongitude) to find the next opposition event. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to begin the search. The maximum elongation event found will always be the first one that occurs after this date and time. |

**Returns:** Either an exception will be thrown, or the function will return a valid value.

<a name="Astronomy.SearchMoonNode"></a>
### Astronomy.SearchMoonNode(startTime) &#8658; [`NodeEventInfo`](#NodeEventInfo)

**Searches for a time when the Moon's center crosses through the ecliptic plane.**

Searches for the first ascending or descending node of the Moon after `startTime`.
An ascending node is when the Moon's center passes through the ecliptic plane
(the plane of the Earth's orbit around the Sun) from south to north.
A descending node is when the Moon's center passes through the ecliptic plane
from north to south. Nodes indicate possible times of solar or lunar eclipses,
if the Moon also happens to be in the correct phase (new or full, respectively).
Call `Astronomy.SearchMoonNode` to find the first of a series of nodes.
Then call [`Astronomy.NextMoonNode`](#Astronomy.NextMoonNode) to find as many more consecutive nodes as desired.

See [`Astronomy.MoonNodesAfter`](#Astronomy.MoonNodesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for an ascending or descending node of the Moon. |

<a name="Astronomy.SearchMoonPhase"></a>
### Astronomy.SearchMoonPhase(targetLon, startTime, limitDays) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time that the Moon reaches a specified phase.**

Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
longitude with respect to the Sun's geocentric ecliptic longitude.
When the Moon and the Sun have the same longitude, that is defined as a new moon.
When their longitudes are 180 degrees apart, that is defined as a full moon.

This function searches for any value of the lunar phase expressed as an
angle in degrees in the range [0, 360).

If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
it is much easier to call the functions [`Astronomy.SearchMoonQuarter`](#Astronomy.SearchMoonQuarter) and [`Astronomy.NextMoonQuarter`](#Astronomy.NextMoonQuarter).
This function is useful for finding general phase angles outside those four quarters.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `targetLon` | The difference in geocentric longitude between the Sun and Moon that specifies the lunar phase being sought. This can be any value in the range [0, 360). Certain values have conventional names: 0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter. |
| [`AstroTime`](#AstroTime) | `startTime` | The beginning of the time window in which to search for the Moon reaching the specified phase. |
| `double` | `limitDays` | The number of days away from `startTime` that limits the time window for the search. If the value is negative, the search is performed into the past from `startTime`. Otherwise, the search is performed into the future from `startTime`. |

**Returns:** If successful, returns the date and time the moon reaches the phase specified by `targetlon`. This function will return `null` if the phase does not occur within `limitDays` of `startTime`; that is, if the search window is too small.

<a name="Astronomy.SearchMoonQuarter"></a>
### Astronomy.SearchMoonQuarter(startTime) &#8658; [`MoonQuarterInfo`](#MoonQuarterInfo)

**Finds the first lunar quarter after the specified date and time.**

A lunar quarter is one of the following four lunar phase events:
new moon, first quarter, full moon, third quarter.
This function finds the lunar quarter that happens soonest
after the specified date and time.

To continue iterating through consecutive lunar quarters, call this function once,
followed by calls to [`Astronomy.NextMoonQuarter`](#Astronomy.NextMoonQuarter) as many times as desired.

See [`Astronomy.MoonQuartersAfter`](#Astronomy.MoonQuartersAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start the search. |

**Returns:** A [`MoonQuarterInfo`](#MoonQuarterInfo) structure reporting the next quarter phase and the time it will occur.

<a name="Astronomy.SearchPeakMagnitude"></a>
### Astronomy.SearchPeakMagnitude(body, startTime) &#8658; [`IllumInfo`](#IllumInfo)

**Searches for the date and time Venus will next appear brightest as seen from the Earth.**

This function searches for the date and time Venus appears brightest as seen from the Earth.
Currently only Venus is supported for the `body` parameter, though this could change in the future.
Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth,
so peak magnitude events have little practical value for that planet.
Planets other than Venus and Mercury reach peak magnitude at opposition, which can
be found using [`Astronomy.SearchRelativeLongitude`](#Astronomy.SearchRelativeLongitude).
The Moon reaches peak magnitude at full moon, which can be found using
[`Astronomy.SearchMoonQuarter`](#Astronomy.SearchMoonQuarter) or [`Astronomy.SearchMoonPhase`](#Astronomy.SearchMoonPhase).
The Sun reaches peak magnitude at perihelion, which occurs each year in January.
However, the difference is minor and has little practical value.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | Currently only `Body.Venus` is allowed. Any other value causes an exception. See remarks above for more details. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time to start searching for the next peak magnitude event. |

**Returns:** See documentation about the return value from [`Astronomy.Illumination`](#Astronomy.Illumination).

<a name="Astronomy.SearchPlanetApsis"></a>
### Astronomy.SearchPlanetApsis(body, startTime) &#8658; [`ApsisInfo`](#ApsisInfo)

**Finds the date and time of a planet's perihelion (closest approach to the Sun)
or aphelion (farthest distance from the Sun) after a given time.**

Given a date and time to start the search in `startTime`, this function finds the
next date and time that the center of the specified planet reaches the closest or farthest point
in its orbit with respect to the center of the Sun, whichever comes first
after `startTime`.

The closest point is called *perihelion* and the farthest point is called *aphelion*.
The word *apsis* refers to either event.

To iterate through consecutive alternating perihelion and aphelion events,
call `Astronomy.SearchPlanetApsis` once, then use the return value to call
[`Astronomy.NextPlanetApsis`](#Astronomy.NextPlanetApsis). After that, keep feeding the previous return value
from `Astronomy.NextPlanetApsis` into another call of `Astronomy.NextPlanetApsis`
as many times as desired.

See [`Astronomy.PlanetApsidesAfter`](#Astronomy.PlanetApsidesAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet for which to find the next perihelion/aphelion event. Not allowed to be `Body.Sun` or `Body.Moon`. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start searching for the next perihelion or aphelion. |

**Returns:** Returns a structure in which `time` holds the date and time of the next planetary apsis, `kind` holds either `ApsisKind.Pericenter` for perihelion or `ApsisKind.Apocenter` for aphelion. and distance values `dist_au` (astronomical units) and `dist_km` (kilometers).

<a name="Astronomy.SearchRelativeLongitude"></a>
### Astronomy.SearchRelativeLongitude(body, targetRelLon, startTime) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time when the Earth and another planet are separated by a specified angle
in ecliptic longitude, as seen from the Sun.**

A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit
(the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored.
If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around
that plane from one body to the other in the direction the planets orbit the Sun, you will get an
angle somewhere between 0 and 360 degrees. This is the relative longitude.

Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
this function searches for the next time that the relative longitude measured from the planet
to the Earth is `targetRelLon`.

Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:

- When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun.
  For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*,
  a time when the other planet becomes very difficult to see because of being lost in the Sun's glare.
  (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing
  between the Earth and the Sun.)

- When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
  this is known as *opposition*.  Opposition is when the planet is closest to the Earth, and
  also when it is visible for most of the night, so it is considered the best time to observe the planet.

- When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun
  from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is
  very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | A planet other than the Earth. If `body` is `Body.Earth`, `Body.Sun`, or `Body.Moon`, this function throws an exception. |
| `double` | `targetRelLon` | The desired relative longitude, expressed in degrees. Must be in the range [0, 360). |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to begin the search. |

**Returns:** The date and time of the relative longitude event.

<a name="Astronomy.SearchRiseSet"></a>
### Astronomy.SearchRiseSet(body, observer, direction, startTime, limitDays, metersAboveGround) &#8658; [`AstroTime`](#AstroTime)

**Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.**

This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
Rise time is when the body first starts to be visible above the horizon.
For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
Set time is the moment when the body appears to vanish below the horizon.
Therefore, this function adjusts for the apparent angular radius of the observed body
(significant only for the Sun and Moon).

This function corrects for a typical value of atmospheric refraction, which causes celestial
bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
Astronomy Engine uses a correction of 34 arcminutes. Real-world refraction varies based
on air temperature, pressure, and humidity; such weather-based conditions are outside
the scope of Astronomy Engine.

Note that rise or set may not occur in every 24 hour period.
For example, near the Earth's poles, there are long periods of time where
the Sun stays below the horizon, never rising.
Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
significant amount during each rotation of the Earth.
Therefore callers must not assume that the function will always succeed.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [`Astronomy.DefineStar`](#Astronomy.DefineStar). |
| [`Observer`](#Observer) | `observer` | The location where observation takes place. |
| [`Direction`](#Direction) | `direction` | Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time at which to start the search. |
| `double` | `limitDays` | Limits how many days to search for a rise or set time, and defines the direction in time to search. When `limitDays` is positive, the search is performed into the future, after `startTime`. When negative, the search is performed into the past, before `startTime`. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |
| `double` | `metersAboveGround` | Usually the observer is located at ground level. Then this parameter should be zero. But if the observer is significantly higher than ground level, for example in an airplane, this parameter should be a positive number indicating how far above the ground the observer is. An exception occurs if `metersAboveGround` is negative. |

**Returns:** On success, returns the date and time of the rise or set time as requested. If the function returns `null`, it means the rise or set event does not occur within `limitDays` days of `startTime`. This is a normal condition, not an error.

<a name="Astronomy.SearchSunLongitude"></a>
### Astronomy.SearchSunLongitude(targetLon, startTime, limitDays) &#8658; [`AstroTime`](#AstroTime)

**Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.**

This function finds the moment in time, if any exists in the given time window,
that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

This function can be used to determine equinoxes and solstices.
However, it is usually more convenient and efficient to call [`Astronomy.Seasons`](#Astronomy.Seasons)
to calculate all equinoxes and solstices for a given calendar year.

The function searches the window of time specified by `startTime` and `startTime+limitDays`.
The search will return `null` if the Sun never reaches the longitude `targetLon` or
if the window is so large that the longitude ranges more than 180 degrees within it.
It is recommended to keep the window smaller than 10 days when possible.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `targetLon` | The desired ecliptic longitude in degrees, relative to the true equinox of date. This may be any value in the range [0, 360), although certain values have conventional meanings: 0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for the desired longitude event. |
| `double` | `limitDays` | The real-valued number of days, which when added to `startTime`, limits the range of time over which the search looks. It is recommended to keep this value between 1 and 10 days. See remarks above for more details. |

**Returns:** The date and time when the Sun reaches the specified apparent ecliptic longitude.

<a name="Astronomy.SearchTransit"></a>
### Astronomy.SearchTransit(body, startTime) &#8658; [`TransitInfo`](#TransitInfo)

**Searches for the first transit of Mercury or Venus after a given date.**

Finds the first transit of Mercury or Venus after a specified date.
A transit is when an inferior planet passes between the Sun and the Earth
so that the silhouette of the planet is visible against the Sun in the background.
To continue the search, pass the `finish` time in the returned structure to
[`Astronomy.NextTransit`](#Astronomy.NextTransit).

See [`Astronomy.TransitsAfter`](#Astronomy.TransitsAfter) for a convenient enumerator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`. |
| [`AstroTime`](#AstroTime) | `startTime` | The date and time for starting the search for a transit. |

<a name="Astronomy.Seasons"></a>
### Astronomy.Seasons(year) &#8658; [`SeasonsInfo`](#SeasonsInfo)

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

**Returns:** A [`SeasonsInfo`](#SeasonsInfo) structure that contains four [`AstroTime`](#AstroTime) values: the March and September equinoxes and the June and December solstices.

<a name="Astronomy.SiderealTime"></a>
### Astronomy.SiderealTime(time) &#8658; `double`

**Calculates Greenwich Apparent Sidereal Time (GAST).**

Given a date and time, this function calculates the rotation of the
Earth, represented by the equatorial angle of the Greenwich prime meridian
with respect to distant stars (not the Sun, which moves relative to background
stars by almost one degree per day).
This angle is called Greenwich Apparent Sidereal Time (GAST).
GAST is measured in sidereal hours in the half-open range [0, 24).
When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
corrected at that time for precession and nutation of the Earth's axis.
In this context, the "equinox" is the direction in space where the Earth's
orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
at the location on the Earth's orbit of the (seasonal) March equinox.
As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
then starts over at 0.
To convert to degrees, multiply the return value by 15.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to find GAST. As an optimization, this function caches the sidereal time value in `time`, unless it has already been cached, in which case the cached value is reused. |

**Returns:** GAST in sidereal hours.

<a name="Astronomy.SphereFromVector"></a>
### Astronomy.SphereFromVector(vector) &#8658; [`Spherical`](#Spherical)

**Converts Cartesian coordinates to spherical coordinates.**

Given a Cartesian vector, returns latitude, longitude, and distance.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `vector` | Cartesian vector to be converted to spherical coordinates. |

**Returns:** Spherical coordinates that are equivalent to the given vector.

<a name="Astronomy.SunPosition"></a>
### Astronomy.SunPosition(time) &#8658; [`Ecliptic`](#Ecliptic)

**Calculates geocentric ecliptic coordinates for the Sun.**

This function calculates the position of the Sun as seen from the Earth.
The returned value includes both Cartesian and spherical coordinates.
The x-coordinate and longitude values in the returned structure are based
on the *true equinox of date*: one of two points in the sky where the instantaneous
plane of the Earth's equator at the given date and time (the *equatorial plane*)
intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
By convention, the apparent location of the Sun at the March equinox is chosen
as the longitude origin and x-axis direction, instead of the one for September.

`SunPosition` corrects for precession and nutation of the Earth's axis
in order to obtain the exact equatorial plane at the given time.

This function can be used for calculating changes of seasons: equinoxes and solstices.
In fact, the function [`Astronomy.Seasons`](#Astronomy.Seasons) does use this function for that purpose.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time for which to calculate the Sun's position. |

**Returns:** The ecliptic coordinates of the Sun using the Earth's true equator of date.

<a name="Astronomy.TransitsAfter"></a>
### Astronomy.TransitsAfter(body, startTime) &#8658; `IEnumerable<`[`TransitInfo`](#TransitInfo)`>`

**Enumerates a series of transits of Mercury or Venus.**

This is a convenience wrapper around the functions
[`Astronomy.SearchTransit`](#Astronomy.SearchTransit) and [`Astronomy.NextTransit`](#Astronomy.NextTransit).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The planet whose transits are to be enumerated. Must be `Body.Mercury` or `Body.Venus`. |
| [`AstroTime`](#AstroTime) | `startTime` | Specifies the time to begin searching for consecutive transits. |

<a name="Astronomy.VectorFromHorizon"></a>
### Astronomy.VectorFromHorizon(sphere, time, refraction) &#8658; [`AstroVector`](#AstroVector)

**Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | A structure that contains apparent horizontal coordinates: `lat` holds the refracted altitude angle, `lon` holds the azimuth in degrees clockwise from north, and `dist` holds the distance from the observer to the object in AU. |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. This is needed because the returned [`AstroVector`](#AstroVector) requires a valid time value when passed to certain other functions. |
| [`Refraction`](#Refraction) | `refraction` | The refraction option used to model atmospheric lensing. See [`Astronomy.RefractionAngle`](#Astronomy.RefractionAngle). This specifies how refraction is to be removed from the altitude stored in `sphere.lat`. |

**Returns:** A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).

<a name="Astronomy.VectorFromSphere"></a>
### Astronomy.VectorFromSphere(sphere, time) &#8658; [`AstroVector`](#AstroVector)

**Converts spherical coordinates to Cartesian coordinates.**

Given spherical coordinates and a time at which they are valid,
returns a vector of Cartesian coordinates. The returned value
includes the time, as required by the type [`AstroVector`](#AstroVector).

| Type | Parameter | Description |
| --- | --- | --- |
| [`Spherical`](#Spherical) | `sphere` | Spherical coordinates to be converted. |
| [`AstroTime`](#AstroTime) | `time` | The time that should be included in the return value. |

**Returns:** The vector form of the supplied spherical coordinates.

<a name="Astronomy.VectorObserver"></a>
### Astronomy.VectorObserver(vector, equdate) &#8658; [`Observer`](#Observer)

**Calculates the geographic location corresponding to an equatorial vector.**

This is the inverse function of [`Astronomy.ObserverVector`](#Astronomy.ObserverVector).
Given a geocentric equatorial vector, it returns the geographic
latitude, longitude, and elevation for that vector.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `vector` | The geocentric equatorial position vector for which to find geographic coordinates. The components are expressed in Astronomical Units (AU). You can calculate AU by dividing kilometers by the constant [`Astronomy.KM_PER_AU`](#Astronomy.KM_PER_AU). The time `vector.t` determines the Earth's rotation. |
| [`EquatorEpoch`](#EquatorEpoch) | `equdate` | Selects the date of the Earth's equator in which `vector` is expressed. The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator at noon UTC on January 1, 2000, in which case this function corrects for precession and nutation of the Earth as it was at the moment specified by `vector.t`. Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `vector.t` as the orientation. |

**Returns:** The geographic latitude, longitude, and elevation above sea level that corresponds to the given equatorial vector.

---

<a name="types"></a>
## Types

---

<a name="Aberration"></a>
## `enum Aberration`

**Aberration calculation options.**

[Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect
causing the apparent direction of an observed body to be shifted due to transverse
movement of the Earth with respect to the rays of light coming from that body.
This angular correction can be anywhere from 0 to about 20 arcseconds,
depending on the position of the observed body relative to the instantaneous
velocity vector of the Earth.

Some Astronomy Engine functions allow optional correction for aberration by
passing in a value of this enumerated type.

Aberration correction is useful to improve accuracy of coordinates of
apparent locations of bodies seen from the Earth.
However, because aberration affects not only the observed body (such as a planet)
but the surrounding stars, aberration may be unhelpful (for example)
for determining exactly when a planet crosses from one constellation to another.

| Value | Description |
| --- | --- |
| `Corrected` | Request correction for aberration. |
| `None` | Do not correct for aberration. |

---

<a name="ApsisInfo"></a>
## `struct ApsisInfo`

**An apsis event: pericenter (closest approach) or apocenter (farthest distance).**

For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
event where the orbiting body reaches its closest or farthest point from the primary body.
The closest approach is called *pericenter* and the farthest point is *apocenter*.

More specific terminology is common for particular orbiting bodies.
The Moon's closest approach to the Earth is called *perigee* and its farthest
point is called *apogee*. The closest approach of a planet to the Sun is called
*perihelion* and the furthest point is called *aphelion*.

This data structure is returned by [`Astronomy.SearchLunarApsis`](#Astronomy.SearchLunarApsis) and [`Astronomy.NextLunarApsis`](#Astronomy.NextLunarApsis)
to iterate through consecutive alternating perigees and apogees.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the apsis. |
| [`ApsisKind`](#ApsisKind) | `kind` | Whether this is a pericenter or apocenter event. |
| `double` | `dist_au` | The distance between the centers of the bodies in astronomical units. |
| `double` | `dist_km` | The distance between the centers of the bodies in kilometers. |

---

<a name="ApsisKind"></a>
## `enum ApsisKind`

**The type of apsis: pericenter (closest approach) or apocenter (farthest distance).**

| Value | Description |
| --- | --- |
| `Pericenter` | The body is at its closest approach to the object it orbits. |
| `Apocenter` | The body is at its farthest distance from the object it orbits. |

---

<a name="AstroTime"></a>
## `class AstroTime`

**A date and time used for astronomical calculations.**

### constructors

### new AstroTime(ut)

**Creates an `AstroTime` object from a Universal Time day value.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` | The number of days after the J2000 epoch. |

### new AstroTime(d)

**Creates an `AstroTime` object from a .NET `DateTime` object.**

| Type | Parameter | Description |
| --- | --- | --- |
| `DateTime` | `d` | The date and time to be converted to AstroTime format. |

### new AstroTime(year, month, day, hour, minute, second)

**Creates an `AstroTime` object from a UTC year, month, day, hour, minute and second.**

| Type | Parameter | Description |
| --- | --- | --- |
| `int` | `year` | The UTC year value. |
| `int` | `month` | The UTC month value 1..12. |
| `int` | `day` | The UTC day of the month 1..31. |
| `int` | `hour` | The UTC hour value 0..23. |
| `int` | `minute` | The UTC minute value 0..59. |
| `double` | `second` | The UTC second value [0, 60). |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `ut` | UT1/UTC number of days since noon on January 1, 2000. |
| `double` | `tt` | Terrestrial Time days since noon on January 1, 2000. |

### properties

| Type | Name | Description |
| --- | --- | --- |
| `double` | `Psi` | Nutation angle `psi`. Intended for unit testing only. |
| `double` | `Eps` | Nutation angle `eps`. Intended for unit testing only. |

### member functions

<a name="AstroTime.AddDays"></a>
### AstroTime.AddDays(days) &#8658; [`AstroTime`](#AstroTime)

**Calculates the sum or difference of an [`AstroTime`](#AstroTime) with a specified floating point number of days.**

Sometimes we need to adjust a given [`AstroTime`](#AstroTime) value by a certain amount of time.
This function adds the given real number of days in `days` to the date and time in this object.

More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
using a best-fit piecewise polynomial model devised by
[Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `days` | A floating point number of days by which to adjust `time`. May be negative, 0, or positive. |

**Returns:** A date and time that is conceptually equal to `time + days`.

<a name="AstroTime.FromTerrestrialTime"></a>
### AstroTime.FromTerrestrialTime(tt) &#8658; [`AstroTime`](#AstroTime)

**Creates an `AstroTime` object from a Terrestrial Time day value.**

This function can be used in rare cases where a time must be based
on Terrestrial Time (TT) rather than Universal Time (UT).
Most developers will want to invoke `new AstroTime(ut)` with a universal time
instead of this function, because usually time is based on civil time adjusted
by leap seconds to match the Earth's rotation, rather than the uniformly
flowing TT used to calculate solar system dynamics. In rare cases
where the caller already knows TT, this function is provided to create
an `AstroTime` value that can be passed to Astronomy Engine functions.

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `tt` | The number of days after the J2000 epoch. |

<a name="AstroTime.ToCalendarDateTime"></a>
### AstroTime.ToCalendarDateTime() &#8658; [`CalendarDateTime`](#CalendarDateTime)

**Converts this object to our custom type [`CalendarDateTime`](#CalendarDateTime).**

The .NET type `DateTime` can only represent years in the range 0000..9999.
However, the Astronomy Engine type [`CalendarDateTime`](#CalendarDateTime) can represent
years in the range -999999..+999999. This is a time span of nearly 2 million years.
This function converts this `AstroTime` object to an equivalent Gregorian calendar representation.

<a name="AstroTime.ToString"></a>
### AstroTime.ToString() &#8658; `string`

**Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.**

**Returns:** Example: "2019-08-30T17:45:22.763Z".

<a name="AstroTime.ToUtcDateTime"></a>
### AstroTime.ToUtcDateTime() &#8658; `DateTime`

**Converts this object to .NET `DateTime` format.**

**Returns:** a UTC `DateTime` object for this `AstroTime` value.

---

<a name="AstroVector"></a>
## `struct AstroVector`

**A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).**

### constructors

### new AstroVector(x, y, z, t)

**Creates an AstroVector.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `x` | A Cartesian x-coordinate expressed in AU. |
| `double` | `y` | A Cartesian y-coordinate expressed in AU. |
| `double` | `z` | A Cartesian z-coordinate expressed in AU. |
| [`AstroTime`](#AstroTime) | `t` | The date and time at which this vector is valid. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `x` | The Cartesian x-coordinate of the vector in AU. |
| `double` | `y` | The Cartesian y-coordinate of the vector in AU. |
| `double` | `z` | The Cartesian z-coordinate of the vector in AU. |
| [`AstroTime`](#AstroTime) | `t` | The date and time at which this vector is valid. |

### member functions

<a name="AstroVector.Length"></a>
### AstroVector.Length() &#8658; `double`

**Calculates the total distance in AU represented by this vector.**

**Returns:** The nonnegative length of the Cartisian vector in AU.

<a name="AstroVector.ToString"></a>
### AstroVector.ToString() &#8658; `string`

**Converts the vector to a string of the format (x, y, z, t).**

---

<a name="AtmosphereInfo"></a>
## `struct AtmosphereInfo`

**Information about idealized atmospheric variables at a given elevation.**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `pressure` | Atmospheric pressure in pascals. |
| `double` | `temperature` | Atmospheric temperature in kelvins. |
| `double` | `density` | Atmospheric density relative to sea level. |

---

<a name="AxisInfo"></a>
## `struct AxisInfo`

**Information about a body's rotation axis at a given time.**

This structure is returned by [`Astronomy.RotationAxis`](#Astronomy.RotationAxis) to report
the orientation of a body's rotation axis at a given moment in time.
The axis is specified by the direction in space that the body's north pole
points, using angular equatorial coordinates in the J2000 system (EQJ).

Thus `ra` is the right ascension, and `dec` is the declination, of the
body's north pole vector at the given moment in time. The north pole
of a body is defined as the pole that lies on the north side of the
[Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
regardless of the body's direction of rotation.

The `spin` field indicates the angular position of a prime meridian
arbitrarily recommended for the body by the International Astronomical
Union (IAU).

The fields `ra`, `dec`, and `spin` correspond to the variables
α0, δ0, and W, respectively, from
[Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).

The field `north` is a unit vector pointing in the direction of the body's north pole.
It is expressed in the J2000 mean equator system (EQJ).

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `ra` | The J2000 right ascension of the body's north pole direction, in sidereal hours. |
| `double` | `dec` | The J2000 declination of the body's north pole direction, in degrees. |
| `double` | `spin` | Rotation angle of the body's prime meridian, in degrees. |
| [`AstroVector`](#AstroVector) | `north` | A J2000 dimensionless unit vector pointing in the direction of the body's north pole. |

---

<a name="Body"></a>
## `enum Body`

**The enumeration of celestial bodies supported by Astronomy Engine.**

| Value | Description |
| --- | --- |
| `Invalid` | A placeholder value representing an invalid or unknown celestial body. |
| `Mercury` | The planet Mercury. |
| `Venus` | The planet Venus. |
| `Earth` | The planet Earth. Some functions that accept a `Body` parameter will fail if passed this value because they assume that an observation is being made from the Earth, and therefore the Earth is not a target of observation. |
| `Mars` | The planet Mars. |
| `Jupiter` | The planet Jupiter. |
| `Saturn` | The planet Saturn. |
| `Uranus` | The planet Uranus. |
| `Neptune` | The planet Neptune. |
| `Pluto` | The planet Pluto. |
| `Sun` | The Sun. |
| `Moon` | The Earth's natural satellite, the Moon. |
| `EMB` | The Earth/Moon Barycenter. |
| `SSB` | The Solar System Barycenter. |
| `Star1` | User-defined star #1. |
| `Star2` | User-defined star #2. |
| `Star3` | User-defined star #3. |
| `Star4` | User-defined star #4. |
| `Star5` | User-defined star #5. |
| `Star6` | User-defined star #6. |
| `Star7` | User-defined star #7. |
| `Star8` | User-defined star #8. |

---

<a name="CalendarDateTime"></a>
## `struct CalendarDateTime`

**Represents a Gregorian calendar date and time within plus or minus 1 million years from the year 0.**

The C# standard type `System.DateTime` only allows years from 0001 to 9999.
However, the [`AstroTime`](#AstroTime) class can represent years in the range -999999 to +999999.
In order to support formatting dates with extreme year values in an extrapolated
Gregorian calendar, the `CalendarDateTime` class breaks out the components of
a date into separate fields.

### constructors

### new CalendarDateTime(ut)

**Convert a J2000 day value to a Gregorian calendar date.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `ut` | The real-valued number of days since the J2000 epoch. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `int` | `year` | The year value in the range -999999 to +999999. |
| `int` | `month` | The calendar month in the range 1..12. |
| `int` | `day` | The day of the month in the reange 1..31. |
| `int` | `hour` | The hour in the range 0..23. |
| `int` | `minute` | The minute in the range 0..59. |
| `double` | `second` | The real-valued second in the half-open range [0, 60). |

### member functions

<a name="CalendarDateTime.ToString"></a>
### CalendarDateTime.ToString() &#8658; `string`

**Converts this `CalendarDateTime` to ISO 8601 format, expressed in UTC with millisecond resolution.**

**Returns:** Example: "2019-08-30T17:45:22.763Z".

---

<a name="ConstellationInfo"></a>
## `struct ConstellationInfo`

**Reports the constellation that a given celestial point lies within.**

The [`Astronomy.Constellation`](#Astronomy.Constellation) function returns this struct
to report which constellation corresponds with a given point in the sky.
Constellations are defined with respect to the B1875 equatorial system
per IAU standard. Although `Astronomy.Constellation` requires J2000 equatorial
coordinates, the struct contains converted B1875 coordinates for reference.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `string` | `Symbol` | 3-character mnemonic symbol for the constellation, e.g. "Ori". |
| `string` | `Name` | Full name of constellation, e.g. "Orion". |
| `double` | `Ra1875` | Right ascension expressed in B1875 coordinates. |
| `double` | `Dec1875` | Declination expressed in B1875 coordinates. |

---

<a name="DeltaTimeFunc"></a>
## `class DeltaTimeFunc`

**Defines a function type for calculating Delta T.**

Delta T is the discrepancy between times measured using an atomic clock
and times based on observations of the Earth's rotation, which is gradually
slowing down over time. Delta T = TT - UT, where
TT = Terrestrial Time, based on atomic time, and
UT = Universal Time, civil time based on the Earth's rotation.
Astronomy Engine defaults to using a Delta T function defined by
Espenak and Meeus in their "Five Millennium Canon of Solar Eclipses".
See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

### member functions

---

<a name="Direction"></a>
## `enum Direction`

**Selects whether to search for a rising event or a setting event for a celestial body.**

| Value | Description |
| --- | --- |
| `Rise` | Indicates a rising event: a celestial body is observed to rise above the horizon by an observer on the Earth. |
| `Set` | Indicates a setting event: a celestial body is observed to sink below the horizon by an observer on the Earth. |

---

<a name="EarthNotAllowedException"></a>
## `class EarthNotAllowedException`

**This exception is thrown by certain Astronomy Engine functions
when an invalid attempt is made to use the Earth as the observed
celestial body. Usually this happens for cases where the Earth itself
is the location of the observer.**

---

<a name="EclipseEvent"></a>
## `struct EclipseEvent`

**Holds a time and the observed altitude of the Sun at that time.**

When reporting a solar eclipse observed at a specific location on the Earth
(a "local" solar eclipse), a series of events occur. In addition
to the time of each event, it is important to know the altitude of the Sun,
because each event may be invisible to the observer if the Sun is below
the horizon.

If `altitude` is negative, the event is theoretical only; it would be
visible if the Earth were transparent, but the observer cannot actually see it.
If `altitude` is positive but less than a few degrees, visibility will be impaired by
atmospheric interference (sunrise or sunset conditions).

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the event. |
| `double` | `altitude` | The angular altitude of the center of the Sun above/below the horizon, at `time`, corrected for atmospheric refraction and expressed in degrees. |

---

<a name="EclipseKind"></a>
## `enum EclipseKind`

**The different kinds of lunar/solar eclipses.**

| Value | Description |
| --- | --- |
| `None` | No eclipse found. |
| `Penumbral` | A penumbral lunar eclipse. (Never used for a solar eclipse.) |
| `Partial` | A partial lunar/solar eclipse. |
| `Annular` | An annular solar eclipse. (Never used for a lunar eclipse.) |
| `Total` | A total lunar/solar eclipse. |

---

<a name="Ecliptic"></a>
## `struct Ecliptic`

**Ecliptic angular and Cartesian coordinates.**

Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `vec` | Cartesian ecliptic vector, with components as follows: x: the direction of the equinox along the ecliptic plane. y: in the ecliptic plane 90 degrees prograde from the equinox. z: perpendicular to the ecliptic plane. Positive is north. |
| `double` | `elat` | Latitude in degrees north (positive) or south (negative) of the ecliptic plane. |
| `double` | `elon` | Longitude in degrees around the ecliptic plane prograde from the equinox. |

---

<a name="ElongationInfo"></a>
## `struct ElongationInfo`

**Contains information about the visibility of a celestial body at a given date and time.
See [`Astronomy.Elongation`](#Astronomy.Elongation) for more detailed information about the members of this structure.
See also [`Astronomy.SearchMaxElongation`](#Astronomy.SearchMaxElongation) for how to search for maximum elongation events.**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |
| [`Visibility`](#Visibility) | `visibility` | Whether the body is best seen in the morning or the evening. |
| `double` | `elongation` | The angle in degrees between the body and the Sun, as seen from the Earth. |
| `double` | `ecliptic_separation` | The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth. |

---

<a name="EquatorEpoch"></a>
## `enum EquatorEpoch`

**Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.**

The Earth's equator is not always in the same plane due to precession and nutation.

Sometimes it is useful to have a fixed plane of reference for equatorial coordinates
across different calendar dates.  In these cases, a fixed *epoch*, or reference time,
is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers
to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.

For some other purposes, it is more helpful to represent coordinates using the Earth's
equator exactly as it is on that date. For example, when calculating rise/set times
or horizontal coordinates, it is most accurate to use the orientation of the Earth's
equator at that same date and time. For these uses, Astronomy Engine allows *of-date*
calculations.

| Value | Description |
| --- | --- |
| `J2000` | Represent equatorial coordinates in the J2000 epoch. |
| `OfDate` | Represent equatorial coordinates using the Earth's equator at the given date and time. |

---

<a name="Equatorial"></a>
## `struct Equatorial`

**Equatorial angular and cartesian coordinates.**

Coordinates of a celestial body as seen from the Earth
(geocentric or topocentric, depending on context),
oriented with respect to the projection of the Earth's equator onto the sky.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `ra` | Right ascension in sidereal hours. |
| `double` | `dec` | Declination in degrees. |
| `double` | `dist` | Distance to the celestial body in AU. |
| [`AstroVector`](#AstroVector) | `vec` | Equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north. |

---

<a name="GlobalSolarEclipseInfo"></a>
## `struct GlobalSolarEclipseInfo`

**Reports the time and geographic location of the peak of a solar eclipse.**

Returned by [`Astronomy.SearchGlobalSolarEclipse`](#Astronomy.SearchGlobalSolarEclipse) or [`Astronomy.NextGlobalSolarEclipse`](#Astronomy.NextGlobalSolarEclipse)
to report information about a solar eclipse event.

The eclipse is classified as partial, annular, or total, depending on the
maximum amount of the Sun's disc obscured, as seen at the peak location
on the surface of the Earth.

The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
observes either a total or annular eclipse.

If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
fields give the geographic coordinates of the center of the Moon's shadow projected
onto the daytime side of the Earth at the instant of the eclipse's peak.
If `kind` has any other value, `latitude` and `longitude` are undefined and should
not be used.

For total or annular eclipses, the `obscuration` field holds the fraction (0, 1]
of the Sun's apparent disc area that is blocked from view by the Moon's silhouette,
as seen by an observer located at the geographic coordinates `latitude`, `longitude`
at the darkest time `peak`. The value will always be 1 for total eclipses, and less than
1 for annular eclipses.
For partial eclipses, `obscuration` is undefined and should not be used.
This is because there is little practical use for an obscuration value of
a partial eclipse without supplying a particular observation location.
Developers who wish to find an obscuration value for partial solar eclipses should therefore use
[`Astronomy.SearchLocalSolarEclipse`](#Astronomy.SearchLocalSolarEclipse) and provide the geographic coordinates of an observer.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`. |
| `double` | `obscuration` | The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only). |
| [`AstroTime`](#AstroTime) | `peak` | The date and time when the solar eclipse is at its darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center. |
| `double` | `distance` | The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers. |
| `double` | `latitude` | The geographic latitude at the center of the peak eclipse shadow. |
| `double` | `longitude` | The geographic longitude at the center of the peak eclipse shadow. |

---

<a name="GravitySimulator"></a>
## `class GravitySimulator`

**A simulation of zero or more small bodies moving through the Solar System.**

This class calculates the movement of arbitrary small bodies,
such as asteroids or comets, that move through the Solar System.
It does so by calculating the gravitational forces on the small bodies
from the Sun and planets. The user of this class supplies an enumeration
of initial positions and velocities for the small bodies.
Then the class can update the positions and velocities over small time steps.

### constructors

### new GravitySimulator(originBody, time, bodyStates)

**Creates a gravity simulation object.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `originBody` | Specifies the origin of the reference frame. All position vectors and velocity vectors will use `originBody` as the origin of the coordinate system. This origin applies to all the input vectors provided in the `bodyStates` parameter of this function, along with all output vectors returned by [`GravitySimulator.Update`](#GravitySimulator.Update). Most callers will want to provide one of the following: `Body.Sun` for heliocentric coordinates, `Body.SSB` for solar system barycentric coordinates, or `Body.Earth` for geocentric coordinates. Note that the gravity simulator does not correct for light travel time; all state vectors are tied to a Newtonian "instantaneous" time. |
| [`AstroTime`](#AstroTime) | `time` | The initial time at which to start the simulation. |
| `IEnumerable<`[`StateVector`](#StateVector)`>` | `bodyStates` | An enumeration of zero or more initial state vectors (positions and velocities) of the small bodies to be simulated. The caller must know the positions and velocities of the small bodies at an initial moment in time. Their positions and velocities are expressed with respect to `originBody`, using equatorial J2000 orientation (EQJ). Positions are expressed in astronomical units (AU). Velocities are expressed in AU/day. All the times embedded within the state vectors must be exactly equal to `time`, or this constructor will throw an exception. If `bodyStates` is null, the gravity simulator will contain zero small bodies. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`Body`](#Body) | `OriginBody` | The origin of the reference frame. See constructor for more info. |

### properties

| Type | Name | Description |
| --- | --- | --- |
| `int` | `NumSmallBodies` | The number of small bodies that are included in this gravity simulation. |
| [`AstroTime`](#AstroTime) | `Time` | The time represented by the current step of the gravity simulation. |

### member functions

<a name="GravitySimulator.SolarSystemBodyState"></a>
### GravitySimulator.SolarSystemBodyState(body) &#8658; [`StateVector`](#StateVector)

**Get the position and velocity of a Solar System body included in the simulation.**

In order to simulate the movement of small bodies through the Solar System,
the simulator needs to calculate the state vectors for the Sun and planets.

If an application wants to know the positions of one or more of the planets
in addition to the small bodies, this function provides a way to obtain
their state vectors. This is provided for the sake of efficiency, to avoid
redundant calculations.

The state vector is returned relative to the position and velocity
of the `originBody` parameter that was passed to this object's constructor.

| Type | Parameter | Description |
| --- | --- | --- |
| [`Body`](#Body) | `body` | The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune. |

<a name="GravitySimulator.Swap"></a>
### GravitySimulator.Swap() &#8658; `void`

**Exchange the current time step with the previous time step.**

Sometimes it is helpful to "explore" various times near a given
simulation time step, while repeatedly returning to the original
time step. For example, when backdating a position for light travel
time, the caller may wish to repeatedly try different amounts of
backdating. When the backdating solver has converged, the caller
wants to leave the simulation in its original state.

This function allows a single "undo" of a simulation, and does so
very efficiently.

Usually this function will be called immediately after a matching
call to [`GravitySimulator.Update`](#GravitySimulator.Update). It has the effect of rolling
back the most recent update. If called twice in a row, it reverts
the swap and thus has no net effect.

The constructor initializes the current state and previous
state to be identical. Both states represent the `time` parameter that was
passed into the constructor. Therefore, `Swap` will
have no effect from the caller's point of view when passed a simulator
that has not yet been updated by a call to [`GravitySimulator.Update`](#GravitySimulator.Update).

<a name="GravitySimulator.Update"></a>
### GravitySimulator.Update(time, bodyStates) &#8658; `void`

**Advances a gravity simulation by a small time step.**

Updates the simulation of the user-supplied small bodies
to the time indicated by the `time` parameter.
Updates the supplied array `bodyStates` of state vectors for the small bodies.
This array must be the same size as the number of bodies supplied
to the constructor of this object.
The positions and velocities in the returned array are referenced
to the `originBody` that was used to construct this simulator.

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | A time that is a small increment away from the current simulation time. It is up to the developer to figure out an appropriate time increment. Depending on the trajectories, a smaller or larger increment may be needed for the desired accuracy. Some experimentation may be needed. Generally, bodies that stay in the outer Solar System and move slowly can use larger time steps. Bodies that pass into the inner Solar System and move faster will need a smaller time step to maintain accuracy. The `time` value may be after or before the current simulation time to move forward or backward in time. |
| [`StateVector`](#StateVector)`[]` | `bodyStates` | If this array is not null, it must contain exactly the same number of elements as the number of small bodies that were added when this simulator was created. The non-null array receives updated state vectors for the simulated small bodies. If `bodyStates` is null, the simulation is updated but without returning the state vectors. |

---

<a name="HourAngleInfo"></a>
## `struct HourAngleInfo`

**Information about a celestial body crossing a specific hour angle.**

Returned by the function [`Astronomy.SearchHourAngle`](#Astronomy.SearchHourAngle) to report information about
a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time when the body crosses the specified hour angle. |
| [`Topocentric`](#Topocentric) | `hor` | Apparent coordinates of the body at the time it crosses the specified hour angle. |

---

<a name="IllumInfo"></a>
## `struct IllumInfo`

**Information about the brightness and illuminated shape of a celestial body.**

Returned by the functions [`Astronomy.Illumination`](#Astronomy.Illumination) and [`Astronomy.SearchPeakMagnitude`](#Astronomy.SearchPeakMagnitude)
to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the observation. |
| `double` | `mag` | The visual magnitude of the body. Smaller values are brighter. |
| `double` | `phase_angle` | The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. |
| `double` | `phase_fraction` | A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth. |
| `double` | `helio_dist` | The distance between the Sun and the body at the observation time. |
| `double` | `ring_tilt` | For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0. |

---

<a name="InternalError"></a>
## `class InternalError`

**This exception indicates an unexpected error occurred inside Astronomy Engine.
Please report any such errors by creating an issue at:
https://github.com/cosinekitty/astronomy/issues**

---

<a name="InvalidBodyException"></a>
## `class InvalidBodyException`

**This exception is thrown by certain Astronomy Engine functions
when a body is specified that is not appropriate for the given operation.**

---

<a name="IPositionFunction"></a>
## `struct IPositionFunction`

**A function for which to solve a light-travel time problem.**

The function [`Astronomy.CorrectLightTravel`](#Astronomy.CorrectLightTravel) solves a generalized
problem of deducing how far in the past light must have left
a target object to be seen by an observer at a specified time.
This interface expresses an arbitrary position vector as
function of time that is passed to [`Astronomy.CorrectLightTravel`](#Astronomy.CorrectLightTravel).

### member functions

<a name="IPositionFunction.Position"></a>
### IPositionFunction.Position(time) &#8658; [`AstroVector`](#AstroVector)

**Returns a relative position vector for a given time.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The time at which to evaluate a relative position vector. |

---

<a name="JupiterMoonsInfo"></a>
## `struct JupiterMoonsInfo`

**Holds the positions and velocities of Jupiter's major 4 moons.**

The [`Astronomy.JupiterMoons`](#Astronomy.JupiterMoons) function returns an object of this type
to report position and velocity vectors for Jupiter's largest 4 moons
Io, Europa, Ganymede, and Callisto. Each position vector is relative
to the center of Jupiter. Both position and velocity are oriented in
the EQJ system (that is, using Earth's equator at the J2000 epoch).
The positions are expressed in astronomical units (AU),
and the velocities in AU/day.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`StateVector`](#StateVector) | `io` | The position and velocity of Jupiter's moon Io. |
| [`StateVector`](#StateVector) | `europa` | The position and velocity of Jupiter's moon Europa. |
| [`StateVector`](#StateVector) | `ganymede` | The position and velocity of Jupiter's moon Ganymede. |
| [`StateVector`](#StateVector) | `callisto` | The position and velocity of Jupiter's moon Callisto. |

---

<a name="LibrationInfo"></a>
## `struct LibrationInfo`

**Lunar libration angles, returned by [`Astronomy.Libration`](#Astronomy.Libration).**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `elat` | Sub-Earth libration ecliptic latitude angle, in degrees. |
| `double` | `elon` | Sub-Earth libration ecliptic longitude angle, in degrees. |
| `double` | `mlat` | Moon's geocentric ecliptic latitude in degrees. |
| `double` | `mlon` | Moon's geocentric ecliptic longitude in degrees. |
| `double` | `dist_km` | Distance between the centers of the Earth and Moon in kilometers. |
| `double` | `diam_deg` | The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth. |

---

<a name="LocalSolarEclipseInfo"></a>
## `struct LocalSolarEclipseInfo`

**Information about a solar eclipse as seen by an observer at a given time and geographic location.**

Returned by [`Astronomy.SearchLocalSolarEclipse`](#Astronomy.SearchLocalSolarEclipse) or [`Astronomy.NextLocalSolarEclipse`](#Astronomy.NextLocalSolarEclipse)
to report information about a solar eclipse as seen at a given geographic location.

When a solar eclipse is found, it is classified as partial, annular, or total.
The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
A partial solar eclipse is when the Moon does not line up directly enough with the Sun
to completely block the Sun's light from reaching the observer.
An annular eclipse occurs when the Moon's disc is completely visible against the Sun
but the Moon is too far away to completely block the Sun's light; this leaves the
Sun with a ring-like appearance.
A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
Sun just right to completely block all sunlight from reaching the observer.

The `obscuration` field reports what fraction of the Sun's disc appears blocked
by the Moon when viewed by the observer at the peak eclipse time.
This is a value that ranges from 0 (no blockage) to 1 (total eclipse).
The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses.
The value will be exactly 1 for total eclipses. Obscuration gives an indication
of how dark the eclipse appears.

There are 5 "event" fields, each of which contains a time and a solar altitude.
Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
The fields `partial_begin` and `partial_end` are always set, and indicate when
the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
`total_begin` and `total_end` indicate when the total/annular phase begins/ends.
When an event field is valid, the caller must also check its `altitude` field to
see whether the Sun is above the horizon at the time indicated by the `time` field.
See [`EclipseEvent`](#EclipseEvent) for more information.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`. |
| `double` | `obscuration` | The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak. |
| [`EclipseEvent`](#EclipseEvent) | `partial_begin` | The time and Sun altitude at the beginning of the eclipse. |
| [`EclipseEvent`](#EclipseEvent) | `total_begin` | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid. |
| [`EclipseEvent`](#EclipseEvent) | `peak` | The time and Sun altitude when the eclipse reaches its peak. |
| [`EclipseEvent`](#EclipseEvent) | `total_end` | If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid. |
| [`EclipseEvent`](#EclipseEvent) | `partial_end` | The time and Sun altitude at the end of the eclipse. |

---

<a name="LunarEclipseInfo"></a>
## `struct LunarEclipseInfo`

**Information about a lunar eclipse.**

Returned by [`Astronomy.SearchLunarEclipse`](#Astronomy.SearchLunarEclipse) or [`Astronomy.NextLunarEclipse`](#Astronomy.NextLunarEclipse)
to report information about a lunar eclipse event.
When a lunar eclipse is found, it is classified as penumbral, partial, or total.
Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
Total eclipses occur when the entire Moon passes into the Earth's umbra.

The `kind` field thus holds `EclipseKind.Penumbral`, `EclipseKind.Partial`,
or `EclipseKind.Total`, depending on the kind of lunar eclipse found.

The `obscuration` field holds a value in the range [0, 1] that indicates what fraction
of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak.
This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration
is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses,
the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.

Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.

Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
of the eclipse, which is half of the amount of time the eclipse spends in each
phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
By converting from minutes to days, and subtracting/adding with `peak`, the caller
may determine the date and time of the beginning/end of each eclipse phase.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`EclipseKind`](#EclipseKind) | `kind` | The type of lunar eclipse found. |
| `double` | `obscuration` | The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra. |
| [`AstroTime`](#AstroTime) | `peak` | The time of the eclipse at its peak. |
| `double` | `sd_penum` | The semi-duration of the penumbral phase in minutes. |
| `double` | `sd_partial` | The semi-duration of the partial phase in minutes, or 0.0 if none. |
| `double` | `sd_total` | The semi-duration of the total phase in minutes, or 0.0 if none. |

---

<a name="MoonQuarterInfo"></a>
## `struct MoonQuarterInfo`

**A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `int` | `quarter` | 0=new moon, 1=first quarter, 2=full moon, 3=third quarter. |
| [`AstroTime`](#AstroTime) | `time` | The date and time of the lunar quarter. |

---

<a name="NodeEventInfo"></a>
## `struct NodeEventInfo`

**Information about an ascending or descending node of a body.**

This structure is returned by [`Astronomy.SearchMoonNode`](#Astronomy.SearchMoonNode) and [`Astronomy.NextMoonNode`](#Astronomy.NextMoonNode)
to report information about the center of the Moon passing through the ecliptic plane.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The time when the body passes through the ecliptic plane. |
| [`NodeEventKind`](#NodeEventKind) | `kind` | Whether the node is ascending (south to north) or descending (north to south). |

---

<a name="NodeEventKind"></a>
## `enum NodeEventKind`

**Indicates whether a crossing through the ecliptic plane is ascending or descending.**

| Value | Description |
| --- | --- |
| `Invalid` | Placeholder value for a missing or invalid node. |
| `Ascending` | The body passes through the ecliptic plane from south to north. |
| `Descending` | The body passes through the ecliptic plane from north to south. |

---

<a name="Observer"></a>
## `struct Observer`

**The location of an observer on (or near) the surface of the Earth.**

This structure is passed to functions that calculate phenomena as observed
from a particular place on the Earth.

### constructors

### new Observer(latitude, longitude, height)

**Creates an Observer object.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `latitude` | Geographic latitude in degrees north (positive) or south (negative) of the equator. |
| `double` | `longitude` | Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. |
| `double` | `height` | The height above (positive) or below (negative) sea level, expressed in meters. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `latitude` | Geographic latitude in degrees north (positive) or south (negative) of the equator. |
| `double` | `longitude` | Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. |
| `double` | `height` | The height above (positive) or below (negative) sea level, expressed in meters. |

### member functions

<a name="Observer.ToString"></a>
### Observer.ToString() &#8658; `string`

**Converts an `Observer` to a string representation like `(N 26.728965, W 093.157562, 1234.567 m)`.**

---

<a name="Refraction"></a>
## `enum Refraction`

**Selects whether to correct for atmospheric refraction, and if so, how.**

| Value | Description |
| --- | --- |
| `None` | No atmospheric refraction correction (airless). |
| `Normal` | Recommended correction for standard atmospheric refraction. |
| `JplHor` | Used only for compatibility testing with JPL Horizons online tool. |

---

<a name="RotationMatrix"></a>
## `struct RotationMatrix`

**A rotation matrix that can be used to transform one coordinate system to another.**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double[3,3]` | `rot` | A normalized 3x3 rotation matrix. |

---

<a name="SearchContext"></a>
## `class SearchContext`

**Represents a function whose ascending root is to be found.
See [`Astronomy.Search`](#Astronomy.Search).**

### member functions

<a name="SearchContext.Eval"></a>
### SearchContext.Eval(time) &#8658; `double`

**Evaluates the function at a given time**

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `time` | The time at which to evaluate the function. |

**Returns:** The floating point value of the function at the specified time.

---

<a name="SeasonsInfo"></a>
## `struct SeasonsInfo`

**The dates and times of changes of season for a given calendar year.
Call [`Astronomy.Seasons`](#Astronomy.Seasons) to calculate this data structure for a given year.**

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `mar_equinox` | The date and time of the March equinox for the specified year. |
| [`AstroTime`](#AstroTime) | `jun_solstice` | The date and time of the June soltice for the specified year. |
| [`AstroTime`](#AstroTime) | `sep_equinox` | The date and time of the September equinox for the specified year. |
| [`AstroTime`](#AstroTime) | `dec_solstice` | The date and time of the December solstice for the specified year. |

---

<a name="Spherical"></a>
## `struct Spherical`

**Spherical coordinates: latitude, longitude, distance.**

### constructors

### new Spherical(lat, lon, dist)

**Creates a set of spherical coordinates.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `lat` | The latitude angle: -90..+90 degrees. |
| `double` | `lon` | The longitude angle: 0..360 degrees. |
| `double` | `dist` | Distance in AU. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `lat` | The latitude angle: -90..+90 degrees. |
| `double` | `lon` | The longitude angle: 0..360 degrees. |
| `double` | `dist` | Distance in AU. |

---

<a name="StateVector"></a>
## `struct StateVector`

**A combination of a position vector and a velocity vector at a given moment in time.**

A state vector represents the dynamic state of a point at a given moment.
It includes the position vector of the point, expressed in Astronomical Units (AU)
along with the velocity vector of the point, expressed in AU/day.

### constructors

### new StateVector(x, y, z, vx, vy, vz, t)

**Creates an AstroVector.**

| Type | Parameter | Description |
| --- | --- | --- |
| `double` | `x` | A position x-coordinate expressed in AU. |
| `double` | `y` | A position y-coordinate expressed in AU. |
| `double` | `z` | A position z-coordinate expressed in AU. |
| `double` | `vx` | A velocity x-component expressed in AU/day. |
| `double` | `vy` | A velocity y-component expressed in AU/day. |
| `double` | `vz` | A velocity z-component expressed in AU/day. |
| [`AstroTime`](#AstroTime) | `t` | The date and time at which this state vector is valid. |

### new StateVector(pos, vel, time)

**Combines a position vector and a velocity vector into a single state vector.**

| Type | Parameter | Description |
| --- | --- | --- |
| [`AstroVector`](#AstroVector) | `pos` | A position vector. |
| [`AstroVector`](#AstroVector) | `vel` | A velocity vector. |
| [`AstroTime`](#AstroTime) | `time` | The common time that represents the given position and velocity. |


### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `x` | The position x-coordinate in AU. |
| `double` | `y` | The position y-coordinate in AU. |
| `double` | `z` | The position z-coordinate in AU. |
| `double` | `vx` | The velocity x-component in AU/day. |
| `double` | `vy` | The velocity y-component in AU/day. |
| `double` | `vz` | The velocity z-component in AU/day. |
| [`AstroTime`](#AstroTime) | `t` | The date and time at which this vector is valid. |

### member functions

<a name="StateVector.Position"></a>
### StateVector.Position() &#8658; [`AstroVector`](#AstroVector)

**Returns the position vector associated with this state vector.**

<a name="StateVector.ToString"></a>
### StateVector.ToString() &#8658; `string`

**Converts the state vector to a string of the format (x, y, z, vx, vy, vz, t).**

<a name="StateVector.Velocity"></a>
### StateVector.Velocity() &#8658; [`AstroVector`](#AstroVector)

**Returns the velocity vector associated with this state vector.**

---

<a name="Topocentric"></a>
## `struct Topocentric`

**Coordinates of a celestial body as seen by a topocentric observer.**

Contains horizontal and equatorial coordinates seen by an observer on or near
the surface of the Earth (a topocentric observer).
Optionally corrected for atmospheric refraction.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| `double` | `azimuth` | Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. |
| `double` | `altitude` | Angle in degrees above (positive) or below (negative) the observer's horizon. |
| `double` | `ra` | Right ascension in sidereal hours. |
| `double` | `dec` | Declination in degrees. |

---

<a name="TransitInfo"></a>
## `struct TransitInfo`

**Information about a transit of Mercury or Venus, as seen from the Earth.**

Returned by [`Astronomy.SearchTransit`](#Astronomy.SearchTransit) or [`Astronomy.NextTransit`](#Astronomy.NextTransit) to report
information about a transit of Mercury or Venus.
A transit is when Mercury or Venus passes between the Sun and Earth so that
the other planet is seen in silhouette against the Sun.

The `start` field reports the moment in time when the planet first becomes
visible against the Sun in its background.
The `peak` field reports when the planet is most aligned with the Sun,
as seen from the Earth.
The `finish` field reports the last moment when the planet is visible
against the Sun in its background.

The calculations are performed from the point of view of a geocentric observer.

### member variables

| Type | Name | Description |
| --- | --- | --- |
| [`AstroTime`](#AstroTime) | `start` | Date and time at the beginning of the transit. |
| [`AstroTime`](#AstroTime) | `peak` | Date and time of the peak of the transit. |
| [`AstroTime`](#AstroTime) | `finish` | Date and time at the end of the transit. |
| `double` | `separation` | Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`. |

---

<a name="Visibility"></a>
## `enum Visibility`

**Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.**

| Value | Description |
| --- | --- |
| `Morning` | The body is best visible in the morning, before sunrise. |
| `Evening` | The body is best visible in the evening, after sunset. |

---

