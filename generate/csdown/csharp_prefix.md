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

