# Astronomy Engine (C/C++)

This is the complete programming reference for the C version of 
[Astronomy Engine](../../). It can be used directly from C++ programs also.
Other programming languages are supported. See the [home page](../../) for more info.

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

