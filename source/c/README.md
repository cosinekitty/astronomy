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
`Astronomy_AddDays`



---

<a name="Astronomy_AngleFromSun"></a>
`Astronomy_AngleFromSun`



---

<a name="Astronomy_BodyCode"></a>
`Astronomy_BodyCode`



---

<a name="Astronomy_BodyName"></a>
`Astronomy_BodyName`



---

<a name="Astronomy_CurrentTime"></a>
`Astronomy_CurrentTime`



---

<a name="Astronomy_Ecliptic"></a>
`Astronomy_Ecliptic`



---

<a name="Astronomy_EclipticLongitude"></a>
`Astronomy_EclipticLongitude`



---

<a name="Astronomy_Elongation"></a>
`Astronomy_Elongation`



---

<a name="Astronomy_Equator"></a>
`Astronomy_Equator`



---

<a name="Astronomy_GeoMoon"></a>
`Astronomy_GeoMoon`



---

<a name="Astronomy_GeoVector"></a>
`Astronomy_GeoVector`



---

<a name="Astronomy_HelioVector"></a>
`Astronomy_HelioVector`



---

<a name="Astronomy_Horizon"></a>
`Astronomy_Horizon`



---

<a name="Astronomy_Illumination"></a>
`Astronomy_Illumination`



---

<a name="Astronomy_LongitudeFromSun"></a>
`Astronomy_LongitudeFromSun`



---

<a name="Astronomy_MakeObserver"></a>
`Astronomy_MakeObserver`



---

<a name="Astronomy_MakeTime"></a>
`Astronomy_MakeTime`



---

<a name="Astronomy_MoonPhase"></a>
`Astronomy_MoonPhase`



---

<a name="Astronomy_NextLunarApsis"></a>
`Astronomy_NextLunarApsis`



---

<a name="Astronomy_NextMoonQuarter"></a>
`Astronomy_NextMoonQuarter`



---

<a name="Astronomy_Search"></a>
`Astronomy_Search`



---

<a name="Astronomy_SearchHourAngle"></a>
`Astronomy_SearchHourAngle`



---

<a name="Astronomy_SearchLunarApsis"></a>
`Astronomy_SearchLunarApsis`



---

<a name="Astronomy_SearchMaxElongation"></a>
`Astronomy_SearchMaxElongation`



---

<a name="Astronomy_SearchMoonPhase"></a>
`Astronomy_SearchMoonPhase`



---

<a name="Astronomy_SearchMoonQuarter"></a>
`Astronomy_SearchMoonQuarter`



---

<a name="Astronomy_SearchPeakMagnitude"></a>
`Astronomy_SearchPeakMagnitude`



---

<a name="Astronomy_SearchRelativeLongitude"></a>
`Astronomy_SearchRelativeLongitude`



---

<a name="Astronomy_SearchRiseSet"></a>
`Astronomy_SearchRiseSet`



---

<a name="Astronomy_SearchSunLongitude"></a>
`Astronomy_SearchSunLongitude`



---

<a name="Astronomy_Seasons"></a>
`Astronomy_Seasons`



---

<a name="Astronomy_SunPosition"></a>
`Astronomy_SunPosition`



---

<a name="Astronomy_TimeFromUtc"></a>
`Astronomy_TimeFromUtc`



---

<a name="Astronomy_UtcFromTime"></a>
`Astronomy_UtcFromTime`



---

<a name="Astronomy_VectorLength"></a>
`Astronomy_VectorLength`

## Enumerated Types



---

<a name="astro_apsis_kind_t"></a>
`astro_apsis_kind_t`



---

<a name="astro_body_t"></a>
`astro_body_t`



---

<a name="astro_refraction_t"></a>
`astro_refraction_t`



---

<a name="astro_status_t"></a>
`astro_status_t`



---

<a name="astro_visibility_t"></a>
`astro_visibility_t`

## Structures



---

<a name="astro_angle_result_t"></a>
`astro_angle_result_t`



---

<a name="astro_apsis_t"></a>
`astro_apsis_t`



---

<a name="astro_ecliptic_t"></a>
`astro_ecliptic_t`



---

<a name="astro_elongation_t"></a>
`astro_elongation_t`



---

<a name="astro_equatorial_t"></a>
`astro_equatorial_t`



---

<a name="astro_func_result_t"></a>
`astro_func_result_t`



---

<a name="astro_horizon_t"></a>
`astro_horizon_t`



---

<a name="astro_hour_angle_t"></a>
`astro_hour_angle_t`



---

<a name="astro_illum_t"></a>
`astro_illum_t`



---

<a name="astro_moon_quarter_t"></a>
`astro_moon_quarter_t`



---

<a name="astro_observer_t"></a>
`astro_observer_t`



---

<a name="astro_search_result_t"></a>
`astro_search_result_t`



---

<a name="astro_seasons_t"></a>
`astro_seasons_t`



---

<a name="astro_time_t"></a>
`astro_time_t`



---

<a name="astro_utc_t"></a>
`astro_utc_t`



---

<a name="astro_vector_t"></a>
`astro_vector_t`

## Type Definitions



---

<a name="astro_search_func_t"></a>
`astro_search_func_t`

