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

| [HelioVector](#Astronomy_HelioVector) | Calculates vector with respect to the center of the Sun. |
| [GeoVector](#Astronomy_GeoVector)     | Calculates vector with respect to the center of the Earth. |
| [Equator](#Astronomy_Equator)         | Calculates right ascension and declination. |
| [Ecliptic](#Astronomy_Ecliptic)       | Calculates ecliptic latitude, longitude, and Cartesian coordinates. |
| [Horizon](#Astronomy_Horizon)         | Calculates horizontal coordinates (azimuth, altitude) for a given observer on the Earth. |

### Rise, set, and culmination times

| [SearchRiseSet](#Astronomy_SearchRiseSet) | Finds time of rise or set for a body as seen by an observer on the Earth. |
| [SearchHourAngle](#Astronomy_SearchHourAngle) | Finds when body reaches a given hour angle for an observer on the Earth. Hour angle = 0 finds culmination, the highest point in the sky. |

### Moon phases

| [MoonPhase](#Astronomy_MoonPhase) | Determines the Moon's phase expressed as an ecliptic longitude. |
| [SearchMoonQuarter](#Astronomy_SearchMoonQuarter) | Find the first quarter moon phase after a given date and time. |
| [NextMoonQuarter](#Astronomy_NextMoonQuarter) | Find the next quarter moon phase after a previous one that has been found. |

### Lunar perigee and apogee

| [SearchLunarApsis](#Astronomy_SearchLunarApsis) | Finds the next perigee or apogee of the Moon after a specified date. |
| [NextLunarApsis](#Astronomy_NextLunarApsis) | Given an already-found apsis, find the next perigee or apogee of the Moon. |

### Visual magnitude and elongation

| [Illumination](#Astronomy_Illumination) | Calculates visual magnitude and phase angle of bodies as seen from the Earth. |
| [SearchPeakMagnitude](#Astronomy_SearchPeakMagnitude) | Searches for the date and time Venus will next appear brightest as seen from the Earth. |
| [AngleFromSun](#Astronomy_AngleFromSun) | Returns full angle seen from Earth between body and Sun. |
| [Elongation](#Astronomy_Elongation) | Calculates ecliptic longitude angle between a body and the Sun, as seen from the Earth. |
| [SearchMaxElongation](#Astronomy_SearchMaxElongation) | Searches for the next maximum elongation event for Mercury or Venus that occurs after the given date. |

### Oppositions and conjunctions

| [SearchRelativeLongitude](#Astronomy_SearchRelativeLongitude) | Find oppositions and conjunctions of planets. |

### Equinoxes and solstices

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


---

<a name="astro_apsis_t"></a>


---

<a name="astro_ecliptic_t"></a>


---

<a name="astro_elongation_t"></a>


---

<a name="astro_equatorial_t"></a>


---

<a name="astro_func_result_t"></a>


---

<a name="astro_horizon_t"></a>


---

<a name="astro_hour_angle_t"></a>


---

<a name="astro_illum_t"></a>


---

<a name="astro_moon_quarter_t"></a>


---

<a name="astro_observer_t"></a>


---

<a name="astro_search_result_t"></a>


---

<a name="astro_seasons_t"></a>


---

<a name="astro_time_t"></a>


---

<a name="astro_utc_t"></a>


---

<a name="astro_vector_t"></a>
## Type Definitions



---

<a name="astro_search_func_t"></a>
