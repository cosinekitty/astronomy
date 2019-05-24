/*
    Astronomy library for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef __ASTRONOMY_H
#define __ASTRONOMY_H

#ifdef __cplusplus
extern "C" {
#endif

/*---------- types ----------*/

typedef enum
{
    ASTRO_SUCCESS,
    ASTRO_NOT_INITIALIZED,
    ASTRO_INVALID_BODY,
    ASTRO_NO_CONVERGE,
    ASTRO_BAD_TIME,
    ASTRO_BAD_VECTOR,
    ASTRO_SEARCH_FAILURE,
    ASTRO_EARTH_NOT_ALLOWED,
    ASTRO_NO_MOON_QUARTER,
    ASTRO_WRONG_MOON_QUARTER,
    ASTRO_INTERNAL_ERROR,
    ASTRO_INVALID_PARAMETER
}
astro_status_t;

typedef struct
{
    double ut;
    double tt;
}
astro_time_t;

typedef struct
{
    astro_status_t status;
    double x;
    double y;
    double z;
    astro_time_t t;
}
astro_vector_t;

typedef struct
{
    astro_status_t status;
    double angle;
}
astro_angle_result_t;

typedef enum
{
    BODY_INVALID = -1,
    BODY_MERCURY,
    BODY_VENUS,
    BODY_EARTH,
    BODY_MARS,
    BODY_JUPITER,
    BODY_SATURN,
    BODY_URANUS,
    BODY_NEPTUNE,
    BODY_PLUTO,
    BODY_SUN,
    BODY_MOON
}
astro_body_t;

#define MIN_BODY    BODY_MERCURY
#define MAX_BODY    BODY_MOON

typedef struct
{
    double latitude;
    double longitude;
    double height;
}
astro_observer_t;

typedef struct
{
    astro_status_t status;
    double ra;
    double dec;
    double dist;
}
astro_equatorial_t;

typedef struct
{
    astro_status_t status;
    double ex;
    double ey;
    double ez;
    double elat;
    double elon;
}
astro_ecliptic_t;

typedef struct
{
    double azimuth;
    double altitude;
    double ra;
    double dec;
}
astro_horizon_t;

typedef enum
{
    REFRACTION_NONE,
    REFRACTION_NORMAL,
    REFRACTION_JPLHOR
}
astro_refraction_t;

typedef struct
{
    astro_status_t  status;
    astro_time_t    time;
    int             iter;
}
astro_search_result_t;

typedef struct
{
    astro_status_t  status;
    astro_time_t    mar_equinox;
    astro_time_t    jun_solstice;
    astro_time_t    sep_equinox;
    astro_time_t    dec_solstice;
}
astro_seasons_t;

typedef struct
{
    astro_status_t  status;
    int             quarter;
    astro_time_t    time;
}
astro_moon_quarter_t;

typedef struct
{
    astro_status_t status;
    double value;    
}
astro_func_result_t;

typedef astro_func_result_t (* astro_search_func_t) (void *context, astro_time_t time);

typedef enum
{
    VISIBLE_MORNING,
    VISIBLE_EVENING
}
astro_visibility_t;

typedef struct
{
    astro_status_t      status;
    astro_time_t        time;
    astro_visibility_t  visibility;
    double              elongation;
    double              relative_longitude;
}
astro_elongation_t;

typedef struct
{
    astro_status_t      status;
    astro_time_t        time;
    astro_horizon_t     hor;
    int                 iter;
}
astro_hour_angle_t;

typedef struct
{
    astro_status_t      status;
    astro_time_t        time;
    double              mag;
    double              phase_angle;
    double              helio_dist;
    double              ring_tilt;
}
astro_illum_t;

typedef enum
{
    APSIS_PERICENTER,
    APSIS_APOCENTER,
    APSIS_INVALID
}
astro_apsis_kind_t;

typedef struct
{
    astro_status_t      status;
    astro_time_t        time;
    astro_apsis_kind_t  kind;
    double              dist_au;
    double              dist_km;
}
astro_apsis_t;

/*---------- functions ----------*/

double Astronomy_VectorLength(astro_vector_t vector);
const char *Astronomy_BodyName(astro_body_t body);
astro_body_t Astronomy_BodyCode(const char *name);
astro_observer_t Astronomy_MakeObserver(double latitude, double longitude, double height);
astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second);
astro_time_t Astronomy_AddDays(astro_time_t time, double days);
astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time);
astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time, int correct_aberration);
astro_vector_t Astronomy_GeoMoon(astro_time_t time);

astro_equatorial_t Astronomy_Equator(
    astro_body_t body, 
    astro_time_t time, 
    astro_observer_t observer,
    int ofdate,
    int aberration
);

astro_ecliptic_t Astronomy_SunPosition(astro_time_t time);
astro_ecliptic_t Astronomy_Ecliptic(astro_vector_t equ);
astro_angle_result_t Astronomy_EclipticLongitude(astro_body_t body, astro_time_t time);

astro_horizon_t Astronomy_Horizon(
    astro_time_t time, 
    astro_observer_t observer, 
    double ra, 
    double dec, 
    astro_refraction_t refraction);

astro_angle_result_t Astronomy_AngleFromSun(astro_body_t body, astro_time_t time);
astro_elongation_t Astronomy_Elongation(astro_body_t body, astro_time_t time);
astro_elongation_t Astronomy_SearchMaxElongation(astro_body_t body, astro_time_t startDate);
astro_angle_result_t Astronomy_LongitudeFromSun(astro_body_t body, astro_time_t time);
astro_search_result_t Astronomy_SearchRelativeLongitude(astro_body_t body, double targetRelLon, astro_time_t startDate);
astro_angle_result_t Astronomy_MoonPhase(astro_time_t time);
astro_search_result_t Astronomy_SearchMoonPhase(double targetLon, astro_time_t dateStart, double limitDays);
astro_moon_quarter_t Astronomy_SearchMoonQuarter(astro_time_t dateStart);
astro_moon_quarter_t Astronomy_NextMoonQuarter(astro_moon_quarter_t mq);

astro_search_result_t Astronomy_Search(
    astro_search_func_t func,
    void *context,
    astro_time_t t1,
    astro_time_t t2,
    double dt_tolerance_seconds);

astro_search_result_t Astronomy_SearchSunLongitude(
    double targetLon, 
    astro_time_t dateStart,
    double limitDays);

astro_hour_angle_t Astronomy_SearchHourAngle(
    astro_body_t body,
    astro_observer_t observer,
    double hourAngle,
    astro_time_t dateStart);

astro_search_result_t Astronomy_SearchRiseSet(
    astro_body_t body,
    astro_observer_t observer,
    int direction,
    astro_time_t dateStart,
    double limitDays);

astro_seasons_t Astronomy_Seasons(int calendar_year);
astro_illum_t Astronomy_Illumination(astro_body_t body, astro_time_t time);
astro_illum_t Astronomy_SearchPeakMagnitude(astro_body_t body, astro_time_t startDate);
astro_apsis_t Astronomy_SearchLunarApsis(astro_time_t startTime);

#ifdef __cplusplus
}
#endif

#endif  /* ifndef __ASTRONOMY_H */
