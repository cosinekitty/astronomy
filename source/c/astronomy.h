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

/**
 * @brief Indicates success/failure of an Astronomy Engine function call.
 */
typedef enum
{
    ASTRO_SUCCESS,                  /**< The operation was successful. */
    ASTRO_NOT_INITIALIZED,          /**< A placeholder that can be used for data that is not yet initialized. */
    ASTRO_INVALID_BODY,             /**< The celestial body was not valid. Different sets of bodies are supported depending on the function. */
    ASTRO_NO_CONVERGE,              /**< A numeric solver failed to converge. This should not happen unless there is a bug in Astronomy Engine. */
    ASTRO_BAD_TIME,                 /**< Cannot calculate Pluto's position outside the year range 1700..2200. */
    ASTRO_BAD_VECTOR,               /**< Vector magnitude is too small to be normalized into a unit vector. */
    ASTRO_SEARCH_FAILURE,           /**< Search was not able to find an ascending root crossing of the function in the specified time interval. */
    ASTRO_EARTH_NOT_ALLOWED,        /**< The Earth cannot be treated as a celestial body seen from an observer on the Earth itself. */
    ASTRO_NO_MOON_QUARTER,          /**< No lunar quarter occurs inside the specified time range. */
    ASTRO_WRONG_MOON_QUARTER,       /**< Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter. */
    ASTRO_INTERNAL_ERROR,           /**< A self-check failed inside the code somewhere, indicating a bug needs to be fixed. */
    ASTRO_INVALID_PARAMETER         /**< A parameter value passed to a function was not valid. */
}
astro_status_t;

/**
 * @brief A date and time used for astronomical calculations.
 */
typedef struct
{
    double ut;      /**< UT1/UTC number of days since noon on January 1, 2000 */
    double tt;      /**< Terrestrial Time days since noon on January 1, 2000 */
}
astro_time_t;

/**
 * @brief A calendar date and time expressed in UTC.
 */
typedef struct
{
    int     year;       /**< The year value, e.g. 2019. */
    int     month;      /**< The month value: 1=January, 2=February, ..., 12=December. */
    int     day;        /**< The day of the month in the range 1..31. */
    int     hour;       /**< The hour of the day in the range 0..23. */
    int     minute;     /**< The minute of the hour in the range 0..59. */
    double  second;     /**< The floating point number of seconds in the range [0,60). */
}
astro_utc_t;

/**
 * @brief A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
 */
typedef struct
{
    astro_status_t status;  /**< ASTRO_SUCCESS if this struct is valid; otherwise an error code. */
    double x;               /**< The Cartesian x-coordinate of the vector in AU. */
    double y;               /**< The Cartesian y-coordinate of the vector in AU. */
    double z;               /**< The Cartesian z-coordinate of the vector in AU. */
    astro_time_t t;         /**< The date and time at which this vector is valid. */
}
astro_vector_t;

/**
 * @brief An angular value expressed in degrees.
 */
typedef struct
{
    astro_status_t status;  /**< ASTRO_SUCCESS if this struct is valid; otherwise an error code. */
    double angle;           /**< An angle expressed in degrees. */
}
astro_angle_result_t;

/**
 * @brief A celestial body.
 */
typedef enum
{
    BODY_INVALID = -1,      /**< An invalid or undefined celestial body. */
    BODY_MERCURY,           /**< Mercury */
    BODY_VENUS,             /**< Venus */
    BODY_EARTH,             /**< Earth */
    BODY_MARS,              /**< Mars */
    BODY_JUPITER,           /**< Jupiter */
    BODY_SATURN,            /**< Saturn */
    BODY_URANUS,            /**< Uranus */
    BODY_NEPTUNE,           /**< Neptune */
    BODY_PLUTO,             /**< Pluto */
    BODY_SUN,               /**< Sun */
    BODY_MOON               /**< Moon */
}
astro_body_t;

#define MIN_BODY    BODY_MERCURY    /**< Minimum valid `astro_body_t` value; useful for iteration. */
#define MAX_BODY    BODY_MOON       /**< Maximum valid astro_body_t value; useful for iteration. */

/**
 * @brief The location of an observer on (or near) the surface of the Earth.
 * 
 * This structure is passed to functions that calculate phenomena as observed
 * from a particular place on the Earth.
 */
typedef struct
{
    double latitude;        /**< Geographic latitude in degrees north (positive) or south (negative) of the equator. */
    double longitude;       /**< Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. */
    double height;          /**< The height above (positive) or below (negative) sea level, expressed in meters. */
}
astro_observer_t;

/**
 * @brief Equatorial angular coordinates.
 * 
 * Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context),
 * oriented with respect to the projection of the Earth's equator onto the sky.
 */
typedef struct
{
    astro_status_t status;  /**< ASTRO_SUCCESS if this struct is valid; otherwise an error code. */
    double ra;              /**< right ascension in sidereal hours. */
    double dec;             /**< declination in degrees */
    double dist;            /**< distance to the celestial body in AU. */
}
astro_equatorial_t;

/**
 * @brief Ecliptic angular and Cartesian coordinates.
 * 
 * Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
 * oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).
 */
typedef struct
{
    astro_status_t status;  /**< ASTRO_SUCCESS if this struct is valid; otherwise an error code. */
    double ex;              /**< Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane. */
    double ey;              /**< Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox. */
    double ez;              /**< Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north. */
    double elat;            /**< Latitude in degrees north (positive) or south (negative) of the ecliptic plane. */
    double elon;            /**< Longitude in degrees around the ecliptic plane prograde from the equinox. */
}
astro_ecliptic_t;

/**
 * @brief Coordinates of a celestial body as seen by a topocentric observer.
 * 
 * Contains horizontal and equatorial coordinates seen by an observer on or near
 * the surface of the Earth (a topocentric observer).
 * Optionally corrected for atmospheric refraction.
 */
typedef struct
{
    double azimuth;     /**< Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. */
    double altitude;    /**< Angle in degrees above (positive) or below (negative) the observer's horizon. */
    double ra;          /**< Right ascension in sidereal hours. */
    double dec;         /**< Declination in degrees. */
}
astro_horizon_t;

/**
 * @brief Selects whether to correct for atmospheric refraction, and if so, how.
 */
typedef enum
{
    REFRACTION_NONE,    /**< No atmospheric refraction corection (airless). */
    REFRACTION_NORMAL,  /**< Recommended correction for standard atmospheric refraction. */
    REFRACTION_JPLHOR   /**< Used only for compatibility testing with JPL Horizons online tool. */
}
astro_refraction_t;

/**
 * @brief The result of a search for an astronomical event.
 */
typedef struct
{
    astro_status_t  status;     /**< ASTRO_SUCCESS if this struct is valid; otherwise an error code. */
    astro_time_t    time;       /**< The time at which a searched-for event occurs. */
    int             iter;       /**< The number of iterations required to numerically solve the search. */
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
astro_time_t Astronomy_CurrentTime(void);
astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second);
astro_time_t Astronomy_TimeFromUtc(astro_utc_t utc);
astro_utc_t  Astronomy_UtcFromTime(astro_time_t time);
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
astro_apsis_t Astronomy_NextLunarApsis(astro_apsis_t apsis);

#ifdef __cplusplus
}
#endif

#endif  /* ifndef __ASTRONOMY_H */
