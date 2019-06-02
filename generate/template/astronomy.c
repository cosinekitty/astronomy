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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "astronomy.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @cond DOXYGEN_SKIP */
#define PI      3.14159265358979323846
/** @endcond */

static const double T0        = 2451545.0;
static const double MJD_BASIS = 2400000.5;
static const double DEG2RAD = 0.017453292519943296;
static const double RAD2DEG = 57.295779513082321;
static const double ASEC360 = 1296000.0;
static const double ASEC2RAD = 4.848136811095359935899141e-6;
static const double PI2 = 2.0 * PI;
static const double ARC = 3600.0 * 180.0 / PI;          /* arcseconds per radian */
static const double C_AUDAY = 173.1446326846693;        /* speed of light in AU/day */
static const double ERAD = 6378136.6;                   /* mean earth radius in meters */
static const double AU = 1.4959787069098932e+11;        /* astronomical unit in meters */
static const double KM_PER_AU = 1.4959787069098932e+8;
static const double ANGVEL = 7.2921150e-5;
static const double SECONDS_PER_DAY = 24.0 * 3600.0;
static const double SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;
static const double MEAN_SYNODIC_MONTH = 29.530588;     /* average number of days for Moon to return to the same phase */
static const double EARTH_ORBITAL_PERIOD = 365.256;
static const double REFRACTION_NEAR_HORIZON = 34.0 / 60.0;   /* degrees of refractive "lift" seen for objects near horizon */
static const double SUN_RADIUS_AU  = 4.6505e-3;
static const double MOON_RADIUS_AU = 1.15717e-5;
static const double ASEC180 = 180.0 * 60.0 * 60.0;        /* arcseconds per 180 degrees (or pi radians) */

/** @cond DOXYGEN_SKIP */
#define ARRAYSIZE(x)    (sizeof(x) / sizeof(x[0]))
#define AU_PER_PARSEC   (ASEC180 / PI)             /* exact definition of how many AU = one parsec */
#define Y2000_IN_MJD    (T0 - MJD_BASIS)
/** @endcond */

static astro_time_t UniversalTime(double ut);
static astro_ecliptic_t RotateEquatorialToEcliptic(const double pos[3], double obliq_radians);
static int QuadInterp(
    double tm, double dt, double fa, double fm, double fb,
    double *x, double *t, double *df_dt);

static double LongitudeOffset(double diff)
{
    double offset = diff;

    while (offset <= -180.0)
        offset += 360.0;
    
    while (offset > 180.0)
        offset -= 360.0;

    return offset;
}

static double NormalizeLongitude(double lon)
{
    while (lon < 0.0)
        lon += 360.0;

    while (lon >= 360.0)
        lon -= 360.0;

    return lon;
}

/**
 * @brief Calculates the length of the given vector. 
 * 
 * Calculates the non-negative length of the given vector.
 * The length is expressed in the same units as the vector's components,
 * usually astronomical units (AU).
 * 
 * @param vector The vector whose length is to be calculated.
 * @return The length of the vector.
 */
double Astronomy_VectorLength(astro_vector_t vector)
{
    return sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
}

/**
 * @brief Finds the name of a celestial body.
 * @param body The celestial body whose name is to be found.
 * @return The English-language name of the celestial body, or "" if the body is not valid.
 */
const char *Astronomy_BodyName(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:  return "Mercury";
    case BODY_VENUS:    return "Venus";
    case BODY_EARTH:    return "Earth";
    case BODY_MARS:     return "Mars";
    case BODY_JUPITER:  return "Jupiter";
    case BODY_SATURN:   return "Saturn";
    case BODY_URANUS:   return "Uranus";
    case BODY_NEPTUNE:  return "Neptune";
    case BODY_PLUTO:    return "Pluto";
    case BODY_SUN:      return "Sun";
    case BODY_MOON:     return "Moon";
    default:            return "";
    }
}

/**
 * @brief Returns the #astro_body_t value corresponding to the given English name.
 * @param name One of the following strings: Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto.
 * @return If `name` is one of the strings (case-sensitive) listed above, the returned value is the corresponding #astro_body_t value, otherwise it is #BODY_INVALID.
 */
astro_body_t Astronomy_BodyCode(const char *name)
{
    if (name != NULL)
    {
        if (!strcmp(name, "Mercury"))   return BODY_MERCURY;
        if (!strcmp(name, "Venus"))     return BODY_VENUS;
        if (!strcmp(name, "Earth"))     return BODY_EARTH;
        if (!strcmp(name, "Mars"))      return BODY_MARS;
        if (!strcmp(name, "Jupiter"))   return BODY_JUPITER;
        if (!strcmp(name, "Saturn"))    return BODY_SATURN;
        if (!strcmp(name, "Uranus"))    return BODY_URANUS;
        if (!strcmp(name, "Neptune"))   return BODY_NEPTUNE;
        if (!strcmp(name, "Pluto"))     return BODY_PLUTO;
        if (!strcmp(name, "Sun"))       return BODY_SUN;
        if (!strcmp(name, "Moon"))      return BODY_MOON;
    }
    return BODY_INVALID;
}

static int IsSuperiorPlanet(astro_body_t body)
{
    switch (body)
    {
    case BODY_MARS:
    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
    case BODY_PLUTO:
        return 1;

    default:
        return 0;
    }
}

static double PlanetOrbitalPeriod(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:  return     87.969;
    case BODY_VENUS:    return    224.701;
    case BODY_EARTH:    return    EARTH_ORBITAL_PERIOD;
    case BODY_MARS:     return    686.980;
    case BODY_JUPITER:  return   4332.589;
    case BODY_SATURN:   return  10759.22;
    case BODY_URANUS:   return  30685.4;
    case BODY_NEPTUNE:  return  60189.0;
    case BODY_PLUTO:    return  90560.0;
    default:            return  0.0;        /* invalid body */
    }
}

static void FatalError(const char *message)
{
    fprintf(stderr, "FATAL: %s\n", message);
    exit(1);
}

static astro_vector_t VecError(astro_status_t status, astro_time_t time)
{
    astro_vector_t vec;
    vec.x = vec.y = vec.z = NAN;
    vec.t = time;
    vec.status = status;
    return vec;
}

static astro_equatorial_t EquError(astro_status_t status)
{
    astro_equatorial_t equ;
    equ.ra = equ.dec = equ.dist = NAN;
    equ.status = status;
    return equ;
}

static astro_ecliptic_t EclError(astro_status_t status)
{
    astro_ecliptic_t ecl;
    ecl.status = status;
    ecl.ex = ecl.ey = ecl.ez = ecl.elat = ecl.elon = NAN;
    return ecl;
}

static astro_angle_result_t AngleError(astro_status_t status)
{
    astro_angle_result_t result;
    result.status = status;
    result.angle = NAN;
    return result;
}

static astro_func_result_t FuncError(astro_status_t status)
{
    astro_func_result_t result;
    result.status = status;
    result.value = NAN;
    return result;
}

static astro_moon_quarter_t MoonQuarterError(astro_status_t status)
{
    astro_moon_quarter_t result;
    result.status = status;
    result.quarter = -1;
    result.time.tt = result.time.ut = NAN;
    return result;
}

static astro_elongation_t ElongError(astro_status_t status)
{
    astro_elongation_t result;

    result.status = status;
    result.elongation = NAN;
    result.relative_longitude = NAN;
    result.time.tt = result.time.ut = NAN;
    result.visibility = (astro_visibility_t)(-1);

    return result;
}

static astro_hour_angle_t HourAngleError(astro_status_t status)
{
    astro_hour_angle_t result;

    result.status = status;
    result.time.tt = result.time.ut = NAN;
    result.hor.altitude = result.hor.azimuth = result.hor.dec = result.hor.ra = NAN;

    return result;
}

static astro_illum_t IllumError(astro_status_t status)
{
    astro_illum_t result;

    result.status = status;
    result.time.tt = result.time.ut = NAN;
    result.mag = NAN;
    result.phase_angle = NAN;
    result.helio_dist = NAN;
    result.ring_tilt = NAN;

    return result;
}

static astro_apsis_t ApsisError(astro_status_t status)
{
    astro_apsis_t result;

    result.status = status;
    result.time.tt = result.time.ut = NAN;
    result.kind = APSIS_INVALID;
    result.dist_km = result.dist_au = NAN;

    return result;
}

static astro_search_result_t SearchError(astro_status_t status)
{
    astro_search_result_t result;
    result.time.tt = result.time.ut = NAN;
    result.status = status;
    return result;
}

static astro_func_result_t SynodicPeriod(astro_body_t body)
{
    static const double Te = 365.256;  /* Earth's orbital period in days */
    double Tp;                         /* planet's orbital period in days */
    astro_func_result_t result;

    /* The Earth does not have a synodic period as seen from itself. */
    if (body == BODY_EARTH)
        return FuncError(ASTRO_EARTH_NOT_ALLOWED);

    if (body == BODY_MOON)
    {
        result.status = ASTRO_SUCCESS;
        result.value = MEAN_SYNODIC_MONTH;
        return result;
    }

    Tp = PlanetOrbitalPeriod(body);
    if (Tp <= 0.0)
        return FuncError(ASTRO_INVALID_BODY);

    result.status = ASTRO_SUCCESS;
    result.value = fabs(Te / (Te/Tp - 1.0));
    return result;
}

static astro_angle_result_t AngleBetween(astro_vector_t a, astro_vector_t b)
{
    double r, dot;
    astro_angle_result_t result;

    r = Astronomy_VectorLength(a) * Astronomy_VectorLength(b);
    if (r < 1.0e-8)
        return AngleError(ASTRO_BAD_VECTOR);

    dot = (a.x*b.x + a.y*b.y + a.z*b.z) / r;

    if (dot <= -1.0)
        result.angle = 180.0;
    else if (dot >= +1.0)
        result.angle = 0.0;
    else
        result.angle = RAD2DEG * acos(dot);
    
    result.status = ASTRO_SUCCESS;
    return result;
}

/** @cond DOXYGEN_SKIP */
typedef struct
{
    double mjd;
    double dt;
}
deltat_entry_t;
/** @endcond */

static const deltat_entry_t DT[] = $ASTRO_DELTA_T();

/** @cond DOXYGEN_SKIP */
#define DT_LENGTH     (sizeof(DT) / sizeof(DT[0]))
/** @endcond */

static double DeltaT(double mjd)
{
    int lo, hi, c;
    double frac;

    if (mjd <= DT[0].mjd)
        return DT[0].dt;

    if (mjd >= DT[DT_LENGTH-1].mjd)
        return DT[DT_LENGTH-1].dt;

    /* Do a binary search to find the pair of indexes this mjd lies between. */

    lo = 0;
    hi = DT_LENGTH-2;   /* make sure there is always an array element after the one we are looking at. */
    for(;;)
    {
        if (lo > hi)
        {
            /* This should never happen unless there is a bug in the binary search. */
            FatalError("DeltaT: could not find delta-t value");
        }

        c = (lo + hi) / 2;
        if (mjd < DT[c].mjd)
            hi = c-1;
        else if (mjd > DT[c+1].mjd)
            lo = c+1;
        else
        {
            frac = (mjd - DT[c].mjd) / (DT[c+1].mjd - DT[c].mjd);
            return DT[c].dt + frac*(DT[c+1].dt - DT[c].dt);
        }
    }
}

static double TerrestrialTime(double ut)
{
    return ut + DeltaT(ut + Y2000_IN_MJD)/86400.0;
}

static astro_time_t UniversalTime(double ut)
{
    astro_time_t  time;
    time.ut = ut;
    time.tt = TerrestrialTime(ut);
    return time;
}

/**
 * @brief Returns the computer's current date and time in the form of an #astro_time_t.
 * 
 * Uses the computer's system clock to find the current UTC date and time with 1-second granularity.
 * Converts that date and time to an #astro_time_t value and returns the result.
 * Callers can pass this value to other Astronomy Engine functions to calculate
 * current observational conditions.
 */
astro_time_t Astronomy_CurrentTime(void)
{
    astro_time_t t;

    /* Get seconds since midnight January 1, 1970, divide to convert to days, */
    /* then subtract to get days since noon on January 1, 2000. */

    t.ut = (time(NULL) / SECONDS_PER_DAY) - 10957.5;
    t.tt = TerrestrialTime(t.ut);
    return t;
}

/**
 * @brief Creates an #astro_time_t value from a given calendar date and time.
 * 
 * Given a UTC calendar date and time, calculates an #astro_time_t value that can
 * be passed to other Astronomy Engine functions for performing various calculations
 * relating to that date and time.
 * 
 * It is the caller's responsibility to ensure that the parameter values are correct.
 * The parameters are not checked for validity, 
 * and this function never returns any indication of an error.
 * Invalid values, for example passing in February 31, may cause unexpected return values.
 * 
 * @param year      The UTC calendar year, e.g. 2019.
 * @param month     The UTC calendar month in the range 1..12.
 * @param day       The UTC calendar day in the range 1..31.
 * @param hour      The UTC hour of the day in the range 0..23.
 * @param minute    The UTC minute in the range 0..59.
 * @param second    The UTC floating-point second in the range [0, 60).
 * 
 * @return  An #astro_time_t value that represents the given calendar date and time.
 */
astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second)
{
    astro_time_t time;
    long int jd12h;
    long int y2000;

    /* This formula is adapted from NOVAS C 3.1 function julian_date() */
    jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
        + ((long) month - 14L) / 12L) / 4L
        + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
        / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
        / 100L) / 4L;    

    y2000 = jd12h - 2451545L;

    time.ut = (double)y2000 - 0.5 + (hour / 24.0) + (minute / (24.0 * 60.0)) + (second / (24.0 * 3600.0));

    time.tt = TerrestrialTime(time.ut);

    return time;
}

/**
 * @brief   Calculates the sum or difference of an #astro_time_t with a specified floating point number of days.
 * 
 * Sometimes we need to adjust a given #astro_time_t value by a certain amount of time.
 * This function adds the given real number of days in `days` to the date and time in `time`.
 * 
 * More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
 * the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time,
 * according to the historical and predictive Delta-T model provided by the
 * [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).
 * 
 * The value stored in `time` will not be modified; it is passed by value.
 * 
 * @param time  A date and time for which to calculate an adjusted date and time.
 * @param days  A floating point number of days by which to adjust `time`. May be negative, 0, or positive.
 * @return  A date and time that is conceptually equal to `time + days`.
 */
astro_time_t Astronomy_AddDays(astro_time_t time, double days)
{
    /* 
        This is slightly wrong, but the error is tiny.
        We really should be adding to TT, not to UT.
        But using TT would require creating an inverse function for DeltaT,
        which would be quite a bit of extra calculation.
        I estimate the error is in practice on the order of 10^(-7)
        times the value of 'days'.
        This is based on a typical drift of 1 second per year between UT and TT.
    */

    astro_time_t sum;

    sum.ut = time.ut + days;
    sum.tt = TerrestrialTime(sum.ut);

    return sum;   
}

/**
 * @brief   Creates an #astro_time_t value from a given calendar date and time.
 * 
 * This function is similar to #Astronomy_MakeTime, only it receives a
 * UTC calendar date and time in the form of an #astro_utc_t structure instead of
 * as separate numeric parameters.  Astronomy_TimeFromUtc is the inverse of
 * #Astronomy_UtcFromTime.
 * 
 * @param utc   The UTC calendar date and time to be converted to #astro_time_t.
 * @return  A value that can be used for astronomical calculations for the given date and time.
 */
astro_time_t Astronomy_TimeFromUtc(astro_utc_t utc)
{
    return Astronomy_MakeTime(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
}

/**
 * @brief Determines the calendar year, month, day, and time from an #astro_time_t value.
 * 
 * After calculating the date and time of an astronomical event in the form of
 * an #astro_time_t value, it is often useful to display the result in a human-readable
 * form. This function converts the linear time scales in the `ut` field of #astro_time_t
 * into a calendar date and time: year, month, day, hours, minutes, and seconds, expressed
 * in UTC.
 * 
 * @param time  The astronomical time value to be converted to calendar date and time.
 * @return  A date and time broken out into conventional year, month, day, hour, minute, and second.
 */
astro_utc_t Astronomy_UtcFromTime(astro_time_t time)
{
    /* Adapted from the NOVAS C 3.1 function cal_date() */
    astro_utc_t utc;
    long jd, k, m, n;
    double djd, x;

    djd = time.ut + 2451545.5;
    jd = (long)djd;

    x = 24.0 * fmod(djd, 1.0);
    utc.hour = (int)x;
    x = 60.0 * fmod(x, 1.0);
    utc.minute = (int)x;
    utc.second = 60.0 * fmod(x, 1.0);

    k = jd + 68569L;
    n = 4L * k / 146097L;
    k = k - (146097L * n + 3L) / 4L;
    m = 4000L * (k + 1L) / 1461001L;
    k = k - 1461L * m / 4L + 31L;

    utc.month = (int) (80L * k / 2447L);
    utc.day = (int) (k - 2447L * (long)utc.month / 80L);
    k = (long) utc.month / 11L;

    utc.month = (int) ((long)utc.month + 2L - 12L * k);
    utc.year = (int) (100L * (n - 49L) + m + k);

    return utc;
}

/**
 * @brief   Creates an observer object that represents a location on or near the surface of the Earth.
 * 
 * Some Astronomy Engine functions calculate values pertaining to an observer on the Earth.
 * These functions require a value of type #astro_observer_t that represents the location
 * of such an observer.
 * 
 * @param latitude      The geographic latitude of the observer in degrees north (positive) or south (negative) of the equator.
 * @param longitude     The geographic longitude of the observer in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.
 * @param height        The height of the observer in meters above mean sea level.
 * @return An observer object that can be passed to astronomy functions that require a geographic location.
 */
astro_observer_t Astronomy_MakeObserver(double latitude, double longitude, double height)
{
    astro_observer_t observer;

    observer.latitude = latitude;
    observer.longitude = longitude;
    observer.height = height;

    return observer;
}

static void iau2000b(astro_time_t time, double *dpsi, double *deps)
{
    /* Adapted from the NOVAS C 3.1 function of the same name. */

    static const short int nals_t[77][5] = 
    {
        { 0,    0,    0,    0,    1},
        { 0,    0,    2,   -2,    2},
        { 0,    0,    2,    0,    2},
        { 0,    0,    0,    0,    2},
        { 0,    1,    0,    0,    0},
        { 0,    1,    2,   -2,    2},
        { 1,    0,    0,    0,    0},
        { 0,    0,    2,    0,    1},
        { 1,    0,    2,    0,    2},
        { 0,   -1,    2,   -2,    2},
        { 0,    0,    2,   -2,    1},
        {-1,    0,    2,    0,    2},
        {-1,    0,    0,    2,    0},
        { 1,    0,    0,    0,    1},
        {-1,    0,    0,    0,    1},
        {-1,    0,    2,    2,    2},
        { 1,    0,    2,    0,    1},
        {-2,    0,    2,    0,    1},
        { 0,    0,    0,    2,    0},
        { 0,    0,    2,    2,    2},
        { 0,   -2,    2,   -2,    2},
        {-2,    0,    0,    2,    0},
        { 2,    0,    2,    0,    2},
        { 1,    0,    2,   -2,    2},
        {-1,    0,    2,    0,    1},
        { 2,    0,    0,    0,    0},
        { 0,    0,    2,    0,    0},
        { 0,    1,    0,    0,    1},
        {-1,    0,    0,    2,    1},
        { 0,    2,    2,   -2,    2},
        { 0,    0,   -2,    2,    0},
        { 1,    0,    0,   -2,    1},
        { 0,   -1,    0,    0,    1},
        {-1,    0,    2,    2,    1},
        { 0,    2,    0,    0,    0},
        { 1,    0,    2,    2,    2},
        {-2,    0,    2,    0,    0},
        { 0,    1,    2,    0,    2},
        { 0,    0,    2,    2,    1},
        { 0,   -1,    2,    0,    2},
        { 0,    0,    0,    2,    1},
        { 1,    0,    2,   -2,    1},
        { 2,    0,    2,   -2,    2},
        {-2,    0,    0,    2,    1},
        { 2,    0,    2,    0,    1},
        { 0,   -1,    2,   -2,    1},
        { 0,    0,    0,   -2,    1},
        {-1,   -1,    0,    2,    0},
        { 2,    0,    0,   -2,    1},
        { 1,    0,    0,    2,    0},
        { 0,    1,    2,   -2,    1},
        { 1,   -1,    0,    0,    0},
        {-2,    0,    2,    0,    2},
        { 3,    0,    2,    0,    2},
        { 0,   -1,    0,    2,    0},
        { 1,   -1,    2,    0,    2},
        { 0,    0,    0,    1,    0},
        {-1,   -1,    2,    2,    2},
        {-1,    0,    2,    0,    0},
        { 0,   -1,    2,    2,    2},
        {-2,    0,    0,    0,    1},
        { 1,    1,    2,    0,    2},
        { 2,    0,    0,    0,    1},
        {-1,    1,    0,    1,    0},
        { 1,    1,    0,    0,    0},
        { 1,    0,    2,    0,    0},
        {-1,    0,    2,   -2,    1},
        { 1,    0,    0,    0,    2},
        {-1,    0,    0,    1,    0},
        { 0,    0,    2,    1,    2},
        {-1,    0,    2,    4,    2},
        {-1,    1,    0,    1,    1},
        { 0,   -2,    2,   -2,    1},
        { 1,    0,    2,    2,    1},
        {-2,    0,    2,    2,    2},
        {-1,    0,    0,    0,    2},
        { 1,    1,    2,   -2,    2}
    };

   static const double cls_t[77][6] = 
   {
        {-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0},
        { -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0},
        {  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0},
        {   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0},
        {   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0},
        {   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0},
        {    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0},
        {   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0},
        {   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0},
        {    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0},
        {    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0},
        {    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0},
        {    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0},
        {     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0},
        {    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0},
        {    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0},
        {    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0},
        {     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0},
        {     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0},
        {    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0},
        {     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0},
        {    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0},
        {    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0},
        {     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0},
        {     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0},
        {     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0},
        {     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0},
        {    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0},
        {     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0},
        {    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0},
        {     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0},
        {    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0},
        {    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0},
        {    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0},
        {     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0},
        {     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0},
        {    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0},
        {      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0},
        {     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0},
        {     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0},
        {     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0},
        {      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0},
        {      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0},
        {     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0},
        {     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0},
        {     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0},
        {     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0},
        {      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0},
        {      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0},
        {      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0},
        {      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0},
        {      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0},
        {     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0},
        {     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0},
        {      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0},
        {     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0},
        {     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0},
        {     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0},
        {     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0},
        {     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0},
        {     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0},
        {      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0},
        {      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0},
        {      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0},
        {     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0},
        {      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0},
        {     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0},
        {     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0},
        {      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0},
        {      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0},
        {     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0},
        {      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0},
        {     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0},
        {     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0},
        {      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0},
        {      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0},
        {      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0}
    };

    double t, el, elp, f, d, om, arg, dp, de, sarg, carg;
    int i;

    t = time.tt / 36525;
    el  = fmod(485868.249036 + t * 1717915923.2178, ASEC360) * ASEC2RAD;
    elp = fmod(1287104.79305 + t * 129596581.0481,  ASEC360) * ASEC2RAD;
    f   = fmod(335779.526232 + t * 1739527262.8478, ASEC360) * ASEC2RAD;
    d   = fmod(1072260.70369 + t * 1602961601.2090, ASEC360) * ASEC2RAD;
    om  = fmod(450160.398036 - t * 6962890.5431,    ASEC360) * ASEC2RAD;
    dp = 0;
    de = 0;
    for (i=76; i >= 0; --i) 
    {
        arg = fmod((nals_t[i][0]*el + nals_t[i][1]*elp + nals_t[i][2]*f + nals_t[i][3]*d + nals_t[i][4]*om), PI2);
        sarg = sin(arg);
        carg = cos(arg);
        dp += (cls_t[i][0] + cls_t[i][1] * t) * sarg + cls_t[i][2] * carg;
        de += (cls_t[i][3] + cls_t[i][4] * t) * carg + cls_t[i][5] * sarg;
    }

    *dpsi = -0.000135 + (dp * 1.0e-7);
    *deps = +0.000388 + (de * 1.0e-7);
}

static double mean_obliq(double tt)
{
    double t = tt / 36525.0;
    double asec = 
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406;

    return asec / 3600.0;
}

/** @cond DOXYGEN_SKIP */
typedef struct 
{
    double tt;
    double dpsi;
    double deps;
    double ee;
    double mobl;
    double tobl;
}
earth_tilt_t;
/** @endcond */

static earth_tilt_t e_tilt(astro_time_t time)
{
    earth_tilt_t et;

    iau2000b(time, &et.dpsi, &et.deps);
    et.mobl = mean_obliq(time.tt);
    et.tobl = et.mobl + (et.deps / 3600.0);
    et.tt = time.tt;
    et.ee = et.dpsi * cos(et.mobl * DEG2RAD) / 15.0;

    return et;
}

static void ecl2equ_vec(astro_time_t time, const double ecl[3], double equ[3])
{
    double obl = mean_obliq(time.tt) * DEG2RAD;
    double cos_obl = cos(obl);
    double sin_obl = sin(obl);

    equ[0] = ecl[0];
    equ[1] = ecl[1]*cos_obl - ecl[2]*sin_obl;
    equ[2] = ecl[1]*sin_obl + ecl[2]*cos_obl;
}

static void precession(double tt1, const double pos1[3], double tt2, double pos2[3])
{
    double xx, yx, zx, xy, yy, zy, xz, yz, zz;
    double t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;
    double eps0 = 84381.406;

    if ((tt1 != 0.0) && (tt2 != 0.0))
        FatalError("precession: one of (tt1, tt2) must be zero.");
    
    t = (tt2 - tt1) / 36525;
    if (tt2 == 0)
        t = -t;

    psia   = (((((-    0.0000000951  * t
                 +    0.000132851 ) * t
                 -    0.00114045  ) * t
                 -    1.0790069   ) * t
                 + 5038.481507    ) * t);

    omegaa = (((((+    0.0000003337  * t
                 -    0.000000467 ) * t
                 -    0.00772503  ) * t
                 +    0.0512623   ) * t
                 -    0.025754    ) * t + eps0);

    chia   = (((((-    0.0000000560  * t
                 +    0.000170663 ) * t
                 -    0.00121197  ) * t
                 -    2.3814292   ) * t
                 +   10.556403    ) * t);

    eps0 = eps0 * ASEC2RAD;
    psia = psia * ASEC2RAD;
    omegaa = omegaa * ASEC2RAD;
    chia = chia * ASEC2RAD;

    sa = sin(eps0);
    ca = cos(eps0);
    sb = sin(-psia);
    cb = cos(-psia);
    sc = sin(-omegaa);
    cc = cos(-omegaa);
    sd = sin(chia);
    cd = cos(chia);

    xx =  cd * cb - sb * sd * cc;
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
    xy = -sd * cb - sb * cd * cc;
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
    xz =  sb * sc;
    yz = -sc * cb * ca - sa * cc;
    zz = -sc * cb * sa + cc * ca;

    if (tt2 == 0.0)
    { 
        /* Perform rotation from other epoch to J2000.0. */
        pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
        pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
        pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
    }
    else
    {
        /* Perform rotation from J2000.0 to other epoch. */
        pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
        pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
        pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
    }
}

static astro_equatorial_t vector2radec(const double pos[3])
{
    astro_equatorial_t equ;
    double xyproj;

    xyproj = pos[0]*pos[0] + pos[1]*pos[1];
    equ.dist = sqrt(xyproj + pos[2]*pos[2]);
    equ.status = ASTRO_SUCCESS;
    if (xyproj == 0.0)
    {
        if (pos[2] == 0.0)
        {
            /* Indeterminate coordinates; pos vector has zero length. */
            equ = EquError(ASTRO_BAD_VECTOR);
        }
        else if (pos[2] < 0)
        {
            equ.ra = 0.0;
            equ.dec = -90.0;
        }
        else
        {
            equ.ra = 0.0;
            equ.dec = +90.0;
        }        
    }
    else
    {
        equ.ra = atan2(pos[1], pos[0]) / (DEG2RAD * 15.0);
        if (equ.ra < 0)
            equ.ra += 24.0;

        equ.dec = RAD2DEG * atan2(pos[2], sqrt(xyproj));
    }    

    return equ;
}

static void nutation(astro_time_t time, int direction, const double inpos[3], double outpos[3])
{
    earth_tilt_t tilt = e_tilt(time);
    double oblm = tilt.mobl * DEG2RAD;
    double oblt = tilt.tobl * DEG2RAD;
    double psi = tilt.dpsi * ASEC2RAD;
    double cobm = cos(oblm);
    double sobm = sin(oblm);
    double cobt = cos(oblt);
    double sobt = sin(oblt);
    double cpsi = cos(psi);
    double spsi = sin(psi);

    double xx = cpsi;
    double yx = -spsi * cobm;
    double zx = -spsi * sobm;
    double xy = spsi * cobt;
    double yy = cpsi * cobm * cobt + sobm * sobt;
    double zy = cpsi * sobm * cobt - cobm * sobt;
    double xz = spsi * sobt;
    double yz = cpsi * cobm * sobt - sobm * cobt;
    double zz = cpsi * sobm * sobt + cobm * cobt; 

    if (direction == 0) 
    {
        /* forward rotation */
        outpos[0] = xx * inpos[0] + yx * inpos[1] + zx * inpos[2];
        outpos[1] = xy * inpos[0] + yy * inpos[1] + zy * inpos[2];
        outpos[2] = xz * inpos[0] + yz * inpos[1] + zz * inpos[2];
    }
    else
    {
        /* inverse rotation */
        outpos[0] = xx * inpos[0] + xy * inpos[1] + xz * inpos[2];
        outpos[1] = yx * inpos[0] + yy * inpos[1] + yz * inpos[2];
        outpos[2] = zx * inpos[0] + zy * inpos[1] + zz * inpos[2];
    }
}

static double era(astro_time_t time)        /* Earth Rotation Angle */
{
    double thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut;
    double thet3 = fmod(time.ut, 1.0);
    double theta = 360.0 * fmod(thet1 + thet3, 1.0);
    if (theta < 0.0)
        theta += 360.0;

    return theta;
}

static double sidereal_time(astro_time_t time)
{
    double t = time.tt / 36525.0;
    double eqeq = 15.0 * e_tilt(time).ee;    /* Replace with eqeq=0 to get GMST instead of GAST (if we ever need it) */
    double theta = era(time);
    double st = (eqeq + 0.014506 +
        (((( -    0.0000000368   * t
            -    0.000029956  ) * t
            -    0.00000044   ) * t
            +    1.3915817    ) * t
            + 4612.156534     ) * t);

    double gst = fmod(st/3600.0 + theta, 360.0) / 15.0;
    if (gst < 0.0) 
        gst += 24.0;

    return gst;
}

static void terra(astro_observer_t observer, double st, double pos[3], double vel[3]) 
{
    double erad_km = ERAD / 1000.0;
    double df = 1.0 - 0.003352819697896;    /* flattening of the Earth */
    double df2 = df * df;
    double phi = observer.latitude * DEG2RAD;
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double c = 1.0 / sqrt(cosphi*cosphi + df2*sinphi*sinphi);
    double s = df2 * c;
    double ht_km = observer.height / 1000.0;
    double ach = erad_km*c + ht_km;
    double ash = erad_km*s + ht_km;
    double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
    double sinst = sin(stlocl);
    double cosst = cos(stlocl);

    pos[0] = ach * cosphi * cosst / KM_PER_AU;
    pos[1] = ach * cosphi * sinst / KM_PER_AU;
    pos[2] = ash * sinphi / KM_PER_AU;

    vel[0] = -ANGVEL * ach * cosphi * sinst * 86400.0;
    vel[1] = +ANGVEL * ach * cosphi * cosst * 86400.0;
    vel[2] = 0.0;
}

static void geo_pos(astro_time_t time, astro_observer_t observer, double outpos[3])
{
    double gast, vel[3], pos1[3], pos2[3];

    gast = sidereal_time(time);
    terra(observer, gast, pos1, vel);
    nutation(time, -1, pos1, pos2);
    precession(time.tt, pos2, 0.0, outpos);
}

static void spin(double angle, const double pos1[3], double vec2[3]) 
{
    double angr = angle * DEG2RAD;
    double cosang = cos(angr);
    double sinang = sin(angr);
    double xx = cosang;
    double yx = sinang;
    double zx = 0;
    double xy = -sinang;
    double yy = cosang;
    double zy = 0;
    double xz = 0;
    double yz = 0;
    double zz = 1;

    vec2[0] = xx*pos1[0] + yx*pos1[1] + zx*pos1[2];
    vec2[1] = xy*pos1[0] + yy*pos1[1] + zy*pos1[2];
    vec2[2] = xz*pos1[0] + yz*pos1[1] + zz*pos1[2];
}

static void ter2cel(astro_time_t time, const double vec1[3], double vec2[3])
{
    double gast = sidereal_time(time);
    spin(-15.0 * gast, vec1, vec2);
}

/*------------------ CalcMoon ------------------*/

/** @cond DOXYGEN_SKIP */

#define DECLARE_PASCAL_ARRAY_1(elemtype,name,xmin,xmax) \
    elemtype name[(xmax)-(xmin)+1]

#define DECLARE_PASCAL_ARRAY_2(elemtype,name,xmin,xmax,ymin,ymax) \
    elemtype name[(xmax)-(xmin)+1][(ymax)-(ymin)+1]

#define ACCESS_PASCAL_ARRAY_1(name,xmin,x) \
    ((name)[(x)-(xmin)])

#define ACCESS_PASCAL_ARRAY_2(name,xmin,ymin,x,y) \
    ((name)[(x)-(xmin)][(y)-(ymin)])

typedef struct 
{
    double t;
    double dgam;
    double dlam, n, gam1c, sinpi;
    double l0, l, ls, f, d, s;
    double dl0, dl, dls, df, dd, ds;
    DECLARE_PASCAL_ARRAY_2(double,co,-6,6,1,4);   /* ARRAY[-6..6,1..4] OF REAL */
    DECLARE_PASCAL_ARRAY_2(double,si,-6,6,1,4);   /* ARRAY[-6..6,1..4] OF REAL */
}
MoonContext;

#define T           (ctx->t)
#define DGAM        (ctx->dgam)
#define DLAM        (ctx->dlam)
#define N           (ctx->n)
#define GAM1C       (ctx->gam1c)
#define SINPI       (ctx->sinpi)
#define L0          (ctx->l0)
#define L           (ctx->l)
#define LS          (ctx->ls)
#define F           (ctx->f)
#define D           (ctx->d)
#define S           (ctx->s)
#define DL0         (ctx->dl0)
#define DL          (ctx->dl)
#define DLS         (ctx->dls)
#define DF          (ctx->df)
#define DD          (ctx->dd)
#define DS          (ctx->ds)
#define CO(x,y)     ACCESS_PASCAL_ARRAY_2(ctx->co,-6,1,x,y)
#define SI(x,y)     ACCESS_PASCAL_ARRAY_2(ctx->si,-6,1,x,y)

static double Frac(double x)
{
    return x - floor(x);
}

static void AddThe(
    double c1, double s1, double c2, double s2,
    double *c, double *s)
{
    *c = c1*c2 - s1*s2;
    *s = s1*c2 + c1*s2;
}

static double Sine(double phi)
{
    /* sine, of phi in revolutions, not radians */
    return sin(PI2 * phi);
}

static void LongPeriodic(MoonContext *ctx)
{
    double S1 = Sine(0.19833+0.05611*T); 
    double S2 = Sine(0.27869+0.04508*T);
    double S3 = Sine(0.16827-0.36903*T); 
    double S4 = Sine(0.34734-5.37261*T);
    double S5 = Sine(0.10498-5.37899*T); 
    double S6 = Sine(0.42681-0.41855*T);
    double S7 = Sine(0.14943-5.37511*T);

    DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
    DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
    DLS =-6.40*S1                                   -1.89*S6;
    DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
    DD  = DL0-DLS;
    DGAM  = -3332E-9 * Sine(0.59734-5.37261*T)
             -539E-9 * Sine(0.35498-5.37899*T)
              -64E-9 * Sine(0.39943-5.37511*T);
}

static void Init(MoonContext *ctx)
{
    int I, J, MAX;
    double T2, ARG, FAC;

    T2 = T*T;
    DLAM = 0; 
    DS = 0; 
    GAM1C = 0; 
    SINPI = 3422.7000;
    LongPeriodic(ctx);
    L0 = PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/ARC;
    L  = PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /ARC;
    LS = PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/ARC;
    F  = PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /ARC;
    D  = PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /ARC;
    for (I=1; I<=4; ++I)
    {
        switch(I)
        {
            case 1: ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2: ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3: ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            case 4: ARG=D;  MAX=6; FAC=1.0;                       break;
        }
        CO(0,1) = 1.0; 
        CO(1,I) = cos(ARG)*FAC;
        SI(0,I) = 0.0; 
        SI(1,I) = sin(ARG)*FAC;
        for (J=2; J<=MAX; ++J)
            AddThe(CO(J-1,I), SI(J-1,I), CO(1,I), SI(1,I), &CO(J,I), &SI(J,I));

        for (J=1; J<=MAX; ++J)
        {
            CO(-J,I) =  CO(J,I);
            SI(-J,I) = -SI(J,I);
        }
    }
}

static void Term(MoonContext *ctx, int p, int q, int r, int s, double *x, double *y)
{
    int k;
    DECLARE_PASCAL_ARRAY_1(int, i, 1, 4);
    #define I(n) ACCESS_PASCAL_ARRAY_1(i,1,n)

    I(1) = p;
    I(2) = q;
    I(3) = r;
    I(4) = s;
    *x = 1.0;
    *y = 0.0;

    for (k=1; k<=4; ++k)
        if (I(k) != 0.0)
            AddThe(*x, *y, CO(I(k), k), SI(I(k), k), x, y);

    #undef I
}

static void AddSol(
    MoonContext *ctx, 
    double coeffl,
    double coeffs,
    double coeffg,
    double coeffp,
    int p,
    int q,
    int r,
    int s)
{
    double x, y;
    Term(ctx, p, q, r, s, &x, &y);
    DLAM += coeffl*y;
    DS += coeffs*y;
    GAM1C += coeffg*x;
    SINPI += coeffp*x;
}

static void Solar1(MoonContext *ctx)
{
    AddSol(ctx,    13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4);
    AddSol(ctx,     0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3);
    AddSol(ctx,  2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2);
    AddSol(ctx,  -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1);
    AddSol(ctx,     1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4);
    AddSol(ctx,   191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2);
    AddSol(ctx,    -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1);
    AddSol(ctx, 22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0);
    AddSol(ctx,    18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1);
    AddSol(ctx, -4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2);
    AddSol(ctx,    +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3);
    AddSol(ctx,   -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4);
    AddSol(ctx,    -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6);
    AddSol(ctx,    -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4);
    AddSol(ctx,   -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2);
    AddSol(ctx,    18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1);
    AddSol(ctx,  -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0);
    AddSol(ctx,     0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1);
    AddSol(ctx,  -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2);
    AddSol(ctx,    -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4);
    AddSol(ctx,     0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4);
    AddSol(ctx,    14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2);
    AddSol(ctx,    -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1);
    AddSol(ctx,   769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0);
    AddSol(ctx,    +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1);
    AddSol(ctx,  -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2);
    AddSol(ctx,    +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3);
    AddSol(ctx,   -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4);
    AddSol(ctx,    -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6);
    AddSol(ctx,    -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2);
    AddSol(ctx,    +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1);
    AddSol(ctx,  -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0);
    AddSol(ctx,  -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2);
    AddSol(ctx,     0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3);
    AddSol(ctx,    -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4);
}

static void Solar2(MoonContext *ctx)
{
    AddSol(ctx,     0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4);
    AddSol(ctx,    14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2);
    AddSol(ctx,   147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0);
    AddSol(ctx,    -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1);
    AddSol(ctx,    28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2);
    AddSol(ctx,    -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3);
    AddSol(ctx,     0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4);
    AddSol(ctx,    -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2);
    AddSol(ctx,    -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0);
    AddSol(ctx,    -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2);
    AddSol(ctx,    -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2);
    AddSol(ctx,     0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1);
    AddSol(ctx,  -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0);
    AddSol(ctx,     0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1);
    AddSol(ctx,   -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2);
    AddSol(ctx,     0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3);
    AddSol(ctx,    +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4);
    AddSol(ctx,     1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2);
    AddSol(ctx,    36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0);
    AddSol(ctx,   -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2);
    AddSol(ctx,    -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4);
    AddSol(ctx,    -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6);
    AddSol(ctx,    -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2);
    AddSol(ctx,    -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0);
    AddSol(ctx,    -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2);
    AddSol(ctx,    -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4);
    AddSol(ctx,     1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2);
    AddSol(ctx,     9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0);
    AddSol(ctx,    -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1);
    AddSol(ctx,    -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2);
    AddSol(ctx,     0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4);
    AddSol(ctx,    -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0);
    AddSol(ctx,    -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2);
    AddSol(ctx,    -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4);
    AddSol(ctx,    +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2);
    AddSol(ctx,    +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0);
    AddSol(ctx,    +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2);
    AddSol(ctx,    -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2);
    AddSol(ctx,    -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2);
    AddSol(ctx,   -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0);
    AddSol(ctx,    -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2);
    AddSol(ctx,    -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4);
    AddSol(ctx,    -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2);
    AddSol(ctx,    39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0);
    AddSol(ctx,     9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2);
    AddSol(ctx,     0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4);
}

static void Solar3(MoonContext *ctx)
{
    AddSol(ctx,     0.415,    0.10, 0.0  ,  0.0013,0, 1, 2, 0);
    AddSol(ctx,    -2.152,   -2.26, 0.0  , -0.0066,0, 1, 2,-2);
    AddSol(ctx,    -1.440,   -1.30, 0.0  , +0.0014,0, 1,-2, 2);
    AddSol(ctx,     0.384,   -0.04, 0.0  ,  0.0   ,0, 1,-2,-2);
    AddSol(ctx,    +1.938,   +3.60,-0.145, +0.0401,4, 0, 0, 0);
    AddSol(ctx,    -0.952,   -1.58,+0.052, -0.0130,4, 0, 0,-2);
    AddSol(ctx,    -0.551,   -0.94,+0.032, -0.0097,3, 1, 0, 0);
    AddSol(ctx,    -0.482,   -0.57,+0.005, -0.0045,3, 1, 0,-2);
    AddSol(ctx,     0.681,    0.96,-0.026,  0.0115,3,-1, 0, 0);
    AddSol(ctx,    -0.297,   -0.27, 0.002, -0.0009,2, 2, 0,-2);
    AddSol(ctx,     0.254,   +0.21,-0.003,  0.0   ,2,-2, 0,-2);
    AddSol(ctx,    -0.250,   -0.22, 0.004,  0.0014,1, 3, 0,-2);
    AddSol(ctx,    -3.996,    0.0 , 0.0  , +0.0004,2, 0, 2, 0);
    AddSol(ctx,     0.557,   -0.75, 0.0  , -0.0090,2, 0, 2,-2);
    AddSol(ctx,    -0.459,   -0.38, 0.0  , -0.0053,2, 0,-2, 2);
    AddSol(ctx,    -1.298,    0.74, 0.0  , +0.0004,2, 0,-2, 0);
    AddSol(ctx,     0.538,    1.14, 0.0  , -0.0141,2, 0,-2,-2);
    AddSol(ctx,     0.263,    0.02, 0.0  ,  0.0   ,1, 1, 2, 0);
    AddSol(ctx,     0.426,   +0.07, 0.0  , -0.0006,1, 1,-2,-2);
    AddSol(ctx,    -0.304,   +0.03, 0.0  , +0.0003,1,-1, 2, 0);
    AddSol(ctx,    -0.372,   -0.19, 0.0  , -0.0027,1,-1,-2, 2);
    AddSol(ctx,    +0.418,    0.0 , 0.0  ,  0.0   ,0, 0, 4, 0);
    AddSol(ctx,    -0.330,   -0.04, 0.0  ,  0.0   ,3, 0, 2, 0);
}

#define ADDN(coeffn,p,q,r,s)    ( Term(ctx, (p),(q),(r),(s),&x,&y), (N += (coeffn)*y) )

static void SolarN(MoonContext *ctx)
{
    double x, y;

    N = 0.0;
    ADDN(-526.069, 0, 0,1,-2); 
    ADDN(  -3.352, 0, 0,1,-4);
    ADDN( +44.297,+1, 0,1,-2); 
    ADDN(  -6.000,+1, 0,1,-4);
    ADDN( +20.599,-1, 0,1, 0); 
    ADDN( -30.598,-1, 0,1,-2);
    ADDN( -24.649,-2, 0,1, 0); 
    ADDN(  -2.000,-2, 0,1,-2);
    ADDN( -22.571, 0,+1,1,-2); 
    ADDN( +10.985, 0,-1,1,-2);
}

static void Planetary(MoonContext *ctx)
{
    DLAM +=
        +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
        +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
        +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
        +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
        +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
        +0.33*Sine(0.3132   +6.3368*T);
}

static void CalcMoon(
    double centuries_since_j2000,
    double *geo_eclip_lon,      /* (LAMBDA) equinox of date */
    double *geo_eclip_lat,      /* (BETA)   equinox of date */
    double *distance_au)        /* (R) */
{
    double lat_seconds;
    MoonContext context;
    MoonContext *ctx = &context;    /* goofy, but makes macros work inside this function */

    context.t = centuries_since_j2000;
    Init(ctx);
    Solar1(ctx);
    Solar2(ctx);
    Solar3(ctx);
    SolarN(ctx);
    Planetary(ctx);
    S = F + DS/ARC;

    lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*sin(S)-6.24*sin(3*S) + N;

    *geo_eclip_lon = PI2 * Frac((L0+DLAM/ARC) / PI2);
    *geo_eclip_lat = lat_seconds * (DEG2RAD / 3600.0);
    *distance_au = (ARC * (ERAD / AU)) / (0.999953253 * SINPI);
}

#undef T
#undef DGAM
#undef DLAM
#undef N
#undef GAM1C
#undef SINPI
#undef L0
#undef L
#undef LS
#undef F
#undef D
#undef S
#undef DL0
#undef DL
#undef DLS
#undef DF
#undef DD
#undef DS
#undef CO
#undef SI

/** @endcond */

/**
 * @brief Calculates the geocentric position of the Moon at a given time.
 * 
 * Given a time of observation, calculates the Moon's position as a vector.
 * The vector gives the location of the Moon's center relative to the Earth's center
 * with x-, y-, and z-components measured in astronomical units.
 * 
 * This algorithm is based on Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
 * which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
 * It is adapted from Turbo Pascal code from the book
 * [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
 * by Montenbruck and Pfleger.
 * 
 * @param time  The date and time for which to calculate the Moon's position.
 * @return The Moon's position as a vector in J2000 Cartesian equatorial coordinates.
 */
astro_vector_t Astronomy_GeoMoon(astro_time_t time)
{
    double geo_eclip_lon, geo_eclip_lat, distance_au;
    double dist_cos_lat;
    astro_vector_t vector;
    double gepos[3];
    double mpos1[3];
    double mpos2[3];

    CalcMoon(time.tt / 36525.0, &geo_eclip_lon, &geo_eclip_lat, &distance_au);

    /* Convert geocentric ecliptic spherical coordinates to Cartesian coordinates. */
    dist_cos_lat = distance_au * cos(geo_eclip_lat);
    gepos[0] = dist_cos_lat * cos(geo_eclip_lon);
    gepos[1] = dist_cos_lat * sin(geo_eclip_lon);
    gepos[2] = distance_au * sin(geo_eclip_lat);

    /* Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date. */
    ecl2equ_vec(time, gepos, mpos1);

    /* Convert from mean equinox of date to J2000. */
    precession(time.tt, mpos1, 0, mpos2);

    vector.status = ASTRO_SUCCESS;
    vector.x = mpos2[0];
    vector.y = mpos2[1];
    vector.z = mpos2[2];
    vector.t = time;
    return vector;
}

/*------------------ VSOP ------------------*/

/** @cond DOXYGEN_SKIP */
typedef struct 
{
    double amplitude;
    double phase;
    double frequency;
}
vsop_term_t;

typedef struct 
{
    int nterms;
    const vsop_term_t *term; 
}
vsop_series_t;

typedef struct
{
    int nseries;
    const vsop_series_t *series;
}
vsop_formula_t;

typedef struct 
{
    const vsop_formula_t formula[3];
}
vsop_model_t;
/** @endcond */

$ASTRO_C_VSOP(Mercury);
$ASTRO_C_VSOP(Venus);
$ASTRO_C_VSOP(Earth);
$ASTRO_C_VSOP(Mars);
$ASTRO_C_VSOP(Jupiter);
$ASTRO_C_VSOP(Saturn);
$ASTRO_C_VSOP(Uranus);
$ASTRO_C_VSOP(Neptune);

/** @cond DOXYGEN_SKIP */
#define VSOPFORMULA(x)    { ARRAYSIZE(x), x }
/** @endcond */

static const vsop_model_t vsop[] = 
{
    { { VSOPFORMULA(vsop_lat_Mercury),  VSOPFORMULA(vsop_lon_Mercury),  VSOPFORMULA(vsop_rad_Mercury) } },
    { { VSOPFORMULA(vsop_lat_Venus),    VSOPFORMULA(vsop_lon_Venus),    VSOPFORMULA(vsop_rad_Venus)   } },
    { { VSOPFORMULA(vsop_lat_Earth),    VSOPFORMULA(vsop_lon_Earth),    VSOPFORMULA(vsop_rad_Earth)   } },
    { { VSOPFORMULA(vsop_lat_Mars),     VSOPFORMULA(vsop_lon_Mars),     VSOPFORMULA(vsop_rad_Mars)    } },
    { { VSOPFORMULA(vsop_lat_Jupiter),  VSOPFORMULA(vsop_lon_Jupiter),  VSOPFORMULA(vsop_rad_Jupiter) } },
    { { VSOPFORMULA(vsop_lat_Saturn),   VSOPFORMULA(vsop_lon_Saturn),   VSOPFORMULA(vsop_rad_Saturn)  } },
    { { VSOPFORMULA(vsop_lat_Uranus),   VSOPFORMULA(vsop_lon_Uranus),   VSOPFORMULA(vsop_rad_Uranus)  } },
    { { VSOPFORMULA(vsop_lat_Neptune),  VSOPFORMULA(vsop_lon_Neptune),  VSOPFORMULA(vsop_rad_Neptune) } }
};

/** @cond DOXYGEN_SKIP */
#define CalcEarth(time)     CalcVsop(&vsop[BODY_EARTH], (time))
/** @endcond */

static astro_vector_t CalcVsop(const vsop_model_t *model, astro_time_t time)
{
    int k, s, i;
    double t = time.tt / 365250;    /* millennia since 2000 */
    double sphere[3];
    double r_coslat;
    double eclip[3];
    astro_vector_t vector;

    /* Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates. */
    for (k=0; k < 3; ++k)
    {
        double tpower = 1.0;
        const vsop_formula_t *formula = &model->formula[k];
        sphere[k] = 0.0;
        for (s=0; s < formula->nseries; ++s)
        {
            double sum = 0.0;
            const vsop_series_t *series = &formula->series[s];
            for (i=0; i < series->nterms; ++i)
            {
                const vsop_term_t *term = &series->term[i];
                sum += term->amplitude * cos(term->phase + (t * term->frequency));
            }
            sphere[k] += tpower * sum;
            tpower *= t;
        }
    }

    /* Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates. */
    r_coslat = sphere[2] * cos(sphere[1]);
    eclip[0] = r_coslat * cos(sphere[0]);
    eclip[1] = r_coslat * sin(sphere[0]);
    eclip[2] = sphere[2] * sin(sphere[1]);

    /* Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates. */
    vector.status = ASTRO_SUCCESS;
    vector.x = eclip[0] + 0.000000440360*eclip[1] - 0.000000190919*eclip[2];
    vector.y = -0.000000479966*eclip[0] + 0.917482137087*eclip[1] - 0.397776982902*eclip[2];
    vector.z = 0.397776982902*eclip[1] + 0.917482137087*eclip[2];
    vector.t = time;

    return vector;
}

/*------------------ Chebyshev model for Pluto ------------------*/

/** @cond DOXYGEN_SKIP */
typedef struct
{
    double data[3];
}
astro_cheb_coeff_t;

typedef struct
{
    double tt;
    double ndays;
    int ncoeff;
    const astro_cheb_coeff_t *coeff;
}
astro_cheb_record_t;
/** @endcond */

$ASTRO_C_CHEBYSHEV(8);

static double ChebScale(double t_min, double t_max, double t) 
{
    return (2*t - (t_max + t_min)) / (t_max - t_min);
}

static astro_vector_t CalcChebyshev(const astro_cheb_record_t model[], int nrecs, astro_time_t time)
{
    int i, d, k;
    double pos[3];
    double p0, p1, p2, sum;
    astro_vector_t vector;

    /* Search for a record that overlaps the given time value. */
    for (i=0; i < nrecs; ++i) 
    {
        double x = ChebScale(model[i].tt, model[i].tt + model[i].ndays, time.tt);
        if (-1.0 <= x && x <= +1.0)
        {
            for (d=0; d < 3; ++d) 
            {
                p0 = 1.0;
                sum = model[i].coeff[0].data[d];
                p1 = x;
                sum += model[i].coeff[1].data[d] * p1;
                for (k=2; k < model[i].ncoeff; ++k) 
                {
                    p2 = (2 * x * p1) - p0;
                    sum += model[i].coeff[k].data[d] * p2;
                    p0 = p1;
                    p1 = p2;
                }
                pos[d] = sum - model[i].coeff[0].data[d] / 2.0;
            }

            /* We found the position of the body. */
            vector.status = ASTRO_SUCCESS;
            vector.t = time;
            vector.x = pos[0];
            vector.y = pos[1];
            vector.z = pos[2];
            return vector;
        }
    }

    /* The Chebyshev model does not cover this time value. */
    return VecError(ASTRO_BAD_TIME, time);
}

/** @cond DOXYGEN_SKIP */
#define CalcPluto(time)    (CalcChebyshev(cheb_8, ARRAYSIZE(cheb_8), (time)))
/** @endcond */

/*------------------ end of generated code ------------------*/

/**
 * @brief Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.
 * 
 * This function calculates the position of the given celestial body as a vector,
 * using the center of the Sun as the origin.  The result is expressed as a Cartesian
 * vector in the J2000 equatorial system: the coordinates are based on the mean equator
 * of the Earth at noon UTC on 1 January 2000.
 * 
 * The position is not corrected for light travel time or aberration. 
 * This is different from the behavior of #Astronomy_GeoVector.
 * 
 * If given an invalid value for `body`, or the body is `BODY_PLUTO` and the `time` is outside
 * the year range 1700..2200, this function will fail. The caller should always check
 * the `status` field inside the returned #astro_vector_t for `ASTRO_SUCCESS` (success)
 * or any other value (failure) before trusting the resulting vector.
 * 
 * @param body  A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
 * @param time  The date and time for which to calculate the position.
 * @return      A heliocentric position vector of the center of the given body.
 */
astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time)
{
    astro_vector_t vector, earth;

    switch (body)
    {
    case BODY_SUN:
        vector.status = ASTRO_SUCCESS;
        vector.x = 0.0;
        vector.y = 0.0;
        vector.z = 0.0;
        vector.t = time;
        return vector;

    case BODY_MERCURY:
    case BODY_VENUS:
    case BODY_EARTH:
    case BODY_MARS:
    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
        return CalcVsop(&vsop[body], time);

    case BODY_PLUTO:
        return CalcPluto(time);
        
    case BODY_MOON:
        vector = Astronomy_GeoMoon(time);
        earth = CalcEarth(time);
        vector.x += earth.x;
        vector.y += earth.y;
        vector.z += earth.z;
        return vector;

    default:
        return VecError(ASTRO_INVALID_BODY, time);
    }
}

/**
 * @brief Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.
 * 
 * This function calculates the position of the given celestial body as a vector,
 * using the center of the Earth as the origin.  The result is expressed as a Cartesian
 * vector in the J2000 equatorial system: the coordinates are based on the mean equator
 * of the Earth at noon UTC on 1 January 2000.
 * 
 * If given an invalid value for `body`, or the body is `BODY_PLUTO` and the `time` is outside
 * the year range 1700..2200, this function will fail. The caller should always check
 * the `status` field inside the returned #astro_vector_t for `ASTRO_SUCCESS` (success)
 * or any other value (failure) before trusting the resulting vector.
 * 
 * Unlike #Astronomy_HelioVector, this function always corrects for light travel time.
 * This means the position of the body is "back-dated" by the amount of time it takes
 * light to travel from that body to an observer on the Earth.
 * 
 * Also, the position can optionally be corrected for
 * [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
 * causing the apparent direction of the body to be shifted due to transverse
 * movement of the Earth with respect to the rays of light coming from that body.
 * 
 * @param body          A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
 * @param time          The date and time for which to calculate the position.
 * @param aberration    `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected.
 * @return              A heliocentric position vector of the center of the given body.
 */
astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time, astro_aberration_t aberration)
{
    astro_vector_t vector;
    astro_vector_t earth;
    astro_time_t ltime;
    astro_time_t ltime2;
    double dt;
    int iter;

    if (aberration != ABERRATION && aberration != NO_ABERRATION)
        return VecError(ASTRO_INVALID_PARAMETER, time);

    switch (body)
    {
    case BODY_EARTH:
        /* The Earth's geocentric coordinates are always (0,0,0). */
        vector.status = ASTRO_SUCCESS;
        vector.x = 0.0;
        vector.y = 0.0;
        vector.z = 0.0;
        break;

    case BODY_SUN:
        /* The Sun's heliocentric coordinates are always (0,0,0). No need for light travel correction. */
        vector = CalcEarth(time);
        vector.x *= -1.0;
        vector.y *= -1.0;
        vector.z *= -1.0;
        break;

    case BODY_MOON:
        vector = Astronomy_GeoMoon(time);
        break;

    default:
        /* For all other bodies, apply light travel time correction. */

        if (aberration == NO_ABERRATION)
        {
            /* No aberration, so calculate Earth's position once, at the time of observation. */
            earth = CalcEarth(time);
            if (earth.status != ASTRO_SUCCESS)
                return earth;
        }

        ltime = time;
        for (iter=0; iter < 10; ++iter)
        {            
            vector = Astronomy_HelioVector(body, ltime);
            if (vector.status != ASTRO_SUCCESS)
                return vector;

            if (aberration == ABERRATION)
            {
                /* 
                    Include aberration, so make a good first-order approximation
                    by backdating the Earth's position also.
                    This is confusing, but it works for objects within the Solar System
                    because the distance the Earth moves in that small amount of light
                    travel time (a few minutes to a few hours) is well approximated
                    by a line segment that substends the angle seen from the remote
                    body viewing Earth. That angle is pretty close to the aberration
                    angle of the moving Earth viewing the remote body.
                    In other words, both of the following approximate the aberration angle:
                        (transverse distance Earth moves) / (distance to body)
                        (transverse speed of Earth) / (speed of light).
                */
                earth = CalcEarth(ltime);
                if (earth.status != ASTRO_SUCCESS)
                    return earth;
            }

            /* Convert heliocentric vector to geocentric vector. */
            vector.x -= earth.x;
            vector.y -= earth.y;
            vector.z -= earth.z;

            ltime2 = Astronomy_AddDays(time, -Astronomy_VectorLength(vector) / C_AUDAY);
            dt = fabs(ltime2.tt - ltime.tt);
            if (dt < 1.0e-9)
                goto finished;  /* Ensures we patch 'vector.t' with current time, not ante-dated time. */

            ltime = ltime2;
        }
        return VecError(ASTRO_NO_CONVERGE, time);   /* light travel time solver did not converge */
    }

finished:
    vector.t = time;
    return vector;
}

/**
 * @brief   Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.
 * 
 * Calculates topocentric equatorial coordinates in one of two different systems:
 * J2000 or true-equator-of-date, depending on the value of the `equdate` parameter.
 * Equatorial coordinates include right ascension, declination, and distance in astronomical units.
 * 
 * This function corrects for light travel time: it adjusts the apparent location
 * of the observed body based on how long it takes for light to travel from the body to the Earth.
 * 
 * This function corrects for *topocentric parallax*, meaning that it adjusts for the
 * angular shift depending on where the observer is located on the Earth. This is most
 * significant for the Moon, because it is so close to the Earth. However, parallax corection 
 * has a small effect on the apparent positions of other bodies.
 * 
 * Correction for aberration is optional, using the `aberration` parameter.
 * 
 * @param body          The celestial body to be observed. Not allowed to be `BODY_EARTH`.
 * @param time          The date and time at which the observation takes place.
 * @param observer      A location on or near the surface of the Earth.
 * @param equdate       Selects the date of the Earth's equator in which to express the equatorial coordinates.
 * @param aberration    Selects whether or not to correct for aberration.
 */
astro_equatorial_t Astronomy_Equator(
    astro_body_t body, 
    astro_time_t time, 
    astro_observer_t observer,
    astro_equator_date_t equdate,
    astro_aberration_t aberration)
{
    astro_equatorial_t equ;
    astro_vector_t gc;
    double gc_observer[3];
    double j2000[3];
    double temp[3];
    double datevect[3];

    geo_pos(time, observer, gc_observer);
    gc = Astronomy_GeoVector(body, time, aberration);
    if (gc.status != ASTRO_SUCCESS)
        return EquError(gc.status);

    j2000[0] = gc.x - gc_observer[0];
    j2000[1] = gc.y - gc_observer[1];
    j2000[2] = gc.z - gc_observer[2];

    switch (equdate)
    {
    case EQUATOR_OF_DATE:
        precession(0.0, j2000, time.tt, temp);
        nutation(time, 0, temp, datevect);
        equ = vector2radec(datevect);
        return equ;

    case EQUATOR_J2000:
        equ = vector2radec(j2000);
        return equ;

    default:
        return EquError(ASTRO_INVALID_PARAMETER);
    }
}

/**
 * @brief Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.
 * 
 * Given a date and time, the geographic location of an observer on the Earth, and
 * equatorial coordinates (right ascension and declination) of a celestial body,
 * this function returns horizontal coordinates (azimuth and altitude angles) for the body
 * relative to the horizon at the geographic location.
 * 
 * The right ascension `ra` and declination `dec` passed in must be *equator of date*
 * coordinates, based on the Earth's true equator at the date and time of the observation.
 * Otherwise the resulting horizontal coordinates will be inaccurate.
 * Equator of date coordinates can be obtained by calling #Astronomy_Equator, passing in 
 * `EQUATOR_OF_DATE` as its `equdate` parameter. It is also recommended to enable 
 * aberration correction by passing in `ABERRATION` as the `aberration` parameter.
 * 
 * This function optionally corrects for atmospheric refraction.
 * For most uses, it is recommended to pass `REFRACTION_NORMAL` in the `refraction` parameter to
 * correct for optical lensing of the Earth's atmosphere that causes objects 
 * to appear somewhat higher above the horizon than they actually are.
 * However, callers may choose to avoid this correction by passing in `REFRACTION_NONE`.
 * If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
 * in the #astro_horizon_t structure returned by this function will all be corrected for refraction.
 * If refraction is disabled, none of these four coordinates will be corrected; in that case, 
 * the right ascension and declination in the returned structure will be numerically identical
 * to the respective `ra` and `dec` values passed in.
 * 
 * @param time          
 *      The date and time of the observation.
 * 
 * @param observer      
 *      The geographic location of the observer.
 * 
 * @param ra            
 *      The right ascension of the body in sidereal hours.
 *      See remarks above for more details.
 * 
 * @param dec           
 *      The declination of the body in degrees. See remarks above for more details.
 * 
 * @param refraction    
 *      Selects whether to correct for atmospheric refraction, and if so, which model to use.
 *      The recommended value for most uses is `REFRACTION_NORMAL`.
 *      See remarks above for more details.
 * 
 * @return  
 *      The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.
 */
astro_horizon_t Astronomy_Horizon(
    astro_time_t time, astro_observer_t observer, double ra, double dec, astro_refraction_t refraction)
{
    astro_horizon_t hor;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double p[3], pz, pn, pw, proj;
    double az, zd;

    double sinlat = sin(observer.latitude * DEG2RAD);
    double coslat = cos(observer.latitude * DEG2RAD);
    double sinlon = sin(observer.longitude * DEG2RAD);
    double coslon = cos(observer.longitude * DEG2RAD);
    double sindc = sin(dec * DEG2RAD);
    double cosdc = cos(dec * DEG2RAD);
    double sinra = sin(ra * 15 * DEG2RAD);
    double cosra = cos(ra * 15 * DEG2RAD);

    uze[0] = coslat * coslon;
    uze[1] = coslat * sinlon;
    uze[2] = sinlat;

    une[0] = -sinlat * coslon;
    une[1] = -sinlat * sinlon;
    une[2] = coslat;

    uwe[0] = sinlon;
    uwe[1] = -coslon;
    uwe[2] = 0.0;

    ter2cel(time, uze, uz);
    ter2cel(time, une, un);
    ter2cel(time, uwe, uw);

    p[0] = cosdc * cosra;
    p[1] = cosdc * sinra;
    p[2] = sindc;

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];

    proj = sqrt(pn*pn + pw*pw);
    az = 0.0;
    if (proj > 0.0) {
        az = -atan2(pw, pn) * RAD2DEG;
        if (az < 0) 
            az += 360;
        if (az >= 360) 
            az -= 360;
    }
    zd = atan2(proj, pz) * RAD2DEG;
    hor.ra = ra;
    hor.dec = dec;

    if (refraction == REFRACTION_NORMAL || refraction == REFRACTION_JPLHOR) 
    {
        double zd0, refr, hd;
        int j;

        zd0 = zd;

        // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        // JPL Horizons says it uses refraction algorithm from 
        // Meeus "Astronomical Algorithms", 1991, p. 101-102.
        // I found the following Go implementation:
        // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        // This is a translation from the function "Saemundsson" there.
        // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        // This is important because the 'refr' formula below goes crazy near hd = -5.11.
        hd = 90.0 - zd;
        if (hd < -1.0)
            hd = -1.0;

        refr = (1.02 / tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

        if (refraction == REFRACTION_NORMAL && zd > 91.0) 
        {
            // In "normal" mode we gradually reduce refraction toward the nadir
            // so that we never get an altitude angle less than -90 degrees.
            // When horizon angle is -1 degrees, zd = 91, and the factor is exactly 1.
            // As zd approaches 180 (the nadir), the fraction approaches 0 linearly.
            refr *= (180.0 - zd) / 89.0;
        }

        zd -= refr;

        if (refr > 0.0 && zd > 3.0e-4) 
        {
            double sinzd = sin(zd * DEG2RAD);
            double coszd = cos(zd * DEG2RAD);
            double sinzd0 = sin(zd0 * DEG2RAD);
            double coszd0 = cos(zd0 * DEG2RAD);
            double pr[3];

            for (j=0; j<3; ++j)
                pr[j] = ((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd;

            proj = sqrt(pr[0]*pr[0] + pr[1]*pr[1]);
            if (proj > 0) 
            {
                hor.ra = atan2(pr[1], pr[0]) * RAD2DEG / 15;
                if (hor.ra < 0)
                    hor.ra += 24;
                if (hor.ra >= 24)
                    hor.ra -= 24;
            } 
            else 
            {
                hor.ra = 0;
            }
            hor.dec = atan2(pr[2], proj) * RAD2DEG;
        }
    }

    hor.azimuth = az;
    hor.altitude = 90.0 - zd;
    return hor;
}

/**
 * @brief Calculates geocentric ecliptic coordinates for the Sun.
 * 
 * This function calculates the position of the Sun as seen from the Earth.
 * The returned value includes both Cartesian and spherical coordinates.
 * The x-coordinate and longitude values in the returned structure are based
 * on the *true equinox of date*: one of two points in the sky where the instantaneous
 * plane of the Earth's equator at the given date and time (the *equatorial plane*) 
 * intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
 * By convention, the apparent location of the Sun at the March equinox is chosen
 * as the longitude origin and x-axis direction, instead of the one for September.
 * 
 * `Astronomy_SunPosition` corrects for precession and nutation of the Earth's axis
 * in order to obtain the exact equatorial plane at the given time.
 * 
 * This function can be used for calculating changes of seasons: equinoxes and solstices.
 * In fact, the function #Astronomy_Seasons does use this function for that purpose.
 * 
 * @param time  
 *      The date and time for which to calculate the Sun's position.
 * 
 * @return
 *      The ecliptic coordinates of the Sun using the Earth's true equator of date.
 */
astro_ecliptic_t Astronomy_SunPosition(astro_time_t time)
{
    astro_time_t adjusted_time;
    astro_vector_t earth2000;
    double sun2000[3];
    double stemp[3];
    double sun_ofdate[3];
    double true_obliq;

    /* Correct for light travel time from the Sun. */
    /* Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes! */
    adjusted_time = Astronomy_AddDays(time, -1.0 / C_AUDAY);

    earth2000 = CalcEarth(adjusted_time);
    if (earth2000.status != ASTRO_SUCCESS)
        return EclError(earth2000.status);

    /* Convert heliocentric location of Earth to geocentric location of Sun. */
    sun2000[0] = -earth2000.x;
    sun2000[1] = -earth2000.y;
    sun2000[2] = -earth2000.z;

    /* Convert to equatorial Cartesian coordinates of date. */
    precession(0.0, sun2000, adjusted_time.tt, stemp);
    nutation(adjusted_time, 0, stemp, sun_ofdate);

    /* Convert equatorial coordinates to ecliptic coordinates. */
    true_obliq = DEG2RAD * e_tilt(adjusted_time).tobl;
    return RotateEquatorialToEcliptic(sun_ofdate, true_obliq);
}

astro_ecliptic_t Astronomy_Ecliptic(astro_vector_t equ)
{
    /* Based on NOVAS functions equ2ecl() and equ2ecl_vec(). */
    static double ob2000;
    double pos[3];

    if (equ.status != ASTRO_SUCCESS)
        return EclError(equ.status);

    if (ob2000 == 0.0)
    {
        /* Lazy-evaluate and keep the mean obliquity of the ecliptic at J2000. */
        /* This way we don't need to crunch the numbers more than once. */
        /* This is not thread safe. */
        ob2000 = DEG2RAD * e_tilt(Astronomy_MakeTime(2000, 1, 1, 12, 0, 0.0)).mobl;
    }

    pos[0] = equ.x;
    pos[1] = equ.y;
    pos[2] = equ.z;

    return RotateEquatorialToEcliptic(pos, ob2000);
}

astro_angle_result_t Astronomy_EclipticLongitude(astro_body_t body, astro_time_t time)
{
    astro_vector_t hv;
    astro_ecliptic_t eclip;
    astro_angle_result_t result;

    if (body == BODY_SUN)
        return AngleError(ASTRO_INVALID_BODY);      /* cannot calculate heliocentric longitude of the Sun */

    hv = Astronomy_HelioVector(body, time);
    eclip = Astronomy_Ecliptic(hv);     /* checks for errors in hv, so we don't have to here */
    if (eclip.status != ASTRO_SUCCESS)
        return AngleError(eclip.status);

    result.angle = eclip.elon;
    result.status = ASTRO_SUCCESS;
    return result;
}

static astro_ecliptic_t RotateEquatorialToEcliptic(const double pos[3], double obliq_radians)
{
    astro_ecliptic_t ecl;
    double cos_ob, sin_ob;
    double xyproj;

    cos_ob = cos(obliq_radians);
    sin_ob = sin(obliq_radians);

    ecl.ex = +pos[0];
    ecl.ey = +pos[1]*cos_ob + pos[2]*sin_ob;
    ecl.ez = -pos[1]*sin_ob + pos[2]*cos_ob;

    xyproj = sqrt(ecl.ex*ecl.ex + ecl.ey*ecl.ey);
    if (xyproj > 0.0)
    {
        ecl.elon = RAD2DEG * atan2(ecl.ey, ecl.ex);
        if (ecl.elon < 0.0)
            ecl.elon += 360.0;
    }
    else
        ecl.elon = 0.0;

    ecl.elat = RAD2DEG * atan2(ecl.ez, xyproj);
    ecl.status = ASTRO_SUCCESS;
    return ecl;
}

static astro_func_result_t sun_offset(void *context, astro_time_t time)
{
    astro_func_result_t result;
    double targetLon = *((double *)context);
    astro_ecliptic_t ecl = Astronomy_SunPosition(time);
    if (ecl.status != ASTRO_SUCCESS)
        return FuncError(ecl.status);
    result.value = LongitudeOffset(ecl.elon - targetLon);
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_search_result_t Astronomy_SearchSunLongitude(
    double targetLon, 
    astro_time_t dateStart,
    double limitDays)
{
    astro_time_t t2 = Astronomy_AddDays(dateStart, limitDays);
    return Astronomy_Search(sun_offset, &targetLon, dateStart, t2, 1.0);
}

/** @cond DOXYGEN_SKIP */
#define CALLFUNC(f,t)  \
    do { \
        funcres = func(context, (t)); \
        if (funcres.status != ASTRO_SUCCESS) return SearchError(funcres.status); \
        (f) = funcres.value; \
    } while(0)
/** @endcond */

astro_search_result_t Astronomy_Search(
    astro_search_func_t func,
    void *context,
    astro_time_t t1,
    astro_time_t t2,
    double dt_tolerance_seconds)
{
    astro_search_result_t result;
    astro_time_t tmid;
    astro_time_t tq;
    astro_func_result_t funcres;
    double f1, f2, fmid, fq, dt_days, dt, dt_guess;
    double q_x, q_ut, q_df_dt;
    const int iter_limit = 20;
    int iter = 0;
    int calc_fmid = 1;

    dt_days = fabs(dt_tolerance_seconds / SECONDS_PER_DAY);
    CALLFUNC(f1, t1);
    CALLFUNC(f2, t2);

    for(;;)
    {
        if (++iter > iter_limit)
            return SearchError(ASTRO_NO_CONVERGE);

        dt = (t2.tt - t1.tt) / 2.0;
        tmid = Astronomy_AddDays(t1, dt);
        if (fabs(dt) < dt_days)
        {
            /* We are close enough to the event to stop the search. */
            result.time = tmid;
            result.status = ASTRO_SUCCESS;
            return result;
        }

        if (calc_fmid)
            CALLFUNC(fmid, tmid);
        else
            calc_fmid = 1;      /* we already have the correct value of fmid from the previous loop */

        /* Quadratic interpolation: */
        /* Try to find a parabola that passes through the 3 points we have sampled: */
        /* (t1,f1), (tmid,fmid), (t2,f2) */

        if (QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2, &q_x, &q_ut, &q_df_dt))
        {
            tq = UniversalTime(q_ut);
            CALLFUNC(fq, tq);
            if (q_df_dt != 0.0)
            {
                if (fabs(fq / q_df_dt) < dt_days)
                {
                    /* The estimated time error is small enough that we can quit now. */
                    result.time = tq;
                    result.status = ASTRO_SUCCESS;
                    return result;
                }

                /* Try guessing a tighter boundary with the interpolated root at the center. */
                dt_guess = 1.2 * fabs(fq / q_df_dt);
                if (dt_guess < dt/10.0)
                {
                    astro_time_t tleft = Astronomy_AddDays(tq, -dt_guess);
                    astro_time_t tright = Astronomy_AddDays(tq, +dt_guess);
                    if ((tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0) 
                    {
                        if ((tright.ut - t1.ut)*(tright.ut - t2.ut) < 0) 
                        {
                            double fleft, fright;
                            CALLFUNC(fleft, tleft);
                            CALLFUNC(fright, tright);
                            if (fleft<0.0 && fright>=0.0)
                            {
                                f1 = fleft;
                                f2 = fright;
                                t1 = tleft;
                                t2 = tright;
                                fmid = fq;
                                calc_fmid = 0;  /* save a little work -- no need to re-calculate fmid next time around the loop */
                                continue;
                            }
                        }
                    }
                }
            }
        }

        /* After quadratic interpolation attempt. */
        /* Now just divide the region in two parts and pick whichever one appears to contain a root. */
        if (f1 < 0.0 && fmid >= 0.0)
        {
            t2 = tmid;
            f2 = fmid;
            continue;
        }

        if (fmid < 0.0 && f2 >= 0.0)
        {
            t1 = tmid;
            f1 = fmid;
            continue;
        }

        /* Either there is no ascending zero-crossing in this range */
        /* or the search window is too wide (more than one zero-crossing). */
        return SearchError(ASTRO_SEARCH_FAILURE);
    }
}

static int QuadInterp(
    double tm, double dt, double fa, double fm, double fb,
    double *out_x, double *out_t, double *out_df_dt)
{
    double Q, R, S;
    double u, ru, x1, x2;

    Q = (fb + fa)/2.0 - fm;
    R = (fb - fa)/2.0;
    S = fm;

    if (Q == 0.0)
    {
        /* This is a line, not a parabola. */
        if (R == 0.0)
            return 0;       /* This is a HORIZONTAL line... can't make progress! */
        *out_x = -S / R;
        if (*out_x < -1.0 || *out_x > +1.0)
            return 0;   /* out of bounds */            
    }
    else
    {
        /* This really is a parabola. Find roots x1, x2. */
        u = R*R - 4*Q*S;
        if (u <= 0.0)
            return 0;   /* can't solve if imaginary, or if vertex of parabola is tangent. */

        ru = sqrt(u);
        x1 = (-R + ru) / (2.0 * Q);
        x2 = (-R - ru) / (2.0 * Q);
        if (-1.0 <= x1 && x1 <= +1.0)
        {
            if (-1.0 <= x2 && x2 <= +1.0)
                return 0;   /* two roots are within bounds; we require a unique zero-crossing. */
            *out_x = x1;
        }
        else if (-1.0 <= x2 && x2 <= +1.0)
            *out_x = x2;
        else
            return 0;   /* neither root is within bounds */
    }

    *out_t = tm + (*out_x)*dt;
    *out_df_dt = (2*Q*(*out_x) + R) / dt;
    return 1;   /* success */
}

static astro_status_t FindSeasonChange(double targetLon, int year, int month, int day, astro_time_t *time)
{
    astro_time_t startDate = Astronomy_MakeTime(year, month, day, 0, 0, 0.0);
    astro_search_result_t result = Astronomy_SearchSunLongitude(targetLon, startDate, 4.0);
    *time = result.time;
    return result.status;
}

astro_seasons_t Astronomy_Seasons(int year)
{
    astro_seasons_t seasons;
    astro_status_t  status;

    seasons.status = ASTRO_SUCCESS;

    status = FindSeasonChange(  0, year,  3, 19, &seasons.mar_equinox);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange( 90, year,  6, 19, &seasons.jun_solstice);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange(180, year,  9, 21, &seasons.sep_equinox);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange(270, year, 12, 20, &seasons.dec_solstice);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    return seasons;
}

astro_angle_result_t Astronomy_AngleFromSun(astro_body_t body, astro_time_t time)
{
    astro_vector_t sv, bv;

    sv = Astronomy_GeoVector(BODY_SUN, time, NO_ABERRATION);    /* FIXFIXFIX: use aberration or not? */
    if (sv.status != ASTRO_SUCCESS)
        return AngleError(sv.status);

    bv = Astronomy_GeoVector(body, time, NO_ABERRATION);        /* FIXFIXFIX: use aberration or not? */
    if (bv.status != ASTRO_SUCCESS)
        return AngleError(bv.status);

    return AngleBetween(sv, bv);
}

astro_elongation_t Astronomy_Elongation(astro_body_t body, astro_time_t time)
{
    astro_elongation_t result;
    astro_angle_result_t angres;

    angres = Astronomy_LongitudeFromSun(body, time);
    if (angres.status != ASTRO_SUCCESS)
        return ElongError(angres.status);

    if (angres.angle > 180.0)
    {
        result.visibility = VISIBLE_MORNING;
        result.relative_longitude = 360.0 - angres.angle;
    }
    else
    {
        result.visibility = VISIBLE_EVENING;
        result.relative_longitude = angres.angle;
    }

    angres = Astronomy_AngleFromSun(body, time);
    if (angres.status != ASTRO_SUCCESS)
        return ElongError(angres.status);

    result.elongation = angres.angle;
    result.time = time;
    result.status = ASTRO_SUCCESS;

    return result;
}

static astro_func_result_t neg_elong_slope(void *context, astro_time_t time)
{
    static const double dt = 0.1;    
    astro_angle_result_t e1, e2;
    astro_func_result_t result;
    astro_body_t body = *((astro_body_t *)context);
    astro_time_t t1 = Astronomy_AddDays(time, -dt/2.0);
    astro_time_t t2 = Astronomy_AddDays(time, +dt/2.0);

    e1 = Astronomy_AngleFromSun(body, t1);
    if (e1.status != ASTRO_SUCCESS)
        return FuncError(e1.status);

    e2 = Astronomy_AngleFromSun(body, t2);
    if (e2.status)
        return FuncError(e2.status);

    result.value = (e1.angle - e2.angle)/dt;
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_elongation_t Astronomy_SearchMaxElongation(astro_body_t body, astro_time_t startDate)
{
    double s1, s2;
    int iter;
    astro_angle_result_t plon, elon;
    astro_time_t t_start;
    double rlon, rlon_lo, rlon_hi, adjust_days;
    astro_func_result_t syn;
    astro_search_result_t search1, search2, searchx;
    astro_time_t t1, t2;
    astro_func_result_t m1, m2;

    /* Determine the range of relative longitudes within which maximum elongation can occur for this planet. */
    switch (body)
    {
    case BODY_MERCURY:
        s1 = 50.0;
        s2 = 85.0;
        break;

    case BODY_VENUS:
        s1 = 40.0;
        s2 = 50.0;
        break;

    default:
        /* SearchMaxElongation works for Mercury and Venus only. */
        return ElongError(ASTRO_INVALID_BODY);
    }

    syn = SynodicPeriod(body);
    if (syn.status != ASTRO_SUCCESS)
        return ElongError(syn.status);

    iter = 0;
    while (++iter <= 2)
    {
        plon = Astronomy_EclipticLongitude(body, startDate);
        if (plon.status != ASTRO_SUCCESS)
            return ElongError(plon.status);

        elon = Astronomy_EclipticLongitude(BODY_EARTH, startDate);
        if (elon.status != ASTRO_SUCCESS)
            return ElongError(elon.status);

        rlon = LongitudeOffset(plon.angle - elon.angle);    /* clamp to (-180, +180] */

        /* The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees */
        /* because there is a cusp there that causes a discontinuity in the derivative. */
        /* So we need to guard against searching near such times. */
        if (rlon >= -s1 && rlon < +s1)
        {
            /* Seek to the window [+s1, +s2]. */
            adjust_days = 0.0;
            /* Search forward for the time t1 when rel lon = +s1. */
            rlon_lo = +s1;
            /* Search forward for the time t2 when rel lon = +s2. */
            rlon_hi = +s2;
        }
        else if (rlon > +s2 || rlon < -s2)
        {
            /* Seek to the next search window at [-s2, -s1]. */
            adjust_days = 0.0;
            /* Search forward for the time t1 when rel lon = -s2. */
            rlon_lo = -s2;
            /* Search forward for the time t2 when rel lon = -s1. */
            rlon_hi = -s1;
        }
        else if (rlon >= 0.0)
        {
            /* rlon must be in the middle of the window [+s1, +s2]. */
            /* Search BACKWARD for the time t1 when rel lon = +s1. */
            adjust_days = -syn.value / 4.0;
            rlon_lo = +s1;
            rlon_hi = +s2;
            /* Search forward from t1 to find t2 such that rel lon = +s2. */
        }
        else
        {
            /* rlon must be in the middle of the window [-s2, -s1]. */
            /* Search BACKWARD for the time t1 when rel lon = -s2. */
            adjust_days = -syn.value / 4.0;
            rlon_lo = -s2;
            /* Search forward from t1 to find t2 such that rel lon = -s1. */
            rlon_hi = -s1;
        }

        t_start = Astronomy_AddDays(startDate, adjust_days);

        search1 = Astronomy_SearchRelativeLongitude(body, rlon_lo, t_start);
        if (search1.status != ASTRO_SUCCESS)
            return ElongError(search1.status);
        t1 = search1.time;

        search2 = Astronomy_SearchRelativeLongitude(body, rlon_hi, t1);
        if (search2.status != ASTRO_SUCCESS)
            return ElongError(search2.status);
        t2 = search2.time;

        /* Now we have a time range [t1,t2] that brackets a maximum elongation event. */
        /* Confirm the bracketing. */
        m1 = neg_elong_slope(&body, t1);
        if (m1.status != ASTRO_SUCCESS)
            return ElongError(m1.status);

        if (m1.value >= 0)
            return ElongError(ASTRO_INTERNAL_ERROR);    /* there is a bug in the bracketing algorithm! */

        m2 = neg_elong_slope(&body, t2);
        if (m2.status != ASTRO_SUCCESS)
            return ElongError(m2.status);

        if (m2.value <= 0)
            return ElongError(ASTRO_INTERNAL_ERROR);    /* there is a bug in the bracketing algorithm! */

        /* Use the generic search algorithm to home in on where the slope crosses from negative to positive. */
        searchx = Astronomy_Search(neg_elong_slope, &body, t1, t2, 10.0);
        if (searchx.status != ASTRO_SUCCESS)
            return ElongError(searchx.status);

        if (searchx.time.tt >= startDate.tt)
            return Astronomy_Elongation(body, searchx.time);

        /* This event is in the past (earlier than startDate). */
        /* We need to search forward from t2 to find the next possible window. */
        /* We never need to search more than twice. */
        startDate = Astronomy_AddDays(t2, 1.0);
    }

    return ElongError(ASTRO_SEARCH_FAILURE);
}

astro_angle_result_t Astronomy_LongitudeFromSun(astro_body_t body, astro_time_t time)
{
    astro_vector_t sv, bv;
    astro_ecliptic_t se, be;
    astro_angle_result_t result;

    if (body == BODY_EARTH)
        return AngleError(ASTRO_EARTH_NOT_ALLOWED);

    sv = Astronomy_GeoVector(BODY_SUN, time, NO_ABERRATION);    /* FIXFIXFIX: use aberration or not? */
    se = Astronomy_Ecliptic(sv);        /* checks for errors in sv */
    if (se.status != ASTRO_SUCCESS)
        return AngleError(se.status);

    bv = Astronomy_GeoVector(body, time, NO_ABERRATION);        /* FIXFIXFIX: use aberration or not? */
    be = Astronomy_Ecliptic(bv);        /* checks for errors in bv */
    if (be.status != ASTRO_SUCCESS)
        return AngleError(be.status);

    result.status = ASTRO_SUCCESS;
    result.angle = NormalizeLongitude(be.elon - se.elon);
    return result;
}

astro_angle_result_t Astronomy_MoonPhase(astro_time_t time)
{
    return Astronomy_LongitudeFromSun(BODY_MOON, time);
}

static astro_func_result_t moon_offset(void *context, astro_time_t time)
{
    astro_func_result_t result;
    double targetLon = *((double *)context);
    astro_angle_result_t angres = Astronomy_MoonPhase(time);
    if (angres.status != ASTRO_SUCCESS)
        return FuncError(angres.status);
    result.value = LongitudeOffset(angres.angle - targetLon);
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_search_result_t Astronomy_SearchMoonPhase(double targetLon, astro_time_t dateStart, double limitDays)
{
    /*
        To avoid discontinuities in the moon_offset function causing problems,
        we need to approximate when that function will next return 0.
        We probe it with the start time and take advantage of the fact
        that every lunar phase repeats roughly every 29.5 days.
        There is a surprising uncertainty in the quarter timing,
        due to the eccentricity of the moon's orbit.
        I have seen up to 0.826 days away from the simple prediction.
        To be safe, we take the predicted time of the event and search
        +/-0.9 days around it (a 1.8-day wide window).
        But we must return null if the final result goes beyond limitDays after dateStart.
    */
    const double uncertainty = 0.9;
    astro_func_result_t funcres;
    double ya, est_dt, dt1, dt2;
    astro_time_t t1, t2;

    funcres = moon_offset(&targetLon, dateStart);
    if (funcres.status != ASTRO_SUCCESS)
        return SearchError(funcres.status);

    ya = funcres.value;
    if (ya > 0.0) ya -= 360.0;  /* force searching forward in time, not backward */
    est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
    dt1 = est_dt - uncertainty;
    if (dt1 > limitDays)
        return SearchError(ASTRO_NO_MOON_QUARTER);    /* not possible for moon phase to occur within specified window (too short) */
    dt2 = est_dt + uncertainty;
    if (limitDays < dt2)
        dt2 = limitDays;
    t1 = Astronomy_AddDays(dateStart, dt1);
    t2 = Astronomy_AddDays(dateStart, dt2);
    return Astronomy_Search(moon_offset, &targetLon, t1, t2, 1.0);
}

astro_moon_quarter_t Astronomy_SearchMoonQuarter(astro_time_t dateStart)
{
    astro_moon_quarter_t mq;
    astro_angle_result_t angres;
    astro_search_result_t srchres;

    /* Determine what the next quarter phase will be. */
    angres = Astronomy_MoonPhase(dateStart);
    if (angres.status != ASTRO_SUCCESS)
        return MoonQuarterError(angres.status);

    mq.quarter = (1 + (int)floor(angres.angle / 90.0)) % 4;
    srchres = Astronomy_SearchMoonPhase(90.0 * mq.quarter, dateStart, 10.0);
    if (srchres.status != ASTRO_SUCCESS)
        return MoonQuarterError(srchres.status);

    mq.status = ASTRO_SUCCESS;
    mq.time = srchres.time;
    return mq;
}

astro_moon_quarter_t Astronomy_NextMoonQuarter(astro_moon_quarter_t mq)
{
    astro_time_t time;
    astro_moon_quarter_t next_mq;

    /* Skip 6 days past the previous found moon quarter to find the next one. */
    /* This is less than the minimum possible increment. */
    /* So far I have seen the interval well contained by the range (6.5, 8.3) days. */

    time = Astronomy_AddDays(mq.time, 6.0);
    next_mq = Astronomy_SearchMoonQuarter(time);
    if (next_mq.status == ASTRO_SUCCESS)
    {
        /* Verify that we found the expected moon quarter. */
        if (next_mq.quarter != (1 + mq.quarter) % 4)
            return MoonQuarterError(ASTRO_WRONG_MOON_QUARTER);  /* internal error! we found the wrong moon quarter */
    }
    return next_mq;
}

static astro_func_result_t rlon_offset(astro_body_t body, astro_time_t time, int direction, double targetRelLon)
{
    astro_func_result_t result;
    astro_angle_result_t plon, elon;
    double diff;

    plon = Astronomy_EclipticLongitude(body, time);
    if (plon.status != ASTRO_SUCCESS)
        return FuncError(plon.status);

    elon = Astronomy_EclipticLongitude(BODY_EARTH, time);
    if (elon.status != ASTRO_SUCCESS)
        return FuncError(elon.status);

    diff = direction * (elon.angle - plon.angle);
    result.value = LongitudeOffset(diff - targetRelLon);
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_search_result_t Astronomy_SearchRelativeLongitude(astro_body_t body, double targetRelLon, astro_time_t startDate)
{
    astro_search_result_t result;
    astro_func_result_t syn;
    astro_func_result_t error_angle;
    double prev_angle;
    astro_time_t time;
    int iter, direction;

    if (body == BODY_EARTH)
        return SearchError(ASTRO_EARTH_NOT_ALLOWED);

    if (body == BODY_MOON)
        return SearchError(ASTRO_INVALID_BODY);

    syn = SynodicPeriod(body);
    if (syn.status != ASTRO_SUCCESS)
        return SearchError(syn.status);

    direction = IsSuperiorPlanet(body) ? +1 : -1;

    /* Iterate until we converge on the desired event. */
    /* Calculate the error angle, which will be a negative number of degrees, */
    /* meaning we are "behind" the target relative longitude. */

    error_angle = rlon_offset(body, startDate, direction, targetRelLon);
    if (error_angle.status != ASTRO_SUCCESS)
        return SearchError(error_angle.status);

    if (error_angle.value > 0) 
        error_angle.value -= 360;    /* force searching forward in time */

    time = startDate;
    for (iter = 0; iter < 100; ++iter)
    {
        /* Estimate how many days in the future (positive) or past (negative) */
        /* we have to go to get closer to the target relative longitude. */
        double day_adjust = (-error_angle.value/360.0) * syn.value;
        time = Astronomy_AddDays(time, day_adjust);
        if (fabs(day_adjust) * SECONDS_PER_DAY < 1.0)
        {
            result.time = time;
            result.status = ASTRO_SUCCESS;
            return result;
        }

        prev_angle = error_angle.value;
        error_angle = rlon_offset(body, time, direction, targetRelLon);
        if (error_angle.status != ASTRO_SUCCESS)
            return SearchError(error_angle.status);

        if (fabs(prev_angle) < 30.0 && (prev_angle != error_angle.value))
        {
            /* Improve convergence for Mercury/Mars (eccentric orbits) */
            /* by adjusting the synodic period to more closely match the */
            /* variable speed of both planets in this part of their respective orbits. */
            double ratio = prev_angle / (prev_angle - error_angle.value);
            if (ratio > 0.5 && ratio < 2.0)
                syn.value *= ratio;
        }
    }

    return SearchError(ASTRO_NO_CONVERGE);
}

astro_hour_angle_t Astronomy_SearchHourAngle(
    astro_body_t body,
    astro_observer_t observer,
    double hourAngle,
    astro_time_t dateStart)
{
    int iter = 0;
    astro_time_t time;
    astro_equatorial_t ofdate;
    astro_hour_angle_t result;
    double delta_sidereal_hours, delta_days, gast;

    time = dateStart;
    for(;;) 
    {
        ++iter;

        /* Calculate Greenwich Apparent Sidereal Time (GAST) at the given time. */
        gast = sidereal_time(time);

        /* Obtain equatorial coordinates of date for the body. */
        ofdate = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, ABERRATION);
        if (ofdate.status != ASTRO_SUCCESS)
            return HourAngleError(ofdate.status);

        /* Calculate the adjustment needed in sidereal time */
        /* to bring the hour angle to the desired value. */

        delta_sidereal_hours = fmod((hourAngle + ofdate.ra - observer.longitude/15) - gast, 24.0);
        if (iter == 1) 
        {
            /* On the first iteration, always search forward in time. */
            if (delta_sidereal_hours < 0)
                delta_sidereal_hours += 24;
        } 
        else 
        {
            /* On subsequent iterations, we make the smallest possible adjustment, */
            /* either forward or backward in time. */
            if (delta_sidereal_hours < -12.0)
                delta_sidereal_hours += 24.0;
            else if (delta_sidereal_hours > +12.0)
                delta_sidereal_hours -= 24.0;
        }

        /* If the error is tolerable (less than 0.1 seconds), the search has succeeded. */
        if (fabs(delta_sidereal_hours) * 3600.0 < 0.1) 
        {
            result.hor = Astronomy_Horizon(time, observer, ofdate.ra, ofdate.dec, REFRACTION_NORMAL);
            result.time = time;
            result.status = ASTRO_SUCCESS;
            return result;
        }

        /* We need to loop another time to get more accuracy. */
        /* Update the terrestrial time (in solar days) adjusting by sidereal time (sidereal hours). */
        delta_days = (delta_sidereal_hours / 24.0) * SOLAR_DAYS_PER_SIDEREAL_DAY;
        time = Astronomy_AddDays(time, delta_days);
    }
}

/** @cond DOXYGEN_SKIP */
typedef struct
{
    astro_body_t        body;
    int                 direction;
    astro_observer_t    observer;
    double              body_radius_au;
}
context_peak_altitude_t;
/** @endcond */

static astro_func_result_t peak_altitude(void *context, astro_time_t time)
{
    astro_func_result_t result;
    astro_equatorial_t ofdate;
    astro_horizon_t hor;
    const context_peak_altitude_t *p = context;

    /*
        Return the angular altitude above or below the horizon
        of the highest part (the peak) of the given object.
        This is defined as the apparent altitude of the center of the body plus
        the body's angular radius.
        The 'direction' parameter controls whether the angle is measured 
        positive above the horizon or positive below the horizon, 
        depending on whether the caller wants rise times or set times, respectively.
    */

    ofdate = Astronomy_Equator(p->body, time, p->observer, EQUATOR_OF_DATE, ABERRATION);
    if (ofdate.status != ASTRO_SUCCESS)
        return FuncError(ofdate.status);

    /* We calculate altitude without refraction, then subtract fixed refraction near the horizon. */
    /* This gives us the time of rise/set without the extra work. */
    hor = Astronomy_Horizon(time, p->observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
    result.value = p->direction * (hor.altitude + RAD2DEG*(p->body_radius_au / ofdate.dist) + REFRACTION_NEAR_HORIZON);
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_search_result_t Astronomy_SearchRiseSet(
    astro_body_t body,
    astro_observer_t observer,
    astro_direction_t direction,
    astro_time_t dateStart,
    double limitDays)
{    
    context_peak_altitude_t context;
    double ha_before, ha_after;
    astro_time_t time_start, time_before;
    astro_func_result_t alt_before, alt_after;
    astro_hour_angle_t evt_before, evt_after;

    switch (direction)
    {
    case DIRECTION_RISE:
        ha_before = 12.0;   /* minimum altitude (bottom) happens BEFORE the body rises. */
        ha_after = 0.0;     /* maximum altitude (culmination) happens AFTER the body rises. */
        break;

    case DIRECTION_SET:
        ha_before = 0.0;    /* culmination happens BEFORE the body sets. */
        ha_after = 12.0;    /* bottom happens AFTER the body sets. */
        break;

    default:
        return SearchError(ASTRO_INVALID_PARAMETER);
    }

    /* Set up the context structure for the search function 'peak_altitude'. */
    context.body = body;
    context.direction = (int)direction;
    context.observer = observer;
    switch (body)
    {
    case BODY_SUN:  context.body_radius_au = SUN_RADIUS_AU;     break;
    case BODY_MOON: context.body_radius_au = MOON_RADIUS_AU;    break;
    default:        context.body_radius_au = 0.0;               break;
    }

    /*
        See if the body is currently above/below the horizon.
        If we are looking for next rise time and the body is below the horizon,
        we use the current time as the lower time bound and the next culmination
        as the upper bound.
        If the body is above the horizon, we search for the next bottom and use it
        as the lower bound and the next culmination after that bottom as the upper bound.
        The same logic applies for finding set times, only we swap the hour angles.
        The peak_altitude() function already considers the 'direction' parameter.
    */

    time_start = dateStart;
    alt_before = peak_altitude(&context, time_start);
    if (alt_before.status != ASTRO_SUCCESS)
        return SearchError(alt_before.status);

    if (alt_before.value > 0.0)
    {
        /* We are past the sought event, so we have to wait for the next "before" event (culm/bottom). */
        evt_before = Astronomy_SearchHourAngle(body, observer, ha_before, time_start);
        if (evt_before.status != ASTRO_SUCCESS)
            return SearchError(evt_before.status);

        time_before = evt_before.time;

        alt_before = peak_altitude(&context, time_before);
        if (alt_before.status != ASTRO_SUCCESS)
            return SearchError(alt_before.status);
    }
    else
    {
        /* We are before or at the sought event, so we find the next "after" event (bottom/culm), */
        /* and use the current time as the "before" event. */
        time_before = time_start;
    }
    
    evt_after = Astronomy_SearchHourAngle(body, observer, ha_after, time_before);
    if (evt_after.status != ASTRO_SUCCESS)
        return SearchError(evt_after.status);

    alt_after = peak_altitude(&context, evt_after.time);
    if (alt_after.status != ASTRO_SUCCESS)
        return SearchError(alt_after.status);

    for(;;)
    {
        if (alt_before.value <= 0.0 && alt_after.value > 0.0)
        {
            /* Search between evt_before and evt_after for the desired event. */
            astro_search_result_t result = Astronomy_Search(peak_altitude, &context, time_before, evt_after.time, 1.0);

            /* ASTRO_SEARCH_FAILURE is a special error that indicates a normal lack of finding a solution. */
            /* If successful, or any other error, return immediately. */
            if (result.status != ASTRO_SEARCH_FAILURE) 
                return result;
        }
        
        /* If we didn't find the desired event, use time_after to find the next before-event. */
        evt_before = Astronomy_SearchHourAngle(body, observer, ha_before, evt_after.time);
        if (evt_before.status != ASTRO_SUCCESS)
            return SearchError(evt_before.status);

        evt_after = Astronomy_SearchHourAngle(body, observer, ha_after, evt_before.time);
        if (evt_after.status != ASTRO_SUCCESS)
            return SearchError(evt_after.status);

        if (evt_before.time.ut >= time_start.ut + limitDays)
            return SearchError(ASTRO_SEARCH_FAILURE);

        time_before = evt_before.time;

        alt_before = peak_altitude(&context, evt_before.time);
        if (alt_before.status != ASTRO_SUCCESS)
            return SearchError(alt_before.status);            

        alt_after = peak_altitude(&context, evt_after.time);
        if (alt_after.status != ASTRO_SUCCESS)
            return SearchError(alt_after.status);        
    }
}

static double MoonMagnitude(double phase, double helio_dist, double geo_dist) 
{
    /* https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve */
    double rad = phase * DEG2RAD;
    double rad2 = rad * rad;
    double rad4 = rad2 * rad2;
    double mag = -12.717 + 1.49*fabs(rad) + 0.0431*rad4;
    double moon_mean_distance_au = 385000.6 / KM_PER_AU;
    double geo_au = geo_dist / moon_mean_distance_au;
    mag += 5*log10(helio_dist * geo_au);
    return mag;
}

static astro_status_t SaturnMagnitude(
    double phase, 
    double helio_dist, 
    double geo_dist, 
    astro_vector_t gc, 
    astro_time_t time,
    double *mag,
    double *ring_tilt) 
{
    astro_ecliptic_t eclip;
    double ir, Nr, lat, lon, tilt, sin_tilt;

    *mag = *ring_tilt = NAN;

    /* Based on formulas by Paul Schlyter found here: */
    /* http://www.stjarnhimlen.se/comp/ppcomp.html#15 */

    /* We must handle Saturn's rings as a major component of its visual magnitude. */
    /* Find geocentric ecliptic coordinates of Saturn. */
    eclip = Astronomy_Ecliptic(gc);
    if (eclip.status != ASTRO_SUCCESS)
        return eclip.status;

    ir = DEG2RAD * 28.06;   /* tilt of Saturn's rings to the ecliptic, in radians */
    Nr = DEG2RAD * (169.51 + (3.82e-5 * time.tt));    /* ascending node of Saturn's rings, in radians */

    /* Find tilt of Saturn's rings, as seen from Earth. */
    lat = DEG2RAD * eclip.elat;
    lon = DEG2RAD * eclip.elon;
    tilt = asin(sin(lat)*cos(ir) - cos(lat)*sin(ir)*sin(lon-Nr));
    sin_tilt = sin(fabs(tilt));

    *mag = -9.0 + 0.044*phase;
    *mag += sin_tilt*(-2.6 + 1.2*sin_tilt);
    *mag += 5.0 * log10(helio_dist * geo_dist);

    *ring_tilt = RAD2DEG * tilt;

    return ASTRO_SUCCESS;
}

static astro_status_t VisualMagnitude(
    astro_body_t body, 
    double phase, 
    double helio_dist, 
    double geo_dist, 
    double *mag)
{
    /* For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212 */
    double c0, c1=0, c2=0, c3=0, x;
    *mag = NAN;
    switch (body) 
    {
    case BODY_MERCURY:  c0 = -0.60, c1 = +4.98, c2 = -4.88, c3 = +3.02; break;
    case BODY_VENUS:
        if (phase < 163.6)
            c0 = -4.47, c1 = +1.03, c2 = +0.57, c3 = +0.13;
        else
            c0 = 0.98, c1 = -1.02;
        break;
    case BODY_MARS:        c0 = -1.52, c1 = +1.60;   break;
    case BODY_JUPITER:     c0 = -9.40, c1 = +0.50;   break;
    case BODY_URANUS:      c0 = -7.19, c1 = +0.25;   break;
    case BODY_NEPTUNE:     c0 = -6.87;               break;
    case BODY_PLUTO:       c0 = -1.00, c1 = +4.00;   break;
    default: return ASTRO_INVALID_BODY;
    }

    x = phase / 100;
    *mag = c0 + x*(c1 + x*(c2 + x*c3));
    *mag += 5.0 * log10(helio_dist * geo_dist);
    return ASTRO_SUCCESS;
}

astro_illum_t Astronomy_Illumination(astro_body_t body, astro_time_t time)
{
    astro_vector_t earth;   /* vector from Sun to Earth */
    astro_vector_t hc;      /* vector from Sun to body */
    astro_vector_t gc;      /* vector from Earth to body */
    double mag;             /* visual magnitude */   
    astro_angle_result_t phase;     /* phase angle in degrees between Earth and Sun as seen from body */
    double helio_dist;      /* distance from Sun to body */
    double geo_dist;        /* distance from Earth to body */
    double ring_tilt = 0.0; /* Saturn's ring tilt (0 for all other bodies) */
    astro_illum_t illum;
    astro_status_t status;

    if (body == BODY_EARTH)
        return IllumError(ASTRO_EARTH_NOT_ALLOWED);

    earth = CalcEarth(time);
    if (earth.status != ASTRO_SUCCESS)
        return IllumError(earth.status);

    if (body == BODY_SUN)
    {
        gc.status = ASTRO_SUCCESS;
        gc.t = time;
        gc.x = -earth.x;
        gc.y = -earth.y;
        gc.z = -earth.z;

        hc.status = ASTRO_SUCCESS;
        hc.t = time;
        hc.x = 0.0;
        hc.y = 0.0;
        hc.z = 0.0;

        /* The Sun emits light instead of reflecting it, */
        /* so we report a placeholder phase angle of 0. */
        phase.status = ASTRO_SUCCESS;
        phase.angle = 0.0;       
    }
    else
    {
        if (body == BODY_MOON)
        {
            /* For extra numeric precision, use geocentric Moon formula directly. */
            gc = Astronomy_GeoMoon(time);
            if (gc.status != ASTRO_SUCCESS)
                return IllumError(gc.status);

            hc.status = ASTRO_SUCCESS;
            hc.t = time;
            hc.x = earth.x + gc.x;
            hc.y = earth.y + gc.y;
            hc.z = earth.z + gc.z;
        }
        else
        {
            /* For planets, the heliocentric vector is more direct to calculate. */
            hc = Astronomy_HelioVector(body, time);
            if (hc.status != ASTRO_SUCCESS)
                return IllumError(hc.status);

            gc.status = ASTRO_SUCCESS;
            gc.t = time;
            gc.x = hc.x - earth.x;
            gc.y = hc.y - earth.y;
            gc.z = hc.z - earth.z;
        }

        phase = AngleBetween(gc, hc);
        if (phase.status != ASTRO_SUCCESS)
            return IllumError(phase.status);
    }

    geo_dist = Astronomy_VectorLength(gc);
    helio_dist = Astronomy_VectorLength(hc);

    switch (body)
    {
    case BODY_SUN:
        mag = -0.17 + 5.0*log10(geo_dist / AU_PER_PARSEC);
        break;

    case BODY_MOON:
        mag = MoonMagnitude(phase.angle, helio_dist, geo_dist);
        break;

    case BODY_SATURN:
        status = SaturnMagnitude(phase.angle, helio_dist, geo_dist, gc, time, &mag, &ring_tilt);
        if (status != ASTRO_SUCCESS)
            return IllumError(status);
        break;

    default:
        status = VisualMagnitude(body, phase.angle, helio_dist, geo_dist, &mag);
        break;
    }

    illum.status = ASTRO_SUCCESS;
    illum.time = time;
    illum.mag = mag;
    illum.phase_angle = phase.angle;
    illum.helio_dist = helio_dist;
    illum.ring_tilt = ring_tilt;

    return illum;
}

static astro_func_result_t mag_slope(void *context, astro_time_t time)
{
    /* 
        The Search() function finds a transition from negative to positive values.
        The derivative of magnitude y with respect to time t (dy/dt)
        is negative as an object gets brighter, because the magnitude numbers
        get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
        dy/dt > 0.
    */
    static const double dt = 0.01;
    astro_illum_t y1, y2;
    astro_body_t body = *((astro_body_t *)context);
    astro_time_t t1 = Astronomy_AddDays(time, -dt/2);
    astro_time_t t2 = Astronomy_AddDays(time, +dt/2);
    astro_func_result_t result;

    y1 = Astronomy_Illumination(body, t1);
    if (y1.status != ASTRO_SUCCESS)
        return FuncError(y1.status);

    y2 = Astronomy_Illumination(body, t2);
    if (y2.status != ASTRO_SUCCESS)
        return FuncError(y2.status);

    result.value = (y2.mag - y1.mag) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_illum_t Astronomy_SearchPeakMagnitude(astro_body_t body, astro_time_t startDate)
{
    /* s1 and s2 are relative longitudes within which peak magnitude of Venus can occur. */
    static const double s1 = 10.0;
    static const double s2 = 30.0;
    int iter;
    astro_angle_result_t plon, elon;
    astro_search_result_t t1, t2, tx;
    astro_func_result_t syn, m1, m2;
    astro_time_t t_start;
    double rlon, rlon_lo, rlon_hi, adjust_days;

    if (body != BODY_VENUS)
        return IllumError(ASTRO_INVALID_BODY);

    iter = 0;
    while (++iter <= 2)
    {
        /* Find current heliocentric relative longitude between the */
        /* inferior planet and the Earth. */
        plon = Astronomy_EclipticLongitude(body, startDate);
        if (plon.status != ASTRO_SUCCESS)
            return IllumError(plon.status);

        elon = Astronomy_EclipticLongitude(BODY_EARTH, startDate);
        if (elon.status != ASTRO_SUCCESS)
            return IllumError(elon.status);

        rlon = LongitudeOffset(plon.angle - elon.angle);    /* clamp to (-180, +180]. */

        /* The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees */
        /* because there is a cusp there that causes a discontinuity in the derivative. */
        /* So we need to guard against searching near such times. */

        if (rlon >= -s1 && rlon < +s1)
        {
            /* Seek to the window [+s1, +s2]. */
            adjust_days = 0.0;
            /* Search forward for the time t1 when rel lon = +s1. */
            rlon_lo = +s1;
            /* Search forward for the time t2 when rel lon = +s2. */
            rlon_hi = +s2;
        }
        else if (rlon >= +s2 || rlon < -s2 ) 
        {
            /* Seek to the next search window at [-s2, -s1]. */
            adjust_days = 0.0;
            /* Search forward for the time t1 when rel lon = -s2. */
            rlon_lo = -s2;
            /* Search forward for the time t2 when rel lon = -s1. */
            rlon_hi = -s1;
        } 
        else if (rlon >= 0) 
        {
            /* rlon must be in the middle of the window [+s1, +s2]. */
            /* Search BACKWARD for the time t1 when rel lon = +s1. */
            syn = SynodicPeriod(body);
            if (syn.status != ASTRO_SUCCESS)
                return IllumError(syn.status);
            adjust_days = -syn.value / 4;
            rlon_lo = +s1;
            /* Search forward from t1 to find t2 such that rel lon = +s2. */
            rlon_hi = +s2;
        } 
        else 
        {
            /* rlon must be in the middle of the window [-s2, -s1]. */
            /* Search BACKWARD for the time t1 when rel lon = -s2. */
            syn = SynodicPeriod(body);
            if (syn.status != ASTRO_SUCCESS)
                return IllumError(syn.status);
            adjust_days = -syn.value / 4;
            rlon_lo = -s2;
            /* Search forward from t1 to find t2 such that rel lon = -s1. */
            rlon_hi = -s1;
        }
        t_start = Astronomy_AddDays(startDate, adjust_days);
        t1 = Astronomy_SearchRelativeLongitude(body, rlon_lo, t_start);
        if (t1.status != ASTRO_SUCCESS)
            return IllumError(t1.status);
        t2 = Astronomy_SearchRelativeLongitude(body, rlon_hi, t1.time);
        if (t2.status != ASTRO_SUCCESS)
            return IllumError(t2.status);

        /* Now we have a time range [t1,t2] that brackets a maximum magnitude event. */
        /* Confirm the bracketing. */
        m1 = mag_slope(&body, t1.time);
        if (m1.status != ASTRO_SUCCESS)
            return IllumError(m1.status);
        if (m1.value >= 0.0) 
            return IllumError(ASTRO_INTERNAL_ERROR);    /* should never happen! */

        m2 = mag_slope(&body, t2.time);
        if (m2.status != ASTRO_SUCCESS)
            return IllumError(m2.status);
        if (m2.value <= 0.0)
            return IllumError(ASTRO_INTERNAL_ERROR);    /* should never happen! */

        /* Use the generic search algorithm to home in on where the slope crosses from negative to positive. */
        tx = Astronomy_Search(mag_slope, &body, t1.time, t2.time, 10.0);
        if (tx.status != ASTRO_SUCCESS) 
            return IllumError(tx.status);

        if (tx.time.tt >= startDate.tt)
            return Astronomy_Illumination(body, tx.time);

        /* This event is in the past (earlier than startDate). */
        /* We need to search forward from t2 to find the next possible window. */
        /* We never need to search more than twice. */
        startDate = Astronomy_AddDays(t2.time, 1.0);
    }

    return IllumError(ASTRO_SEARCH_FAILURE);
}

static double MoonDistance(astro_time_t t)
{
    double lon, lat, dist;
    CalcMoon(t.tt / 36525.0, &lon, &lat, &dist);
    return dist;
}

static astro_func_result_t distance_slope(void *context, astro_time_t time)
{
    static const double dt = 0.001;
    astro_time_t t1 = Astronomy_AddDays(time, -dt/2.0);
    astro_time_t t2 = Astronomy_AddDays(time, +dt/2.0);
    double dist1, dist2;
    int direction = *((int *)context);
    astro_func_result_t result;

    dist1 = MoonDistance(t1);
    dist2 = MoonDistance(t2);
    result.value = direction * (dist2 - dist1) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}

astro_apsis_t Astronomy_SearchLunarApsis(astro_time_t startTime)
{
    astro_time_t t1, t2;
    astro_search_result_t search;
    astro_func_result_t m1, m2;
    int positive_direction = +1;
    int negative_direction = -1;
    const double increment = 5.0;   /* number of days to skip in each iteration */
    astro_apsis_t result;
    int iter;

    /*
        Check the rate of change of the distance dr/dt at the start time.
        If it is positive, the Moon is currently getting farther away,
        so start looking for apogee.
        Conversely, if dr/dt < 0, start looking for perigee.
        Either way, the polarity of the slope will change, so the product will be negative.
        Handle the crazy corner case of exactly touching zero by checking for m1*m2 <= 0.
    */

    t1 = startTime;
    m1 = distance_slope(&positive_direction, t1);
    if (m1.status != ASTRO_SUCCESS)
        return ApsisError(m1.status);

    for (iter=0; iter * increment < 2.0 * MEAN_SYNODIC_MONTH; ++iter)
    {
        t2 = Astronomy_AddDays(t1, increment);
        m2 = distance_slope(&positive_direction, t2);
        if (m2.status != ASTRO_SUCCESS)
            return ApsisError(m2.status);

        if (m1.value * m2.value <= 0.0)
        {
            /* There is a change of slope polarity within the time range [t1, t2]. */
            /* Therefore this time range contains an apsis. */
            /* Figure out whether it is perigee or apogee. */

            if (m1.value < 0.0 || m2.value > 0.0)
            {
                /* We found a minimum-distance event: perigee. */
                /* Search the time range for the time when the slope goes from negative to positive. */
                search = Astronomy_Search(distance_slope, &positive_direction, t1, t2, 1.0);
                result.kind = APSIS_PERICENTER;
            }
            else if (m1.value > 0.0 || m2.value < 0.0)
            {
                /* We found a maximum-distance event: apogee. */
                /* Search the time range for the time when the slope goes from positive to negative. */
                search = Astronomy_Search(distance_slope, &negative_direction, t1, t2, 1.0);
                result.kind = APSIS_APOCENTER;
            }
            else
            {
                /* This should never happen. It should not be possible for both slopes to be zero. */
                return ApsisError(ASTRO_INTERNAL_ERROR);
            }

            if (search.status != ASTRO_SUCCESS)
                return ApsisError(search.status);

            result.status = ASTRO_SUCCESS;
            result.time = search.time;
            result.dist_au = MoonDistance(search.time);
            result.dist_km = result.dist_au * KM_PER_AU;
            return result;
        }

        /* We have not yet found a slope polarity change. Keep searching. */
        t1 = t2;
        m1 = m2;
    }

    /* It should not be possible to fail to find an apsis within 2 synodic months. */
    return ApsisError(ASTRO_INTERNAL_ERROR);
}

astro_apsis_t Astronomy_NextLunarApsis(astro_apsis_t apsis)
{
    static const double skip = 11.0;    /* number of days to skip to start looking for next apsis event */
    astro_apsis_t next;
    astro_time_t time;
    astro_apsis_kind_t expected = APSIS_INVALID;

    if (apsis.status != ASTRO_SUCCESS)
        return ApsisError(ASTRO_INVALID_PARAMETER);

    time = Astronomy_AddDays(apsis.time, skip);
    next = Astronomy_SearchLunarApsis(time);
    if (next.status == ASTRO_SUCCESS)
    {
        /* Verify that we found the opposite apsis from the previous one. */
        switch (apsis.kind)
        {
        case APSIS_APOCENTER:
            expected = APSIS_PERICENTER;
            break;

        case APSIS_PERICENTER:
            expected = APSIS_APOCENTER;
            break;

        default:
            /* The apsis passed in has an invalid kind. */
            return ApsisError(ASTRO_INVALID_PARAMETER);
        }

        if (next.kind != expected)
            return ApsisError(ASTRO_INTERNAL_ERROR);
    }
    return next;
}

#ifdef __cplusplus
}
#endif
