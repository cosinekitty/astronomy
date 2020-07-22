/*
    Astronomy Engine for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2020 Don Cross <cosinekitty@gmail.com>

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

static const double DAYS_PER_TROPICAL_YEAR = 365.24217;
static const double DEG2RAD = 0.017453292519943296;
static const double RAD2DEG = 57.295779513082321;
static const double ASEC360 = 1296000.0;
static const double ASEC2RAD = 4.848136811095359935899141e-6;
static const double PI2 = 2.0 * PI;
static const double ARC = 3600.0 * 180.0 / PI;          /* arcseconds per radian */
static const double C_AUDAY = 173.1446326846693;        /* speed of light in AU/day */
static const double KM_PER_AU = 1.4959787069098932e+8;
static const double SECONDS_PER_DAY = 24.0 * 3600.0;
static const double SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;
static const double MEAN_SYNODIC_MONTH = 29.530588;     /* average number of days for Moon to return to the same phase */
static const double EARTH_ORBITAL_PERIOD = 365.256;
static const double NEPTUNE_ORBITAL_PERIOD = 60189.0;
static const double REFRACTION_NEAR_HORIZON = 34.0 / 60.0;   /* degrees of refractive "lift" seen for objects near horizon */

static const double SUN_RADIUS_KM  = 695700.0;
#define             SUN_RADIUS_AU  (SUN_RADIUS_KM / KM_PER_AU)

#define EARTH_FLATTENING            0.996647180302104
#define EARTH_EQUATORIAL_RADIUS_KM  6378.1366
#define EARTH_EQUATORIAL_RADIUS_AU  (EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU)
#define EARTH_MEAN_RADIUS_KM        6371.0            /* mean radius of the Earth's geoid, without atmosphere */
#define EARTH_ATMOSPHERE_KM           88.0            /* effective atmosphere thickness for lunar eclipses */
#define EARTH_ECLIPSE_RADIUS_KM     (EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM)
/* Note: if we ever need Earth's polar radius, it is (EARTH_FLATTENING * EARTH_EQUATORIAL_RADIUS_KM) */

#define MOON_EQUATORIAL_RADIUS_KM   1738.1
#define MOON_MEAN_RADIUS_KM         1737.4
#define MOON_POLAR_RADIUS_KM        1736.0
#define MOON_EQUATORIAL_RADIUS_AU   (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU)

static const double ASEC180 = 180.0 * 60.0 * 60.0;      /* arcseconds per 180 degrees (or pi radians) */
static const double EARTH_MOON_MASS_RATIO = 81.30056;
static const double SUN_MASS     = 333054.25318;        /* Sun's mass relative to Earth. */
static const double JUPITER_MASS =    317.84997;        /* Jupiter's mass relative to Earth. */
static const double SATURN_MASS  =     95.16745;        /* Saturn's mass relative to Earth. */
static const double URANUS_MASS  =     14.53617;        /* Uranus's mass relative to Earth. */
static const double NEPTUNE_MASS =     17.14886;        /* Neptune's mass relative to Earth. */

/** @cond DOXYGEN_SKIP */
#define ARRAYSIZE(x)    (sizeof(x) / sizeof(x[0]))
#define AU_PER_PARSEC   (ASEC180 / PI)             /* exact definition of how many AU = one parsec */
#define Y2000_IN_MJD    (T0 - MJD_BASIS)
/** @endcond */

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
    case BODY_EMB:      return "EMB";
    case BODY_SSB:      return "SSB";
    default:            return "";
    }
}

/**
 * @brief Returns the #astro_body_t value corresponding to the given English name.
 * @param name One of the following strings: Sun, Moon, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, EMB, SSB.
 * @return If `name` is one of the listed strings (case-sensitive), the returned value is the corresponding #astro_body_t value, otherwise it is `BODY_INVALID`.
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
        if (!strcmp(name, "EMB"))       return BODY_EMB;
        if (!strcmp(name, "SSB"))       return BODY_SSB;
    }
    return BODY_INVALID;
}

/**
 * @brief Returns 1 for planets that are farther from the Sun than the Earth is, 0 otherwise.
 */
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

/**
 * @brief Returns the number of days it takes for a planet to orbit the Sun.
 */
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
    case BODY_NEPTUNE:  return  NEPTUNE_ORBITAL_PERIOD;
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

static astro_spherical_t SphereError(astro_status_t status)
{
    astro_spherical_t sphere;
    sphere.status = status;
    sphere.dist = sphere.lat = sphere.lon = NAN;
    return sphere;
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

static astro_time_t TimeError(void)
{
    astro_time_t time;
    time.tt = time.ut = time.eps = time.psi = NAN;
    return time;
}

static astro_rotation_t RotationErr(astro_status_t status)
{
    astro_rotation_t rotation;
    int i, j;

    rotation.status = status;
    for (i=0; i<3; ++i)
        for (j=0; j<3; ++j)
            rotation.rot[i][j] = NAN;

    return rotation;
}

static astro_moon_quarter_t MoonQuarterError(astro_status_t status)
{
    astro_moon_quarter_t result;
    result.status = status;
    result.quarter = -1;
    result.time = TimeError();
    return result;
}

static astro_elongation_t ElongError(astro_status_t status)
{
    astro_elongation_t result;

    result.status = status;
    result.elongation = NAN;
    result.ecliptic_separation = NAN;
    result.time = TimeError();
    result.visibility = (astro_visibility_t)(-1);

    return result;
}

static astro_hour_angle_t HourAngleError(astro_status_t status)
{
    astro_hour_angle_t result;

    result.status = status;
    result.time = TimeError();
    result.hor.altitude = result.hor.azimuth = result.hor.dec = result.hor.ra = NAN;

    return result;
}

static astro_illum_t IllumError(astro_status_t status)
{
    astro_illum_t result;

    result.status = status;
    result.time = TimeError();
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
    result.time = TimeError();
    result.kind = APSIS_INVALID;
    result.dist_km = result.dist_au = NAN;

    return result;
}

static astro_search_result_t SearchError(astro_status_t status)
{
    astro_search_result_t result;
    result.time = TimeError();
    result.status = status;
    return result;
}

static astro_constellation_t ConstelErr(astro_status_t status)
{
    astro_constellation_t constel;
    constel.status = status;
    constel.symbol = constel.name = NULL;
    constel.ra_1875 = constel.dec_1875 = NAN;
    return constel;
}

static astro_transit_t TransitErr(astro_status_t status)
{
    astro_transit_t transit;
    transit.status = status;
    transit.start = transit.peak = transit.finish = TimeError();
    transit.separation = NAN;
    return transit;
}

static astro_func_result_t SynodicPeriod(astro_body_t body)
{
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
    result.value = fabs(EARTH_ORBITAL_PERIOD / (EARTH_ORBITAL_PERIOD/Tp - 1.0));
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

/**
 * @brief The default Delta T function used by Astronomy Engine.
 *
 * Espenak and Meeus use a series of piecewise polynomials to
 * approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
 * See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
 * This is the default Delta T function used by Astronomy Engine.
 *
 * @param ut
 *      The floating point number of days since noon UTC on January 1, 2000.
 *
 * @returns
 *      The estimated difference TT-UT on the given date, expressed in seconds.
 */
double Astronomy_DeltaT_EspenakMeeus(double ut)
{
    double y, u, u2, u3, u4, u5, u6, u7;

    /*
        Fred Espenak writes about Delta-T generically here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
        https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html

        He provides polynomial approximations for distant years here:
        https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

        They start with a year value 'y' such that y=2000 corresponds
        to the UTC Date 15-January-2000. Convert difference in days
        to mean tropical years.
    */

    y = 2000 + ((ut - 14) / DAYS_PER_TROPICAL_YEAR);

    if (y < -500)
    {
        u = (y - 1820) / 100;
        return -20 + (32 * u*u);
    }
    if (y < 500)
    {
        u = y / 100;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
        return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6;
    }
    if (y < 1600)
    {
        u = (y - 1000) / 100;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
        return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6;
    }
    if (y < 1700)
    {
        u = y - 1600;
        u2 = u*u; u3 = u*u2;
        return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0;
    }
    if (y < 1800)
    {
        u = y - 1700;
        u2 = u*u; u3 = u*u2; u4 = u2*u2;
        return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000;
    }
    if (y < 1860)
    {
        u = y - 1800;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3; u7 = u3*u4;
        return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7;
    }
    if (y < 1900)
    {
        u = y - 1860;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
        return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174;
    }
    if (y < 1920)
    {
        u = y - 1900;
        u2 = u*u; u3 = u*u2; u4 = u2*u2;
        return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4;
    }
    if (y < 1941)
    {
        u = y - 1920;
        u2 = u*u; u3 = u*u2;
        return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3;
    }
    if (y < 1961)
    {
        u = y - 1950;
        u2 = u*u; u3 = u*u2;
        return 29.07 + 0.407*u - u2/233 + u3/2547;
    }
    if (y < 1986)
    {
        u = y - 1975;
        u2 = u*u; u3 = u*u2;
        return 45.45 + 1.067*u - u2/260 - u3/718;
    }
    if (y < 2005)
    {
        u = y - 2000;
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
        return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5;
    }
    if (y < 2050)
    {
        u = y - 2000;
        return 62.92 + 0.32217*u + 0.005589*u*u;
    }
    if (y < 2150)
    {
        u = (y-1820)/100;
        return -20 + 32*u*u - 0.5628*(2150 - y);
    }

    /* all years after 2150 */
    u = (y - 1820) / 100;
    return -20 + (32 * u*u);
}

/**
 * @brief A Delta T function that approximates the one used by the JPL Horizons tool.
 *
 * In order to support unit tests based on data generated by the JPL Horizons online
 * tool, I had to reverse engineer their Delta T function by generating a table that
 * contained it. The main difference between their tool and the Espenak/Meeus function
 * is that they stop extrapolating the Earth's deceleration after the year 2017.
 *
 * @param ut
 *      The floating point number of days since noon UTC on January 1, 2000.
 *
 * @returns
 *      The estimated difference TT-UT on the given date, expressed in seconds.
 */
double Astronomy_DeltaT_JplHorizons(double ut)
{
    if (ut > 17.0 * DAYS_PER_TROPICAL_YEAR)
        ut = 17.0 * DAYS_PER_TROPICAL_YEAR;

    return Astronomy_DeltaT_EspenakMeeus(ut);
}

static astro_deltat_func DeltaTFunc = Astronomy_DeltaT_EspenakMeeus;

/**
 * @brief Changes the function Astronomy Engine uses to calculate Delta T.
 *
 * Most programs should not call this function. It is for advanced use cases only.
 * By default, Astronomy Engine uses the function #Astronomy_DeltaT_EspenakMeeus
 * to estimate changes in the Earth's rotation rate over time.
 * However, for the sake of unit tests that compare calculations against
 * external data sources that use alternative models for Delta T,
 * it is sometimes useful to replace the Delta T model to match.
 * This function allows replacing the Delta T model with any other
 * desired model.
 *
 * @param func
 *      A pointer to a function to convert UT values to DeltaT values.
 */
void Astronomy_SetDeltaTFunction(astro_deltat_func func)
{
    DeltaTFunc = func;
}

static double TerrestrialTime(double ut)
{
    return ut + DeltaTFunc(ut)/86400.0;
}

/**
 * @brief
 *      Converts a J2000 day value to an #astro_time_t value.
 *
 * This function can be useful for reproducing an #astro_time_t structure
 * from its `ut` field only.
 *
 * @param ut
 *      The floating point number of days since noon UTC on January 1, 2000.
 *
 * @returns
 *      An #astro_time_t value for the given `ut` value.
 */
astro_time_t Astronomy_TimeFromDays(double ut)
{
    astro_time_t  time;
    time.ut = ut;
    time.tt = TerrestrialTime(ut);
    time.psi = time.eps = NAN;
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
    t.psi = t.eps = NAN;
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
    time.psi = time.eps = NAN;

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
    sum.eps = sum.psi = NAN;

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
 * @brief Formats an #astro_time_t value as an ISO 8601 string.
 *
 * Given an #astro_time_t value `time`, formats it as an ISO 8601
 * string to the resolution specified by the `format` parameter.
 * The result is stored in the `text` buffer whose capacity in bytes
 * is specified by `size`.
 *
 * @param time
 *      The date and time whose civil time `time.ut` is to be formatted as an ISO 8601 string.
 *      If the civil time is outside the year range 0000 to 9999, the function fails
 *      and returns `ASTRO_BAD_TIME`. Years prior to 1583 are treated as if they are
 *      using the modern Gregorian calendar, even when the Julian calendar was actually in effect.
 *
 * @param format
 *      Specifies the resolution to which the date and time should be formatted,
 *      as explained at #astro_time_format_t.
 *      If the value of `format` is not recognized, the function fails and
 *      returns `ASTRO_INVALID_PARAMETER`.
 *
 * @param text
 *      A pointer to a text buffer to receive the output.
 *      If `text` is `NULL`, this function returns `ASTRO_INVALID_PARAMETER`.
 *      If the function fails for any reason, and `text` is not `NULL`,
 *      and `size` is greater than 0, the `text` buffer is set to an empty string.
 *
 * @param size
 *      The size in bytes of the buffer pointed to by `text`. The buffer must
 *      be large enough to accomodate the output format selected by the
 *      `format` parameter, as specified at #astro_time_format_t.
 *      If `size` is too small to hold the string as specified by `format`,
 *      the `text` buffer is set to `""` (if possible)
 *      and the function returns `ASTRO_BUFFER_TOO_SMALL`.
 *      A buffer that is `TIME_TEXT_BYTES` (25) bytes or larger is always large enough for this function.
 *
 * @return `ASTRO_SUCCESS` on success; otherwise an error as described in the parameter notes.
 */
astro_status_t Astronomy_FormatTime(
    astro_time_t time,
    astro_time_format_t format,
    char *text,
    size_t size)
{
    int nprinted;
    double rounding;
    size_t min_size;
    astro_utc_t utc;

    if (text == NULL)
        return ASTRO_INVALID_PARAMETER;

    if (size == 0)
        return ASTRO_BUFFER_TOO_SMALL;

    text[0] = '\0';     /* initialize to empty string, in case an error occurs */

    /* Validate 'size' parameter and perform date/time rounding. */
    switch (format)
    {
    case TIME_FORMAT_DAY:
        min_size = 11;                          /* "2020-12-31" */
        rounding = 0.0;                         /* no rounding */
        break;

    case TIME_FORMAT_MINUTE:
        min_size = 18;                          /* "2020-12-31T15:47Z" */
        rounding = 0.5 / (24.0 * 60.0);         /* round to nearest minute */
        break;

    case TIME_FORMAT_SECOND:
        min_size = 21;                          /* "2020-12-31T15:47:59Z" */
        rounding = 0.5 / (24.0 * 3600.0);       /* round to nearest second */
        break;

    case TIME_FORMAT_MILLI:
        min_size = 25;                          /* "2020-12-31T15:47:59.123Z" */
        rounding = 0.5 / (24.0 * 3600000.0);    /* round to nearest millisecond */
        break;

    default:
        return ASTRO_INVALID_PARAMETER;
    }

    /* Check for insufficient buffer size. */
    if (size < min_size)
        return ASTRO_BUFFER_TOO_SMALL;

    /* Perform rounding. */
    time.ut += rounding;

    /* Convert linear J2000 days to Gregorian UTC date/time. */
    utc = Astronomy_UtcFromTime(time);

    /* We require the year to be formatted as a 4-digit non-negative integer. */
    if (utc.year < 0 || utc.year > 9999)
        return ASTRO_BAD_TIME;

    /* Format the string. */
    switch (format)
    {
    case TIME_FORMAT_DAY:
        nprinted = snprintf(text, size, "%04d-%02d-%02d",
            utc.year, utc.month, utc.day);
        break;

    case TIME_FORMAT_MINUTE:
        nprinted = snprintf(text, size, "%04d-%02d-%02dT%02d:%02dZ",
            utc.year, utc.month, utc.day,
            utc.hour, utc.minute);
        break;

    case TIME_FORMAT_SECOND:
        nprinted = snprintf(text, size, "%04d-%02d-%02dT%02d:%02d:%02.0lfZ",
            utc.year, utc.month, utc.day,
            utc.hour, utc.minute, floor(utc.second));
        break;

    case TIME_FORMAT_MILLI:
        nprinted = snprintf(text, size, "%04d-%02d-%02dT%02d:%02d:%06.3lfZ",
            utc.year, utc.month, utc.day,
            utc.hour, utc.minute, floor(1000.0 * utc.second) / 1000.0);
        break;

    default:
        /* We should have already failed for any unknown 'format' value. */
        return ASTRO_INTERNAL_ERROR;
    }

    if (nprinted < 0)
        return ASTRO_INTERNAL_ERROR;    /* should not be possible for snprintf to return a negative number */

    if (1+(int)nprinted != min_size)
        return ASTRO_INTERNAL_ERROR;    /* there must be a bug calculating min_size or formatting the string */

    return ASTRO_SUCCESS;
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

static void iau2000b(astro_time_t *time)
{
    /* Adapted from the NOVAS C 3.1 function of the same name. */

    struct row_t
    {
        int nals[5];
        double cls[6];
    };

    static const struct row_t row[77] =
    {

        { {  0,  0,  0,  0,  1 }, {   -172064161,      -174666,        33386,     92052331,         9086,        15377 } },
        { {  0,  0,  2, -2,  2 }, {    -13170906,        -1675,       -13696,      5730336,        -3015,        -4587 } },
        { {  0,  0,  2,  0,  2 }, {     -2276413,         -234,         2796,       978459,         -485,         1374 } },
        { {  0,  0,  0,  0,  2 }, {      2074554,          207,         -698,      -897492,          470,         -291 } },
        { {  0,  1,  0,  0,  0 }, {      1475877,        -3633,        11817,        73871,         -184,        -1924 } },
        { {  0,  1,  2, -2,  2 }, {      -516821,         1226,         -524,       224386,         -677,         -174 } },
        { {  1,  0,  0,  0,  0 }, {       711159,           73,         -872,        -6750,            0,          358 } },
        { {  0,  0,  2,  0,  1 }, {      -387298,         -367,          380,       200728,           18,          318 } },
        { {  1,  0,  2,  0,  2 }, {      -301461,          -36,          816,       129025,          -63,          367 } },
        { {  0, -1,  2, -2,  2 }, {       215829,         -494,          111,       -95929,          299,          132 } },
        { {  0,  0,  2, -2,  1 }, {       128227,          137,          181,       -68982,           -9,           39 } },
        { { -1,  0,  2,  0,  2 }, {       123457,           11,           19,       -53311,           32,           -4 } },
        { { -1,  0,  0,  2,  0 }, {       156994,           10,         -168,        -1235,            0,           82 } },
        { {  1,  0,  0,  0,  1 }, {        63110,           63,           27,       -33228,            0,           -9 } },
        { { -1,  0,  0,  0,  1 }, {       -57976,          -63,         -189,        31429,            0,          -75 } },
        { { -1,  0,  2,  2,  2 }, {       -59641,          -11,          149,        25543,          -11,           66 } },
        { {  1,  0,  2,  0,  1 }, {       -51613,          -42,          129,        26366,            0,           78 } },
        { { -2,  0,  2,  0,  1 }, {        45893,           50,           31,       -24236,          -10,           20 } },
        { {  0,  0,  0,  2,  0 }, {        63384,           11,         -150,        -1220,            0,           29 } },
        { {  0,  0,  2,  2,  2 }, {       -38571,           -1,          158,        16452,          -11,           68 } },
        { {  0, -2,  2, -2,  2 }, {        32481,            0,            0,       -13870,            0,            0 } },
        { { -2,  0,  0,  2,  0 }, {       -47722,            0,          -18,          477,            0,          -25 } },
        { {  2,  0,  2,  0,  2 }, {       -31046,           -1,          131,        13238,          -11,           59 } },
        { {  1,  0,  2, -2,  2 }, {        28593,            0,           -1,       -12338,           10,           -3 } },
        { { -1,  0,  2,  0,  1 }, {        20441,           21,           10,       -10758,            0,           -3 } },
        { {  2,  0,  0,  0,  0 }, {        29243,            0,          -74,         -609,            0,           13 } },
        { {  0,  0,  2,  0,  0 }, {        25887,            0,          -66,         -550,            0,           11 } },
        { {  0,  1,  0,  0,  1 }, {       -14053,          -25,           79,         8551,           -2,          -45 } },
        { { -1,  0,  0,  2,  1 }, {        15164,           10,           11,        -8001,            0,           -1 } },
        { {  0,  2,  2, -2,  2 }, {       -15794,           72,          -16,         6850,          -42,           -5 } },
        { {  0,  0, -2,  2,  0 }, {        21783,            0,           13,         -167,            0,           13 } },
        { {  1,  0,  0, -2,  1 }, {       -12873,          -10,          -37,         6953,            0,          -14 } },
        { {  0, -1,  0,  0,  1 }, {       -12654,           11,           63,         6415,            0,           26 } },
        { { -1,  0,  2,  2,  1 }, {       -10204,            0,           25,         5222,            0,           15 } },
        { {  0,  2,  0,  0,  0 }, {        16707,          -85,          -10,          168,           -1,           10 } },
        { {  1,  0,  2,  2,  2 }, {        -7691,            0,           44,         3268,            0,           19 } },
        { { -2,  0,  2,  0,  0 }, {       -11024,            0,          -14,          104,            0,            2 } },
        { {  0,  1,  2,  0,  2 }, {         7566,          -21,          -11,        -3250,            0,           -5 } },
        { {  0,  0,  2,  2,  1 }, {        -6637,          -11,           25,         3353,            0,           14 } },
        { {  0, -1,  2,  0,  2 }, {        -7141,           21,            8,         3070,            0,            4 } },
        { {  0,  0,  0,  2,  1 }, {        -6302,          -11,            2,         3272,            0,            4 } },
        { {  1,  0,  2, -2,  1 }, {         5800,           10,            2,        -3045,            0,           -1 } },
        { {  2,  0,  2, -2,  2 }, {         6443,            0,           -7,        -2768,            0,           -4 } },
        { { -2,  0,  0,  2,  1 }, {        -5774,          -11,          -15,         3041,            0,           -5 } },
        { {  2,  0,  2,  0,  1 }, {        -5350,            0,           21,         2695,            0,           12 } },
        { {  0, -1,  2, -2,  1 }, {        -4752,          -11,           -3,         2719,            0,           -3 } },
        { {  0,  0,  0, -2,  1 }, {        -4940,          -11,          -21,         2720,            0,           -9 } },
        { { -1, -1,  0,  2,  0 }, {         7350,            0,           -8,          -51,            0,            4 } },
        { {  2,  0,  0, -2,  1 }, {         4065,            0,            6,        -2206,            0,            1 } },
        { {  1,  0,  0,  2,  0 }, {         6579,            0,          -24,         -199,            0,            2 } },
        { {  0,  1,  2, -2,  1 }, {         3579,            0,            5,        -1900,            0,            1 } },
        { {  1, -1,  0,  0,  0 }, {         4725,            0,           -6,          -41,            0,            3 } },
        { { -2,  0,  2,  0,  2 }, {        -3075,            0,           -2,         1313,            0,           -1 } },
        { {  3,  0,  2,  0,  2 }, {        -2904,            0,           15,         1233,            0,            7 } },
        { {  0, -1,  0,  2,  0 }, {         4348,            0,          -10,          -81,            0,            2 } },
        { {  1, -1,  2,  0,  2 }, {        -2878,            0,            8,         1232,            0,            4 } },
        { {  0,  0,  0,  1,  0 }, {        -4230,            0,            5,          -20,            0,           -2 } },
        { { -1, -1,  2,  2,  2 }, {        -2819,            0,            7,         1207,            0,            3 } },
        { { -1,  0,  2,  0,  0 }, {        -4056,            0,            5,           40,            0,           -2 } },
        { {  0, -1,  2,  2,  2 }, {        -2647,            0,           11,         1129,            0,            5 } },
        { { -2,  0,  0,  0,  1 }, {        -2294,            0,          -10,         1266,            0,           -4 } },
        { {  1,  1,  2,  0,  2 }, {         2481,            0,           -7,        -1062,            0,           -3 } },
        { {  2,  0,  0,  0,  1 }, {         2179,            0,           -2,        -1129,            0,           -2 } },
        { { -1,  1,  0,  1,  0 }, {         3276,            0,            1,           -9,            0,            0 } },
        { {  1,  1,  0,  0,  0 }, {        -3389,            0,            5,           35,            0,           -2 } },
        { {  1,  0,  2,  0,  0 }, {         3339,            0,          -13,         -107,            0,            1 } },
        { { -1,  0,  2, -2,  1 }, {        -1987,            0,           -6,         1073,            0,           -2 } },
        { {  1,  0,  0,  0,  2 }, {        -1981,            0,            0,          854,            0,            0 } },
        { { -1,  0,  0,  1,  0 }, {         4026,            0,         -353,         -553,            0,         -139 } },
        { {  0,  0,  2,  1,  2 }, {         1660,            0,           -5,         -710,            0,           -2 } },
        { { -1,  0,  2,  4,  2 }, {        -1521,            0,            9,          647,            0,            4 } },
        { { -1,  1,  0,  1,  1 }, {         1314,            0,            0,         -700,            0,            0 } },
        { {  0, -2,  2, -2,  1 }, {        -1283,            0,            0,          672,            0,            0 } },
        { {  1,  0,  2,  2,  1 }, {        -1331,            0,            8,          663,            0,            4 } },
        { { -2,  0,  2,  2,  2 }, {         1383,            0,           -2,         -594,            0,           -2 } },
        { { -1,  0,  0,  0,  2 }, {         1405,            0,            4,         -610,            0,            2 } },
        { {  1,  1,  2, -2,  2 }, {         1290,            0,            0,         -556,            0,            0 } }

    };

    double t, el, elp, f, d, om, arg, dp, de, sarg, carg;
    int i;

    if (isnan(time->psi))
    {
        t = time->tt / 36525;
        el  = fmod(485868.249036 + t * 1717915923.2178, ASEC360) * ASEC2RAD;
        elp = fmod(1287104.79305 + t * 129596581.0481,  ASEC360) * ASEC2RAD;
        f   = fmod(335779.526232 + t * 1739527262.8478, ASEC360) * ASEC2RAD;
        d   = fmod(1072260.70369 + t * 1602961601.2090, ASEC360) * ASEC2RAD;
        om  = fmod(450160.398036 - t * 6962890.5431,    ASEC360) * ASEC2RAD;
        dp = 0;
        de = 0;
        for (i=76; i >= 0; --i)
        {
            arg = fmod((row[i].nals[0]*el + row[i].nals[1]*elp + row[i].nals[2]*f + row[i].nals[3]*d + row[i].nals[4]*om), PI2);
            sarg = sin(arg);
            carg = cos(arg);
            dp += (row[i].cls[0] + row[i].cls[1]*t) * sarg + row[i].cls[2]*carg;
            de += (row[i].cls[3] + row[i].cls[4]*t) * carg + row[i].cls[5]*sarg;
        }

        time->psi = -0.000135 + (dp * 1.0e-7);
        time->eps = +0.000388 + (de * 1.0e-7);
    }
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

static earth_tilt_t e_tilt(astro_time_t *time)
{
    earth_tilt_t et;

    iau2000b(time);
    et.dpsi = time->psi;
    et.deps = time->eps;
    et.mobl = mean_obliq(time->tt);
    et.tobl = et.mobl + (et.deps / 3600.0);
    et.tt = time->tt;
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


static astro_rotation_t precession_rot(double tt1, double tt2)
{
    astro_rotation_t rotation;
    double xx, yx, zx, xy, yy, zy, xz, yz, zz;
    double t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;
    double eps0 = 84381.406;

    if ((tt1 != 0.0) && (tt2 != 0.0))
        FatalError("precession_rot: one of (tt1, tt2) must be zero.");

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
        rotation.rot[0][0] = xx;
        rotation.rot[0][1] = yx;
        rotation.rot[0][2] = zx;
        rotation.rot[1][0] = xy;
        rotation.rot[1][1] = yy;
        rotation.rot[1][2] = zy;
        rotation.rot[2][0] = xz;
        rotation.rot[2][1] = yz;
        rotation.rot[2][2] = zz;
    }
    else
    {
        /* Perform rotation from J2000.0 to other epoch. */
        rotation.rot[0][0] = xx;
        rotation.rot[0][1] = xy;
        rotation.rot[0][2] = xz;
        rotation.rot[1][0] = yx;
        rotation.rot[1][1] = yy;
        rotation.rot[1][2] = yz;
        rotation.rot[2][0] = zx;
        rotation.rot[2][1] = zy;
        rotation.rot[2][2] = zz;
    }

    rotation.status = ASTRO_SUCCESS;
    return rotation;
}


static void precession(double tt1, const double pos1[3], double tt2, double pos2[3])
{
    astro_rotation_t r = precession_rot(tt1, tt2);
    pos2[0] = r.rot[0][0]*pos1[0] + r.rot[1][0]*pos1[1] + r.rot[2][0]*pos1[2];
    pos2[1] = r.rot[0][1]*pos1[0] + r.rot[1][1]*pos1[1] + r.rot[2][1]*pos1[2];
    pos2[2] = r.rot[0][2]*pos1[0] + r.rot[1][2]*pos1[1] + r.rot[2][2]*pos1[2];
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


static astro_rotation_t nutation_rot(astro_time_t *time, int direction)
{
    astro_rotation_t rotation;
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
        rotation.rot[0][0] = xx;
        rotation.rot[0][1] = xy;
        rotation.rot[0][2] = xz;
        rotation.rot[1][0] = yx;
        rotation.rot[1][1] = yy;
        rotation.rot[1][2] = yz;
        rotation.rot[2][0] = zx;
        rotation.rot[2][1] = zy;
        rotation.rot[2][2] = zz;
    }
    else
    {
        /* inverse rotation */
        rotation.rot[0][0] = xx;
        rotation.rot[0][1] = yx;
        rotation.rot[0][2] = zx;
        rotation.rot[1][0] = xy;
        rotation.rot[1][1] = yy;
        rotation.rot[1][2] = zy;
        rotation.rot[2][0] = xz;
        rotation.rot[2][1] = yz;
        rotation.rot[2][2] = zz;
    }

    rotation.status = ASTRO_SUCCESS;
    return rotation;
}

static void nutation(astro_time_t *time, int direction, const double inpos[3], double outpos[3])
{
    astro_rotation_t r = nutation_rot(time, direction);
    outpos[0] = r.rot[0][0]*inpos[0] + r.rot[1][0]*inpos[1] + r.rot[2][0]*inpos[2];
    outpos[1] = r.rot[0][1]*inpos[0] + r.rot[1][1]*inpos[1] + r.rot[2][1]*inpos[2];
    outpos[2] = r.rot[0][2]*inpos[0] + r.rot[1][2]*inpos[1] + r.rot[2][2]*inpos[2];
}

static double era(double ut)        /* Earth Rotation Angle */
{
    double thet1 = 0.7790572732640 + 0.00273781191135448 * ut;
    double thet3 = fmod(ut, 1.0);
    double theta = 360.0 * fmod(thet1 + thet3, 1.0);
    if (theta < 0.0)
        theta += 360.0;

    return theta;
}

static double sidereal_time(astro_time_t *time)
{
    double t = time->tt / 36525.0;
    double eqeq = 15.0 * e_tilt(time).ee;    /* Replace with eqeq=0 to get GMST instead of GAST (if we ever need it) */
    double theta = era(time->ut);
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

static void terra(astro_observer_t observer, double st, double pos[3])
{
    double df2 = EARTH_FLATTENING * EARTH_FLATTENING;
    double phi = observer.latitude * DEG2RAD;
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double c = 1.0 / sqrt(cosphi*cosphi + df2*sinphi*sinphi);
    double s = df2 * c;
    double ht_km = observer.height / 1000.0;
    double ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
    double ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
    double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
    double sinst = sin(stlocl);
    double cosst = cos(stlocl);

    pos[0] = ach * cosphi * cosst / KM_PER_AU;
    pos[1] = ach * cosphi * sinst / KM_PER_AU;
    pos[2] = ash * sinphi / KM_PER_AU;

#if 0
    /* If we ever need to calculate the observer's velocity vector, here is how NOVAS C 3.1 does it... */
    static const double ANGVEL = 7.2921150e-5;
    vel[0] = -ANGVEL * ach * cosphi * sinst * 86400.0;
    vel[1] = +ANGVEL * ach * cosphi * cosst * 86400.0;
    vel[2] = 0.0;
#endif
}

static void geo_pos(astro_time_t *time, astro_observer_t observer, double outpos[3])
{
    double gast, pos1[3], pos2[3];

    gast = sidereal_time(time);
    terra(observer, gast, pos1);
    nutation(time, -1, pos1, pos2);
    precession(time->tt, pos2, 0.0, outpos);
}

static void spin(double angle, const double pos1[3], double vec2[3])
{
    double angr = angle * DEG2RAD;
    double cosang = cos(angr);
    double sinang = sin(angr);
    vec2[0] = +cosang*pos1[0] + sinang*pos1[1];
    vec2[1] = -sinang*pos1[0] + cosang*pos1[1];
    vec2[2] = pos1[2];
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
            case 1:  ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2:  ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3:  ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            default: ARG=D;  MAX=6; FAC=1.0;                       break;
        }
        CO(0,I) = 1.0;
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

int _CalcMoonCount;     /* Undocumented global for performance tuning. */

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

    AddSol(ctx,    13.9020,    14.0600,    -0.0010,     0.2607, 0, 0, 0, 4);
    AddSol(ctx,     0.4030,    -4.0100,     0.3940,     0.0023, 0, 0, 0, 3);
    AddSol(ctx,  2369.9120,  2373.3600,     0.6010,    28.2333, 0, 0, 0, 2);
    AddSol(ctx,  -125.1540,  -112.7900,    -0.7250,    -0.9781, 0, 0, 0, 1);
    AddSol(ctx,     1.9790,     6.9800,    -0.4450,     0.0433, 1, 0, 0, 4);
    AddSol(ctx,   191.9530,   192.7200,     0.0290,     3.0861, 1, 0, 0, 2);
    AddSol(ctx,    -8.4660,   -13.5100,     0.4550,    -0.1093, 1, 0, 0, 1);
    AddSol(ctx, 22639.5000, 22609.0700,     0.0790,   186.5398, 1, 0, 0, 0);
    AddSol(ctx,    18.6090,     3.5900,    -0.0940,     0.0118, 1, 0, 0,-1);
    AddSol(ctx, -4586.4650, -4578.1300,    -0.0770,    34.3117, 1, 0, 0,-2);
    AddSol(ctx,     3.2150,     5.4400,     0.1920,    -0.0386, 1, 0, 0,-3);
    AddSol(ctx,   -38.4280,   -38.6400,     0.0010,     0.6008, 1, 0, 0,-4);
    AddSol(ctx,    -0.3930,    -1.4300,    -0.0920,     0.0086, 1, 0, 0,-6);
    AddSol(ctx,    -0.2890,    -1.5900,     0.1230,    -0.0053, 0, 1, 0, 4);
    AddSol(ctx,   -24.4200,   -25.1000,     0.0400,    -0.3000, 0, 1, 0, 2);
    AddSol(ctx,    18.0230,    17.9300,     0.0070,     0.1494, 0, 1, 0, 1);
    AddSol(ctx,  -668.1460,  -126.9800,    -1.3020,    -0.3997, 0, 1, 0, 0);
    AddSol(ctx,     0.5600,     0.3200,    -0.0010,    -0.0037, 0, 1, 0,-1);
    AddSol(ctx,  -165.1450,  -165.0600,     0.0540,     1.9178, 0, 1, 0,-2);
    AddSol(ctx,    -1.8770,    -6.4600,    -0.4160,     0.0339, 0, 1, 0,-4);
    AddSol(ctx,     0.2130,     1.0200,    -0.0740,     0.0054, 2, 0, 0, 4);
    AddSol(ctx,    14.3870,    14.7800,    -0.0170,     0.2833, 2, 0, 0, 2);
    AddSol(ctx,    -0.5860,    -1.2000,     0.0540,    -0.0100, 2, 0, 0, 1);
    AddSol(ctx,   769.0160,   767.9600,     0.1070,    10.1657, 2, 0, 0, 0);
    AddSol(ctx,     1.7500,     2.0100,    -0.0180,     0.0155, 2, 0, 0,-1);
    AddSol(ctx,  -211.6560,  -152.5300,     5.6790,    -0.3039, 2, 0, 0,-2);
    AddSol(ctx,     1.2250,     0.9100,    -0.0300,    -0.0088, 2, 0, 0,-3);
    AddSol(ctx,   -30.7730,   -34.0700,    -0.3080,     0.3722, 2, 0, 0,-4);
    AddSol(ctx,    -0.5700,    -1.4000,    -0.0740,     0.0109, 2, 0, 0,-6);
    AddSol(ctx,    -2.9210,   -11.7500,     0.7870,    -0.0484, 1, 1, 0, 2);
    AddSol(ctx,     1.2670,     1.5200,    -0.0220,     0.0164, 1, 1, 0, 1);
    AddSol(ctx,  -109.6730,  -115.1800,     0.4610,    -0.9490, 1, 1, 0, 0);
    AddSol(ctx,  -205.9620,  -182.3600,     2.0560,     1.4437, 1, 1, 0,-2);
    AddSol(ctx,     0.2330,     0.3600,     0.0120,    -0.0025, 1, 1, 0,-3);
    AddSol(ctx,    -4.3910,    -9.6600,    -0.4710,     0.0673, 1, 1, 0,-4);
    AddSol(ctx,     0.2830,     1.5300,    -0.1110,     0.0060, 1,-1, 0, 4);
    AddSol(ctx,    14.5770,    31.7000,    -1.5400,     0.2302, 1,-1, 0, 2);
    AddSol(ctx,   147.6870,   138.7600,     0.6790,     1.1528, 1,-1, 0, 0);
    AddSol(ctx,    -1.0890,     0.5500,     0.0210,     0.0000, 1,-1, 0,-1);
    AddSol(ctx,    28.4750,    23.5900,    -0.4430,    -0.2257, 1,-1, 0,-2);
    AddSol(ctx,    -0.2760,    -0.3800,    -0.0060,    -0.0036, 1,-1, 0,-3);
    AddSol(ctx,     0.6360,     2.2700,     0.1460,    -0.0102, 1,-1, 0,-4);
    AddSol(ctx,    -0.1890,    -1.6800,     0.1310,    -0.0028, 0, 2, 0, 2);
    AddSol(ctx,    -7.4860,    -0.6600,    -0.0370,    -0.0086, 0, 2, 0, 0);
    AddSol(ctx,    -8.0960,   -16.3500,    -0.7400,     0.0918, 0, 2, 0,-2);
    AddSol(ctx,    -5.7410,    -0.0400,     0.0000,    -0.0009, 0, 0, 2, 2);
    AddSol(ctx,     0.2550,     0.0000,     0.0000,     0.0000, 0, 0, 2, 1);
    AddSol(ctx,  -411.6080,    -0.2000,     0.0000,    -0.0124, 0, 0, 2, 0);
    AddSol(ctx,     0.5840,     0.8400,     0.0000,     0.0071, 0, 0, 2,-1);
    AddSol(ctx,   -55.1730,   -52.1400,     0.0000,    -0.1052, 0, 0, 2,-2);
    AddSol(ctx,     0.2540,     0.2500,     0.0000,    -0.0017, 0, 0, 2,-3);
    AddSol(ctx,     0.0250,    -1.6700,     0.0000,     0.0031, 0, 0, 2,-4);
    AddSol(ctx,     1.0600,     2.9600,    -0.1660,     0.0243, 3, 0, 0, 2);
    AddSol(ctx,    36.1240,    50.6400,    -1.3000,     0.6215, 3, 0, 0, 0);
    AddSol(ctx,   -13.1930,   -16.4000,     0.2580,    -0.1187, 3, 0, 0,-2);
    AddSol(ctx,    -1.1870,    -0.7400,     0.0420,     0.0074, 3, 0, 0,-4);
    AddSol(ctx,    -0.2930,    -0.3100,    -0.0020,     0.0046, 3, 0, 0,-6);
    AddSol(ctx,    -0.2900,    -1.4500,     0.1160,    -0.0051, 2, 1, 0, 2);
    AddSol(ctx,    -7.6490,   -10.5600,     0.2590,    -0.1038, 2, 1, 0, 0);
    AddSol(ctx,    -8.6270,    -7.5900,     0.0780,    -0.0192, 2, 1, 0,-2);
    AddSol(ctx,    -2.7400,    -2.5400,     0.0220,     0.0324, 2, 1, 0,-4);
    AddSol(ctx,     1.1810,     3.3200,    -0.2120,     0.0213, 2,-1, 0, 2);
    AddSol(ctx,     9.7030,    11.6700,    -0.1510,     0.1268, 2,-1, 0, 0);
    AddSol(ctx,    -0.3520,    -0.3700,     0.0010,    -0.0028, 2,-1, 0,-1);
    AddSol(ctx,    -2.4940,    -1.1700,    -0.0030,    -0.0017, 2,-1, 0,-2);
    AddSol(ctx,     0.3600,     0.2000,    -0.0120,    -0.0043, 2,-1, 0,-4);
    AddSol(ctx,    -1.1670,    -1.2500,     0.0080,    -0.0106, 1, 2, 0, 0);
    AddSol(ctx,    -7.4120,    -6.1200,     0.1170,     0.0484, 1, 2, 0,-2);
    AddSol(ctx,    -0.3110,    -0.6500,    -0.0320,     0.0044, 1, 2, 0,-4);
    AddSol(ctx,     0.7570,     1.8200,    -0.1050,     0.0112, 1,-2, 0, 2);
    AddSol(ctx,     2.5800,     2.3200,     0.0270,     0.0196, 1,-2, 0, 0);
    AddSol(ctx,     2.5330,     2.4000,    -0.0140,    -0.0212, 1,-2, 0,-2);
    AddSol(ctx,    -0.3440,    -0.5700,    -0.0250,     0.0036, 0, 3, 0,-2);
    AddSol(ctx,    -0.9920,    -0.0200,     0.0000,     0.0000, 1, 0, 2, 2);
    AddSol(ctx,   -45.0990,    -0.0200,     0.0000,    -0.0010, 1, 0, 2, 0);
    AddSol(ctx,    -0.1790,    -9.5200,     0.0000,    -0.0833, 1, 0, 2,-2);
    AddSol(ctx,    -0.3010,    -0.3300,     0.0000,     0.0014, 1, 0, 2,-4);
    AddSol(ctx,    -6.3820,    -3.3700,     0.0000,    -0.0481, 1, 0,-2, 2);
    AddSol(ctx,    39.5280,    85.1300,     0.0000,    -0.7136, 1, 0,-2, 0);
    AddSol(ctx,     9.3660,     0.7100,     0.0000,    -0.0112, 1, 0,-2,-2);
    AddSol(ctx,     0.2020,     0.0200,     0.0000,     0.0000, 1, 0,-2,-4);
    AddSol(ctx,     0.4150,     0.1000,     0.0000,     0.0013, 0, 1, 2, 0);
    AddSol(ctx,    -2.1520,    -2.2600,     0.0000,    -0.0066, 0, 1, 2,-2);
    AddSol(ctx,    -1.4400,    -1.3000,     0.0000,     0.0014, 0, 1,-2, 2);
    AddSol(ctx,     0.3840,    -0.0400,     0.0000,     0.0000, 0, 1,-2,-2);
    AddSol(ctx,     1.9380,     3.6000,    -0.1450,     0.0401, 4, 0, 0, 0);
    AddSol(ctx,    -0.9520,    -1.5800,     0.0520,    -0.0130, 4, 0, 0,-2);
    AddSol(ctx,    -0.5510,    -0.9400,     0.0320,    -0.0097, 3, 1, 0, 0);
    AddSol(ctx,    -0.4820,    -0.5700,     0.0050,    -0.0045, 3, 1, 0,-2);
    AddSol(ctx,     0.6810,     0.9600,    -0.0260,     0.0115, 3,-1, 0, 0);
    AddSol(ctx,    -0.2970,    -0.2700,     0.0020,    -0.0009, 2, 2, 0,-2);
    AddSol(ctx,     0.2540,     0.2100,    -0.0030,     0.0000, 2,-2, 0,-2);
    AddSol(ctx,    -0.2500,    -0.2200,     0.0040,     0.0014, 1, 3, 0,-2);
    AddSol(ctx,    -3.9960,     0.0000,     0.0000,     0.0004, 2, 0, 2, 0);
    AddSol(ctx,     0.5570,    -0.7500,     0.0000,    -0.0090, 2, 0, 2,-2);
    AddSol(ctx,    -0.4590,    -0.3800,     0.0000,    -0.0053, 2, 0,-2, 2);
    AddSol(ctx,    -1.2980,     0.7400,     0.0000,     0.0004, 2, 0,-2, 0);
    AddSol(ctx,     0.5380,     1.1400,     0.0000,    -0.0141, 2, 0,-2,-2);
    AddSol(ctx,     0.2630,     0.0200,     0.0000,     0.0000, 1, 1, 2, 0);
    AddSol(ctx,     0.4260,     0.0700,     0.0000,    -0.0006, 1, 1,-2,-2);
    AddSol(ctx,    -0.3040,     0.0300,     0.0000,     0.0003, 1,-1, 2, 0);
    AddSol(ctx,    -0.3720,    -0.1900,     0.0000,    -0.0027, 1,-1,-2, 2);
    AddSol(ctx,     0.4180,     0.0000,     0.0000,     0.0000, 0, 0, 4, 0);
    AddSol(ctx,    -0.3300,    -0.0400,     0.0000,     0.0000, 3, 0, 2, 0);

    SolarN(ctx);
    Planetary(ctx);
    S = F + DS/ARC;

    lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*sin(S)-6.24*sin(3*S) + N;

    *geo_eclip_lon = PI2 * Frac((L0+DLAM/ARC) / PI2);
    *geo_eclip_lat = lat_seconds * (DEG2RAD / 3600.0);
    *distance_au = (ARC * EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI);
    ++_CalcMoonCount;
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
 * This algorithm is based on Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
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

static const vsop_term_t vsop_lat_Mercury_0[] =
{
    { 4.40250710144, 0.00000000000, 0.00000000000 },
    { 0.40989414977, 1.48302034195, 26087.90314157420 },
    { 0.05046294200, 4.47785489551, 52175.80628314840 },
    { 0.00855346844, 1.16520322459, 78263.70942472259 },
    { 0.00165590362, 4.11969163423, 104351.61256629678 },
    { 0.00034561897, 0.77930768443, 130439.51570787099 },
    { 0.00007583476, 3.71348404924, 156527.41884944518 }
};

static const vsop_term_t vsop_lat_Mercury_1[] =
{
    { 26087.90313685529, 0.00000000000, 0.00000000000 },
    { 0.01131199811, 6.21874197797, 26087.90314157420 },
    { 0.00292242298, 3.04449355541, 52175.80628314840 },
    { 0.00075775081, 6.08568821653, 78263.70942472259 },
    { 0.00019676525, 2.80965111777, 104351.61256629678 }
};

static const vsop_series_t vsop_lat_Mercury[] =
{
    { 7, vsop_lat_Mercury_0 },
    { 5, vsop_lat_Mercury_1 }
};

static const vsop_term_t vsop_lon_Mercury_0[] =
{
    { 0.11737528961, 1.98357498767, 26087.90314157420 },
    { 0.02388076996, 5.03738959686, 52175.80628314840 },
    { 0.01222839532, 3.14159265359, 0.00000000000 },
    { 0.00543251810, 1.79644363964, 78263.70942472259 },
    { 0.00129778770, 4.83232503958, 104351.61256629678 },
    { 0.00031866927, 1.58088495658, 130439.51570787099 },
    { 0.00007963301, 4.60972126127, 156527.41884944518 }
};

static const vsop_term_t vsop_lon_Mercury_1[] =
{
    { 0.00274646065, 3.95008450011, 26087.90314157420 },
    { 0.00099737713, 3.14159265359, 0.00000000000 }
};

static const vsop_series_t vsop_lon_Mercury[] =
{
    { 7, vsop_lon_Mercury_0 },
    { 2, vsop_lon_Mercury_1 }
};

static const vsop_term_t vsop_rad_Mercury_0[] =
{
    { 0.39528271651, 0.00000000000, 0.00000000000 },
    { 0.07834131818, 6.19233722598, 26087.90314157420 },
    { 0.00795525558, 2.95989690104, 52175.80628314840 },
    { 0.00121281764, 6.01064153797, 78263.70942472259 },
    { 0.00021921969, 2.77820093972, 104351.61256629678 },
    { 0.00004354065, 5.82894543774, 130439.51570787099 }
};

static const vsop_term_t vsop_rad_Mercury_1[] =
{
    { 0.00217347740, 4.65617158665, 26087.90314157420 },
    { 0.00044141826, 1.42385544001, 52175.80628314840 }
};

static const vsop_series_t vsop_rad_Mercury[] =
{
    { 6, vsop_rad_Mercury_0 },
    { 2, vsop_rad_Mercury_1 }
};

;
static const vsop_term_t vsop_lat_Venus_0[] =
{
    { 3.17614666774, 0.00000000000, 0.00000000000 },
    { 0.01353968419, 5.59313319619, 10213.28554621100 },
    { 0.00089891645, 5.30650047764, 20426.57109242200 },
    { 0.00005477194, 4.41630661466, 7860.41939243920 },
    { 0.00003455741, 2.69964447820, 11790.62908865880 },
    { 0.00002372061, 2.99377542079, 3930.20969621960 },
    { 0.00001317168, 5.18668228402, 26.29831979980 },
    { 0.00001664146, 4.25018630147, 1577.34354244780 },
    { 0.00001438387, 4.15745084182, 9683.59458111640 },
    { 0.00001200521, 6.15357116043, 30639.85663863300 }
};

static const vsop_term_t vsop_lat_Venus_1[] =
{
    { 10213.28554621638, 0.00000000000, 0.00000000000 },
    { 0.00095617813, 2.46406511110, 10213.28554621100 },
    { 0.00007787201, 0.62478482220, 20426.57109242200 }
};

static const vsop_series_t vsop_lat_Venus[] =
{
    { 10, vsop_lat_Venus_0 },
    { 3, vsop_lat_Venus_1 }
};

static const vsop_term_t vsop_lon_Venus_0[] =
{
    { 0.05923638472, 0.26702775812, 10213.28554621100 },
    { 0.00040107978, 1.14737178112, 20426.57109242200 },
    { 0.00032814918, 3.14159265359, 0.00000000000 }
};

static const vsop_term_t vsop_lon_Venus_1[] =
{
    { 0.00287821243, 1.88964962838, 10213.28554621100 }
};

static const vsop_series_t vsop_lon_Venus[] =
{
    { 3, vsop_lon_Venus_0 },
    { 1, vsop_lon_Venus_1 }
};

static const vsop_term_t vsop_rad_Venus_0[] =
{
    { 0.72334820891, 0.00000000000, 0.00000000000 },
    { 0.00489824182, 4.02151831717, 10213.28554621100 },
    { 0.00001658058, 4.90206728031, 20426.57109242200 },
    { 0.00001378043, 1.12846591367, 11790.62908865880 },
    { 0.00001632096, 2.84548795207, 7860.41939243920 },
    { 0.00000498395, 2.58682193892, 9683.59458111640 },
    { 0.00000221985, 2.01346696541, 19367.18916223280 },
    { 0.00000237454, 2.55136053886, 15720.83878487840 }
};

static const vsop_term_t vsop_rad_Venus_1[] =
{
    { 0.00034551041, 0.89198706276, 10213.28554621100 }
};

static const vsop_series_t vsop_rad_Venus[] =
{
    { 8, vsop_rad_Venus_0 },
    { 1, vsop_rad_Venus_1 }
};

;
static const vsop_term_t vsop_lat_Earth_0[] =
{
    { 1.75347045673, 0.00000000000, 0.00000000000 },
    { 0.03341656453, 4.66925680415, 6283.07584999140 },
    { 0.00034894275, 4.62610242189, 12566.15169998280 },
    { 0.00003417572, 2.82886579754, 3.52311834900 },
    { 0.00003497056, 2.74411783405, 5753.38488489680 },
    { 0.00003135899, 3.62767041756, 77713.77146812050 },
    { 0.00002676218, 4.41808345438, 7860.41939243920 },
    { 0.00002342691, 6.13516214446, 3930.20969621960 },
    { 0.00001273165, 2.03709657878, 529.69096509460 },
    { 0.00001324294, 0.74246341673, 11506.76976979360 },
    { 0.00000901854, 2.04505446477, 26.29831979980 },
    { 0.00001199167, 1.10962946234, 1577.34354244780 },
    { 0.00000857223, 3.50849152283, 398.14900340820 },
    { 0.00000779786, 1.17882681962, 5223.69391980220 },
    { 0.00000990250, 5.23268072088, 5884.92684658320 },
    { 0.00000753141, 2.53339052847, 5507.55323866740 },
    { 0.00000505267, 4.58292599973, 18849.22754997420 },
    { 0.00000492392, 4.20505711826, 775.52261132400 },
    { 0.00000356672, 2.91954114478, 0.06731030280 },
    { 0.00000284125, 1.89869240932, 796.29800681640 },
    { 0.00000242879, 0.34481445893, 5486.77784317500 },
    { 0.00000317087, 5.84901948512, 11790.62908865880 },
    { 0.00000271112, 0.31486255375, 10977.07880469900 },
    { 0.00000206217, 4.80646631478, 2544.31441988340 },
    { 0.00000205478, 1.86953770281, 5573.14280143310 },
    { 0.00000202318, 2.45767790232, 6069.77675455340 },
    { 0.00000126225, 1.08295459501, 20.77539549240 },
    { 0.00000155516, 0.83306084617, 213.29909543800 }
};

static const vsop_term_t vsop_lat_Earth_1[] =
{
    { 6283.07584999140, 0.00000000000, 0.00000000000 },
    { 0.00206058863, 2.67823455808, 6283.07584999140 },
    { 0.00004303419, 2.63512233481, 12566.15169998280 }
};

static const vsop_term_t vsop_lat_Earth_2[] =
{
    { 0.00008721859, 1.07253635559, 6283.07584999140 }
};

static const vsop_series_t vsop_lat_Earth[] =
{
    { 28, vsop_lat_Earth_0 },
    { 3, vsop_lat_Earth_1 },
    { 1, vsop_lat_Earth_2 }
};

static const vsop_term_t vsop_lon_Earth_1[] =
{
    { 0.00227777722, 3.41376620530, 6283.07584999140 },
    { 0.00003805678, 3.37063423795, 12566.15169998280 }
};

static const vsop_series_t vsop_lon_Earth[] =
{
    { 0, NULL },
    { 2, vsop_lon_Earth_1 }
};

static const vsop_term_t vsop_rad_Earth_0[] =
{
    { 1.00013988784, 0.00000000000, 0.00000000000 },
    { 0.01670699632, 3.09846350258, 6283.07584999140 },
    { 0.00013956024, 3.05524609456, 12566.15169998280 },
    { 0.00003083720, 5.19846674381, 77713.77146812050 },
    { 0.00001628463, 1.17387558054, 5753.38488489680 },
    { 0.00001575572, 2.84685214877, 7860.41939243920 },
    { 0.00000924799, 5.45292236722, 11506.76976979360 },
    { 0.00000542439, 4.56409151453, 3930.20969621960 },
    { 0.00000472110, 3.66100022149, 5884.92684658320 },
    { 0.00000085831, 1.27079125277, 161000.68573767410 },
    { 0.00000057056, 2.01374292245, 83996.84731811189 },
    { 0.00000055736, 5.24159799170, 71430.69561812909 },
    { 0.00000174844, 3.01193636733, 18849.22754997420 },
    { 0.00000243181, 4.27349530790, 11790.62908865880 }
};

static const vsop_term_t vsop_rad_Earth_1[] =
{
    { 0.00103018607, 1.10748968172, 6283.07584999140 },
    { 0.00001721238, 1.06442300386, 12566.15169998280 }
};

static const vsop_term_t vsop_rad_Earth_2[] =
{
    { 0.00004359385, 5.78455133808, 6283.07584999140 }
};

static const vsop_series_t vsop_rad_Earth[] =
{
    { 14, vsop_rad_Earth_0 },
    { 2, vsop_rad_Earth_1 },
    { 1, vsop_rad_Earth_2 }
};

;
static const vsop_term_t vsop_lat_Mars_0[] =
{
    { 6.20347711581, 0.00000000000, 0.00000000000 },
    { 0.18656368093, 5.05037100270, 3340.61242669980 },
    { 0.01108216816, 5.40099836344, 6681.22485339960 },
    { 0.00091798406, 5.75478744667, 10021.83728009940 },
    { 0.00027744987, 5.97049513147, 3.52311834900 },
    { 0.00010610235, 2.93958560338, 2281.23049651060 },
    { 0.00012315897, 0.84956094002, 2810.92146160520 },
    { 0.00008926784, 4.15697846427, 0.01725365220 },
    { 0.00008715691, 6.11005153139, 13362.44970679920 },
    { 0.00006797556, 0.36462229657, 398.14900340820 },
    { 0.00007774872, 3.33968761376, 5621.84292321040 },
    { 0.00003575078, 1.66186505710, 2544.31441988340 },
    { 0.00004161108, 0.22814971327, 2942.46342329160 },
    { 0.00003075252, 0.85696614132, 191.44826611160 },
    { 0.00002628117, 0.64806124465, 3337.08930835080 },
    { 0.00002937546, 6.07893711402, 0.06731030280 },
    { 0.00002389414, 5.03896442664, 796.29800681640 },
    { 0.00002579844, 0.02996736156, 3344.13554504880 },
    { 0.00001528141, 1.14979301996, 6151.53388830500 },
    { 0.00001798806, 0.65634057445, 529.69096509460 },
    { 0.00001264357, 3.62275122593, 5092.15195811580 },
    { 0.00001286228, 3.06796065034, 2146.16541647520 },
    { 0.00001546404, 2.91579701718, 1751.53953141600 },
    { 0.00001024902, 3.69334099279, 8962.45534991020 },
    { 0.00000891566, 0.18293837498, 16703.06213349900 },
    { 0.00000858759, 2.40093811940, 2914.01423582380 },
    { 0.00000832715, 2.46418619474, 3340.59517304760 },
    { 0.00000832720, 4.49495782139, 3340.62968035200 },
    { 0.00000712902, 3.66335473479, 1059.38193018920 },
    { 0.00000748723, 3.82248614017, 155.42039943420 },
    { 0.00000723861, 0.67497311481, 3738.76143010800 },
    { 0.00000635548, 2.92182225127, 8432.76438481560 },
    { 0.00000655162, 0.48864064125, 3127.31333126180 },
    { 0.00000550474, 3.81001042328, 0.98032106820 },
    { 0.00000552750, 4.47479317037, 1748.01641306700 },
    { 0.00000425966, 0.55364317304, 6283.07584999140 },
    { 0.00000415131, 0.49662285038, 213.29909543800 },
    { 0.00000472167, 3.62547124025, 1194.44701022460 },
    { 0.00000306551, 0.38052848348, 6684.74797174860 },
    { 0.00000312141, 0.99853944405, 6677.70173505060 },
    { 0.00000293198, 4.22131299634, 20.77539549240 },
    { 0.00000302375, 4.48618007156, 3532.06069281140 },
    { 0.00000274027, 0.54222167059, 3340.54511639700 },
    { 0.00000281079, 5.88163521788, 1349.86740965880 },
    { 0.00000231183, 1.28242156993, 3870.30339179440 },
    { 0.00000283602, 5.76885434940, 3149.16416058820 },
    { 0.00000236117, 5.75503217933, 3333.49887969900 },
    { 0.00000274033, 0.13372524985, 3340.67973700260 },
    { 0.00000299395, 2.78323740866, 6254.62666252360 }
};

static const vsop_term_t vsop_lat_Mars_1[] =
{
    { 3340.61242700512, 0.00000000000, 0.00000000000 },
    { 0.01457554523, 3.60433733236, 3340.61242669980 },
    { 0.00168414711, 3.92318567804, 6681.22485339960 },
    { 0.00020622975, 4.26108844583, 10021.83728009940 },
    { 0.00003452392, 4.73210393190, 3.52311834900 },
    { 0.00002586332, 4.60670058555, 13362.44970679920 },
    { 0.00000841535, 4.45864030426, 2281.23049651060 }
};

static const vsop_term_t vsop_lat_Mars_2[] =
{
    { 0.00058152577, 2.04961712429, 3340.61242669980 },
    { 0.00013459579, 2.45738706163, 6681.22485339960 }
};

static const vsop_series_t vsop_lat_Mars[] =
{
    { 49, vsop_lat_Mars_0 },
    { 7, vsop_lat_Mars_1 },
    { 2, vsop_lat_Mars_2 }
};

static const vsop_term_t vsop_lon_Mars_0[] =
{
    { 0.03197134986, 3.76832042431, 3340.61242669980 },
    { 0.00298033234, 4.10616996305, 6681.22485339960 },
    { 0.00289104742, 0.00000000000, 0.00000000000 },
    { 0.00031365539, 4.44651053090, 10021.83728009940 },
    { 0.00003484100, 4.78812549260, 13362.44970679920 }
};

static const vsop_term_t vsop_lon_Mars_1[] =
{
    { 0.00217310991, 6.04472194776, 3340.61242669980 },
    { 0.00020976948, 3.14159265359, 0.00000000000 },
    { 0.00012834709, 1.60810667915, 6681.22485339960 }
};

static const vsop_series_t vsop_lon_Mars[] =
{
    { 5, vsop_lon_Mars_0 },
    { 3, vsop_lon_Mars_1 }
};

static const vsop_term_t vsop_rad_Mars_0[] =
{
    { 1.53033488271, 0.00000000000, 0.00000000000 },
    { 0.14184953160, 3.47971283528, 3340.61242669980 },
    { 0.00660776362, 3.81783443019, 6681.22485339960 },
    { 0.00046179117, 4.15595316782, 10021.83728009940 },
    { 0.00008109733, 5.55958416318, 2810.92146160520 },
    { 0.00007485318, 1.77239078402, 5621.84292321040 },
    { 0.00005523191, 1.36436303770, 2281.23049651060 },
    { 0.00003825160, 4.49407183687, 13362.44970679920 },
    { 0.00002306537, 0.09081579001, 2544.31441988340 },
    { 0.00001999396, 5.36059617709, 3337.08930835080 },
    { 0.00002484394, 4.92545639920, 2942.46342329160 },
    { 0.00001960195, 4.74249437639, 3344.13554504880 },
    { 0.00001167119, 2.11260868341, 5092.15195811580 },
    { 0.00001102816, 5.00908403998, 398.14900340820 },
    { 0.00000899066, 4.40791133207, 529.69096509460 },
    { 0.00000992252, 5.83861961952, 6151.53388830500 },
    { 0.00000807354, 2.10217065501, 1059.38193018920 },
    { 0.00000797915, 3.44839203899, 796.29800681640 },
    { 0.00000740975, 1.49906336885, 2146.16541647520 }
};

static const vsop_term_t vsop_rad_Mars_1[] =
{
    { 0.01107433345, 2.03250524857, 3340.61242669980 },
    { 0.00103175887, 2.37071847807, 6681.22485339960 },
    { 0.00012877200, 0.00000000000, 0.00000000000 },
    { 0.00010815880, 2.70888095665, 10021.83728009940 }
};

static const vsop_term_t vsop_rad_Mars_2[] =
{
    { 0.00044242249, 0.47930604954, 3340.61242669980 },
    { 0.00008138042, 0.86998389204, 6681.22485339960 }
};

static const vsop_series_t vsop_rad_Mars[] =
{
    { 19, vsop_rad_Mars_0 },
    { 4, vsop_rad_Mars_1 },
    { 2, vsop_rad_Mars_2 }
};

;
static const vsop_term_t vsop_lat_Jupiter_0[] =
{
    { 0.59954691494, 0.00000000000, 0.00000000000 },
    { 0.09695898719, 5.06191793158, 529.69096509460 },
    { 0.00573610142, 1.44406205629, 7.11354700080 },
    { 0.00306389205, 5.41734730184, 1059.38193018920 },
    { 0.00097178296, 4.14264726552, 632.78373931320 },
    { 0.00072903078, 3.64042916389, 522.57741809380 },
    { 0.00064263975, 3.41145165351, 103.09277421860 },
    { 0.00039806064, 2.29376740788, 419.48464387520 },
    { 0.00038857767, 1.27231755835, 316.39186965660 },
    { 0.00027964629, 1.78454591820, 536.80451209540 },
    { 0.00013589730, 5.77481040790, 1589.07289528380 },
    { 0.00008246349, 3.58227925840, 206.18554843720 },
    { 0.00008768704, 3.63000308199, 949.17560896980 },
    { 0.00007368042, 5.08101194270, 735.87651353180 },
    { 0.00006263150, 0.02497628807, 213.29909543800 },
    { 0.00006114062, 4.51319998626, 1162.47470440780 },
    { 0.00004905396, 1.32084470588, 110.20632121940 },
    { 0.00005305285, 1.30671216791, 14.22709400160 },
    { 0.00005305441, 4.18625634012, 1052.26838318840 },
    { 0.00004647248, 4.69958103684, 3.93215326310 },
    { 0.00003045023, 4.31676431084, 426.59819087600 },
    { 0.00002609999, 1.56667394063, 846.08283475120 },
    { 0.00002028191, 1.06376530715, 3.18139373770 },
    { 0.00001764763, 2.14148655117, 1066.49547719000 },
    { 0.00001722972, 3.88036268267, 1265.56747862640 },
    { 0.00001920945, 0.97168196472, 639.89728631400 },
    { 0.00001633223, 3.58201833555, 515.46387109300 },
    { 0.00001431999, 4.29685556046, 625.67019231240 },
    { 0.00000973272, 4.09764549134, 95.97922721780 }
};

static const vsop_term_t vsop_lat_Jupiter_1[] =
{
    { 529.69096508814, 0.00000000000, 0.00000000000 },
    { 0.00489503243, 4.22082939470, 529.69096509460 },
    { 0.00228917222, 6.02646855621, 7.11354700080 },
    { 0.00030099479, 4.54540782858, 1059.38193018920 },
    { 0.00020720920, 5.45943156902, 522.57741809380 },
    { 0.00012103653, 0.16994816098, 536.80451209540 },
    { 0.00006067987, 4.42422292017, 103.09277421860 },
    { 0.00005433968, 3.98480737746, 419.48464387520 },
    { 0.00004237744, 5.89008707199, 14.22709400160 }
};

static const vsop_term_t vsop_lat_Jupiter_2[] =
{
    { 0.00047233601, 4.32148536482, 7.11354700080 },
    { 0.00030649436, 2.92977788700, 529.69096509460 },
    { 0.00014837605, 3.14159265359, 0.00000000000 }
};

static const vsop_series_t vsop_lat_Jupiter[] =
{
    { 29, vsop_lat_Jupiter_0 },
    { 9, vsop_lat_Jupiter_1 },
    { 3, vsop_lat_Jupiter_2 }
};

static const vsop_term_t vsop_lon_Jupiter_0[] =
{
    { 0.02268615702, 3.55852606721, 529.69096509460 },
    { 0.00109971634, 3.90809347197, 1059.38193018920 },
    { 0.00110090358, 0.00000000000, 0.00000000000 },
    { 0.00008101428, 3.60509572885, 522.57741809380 },
    { 0.00006043996, 4.25883108339, 1589.07289528380 },
    { 0.00006437782, 0.30627119215, 536.80451209540 }
};

static const vsop_term_t vsop_lon_Jupiter_1[] =
{
    { 0.00078203446, 1.52377859742, 529.69096509460 }
};

static const vsop_series_t vsop_lon_Jupiter[] =
{
    { 6, vsop_lon_Jupiter_0 },
    { 1, vsop_lon_Jupiter_1 }
};

static const vsop_term_t vsop_rad_Jupiter_0[] =
{
    { 5.20887429326, 0.00000000000, 0.00000000000 },
    { 0.25209327119, 3.49108639871, 529.69096509460 },
    { 0.00610599976, 3.84115365948, 1059.38193018920 },
    { 0.00282029458, 2.57419881293, 632.78373931320 },
    { 0.00187647346, 2.07590383214, 522.57741809380 },
    { 0.00086792905, 0.71001145545, 419.48464387520 },
    { 0.00072062974, 0.21465724607, 536.80451209540 },
    { 0.00065517248, 5.97995884790, 316.39186965660 },
    { 0.00029134542, 1.67759379655, 103.09277421860 },
    { 0.00030135335, 2.16132003734, 949.17560896980 },
    { 0.00023453271, 3.54023522184, 735.87651353180 },
    { 0.00022283743, 4.19362594399, 1589.07289528380 },
    { 0.00023947298, 0.27458037480, 7.11354700080 },
    { 0.00013032614, 2.96042965363, 1162.47470440780 },
    { 0.00009703360, 1.90669633585, 206.18554843720 },
    { 0.00012749023, 2.71550286592, 1052.26838318840 },
    { 0.00007057931, 2.18184839926, 1265.56747862640 },
    { 0.00006137703, 6.26418240033, 846.08283475120 },
    { 0.00002616976, 2.00994012876, 1581.95934828300 }
};

static const vsop_term_t vsop_rad_Jupiter_1[] =
{
    { 0.01271801520, 2.64937512894, 529.69096509460 },
    { 0.00061661816, 3.00076460387, 1059.38193018920 },
    { 0.00053443713, 3.89717383175, 522.57741809380 },
    { 0.00031185171, 4.88276958012, 536.80451209540 },
    { 0.00041390269, 0.00000000000, 0.00000000000 }
};

static const vsop_series_t vsop_rad_Jupiter[] =
{
    { 19, vsop_rad_Jupiter_0 },
    { 5, vsop_rad_Jupiter_1 }
};

;
static const vsop_term_t vsop_lat_Saturn_0[] =
{
    { 0.87401354025, 0.00000000000, 0.00000000000 },
    { 0.11107659762, 3.96205090159, 213.29909543800 },
    { 0.01414150957, 4.58581516874, 7.11354700080 },
    { 0.00398379389, 0.52112032699, 206.18554843720 },
    { 0.00350769243, 3.30329907896, 426.59819087600 },
    { 0.00206816305, 0.24658372002, 103.09277421860 },
    { 0.00079271300, 3.84007056878, 220.41264243880 },
    { 0.00023990355, 4.66976924553, 110.20632121940 },
    { 0.00016573588, 0.43719228296, 419.48464387520 },
    { 0.00014906995, 5.76903183869, 316.39186965660 },
    { 0.00015820290, 0.93809155235, 632.78373931320 },
    { 0.00014609559, 1.56518472000, 3.93215326310 },
    { 0.00013160301, 4.44891291899, 14.22709400160 },
    { 0.00015053543, 2.71669915667, 639.89728631400 },
    { 0.00013005299, 5.98119023644, 11.04570026390 },
    { 0.00010725067, 3.12939523827, 202.25339517410 },
    { 0.00005863206, 0.23656938524, 529.69096509460 },
    { 0.00005227757, 4.20783365759, 3.18139373770 },
    { 0.00006126317, 1.76328667907, 277.03499374140 },
    { 0.00005019687, 3.17787728405, 433.71173787680 },
    { 0.00004592550, 0.61977744975, 199.07200143640 },
    { 0.00004005867, 2.24479718502, 63.73589830340 },
    { 0.00002953796, 0.98280366998, 95.97922721780 },
    { 0.00003873670, 3.22283226966, 138.51749687070 },
    { 0.00002461186, 2.03163875071, 735.87651353180 },
    { 0.00003269484, 0.77492638211, 949.17560896980 },
    { 0.00001758145, 3.26580109940, 522.57741809380 },
    { 0.00001640172, 5.50504453050, 846.08283475120 },
    { 0.00001391327, 4.02333150505, 323.50541665740 },
    { 0.00001580648, 4.37265307169, 309.27832265580 },
    { 0.00001123498, 2.83726798446, 415.55249061210 },
    { 0.00001017275, 3.71700135395, 227.52618943960 },
    { 0.00000848642, 3.19150170830, 209.36694217490 }
};

static const vsop_term_t vsop_lat_Saturn_1[] =
{
    { 213.29909521690, 0.00000000000, 0.00000000000 },
    { 0.01297370862, 1.82834923978, 213.29909543800 },
    { 0.00564345393, 2.88499717272, 7.11354700080 },
    { 0.00093734369, 1.06311793502, 426.59819087600 },
    { 0.00107674962, 2.27769131009, 206.18554843720 },
    { 0.00040244455, 2.04108104671, 220.41264243880 },
    { 0.00019941774, 1.27954390470, 103.09277421860 },
    { 0.00010511678, 2.74880342130, 14.22709400160 },
    { 0.00006416106, 0.38238295041, 639.89728631400 },
    { 0.00004848994, 2.43037610229, 419.48464387520 },
    { 0.00004056892, 2.92133209468, 110.20632121940 },
    { 0.00003768635, 3.64965330780, 3.93215326310 }
};

static const vsop_term_t vsop_lat_Saturn_2[] =
{
    { 0.00116441330, 1.17988132879, 7.11354700080 },
    { 0.00091841837, 0.07325195840, 213.29909543800 },
    { 0.00036661728, 0.00000000000, 0.00000000000 },
    { 0.00015274496, 4.06493179167, 206.18554843720 }
};

static const vsop_series_t vsop_lat_Saturn[] =
{
    { 33, vsop_lat_Saturn_0 },
    { 12, vsop_lat_Saturn_1 },
    { 4, vsop_lat_Saturn_2 }
};

static const vsop_term_t vsop_lon_Saturn_0[] =
{
    { 0.04330678039, 3.60284428399, 213.29909543800 },
    { 0.00240348302, 2.85238489373, 426.59819087600 },
    { 0.00084745939, 0.00000000000, 0.00000000000 },
    { 0.00030863357, 3.48441504555, 220.41264243880 },
    { 0.00034116062, 0.57297307557, 206.18554843720 },
    { 0.00014734070, 2.11846596715, 639.89728631400 },
    { 0.00009916667, 5.79003188904, 419.48464387520 },
    { 0.00006993564, 4.73604689720, 7.11354700080 },
    { 0.00004807588, 5.43305312061, 316.39186965660 }
};

static const vsop_term_t vsop_lon_Saturn_1[] =
{
    { 0.00198927992, 4.93901017903, 213.29909543800 },
    { 0.00036947916, 3.14159265359, 0.00000000000 },
    { 0.00017966989, 0.51979431110, 426.59819087600 }
};

static const vsop_series_t vsop_lon_Saturn[] =
{
    { 9, vsop_lon_Saturn_0 },
    { 3, vsop_lon_Saturn_1 }
};

static const vsop_term_t vsop_rad_Saturn_0[] =
{
    { 9.55758135486, 0.00000000000, 0.00000000000 },
    { 0.52921382865, 2.39226219573, 213.29909543800 },
    { 0.01873679867, 5.23549604660, 206.18554843720 },
    { 0.01464663929, 1.64763042902, 426.59819087600 },
    { 0.00821891141, 5.93520042303, 316.39186965660 },
    { 0.00547506923, 5.01532618980, 103.09277421860 },
    { 0.00371684650, 2.27114821115, 220.41264243880 },
    { 0.00361778765, 3.13904301847, 7.11354700080 },
    { 0.00140617506, 5.70406606781, 632.78373931320 },
    { 0.00108974848, 3.29313390175, 110.20632121940 },
    { 0.00069006962, 5.94099540992, 419.48464387520 },
    { 0.00061053367, 0.94037691801, 639.89728631400 },
    { 0.00048913294, 1.55733638681, 202.25339517410 },
    { 0.00034143772, 0.19519102597, 277.03499374140 },
    { 0.00032401773, 5.47084567016, 949.17560896980 },
    { 0.00020936596, 0.46349251129, 735.87651353180 },
    { 0.00009796004, 5.20477537945, 1265.56747862640 },
    { 0.00011993338, 5.98050967385, 846.08283475120 },
    { 0.00020839300, 1.52102476129, 433.71173787680 },
    { 0.00015298404, 3.05943814940, 529.69096509460 },
    { 0.00006465823, 0.17732249942, 1052.26838318840 },
    { 0.00011380257, 1.73105427040, 522.57741809380 },
    { 0.00003419618, 4.94550542171, 1581.95934828300 }
};

static const vsop_term_t vsop_rad_Saturn_1[] =
{
    { 0.06182981340, 0.25843511480, 213.29909543800 },
    { 0.00506577242, 0.71114625261, 206.18554843720 },
    { 0.00341394029, 5.79635741658, 426.59819087600 },
    { 0.00188491195, 0.47215589652, 220.41264243880 },
    { 0.00186261486, 3.14159265359, 0.00000000000 },
    { 0.00143891146, 1.40744822888, 7.11354700080 }
};

static const vsop_term_t vsop_rad_Saturn_2[] =
{
    { 0.00436902572, 4.78671677509, 213.29909543800 }
};

static const vsop_series_t vsop_rad_Saturn[] =
{
    { 23, vsop_rad_Saturn_0 },
    { 6, vsop_rad_Saturn_1 },
    { 1, vsop_rad_Saturn_2 }
};

;
static const vsop_term_t vsop_lat_Uranus_0[] =
{
    { 5.48129294297, 0.00000000000, 0.00000000000 },
    { 0.09260408234, 0.89106421507, 74.78159856730 },
    { 0.01504247898, 3.62719260920, 1.48447270830 },
    { 0.00365981674, 1.89962179044, 73.29712585900 },
    { 0.00272328168, 3.35823706307, 149.56319713460 },
    { 0.00070328461, 5.39254450063, 63.73589830340 },
    { 0.00068892678, 6.09292483287, 76.26607127560 },
    { 0.00061998615, 2.26952066061, 2.96894541660 },
    { 0.00061950719, 2.85098872691, 11.04570026390 },
    { 0.00026468770, 3.14152083966, 71.81265315070 },
    { 0.00025710476, 6.11379840493, 454.90936652730 },
    { 0.00021078850, 4.36059339067, 148.07872442630 },
    { 0.00017818647, 1.74436930289, 36.64856292950 },
    { 0.00014613507, 4.73732166022, 3.93215326310 },
    { 0.00011162509, 5.82681796350, 224.34479570190 },
    { 0.00010997910, 0.48865004018, 138.51749687070 },
    { 0.00009527478, 2.95516862826, 35.16409022120 },
    { 0.00007545601, 5.23626582400, 109.94568878850 },
    { 0.00004220241, 3.23328220918, 70.84944530420 },
    { 0.00004051900, 2.27755017300, 151.04766984290 },
    { 0.00003354596, 1.06549007380, 4.45341812490 },
    { 0.00002926718, 4.62903718891, 9.56122755560 },
    { 0.00003490340, 5.48306144511, 146.59425171800 },
    { 0.00003144069, 4.75199570434, 77.75054398390 },
    { 0.00002922333, 5.35235361027, 85.82729883120 },
    { 0.00002272788, 4.36600400036, 70.32818044240 },
    { 0.00002051219, 1.51773566586, 0.11187458460 },
    { 0.00002148602, 0.60745949945, 38.13303563780 },
    { 0.00001991643, 4.92437588682, 277.03499374140 },
    { 0.00001376226, 2.04283539351, 65.22037101170 },
    { 0.00001666902, 3.62744066769, 380.12776796000 },
    { 0.00001284107, 3.11347961505, 202.25339517410 },
    { 0.00001150429, 0.93343589092, 3.18139373770 },
    { 0.00001533221, 2.58594681212, 52.69019803950 },
    { 0.00001281604, 0.54271272721, 222.86032299360 },
    { 0.00001372139, 4.19641530878, 111.43016149680 },
    { 0.00001221029, 0.19900650030, 108.46121608020 },
    { 0.00000946181, 1.19253165736, 127.47179660680 },
    { 0.00001150989, 4.17898916639, 33.67961751290 }
};

static const vsop_term_t vsop_lat_Uranus_1[] =
{
    { 74.78159860910, 0.00000000000, 0.00000000000 },
    { 0.00154332863, 5.24158770553, 74.78159856730 },
    { 0.00024456474, 1.71260334156, 1.48447270830 },
    { 0.00009258442, 0.42829732350, 11.04570026390 },
    { 0.00008265977, 1.50218091379, 63.73589830340 },
    { 0.00009150160, 1.41213765216, 149.56319713460 }
};

static const vsop_series_t vsop_lat_Uranus[] =
{
    { 39, vsop_lat_Uranus_0 },
    { 6, vsop_lat_Uranus_1 }
};

static const vsop_term_t vsop_lon_Uranus_0[] =
{
    { 0.01346277648, 2.61877810547, 74.78159856730 },
    { 0.00062341400, 5.08111189648, 149.56319713460 },
    { 0.00061601196, 3.14159265359, 0.00000000000 },
    { 0.00009963722, 1.61603805646, 76.26607127560 },
    { 0.00009926160, 0.57630380333, 73.29712585900 }
};

static const vsop_term_t vsop_lon_Uranus_1[] =
{
    { 0.00034101978, 0.01321929936, 74.78159856730 }
};

static const vsop_series_t vsop_lon_Uranus[] =
{
    { 5, vsop_lon_Uranus_0 },
    { 1, vsop_lon_Uranus_1 }
};

static const vsop_term_t vsop_rad_Uranus_0[] =
{
    { 19.21264847206, 0.00000000000, 0.00000000000 },
    { 0.88784984413, 5.60377527014, 74.78159856730 },
    { 0.03440836062, 0.32836099706, 73.29712585900 },
    { 0.02055653860, 1.78295159330, 149.56319713460 },
    { 0.00649322410, 4.52247285911, 76.26607127560 },
    { 0.00602247865, 3.86003823674, 63.73589830340 },
    { 0.00496404167, 1.40139935333, 454.90936652730 },
    { 0.00338525369, 1.58002770318, 138.51749687070 },
    { 0.00243509114, 1.57086606044, 71.81265315070 },
    { 0.00190522303, 1.99809394714, 1.48447270830 },
    { 0.00161858838, 2.79137786799, 148.07872442630 },
    { 0.00143706183, 1.38368544947, 11.04570026390 },
    { 0.00093192405, 0.17437220467, 36.64856292950 },
    { 0.00071424548, 4.24509236074, 224.34479570190 },
    { 0.00089806014, 3.66105364565, 109.94568878850 },
    { 0.00039009723, 1.66971401684, 70.84944530420 },
    { 0.00046677296, 1.39976401694, 35.16409022120 },
    { 0.00039025624, 3.36234773834, 277.03499374140 },
    { 0.00036755274, 3.88649278513, 146.59425171800 },
    { 0.00030348723, 0.70100838798, 151.04766984290 },
    { 0.00029156413, 3.18056336700, 77.75054398390 },
    { 0.00022637073, 0.72518687029, 529.69096509460 },
    { 0.00011959076, 1.75043392140, 984.60033162190 },
    { 0.00025620756, 5.25656086672, 380.12776796000 }
};

static const vsop_term_t vsop_rad_Uranus_1[] =
{
    { 0.01479896629, 3.67205697578, 74.78159856730 }
};

static const vsop_series_t vsop_rad_Uranus[] =
{
    { 24, vsop_rad_Uranus_0 },
    { 1, vsop_rad_Uranus_1 }
};

;
static const vsop_term_t vsop_lat_Neptune_0[] =
{
    { 5.31188633046, 0.00000000000, 0.00000000000 },
    { 0.01798475530, 2.90101273890, 38.13303563780 },
    { 0.01019727652, 0.48580922867, 1.48447270830 },
    { 0.00124531845, 4.83008090676, 36.64856292950 },
    { 0.00042064466, 5.41054993053, 2.96894541660 },
    { 0.00037714584, 6.09221808686, 35.16409022120 },
    { 0.00033784738, 1.24488874087, 76.26607127560 },
    { 0.00016482741, 0.00007727998, 491.55792945680 },
    { 0.00009198584, 4.93747051954, 39.61750834610 },
    { 0.00008994250, 0.27462171806, 175.16605980020 }
};

static const vsop_term_t vsop_lat_Neptune_1[] =
{
    { 38.13303563957, 0.00000000000, 0.00000000000 },
    { 0.00016604172, 4.86323329249, 1.48447270830 },
    { 0.00015744045, 2.27887427527, 38.13303563780 }
};

static const vsop_series_t vsop_lat_Neptune[] =
{
    { 10, vsop_lat_Neptune_0 },
    { 3, vsop_lat_Neptune_1 }
};

static const vsop_term_t vsop_lon_Neptune_0[] =
{
    { 0.03088622933, 1.44104372644, 38.13303563780 },
    { 0.00027780087, 5.91271884599, 76.26607127560 },
    { 0.00027623609, 0.00000000000, 0.00000000000 },
    { 0.00015355489, 2.52123799551, 36.64856292950 },
    { 0.00015448133, 3.50877079215, 39.61750834610 }
};

static const vsop_series_t vsop_lon_Neptune[] =
{
    { 5, vsop_lon_Neptune_0 }
};

static const vsop_term_t vsop_rad_Neptune_0[] =
{
    { 30.07013205828, 0.00000000000, 0.00000000000 },
    { 0.27062259632, 1.32999459377, 38.13303563780 },
    { 0.01691764014, 3.25186135653, 36.64856292950 },
    { 0.00807830553, 5.18592878704, 1.48447270830 },
    { 0.00537760510, 4.52113935896, 35.16409022120 },
    { 0.00495725141, 1.57105641650, 491.55792945680 },
    { 0.00274571975, 1.84552258866, 175.16605980020 },
    { 0.00012012320, 1.92059384991, 1021.24889455140 },
    { 0.00121801746, 5.79754470298, 76.26607127560 },
    { 0.00100896068, 0.37702724930, 73.29712585900 },
    { 0.00135134092, 3.37220609835, 39.61750834610 },
    { 0.00007571796, 1.07149207335, 388.46515523820 }
};

static const vsop_series_t vsop_rad_Neptune[] =
{
    { 12, vsop_rad_Neptune_0 }
};

;

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

static double VsopHelioDistance(const vsop_model_t *model, astro_time_t time)
{
    int s, i;
    double t = time.tt / 365250;    /* millennia since 2000 */
    double distance = 0.0;
    double tpower = 1.0;
    const vsop_formula_t *formula = &model->formula[2];     /* [2] is the distance part of the formula */

    /*
        The caller only wants to know the distance between the planet and the Sun.
        So we only need to calculate the radial component of the spherical coordinates.
    */

    for (s=0; s < formula->nseries; ++s)
    {
        double sum = 0.0;
        const vsop_series_t *series = &formula->series[s];
        for (i=0; i < series->nterms; ++i)
        {
            const vsop_term_t *term = &series->term[i];
            sum += term->amplitude * cos(term->phase + (t * term->frequency));
        }
        distance += tpower * sum;
        tpower *= t;
    }

    return distance;
}

/*------------------ TOP2013 model for Pluto ------------------*/

/** @cond DOXYGEN_SKIP */

#define TOP_NCOORDS  6
#define TOP_NSERIES 13

typedef struct
{
    double k;
    double c;
    double s;
}
astro_top_term_t;

typedef struct
{
    int nterms;
    const astro_top_term_t *terms;
}
astro_top_series_t;

typedef struct
{
    int nseries;
    const astro_top_series_t *series;
}
astro_top_model_t;

typedef struct
{
    double a;           /* AU */
    double lambda;      /* rad */
    double k;           /* 1 */
    double h;           /* 1 */
    double q;           /* 1 */
    double p;           /* 1 */
}
top_elliptical_t;
/** @endcond */

static const astro_top_term_t topterms_8_0_0[] =
{
    {        0,  3.9544617144029999e+01,  0.0000000000000000e+00}  /*      0 */
,   {     1402, -1.8891373533434089e-01, -8.5258197635470073e-02}  /*      1 */
,   {     1331, -4.1495877833812339e-02, -3.3387415274886263e-02}  /*      2 */
,   {      522, -4.8502474919249819e-02, -7.3455412272547278e-03}  /*      3 */
,   {       71,  2.8102132918948221e-02, -5.3468437660152152e-03}  /*      4 */
,   {     1261, -9.1600251304608978e-03, -1.2309204554431390e-02}  /*      5 */
,   {      452, -1.2202161344831991e-02, -5.2519856329982890e-03}  /*      6 */
,   {     2875, -9.7845229475467185e-03, -7.2968560644718955e-04}  /*      7 */
,   {       35,  4.8494518209585983e-03, -6.8918970374425084e-03}  /*      8 */
,   {      141,  5.8271375488234932e-03, -3.1946778653436391e-03}  /*      9 */
,   {      137,  1.5300059509576150e-03, -6.0327729791525954e-03}  /*     10 */
,   {        4, -1.9494897408412360e-03,  4.5717130739994704e-03}  /*     11 */
,   {     1190, -1.7524220672664240e-03, -4.2980683502911454e-03}  /*     12 */
,   {      381, -3.1062775803702681e-03, -2.4728667551542258e-03}  /*     13 */
,   {        8, -1.7188663433411050e-03,  2.9756270077158122e-03}  /*     14 */
,   {     1543, -7.7472653128184826e-04,  2.6514626782777680e-03}  /*     15 */
,   {     1115, -1.5405722125111840e-03, -2.0778390548994150e-03}  /*     16 */
,   {     2804, -2.0048397209869230e-03, -4.1957951179189120e-04}  /*     17 */
,   {       67,  9.6850762192148931e-04, -1.5811913714969829e-03}  /*     18 */
,   {      212,  1.5466715821083480e-03, -9.6836654994834209e-04}  /*     19 */
,   {     1119, -1.8820367463891121e-04, -1.4293834479379090e-03}  /*     20 */
,   {    17405, -1.0738845199599739e-03,  9.5985010997943349e-04}  /*     21 */
,   {    28337,  7.5821211083786067e-04,  1.1416213940389449e-03}  /*     22 */
,   {      310, -6.9650370158153983e-04, -1.0024667762364200e-03}  /*     23 */
,   {     1044, -3.2801771454650589e-04, -6.4947140155397116e-04}  /*     24 */
,   {       63,  4.0424767075158291e-04,  5.7315886355325109e-04}  /*     25 */
,   {     1614, -1.3238448498587051e-04,  6.7949492229369074e-04}  /*     26 */
,   {       12, -6.5048480404874143e-04, -1.0430021074697129e-04}  /*     27 */
,   {      283,  3.3929768009878552e-04, -5.1573678840263412e-04}  /*     28 */
,   {      133,  4.0138197254613279e-04,  4.5284627712770291e-04}  /*     29 */
,   {     1421,  5.0823717785117468e-04,  3.1010622577270548e-04}  /*     30 */
,   {     1383, -5.6813906164891585e-04, -1.7225147327178090e-04}  /*     31 */
,   {     4348, -5.1094469101100064e-04,  1.4416513132178369e-04}  /*     32 */
,   {     2733, -4.8192307867155672e-04, -1.8631709444892481e-04}  /*     33 */
,   {      664,  6.3695948832436563e-05,  5.0429756322537705e-04}  /*     34 */
,   {      177, -3.0532744995821309e-04,  4.0226535349675440e-04}  /*     35 */
,   {      204,  3.4677834246970749e-04,  3.2588064363496143e-04}  /*     36 */
,   {     1048,  5.6219036569383140e-05, -4.5145044715373130e-04}  /*     37 */
,   {     1406,  3.6764565360298558e-04, -2.4793326161876619e-04}  /*     38 */
,   {     1398, -6.1399335668996843e-05, -4.0368084856225791e-04}  /*     39 */
,   {      880,  3.8206123841964649e-04,  1.0727867737651879e-04}  /*     40 */
,   {      239, -9.2482984334176681e-05, -3.6357760642322443e-04}  /*     41 */
,   {      541, -3.5492935880988093e-04, -7.2696664262252120e-05}  /*     42 */
,   {    17334, -3.1608420086033548e-04,  1.6282567766535830e-04}  /*     43 */
,   {      503,  3.4204624270103872e-04,  6.9800261361389600e-06}  /*     44 */
,   {    28266,  1.1026801272880990e-04,  3.1928911807394130e-04}  /*     45 */
,   {      974, -1.5403436032839399e-04, -2.6908883645786191e-04}  /*     46 */
,   {      345, -2.1738414747494940e-04,  1.5563797862269319e-04}  /*     47 */
,   {     1924,  2.1398279518720060e-04,  1.4058024849360011e-04}  /*     48 */
,   {      275,  1.7682640345886169e-04,  1.7533174872185701e-04}  /*     49 */
,   {      106,  2.0541555697206479e-04, -1.2965626091678131e-04}  /*     50 */
,   {     1335,  2.1493444414672229e-04, -1.0499061765897920e-04}  /*     51 */
,   {      271,  2.1059070752496800e-04, -9.6647526271973576e-05}  /*     52 */
};

static const astro_top_term_t topterms_8_0_1[] =
{
    {     1402, -2.4541007710854022e-02,  5.3822529675651168e-02}  /*      0 */
,   {        0,  3.7890000000000000e-02,  0.0000000000000000e+00}  /*      1 */
,   {     1331, -1.5981733929364941e-02,  1.9623571765255570e-02}  /*      2 */
,   {      522, -2.1667846285401051e-03,  1.3851523320136261e-02}  /*      3 */
,   {       71,  8.6723020043621556e-04,  4.7241967094509424e-03}  /*      4 */
,   {     1261, -3.8008804453790821e-03,  2.7392780267271222e-03}  /*      5 */
,   {     2875, -5.9015441259769620e-04,  3.3190792202711962e-03}  /*      6 */
,   {       35, -2.0646844712970550e-03, -2.6221746776036530e-03}  /*      7 */
,   {        4, -2.7341177070724500e-03,  8.4711655196688536e-04}  /*      8 */
,   {     1190, -2.1427361371102339e-03,  8.3829844178717706e-04}  /*      9 */
,   {      452, -6.5191457442425422e-04,  1.3962477909825410e-03}  /*     10 */
,   {        8, -1.4174353266390050e-03, -4.8800609526062797e-04}  /*     11 */
,   {      381, -7.7874017402358557e-04,  9.3105764052437356e-04}  /*     12 */
,   {      137, -1.1495119655719700e-03, -2.9568842990989599e-04}  /*     13 */
,   {     2804, -2.9320339800219448e-04,  1.0245691529492771e-03}  /*     14 */
,   {     1119, -9.8344008047702861e-04,  1.2122440012927610e-04}  /*     15 */
,   {     1543,  7.2662531437559094e-04,  2.2380632879329121e-04}  /*     16 */
,   {     1115, -4.9059020449034973e-04,  5.4125269506251448e-04}  /*     17 */
,   {      310, -5.0858314998190228e-04,  3.3071892812221339e-04}  /*     18 */
};

static const astro_top_term_t topterms_8_0_2[] =
{
    {     1402,  6.1564734183986881e-03,  6.9727544398706350e-03}  /*      0 */
,   {     1331,  3.4473718686991932e-03,  5.3547458162748890e-03}  /*      1 */
,   {        0, -6.0199059744408881e-03,  0.0000000000000000e+00}  /*      2 */
,   {      522,  1.8506017571807240e-03,  1.2255527089933200e-03}  /*      3 */
,   {        4, -9.0791196955079877e-04, -1.0597464135497440e-03}  /*      4 */
,   {     1261, -2.5360705302150178e-04,  1.1006253876392330e-03}  /*      5 */
};

static const astro_top_series_t topseries_8_0[] =
{
    {     53, topterms_8_0_0 }
,   {     19, topterms_8_0_1 }
,   {      6, topterms_8_0_2 }
};

static const astro_top_term_t topterms_8_1_0[] =
{
    {        0,  4.1654711248260003e+00,  0.0000000000000000e+00}  /*      0 */
,   {     1402,  2.0610516934972331e-03, -4.5672163937433728e-03}  /*      1 */
,   {        4,  3.4901370659179681e-03,  2.2214931280208532e-03}  /*      2 */
,   {     1473, -2.5274121559055551e-04,  1.4438005411011950e-03}  /*      3 */
,   {        8,  1.1811406701290150e-03,  5.6294240646172070e-04}  /*      4 */
,   {      522,  1.6107717821433110e-04, -1.0623010008063969e-03}  /*      5 */
,   {     1331,  2.7651217171242138e-04, -3.4287860661004818e-04}  /*      6 */
,   {      593,  3.3064342491574279e-05,  3.2281574175741223e-04}  /*      7 */
,   {     2875,  1.7391859631880900e-05, -2.4371046197993790e-04}  /*      8 */
,   {       71,  4.3540142673103847e-05, -2.3324350638665871e-04}  /*      9 */
,   {       12, -1.4733740390542089e-05,  1.6458293257830351e-04}  /*     10 */
,   {       35, -7.7366158840780485e-05, -1.2373453663453069e-04}  /*     11 */
,   {      137,  1.0874981420582371e-04,  2.9642106239502481e-05}  /*     12 */
,   {      452,  3.2837879995135691e-05, -8.0083302669255936e-05}  /*     13 */
,   {      106, -7.1058821662524569e-05, -2.9131217624735831e-05}  /*     14 */
,   {     2945,  1.2419927037655510e-05,  6.9699905790296052e-05}  /*     15 */
,   {     1115,  4.9015184839992363e-05, -3.7626776449287022e-05}  /*     16 */
,   {     1261,  3.7725355301707661e-05, -2.7730601457362580e-05}  /*     17 */
,   {       63,  3.4620953450636643e-05, -2.5727054077819990e-05}  /*     18 */
,   {    17405, -2.4737834776937812e-05, -2.7668964560844859e-05}  /*     19 */
,   {      141, -1.9152627579774910e-05,  3.0608840803230902e-05}  /*     20 */
,   {     1543,  3.4403440209432860e-05,  8.2302490865683481e-06}  /*     21 */
,   {    28337, -2.9450824049753070e-05,  1.9509947544013829e-05}  /*     22 */
,   {      208, -3.3402595486149828e-05,  6.9721024372198099e-07}  /*     23 */
,   {       16, -2.5708226560219951e-05,  5.4690438392817274e-06}  /*     24 */
,   {     2804,  8.2434428446966487e-06, -2.2056820101697141e-05}  /*     25 */
,   {      133,  1.7535496693018641e-05, -1.5317433153649361e-05}  /*     26 */
,   {      177,  9.4209450622295435e-06,  1.8609791865429379e-05}  /*     27 */
,   {      204,  1.3418024861975440e-05, -1.5519809537438971e-05}  /*     28 */
,   {       67,  1.9542040889474939e-05, -2.0438794816549699e-07}  /*     29 */
,   {     1186, -1.1497747150912189e-05,  1.4072246679624251e-05}  /*     30 */
,   {     1421, -7.4843280005381670e-06,  1.2303197921178691e-05}  /*     31 */
,   {     1383,  4.1628765439336240e-06, -1.3719681259667719e-05}  /*     32 */
,   {      275,  7.1803653660818622e-06, -1.2018565522581090e-05}  /*     33 */
,   {      212,  1.3642921909731060e-05,  2.1147480912091179e-06}  /*     34 */
,   {       59, -2.8996742569855378e-06, -1.3093850574478501e-05}  /*     35 */
,   {     4348, -3.6912767653086940e-06, -1.2768209379424250e-05}  /*     36 */
,   {    17476,  8.8159150487736961e-06,  5.9364512059773267e-06}  /*     37 */
,   {       27,  6.2980935945243974e-06, -8.4472797839010716e-06}  /*     38 */
,   {     1406,  5.7734139509693619e-06,  8.7248712048903149e-06}  /*     39 */
,   {     1398,  1.0115989678361140e-05, -1.3220473394572090e-06}  /*     40 */
,   {    28407,  6.7642022362205018e-06, -7.5210649765782616e-06}  /*     41 */
,   {      664,  7.9736063591689162e-06, -2.2024063975673750e-06}  /*     42 */
,   {      381,  3.1761072779613261e-06, -7.3448935507901516e-06}  /*     43 */
,   {      541,  2.5898481389558310e-06, -7.0323742443685127e-06}  /*     44 */
,   {      503, -1.5979599372675240e-07,  7.4586091247535634e-06}  /*     45 */
,   {      271, -3.1728511373578081e-06, -6.6603677156048786e-06}  /*     46 */
,   {      129, -1.8397756711585990e-06, -7.0451043684758007e-06}  /*     47 */
,   {      247,  3.1533300480278400e-07, -7.1604266714610849e-06}  /*     48 */
,   {      200, -2.2707299296152400e-06, -6.6783936815674601e-06}  /*     49 */
,   {       31, -4.2146281840702766e-06, -5.6388188974303981e-06}  /*     50 */
,   {      341, -3.9640903174182786e-06, -5.5968952746387096e-06}  /*     51 */
,   {       20, -2.1086487296606532e-06, -4.2993743513743323e-06}  /*     52 */
,   {      974, -1.9337296464942458e-06,  4.1280092704873532e-06}  /*     53 */
,   {     1492,  1.2386884071711840e-06, -4.0279260453190260e-06}  /*     54 */
,   {     1454, -1.8674475925452479e-07,  4.1994903608977473e-06}  /*     55 */
,   {       55, -4.1054351042334211e-06, -5.3592659852261218e-07}  /*     56 */
,   {     1044,  3.5262751998119361e-06,  1.9367410649084749e-06}  /*     57 */
,   {       39, -7.0347175355875052e-07,  3.9228440423731852e-06}  /*     58 */
,   {      345, -2.7714162211172452e-06, -2.8012028299863999e-06}  /*     59 */
,   {      283, -3.8855399506750791e-06,  4.2353837211089759e-07}  /*     60 */
,   {     4418,  1.9380717555480920e-06,  3.3362906758656841e-06}  /*     61 */
,   {     1708,  3.8363066095252423e-06,  4.0071146882428590e-07}  /*     62 */
,   {      903, -2.4535489999476891e-06,  2.7713257447896069e-06}  /*     63 */
,   {      412, -2.4194398659193562e-06, -2.6818837171571071e-06}  /*     64 */
,   {    17334, -1.5878732753673419e-06, -3.0847459710695832e-06}  /*     65 */
,   {     1410,  2.1528516528967772e-06,  2.5955495586096718e-06}  /*     66 */
,   {      318, -1.6007757166478231e-06,  2.9567402302706461e-06}  /*     67 */
,   {    28266, -3.1282886569147322e-06,  1.0733981547714400e-06}  /*     68 */
,   {     1394,  3.2570224915190650e-06, -5.9986669416035846e-08}  /*     69 */
,   {    72490,  7.8958524750142168e-07,  3.1029527324934821e-06}  /*     70 */
,   {     9220,  2.8325155820290858e-06, -1.4397983817771050e-06}  /*     71 */
,   {     1614,  3.0748937595146630e-06,  5.0326493238056022e-07}  /*     72 */
,   {      408, -3.0722440820616450e-06,  4.9345659005564487e-07}  /*     73 */
,   {      337, -2.9632105019568382e-06,  9.3548023528643972e-08}  /*     74 */
,   {       79,  2.7845485358443480e-06,  7.4025133333557920e-07}  /*     75 */
,   {      354,  1.9360680354843408e-06,  2.0552415375649160e-06}  /*     76 */
,   {      612,  5.3703132984668158e-07,  2.6669817089723940e-06}  /*     77 */
,   {      279, -2.5413275396432281e-07,  2.6502544396045040e-06}  /*     78 */
,   {       75, -2.8830829576295098e-07,  2.6106707165907741e-06}  /*     79 */
,   {      479, -2.4319941706034390e-06,  7.6942886776602959e-07}  /*     80 */
,   {      416,  2.0096007406999570e-06,  1.5520759224768389e-06}  /*     81 */
,   {     1096, -2.3257110689311771e-06,  9.7831153778294669e-07}  /*     82 */
,   {      267, -2.5074175872274202e-06, -1.7701423128546709e-07}  /*     83 */
,   {     3162,  2.0698330022272070e-06, -1.3093627746428469e-06}  /*     84 */
,   {      832, -1.9285373422971391e-06,  1.4160461431375550e-06}  /*     85 */
,   {     1190,  2.2341390259691221e-06, -7.2482790590962755e-07}  /*     86 */
,   {      526,  1.7171371907181510e-06,  1.4529625626913769e-06}  /*     87 */
,   {     2733,  1.0269773239622380e-06, -2.0007533912028411e-06}  /*     88 */
,   {      574, -5.1232176836951890e-07, -2.1792410718247190e-06}  /*     89 */
};

static const astro_top_term_t topterms_8_1_1[] =
{
    {        0,  2.5335660204370001e+01,  0.0000000000000000e+00}  /*      0 */
,   {        4, -2.1897916824442529e-04,  1.7741955864290151e-03}  /*      1 */
,   {     1402, -1.3029527067277161e-03, -5.8987171410210541e-04}  /*      2 */
,   {        8, -5.7857466108113462e-05,  6.2766572590063064e-04}  /*      3 */
,   {      522, -3.0338218770459020e-04, -4.6787236333551149e-05}  /*      4 */
,   {     1331, -1.6234312756262061e-04, -1.3226650376835280e-04}  /*      5 */
,   {       35, -1.6846578574212679e-04,  5.1474337733700072e-05}  /*      6 */
,   {     1473,  1.3716672333829441e-04,  2.7605165443775319e-05}  /*      7 */
,   {       12, -1.1426985570984979e-04,  1.7062621376287921e-05}  /*      8 */
,   {     2875, -8.2654240464229338e-05, -1.4230432684507880e-05}  /*      9 */
,   {     2945,  3.6084370654255798e-05, -3.8286603564628976e-06}  /*     10 */
,   {      593,  3.1050555447062818e-05, -2.3120488185252789e-06}  /*     11 */
,   {       16, -8.7245446550032187e-06, -2.2703299980076952e-05}  /*     12 */
,   {      137,  5.5783261240786444e-06, -2.0642718216052691e-05}  /*     13 */
,   {     1115, -1.3513503169043030e-05, -1.1380785261757990e-05}  /*     14 */
,   {       71,  1.0574660674287450e-05,  1.3281184472910570e-05}  /*     15 */
,   {     1261, -8.2709801392281214e-06, -1.1688074277507450e-05}  /*     16 */
,   {     2804, -1.1505534938221521e-05, -5.1883886955282330e-06}  /*     17 */
,   {       27, -9.0174735475683120e-06, -5.2285448475511598e-06}  /*     18 */
,   {      452, -9.3099342968225354e-06, -3.8372164516008766e-06}  /*     19 */
,   {     1543,  2.4828575108302229e-06, -9.6169483680059550e-06}  /*     20 */
,   {       63, -6.5055787623425069e-06, -7.4385376578198243e-06}  /*     21 */
,   {      106, -4.5512156458035857e-06,  8.6658259433823040e-06}  /*     22 */
,   {      133, -6.3072609842624613e-06, -6.9148116074105151e-06}  /*     23 */
,   {    28337, -4.4530735415054223e-06, -6.7124377273408592e-06}  /*     24 */
,   {      141,  6.9466051844107804e-06,  4.0449257224296616e-06}  /*     25 */
,   {     1398,  1.0522658062293921e-06, -7.0563446922276261e-06}  /*     26 */
,   {       59, -5.7183266993529056e-06,  1.6273101086020201e-06}  /*     27 */
,   {     1383, -5.6807273131313079e-06, -9.5798703870467336e-07}  /*     28 */
,   {       20,  4.8130148431204634e-06, -2.5299767388566989e-06}  /*     29 */
,   {     4348, -5.3269871823774314e-06,  6.0216192993201565e-07}  /*     30 */
,   {      129, -4.2816786763917080e-06,  1.2524207693334740e-06}  /*     31 */
};

static const astro_top_term_t topterms_8_1_2[] =
{
    {        0, -1.8272218839163919e-02,  0.0000000000000000e+00}  /*      0 */
,   {        4, -4.2382205514535369e-04,  6.0951225139627217e-05}  /*      1 */
,   {        8, -2.4214950161520300e-04,  1.3555498866999700e-04}  /*      2 */
,   {     1402, -1.6731953542354990e-04,  1.4877590360992571e-04}  /*      3 */
,   {       12, -4.0514580902466651e-05, -4.6859429953625132e-05}  /*      4 */
,   {     1331, -4.4255479363996823e-05,  2.8523966290637259e-05}  /*      5 */
,   {      522, -2.6626163791203940e-05,  4.0414506964273448e-05}  /*      6 */
,   {       35, -2.3510812229573040e-06,  1.7218835328103482e-05}  /*      7 */
,   {     2875, -7.9375576768093987e-06,  1.3976740463566499e-05}  /*      8 */
,   {       16,  9.7305978044426803e-06, -1.0841117038814680e-05}  /*      9 */
,   {     2945, -4.5741854272088092e-07, -9.4120776625816716e-06}  /*     10 */
};

static const astro_top_term_t topterms_8_1_3[] =
{
    {        0,  1.9409931667071581e-03,  0.0000000000000000e+00}  /*      0 */
,   {        8, -5.2528177259165953e-05, -6.1055439661395280e-05}  /*      1 */
,   {        4, -6.0738316603187738e-05,  3.5387155556321078e-05}  /*      2 */
,   {     1402,  1.5967276563962451e-05,  3.5538891371840628e-05}  /*      3 */
,   {       12,  1.5332065032938759e-05, -2.1295985440213221e-05}  /*      4 */
};

static const astro_top_term_t topterms_8_1_4[] =
{
    {        0,  8.6099959150566781e-05,  0.0000000000000000e+00}  /*      0 */
,   {        4, -5.3544872037241911e-05, -3.2079251781364850e-05}  /*      1 */
};

static const astro_top_series_t topseries_8_1[] =
{
    {     90, topterms_8_1_0 }
,   {     32, topterms_8_1_1 }
,   {     11, topterms_8_1_2 }
,   {      5, topterms_8_1_3 }
,   {      2, topterms_8_1_4 }
};

static const astro_top_term_t topterms_8_2_0[] =
{
    {        0, -1.7873895940349999e-01,  0.0000000000000000e+00}  /*      0 */
,   {     1473,  3.1629832749992988e-03, -2.0985870294942081e-03}  /*      1 */
,   {       71, -6.8180357663860398e-04,  1.1113519163940930e-03}  /*      2 */
,   {     1331,  1.6911781266743710e-04,  1.1885125258969610e-03}  /*      3 */
,   {      593,  5.4774166542017916e-04, -6.3844158559171749e-04}  /*      4 */
,   {     1402, -2.1438473609329679e-05, -5.2026929830140383e-04}  /*      5 */
,   {     1261, -5.7368849879358177e-05,  4.6888445010058260e-04}  /*      6 */
,   {      141, -9.6040157452198532e-05,  3.0415286563684890e-04}  /*      7 */
,   {      452,  1.2124121714193030e-04,  2.7470351829330812e-04}  /*      8 */
,   {     2945,  1.1035693239374779e-04, -1.4867545478165761e-04}  /*      9 */
,   {     1190, -5.9767419503565077e-05,  1.5052450817429230e-04}  /*     10 */
,   {     1543, -1.0081841648710751e-04,  1.1874253775582299e-04}  /*     11 */
,   {      522, -3.7998691236472437e-05, -1.1710923038633279e-04}  /*     12 */
,   {      381,  1.7939534333569889e-05,  1.2165766508553840e-04}  /*     13 */
,   {        8, -5.1065826531599063e-05, -1.0380178618772450e-04}  /*     14 */
,   {        4, -9.6101012457693182e-05, -3.7409967137145340e-05}  /*     15 */
,   {      212, -1.0011374095546651e-05,  9.1861092027193070e-05}  /*     16 */
,   {      208,  6.3118700287276614e-05,  6.6551813665925998e-05}  /*     17 */
,   {      106,  4.6362839297978442e-05,  7.1895265631125496e-05}  /*     18 */
,   {     2804,  2.5951137645879960e-05,  4.8810749391069221e-05}  /*     19 */
,   {     1119, -3.1798344550251363e-05,  4.3442443246298743e-05}  /*     20 */
,   {       35, -2.6365309213148771e-05, -3.7318748312868848e-05}  /*     21 */
,   {     1186,  4.5373467873211630e-05, -3.9315678706168707e-06}  /*     22 */
,   {      310, -5.7440916735318752e-06,  4.2608008327900133e-05}  /*     23 */
,   {       67, -3.2628206231552743e-05,  1.4188485967101389e-05}  /*     24 */
,   {      664, -1.3138986439777871e-05,  2.8745825031870391e-05}  /*     25 */
,   {    17476, -4.5717070878285386e-06, -2.7229104695341480e-05}  /*     26 */
,   {      283,  4.1649905603620409e-06,  2.5778764278582031e-05}  /*     27 */
,   {    28407, -2.6086891392275169e-05,  6.5019643042548361e-07}  /*     28 */
,   {     2875, -1.2604181155312930e-05, -2.2155626126791911e-05}  /*     29 */
,   {     2733,  4.6350251894966598e-06,  2.1642396781861640e-05}  /*     30 */
,   {       12,  1.3787761564064851e-05, -1.6882625493196851e-05}  /*     31 */
,   {       63,  7.3243482301750183e-06, -2.0430415756143061e-05}  /*     32 */
,   {      137,  1.7688914790024399e-05, -2.6466558968384412e-06}  /*     33 */
,   {     1048, -1.3743478634746570e-05,  1.1178051438197709e-05}  /*     34 */
,   {     1614, -1.4527685432017941e-05, -2.5021418914365550e-06}  /*     35 */
,   {     1044, -5.5188522930647246e-06,  1.3429997819426629e-05}  /*     36 */
,   {      239, -6.3224380303283023e-06,  1.2677314304305631e-05}  /*     37 */
,   {      133,  2.9376054111895790e-06, -1.3077911011252820e-05}  /*     38 */
,   {     1492, -9.7624713292315362e-06,  4.8639192065138997e-06}  /*     39 */
,   {     1454,  8.1915641060518482e-06, -7.1484181790528786e-06}  /*     40 */
,   {     4418,  2.9126133839920300e-06, -9.7148097661500440e-06}  /*     41 */
,   {      177, -8.1108828396444682e-06, -5.8479936691072763e-06}  /*     42 */
,   {     2663, -7.4109426380852051e-07,  8.1320485449777430e-06}  /*     43 */
,   {    17334,  7.7062782105104396e-06,  2.1832991677727231e-06}  /*     44 */
,   {     3016, -2.7097054305923721e-06,  7.2985152622098958e-06}  /*     45 */
,   {    28266,  3.0978305562474201e-06, -6.9523478329398256e-06}  /*     46 */
,   {       75, -7.5005963705916821e-06,  5.8552558907898478e-07}  /*     47 */
,   {      354,  2.2767150385710288e-06,  7.1310390742790356e-06}  /*     48 */
,   {      974, -2.8814074492595782e-06,  6.7739175425797296e-06}  /*     49 */
,   {       59, -4.8255490447201409e-06, -5.3633858025303309e-06}  /*     50 */
,   {     1115,  2.4165363833927252e-06, -5.9760260615249186e-06}  /*     51 */
,   {      204, -1.5420460797033710e-07, -6.3645330720578637e-06}  /*     52 */
,   {      612,  4.7791204064686441e-06, -4.1989180399775663e-06}  /*     53 */
,   {      574, -3.2533619740383079e-06,  4.8942036839749464e-06}  /*     54 */
,   {      978, -5.2816161970547299e-06,  2.4194318535884879e-06}  /*     55 */
,   {      129, -4.1704589160563732e-06, -3.8074741128128591e-06}  /*     56 */
,   {      200, -4.3860533690857434e-06, -2.8426817721421870e-06}  /*     57 */
,   {     1685, -4.4937857977759297e-06, -2.6664132359831861e-06}  /*     58 */
,   {     1327, -4.0394363815867892e-06,  3.0479150446235980e-06}  /*     59 */
,   {     1257, -4.4570496522011644e-06,  1.8417244601604851e-06}  /*     60 */
};

static const astro_top_term_t topterms_8_2_1[] =
{
    {        0, -6.1339663802007860e-04,  0.0000000000000000e+00}  /*      0 */
,   {     1331,  5.6664110172313892e-04, -8.1133128701713292e-05}  /*      1 */
,   {     1473, -1.9632805654171590e-04, -2.9843120844069362e-04}  /*      2 */
,   {       71, -1.9721831884265259e-04, -1.2980246221969450e-04}  /*      3 */
,   {     1402, -1.4909525632022481e-04,  5.5470899732068926e-06}  /*      4 */
,   {     1261,  1.4379398981882231e-04,  1.8301247697520619e-05}  /*      5 */
,   {     2945, -7.2143717725975376e-05, -6.1374574108146777e-05}  /*      6 */
,   {     1190,  7.4611359353024315e-05,  3.0289539777815409e-05}  /*      7 */
,   {      593, -5.9939495914062063e-05, -5.1836225465555612e-05}  /*      8 */
,   {        8,  3.9984073707170713e-05, -3.0799116570219068e-05}  /*      9 */
,   {     1543,  3.1094179263405468e-05,  2.6861784078602320e-05}  /*     10 */
,   {      381,  3.7436713283265682e-05, -5.1576484771147526e-06}  /*     11 */
,   {     1119,  2.9759191459225299e-05,  2.2175266912892650e-05}  /*     12 */
,   {      452,  3.2286698454013622e-05, -1.4330028054859710e-05}  /*     13 */
,   {      522, -3.3610639556594409e-05,  1.0732560764842151e-05}  /*     14 */
,   {     2804,  2.6892727962532569e-05, -1.2202088159044841e-05}  /*     15 */
,   {        4,  1.6145362682689730e-05, -2.1535180644982459e-05}  /*     16 */
,   {       35,  1.1619842301292181e-05,  2.1259793782462840e-05}  /*     17 */
,   {      310,  2.1159419772930708e-05,  3.2365387780938209e-06}  /*     18 */
,   {      212, -1.7383881837993001e-05, -1.2570139978872199e-07}  /*     19 */
,   {     2733,  1.5857546421435610e-05, -2.6232455173985150e-06}  /*     20 */
,   {     1048,  9.8372722970562282e-06,  1.2182701032807050e-05}  /*     21 */
,   {       12,  1.1952587407675961e-05,  7.5377841861106607e-06}  /*     22 */
,   {      283, -1.0224581925761621e-05,  1.2303488561812810e-06}  /*     23 */
,   {      239,  8.6293198115122691e-06,  4.6174857696641109e-06}  /*     24 */
};

static const astro_top_term_t topterms_8_2_2[] =
{
    {     1331,  2.3820657956875651e-05, -1.4117217046583029e-04}  /*      0 */
,   {        0,  5.7169784317623151e-05,  0.0000000000000000e+00}  /*      1 */
,   {     1261,  2.8584658843704731e-05, -1.9153125455819128e-05}  /*      2 */
,   {       71, -6.3863669075967202e-06, -3.0182516098252429e-05}  /*      3 */
,   {     2945, -1.6831473476249729e-05,  1.7778057659448761e-05}  /*      4 */
,   {     1402, -8.7282659443619542e-06,  2.1969730180690399e-05}  /*      5 */
,   {     1190,  1.8655664239774939e-05, -1.4275349333837280e-05}  /*      6 */
,   {        8,  2.0572423386833478e-05,  2.9970650695923720e-06}  /*      7 */
};

static const astro_top_series_t topseries_8_2[] =
{
    {     61, topterms_8_2_0 }
,   {     25, topterms_8_2_1 }
,   {      8, topterms_8_2_2 }
};

static const astro_top_term_t topterms_8_3_0[] =
{
    {        0, -1.7340471864230000e-01,  0.0000000000000000e+00}  /*      0 */
,   {     1473,  2.1234813115076170e-03,  3.0130058621335031e-03}  /*      1 */
,   {       71, -1.0505499200359431e-03, -7.0849649152681434e-04}  /*      2 */
,   {     1331,  1.1926139685020671e-03, -1.2467514555563171e-04}  /*      3 */
,   {      593,  6.3328818822552336e-04,  5.1838245948530805e-04}  /*      4 */
,   {     1402, -4.0198902705438299e-04,  3.4879411820834029e-04}  /*      5 */
,   {     1261,  4.6754614871725769e-04,  6.5799333440031267e-05}  /*      6 */
,   {      141, -3.1574545029112658e-04, -1.1057917696717649e-04}  /*      7 */
,   {      452,  2.7842572840314938e-04, -1.1076531434037220e-04}  /*      8 */
,   {     2945,  1.4727652777389929e-04,  1.0316600256389621e-04}  /*      9 */
,   {     1190,  1.4974726446281261e-04,  6.1544751468406387e-05}  /*     10 */
,   {     1543, -1.1418008898386521e-04, -9.8159556344747838e-05}  /*     11 */
,   {      522, -6.9432011255429296e-05,  1.0500648997119170e-04}  /*     12 */
,   {      381,  1.2168346615540400e-04, -1.5724958971846889e-05}  /*     13 */
,   {        4,  3.5694868279606838e-05, -1.1563400718698000e-04}  /*     14 */
,   {        8,  1.0601495232676501e-04, -5.3792770194159239e-05}  /*     15 */
,   {      212, -8.8206829660081605e-05, -5.6186044809908770e-06}  /*     16 */
,   {      208, -6.2116275994912254e-05,  6.2681468828265711e-05}  /*     17 */
,   {      106, -5.9881700950469839e-05,  1.8173312312951261e-05}  /*     18 */
,   {     2804,  4.9734881336053070e-05, -2.4663551556821301e-05}  /*     19 */
,   {     1119,  4.3142176920002318e-05,  3.2166118196149167e-05}  /*     20 */
,   {     1186,  5.3160131153765831e-06,  4.6681999520445390e-05}  /*     21 */
,   {      310,  4.2381917128956993e-05,  6.2114355060802278e-06}  /*     22 */
,   {       67,  9.7801450369382322e-06,  3.8008739374555199e-05}  /*     23 */
,   {       35,  3.4777523130477253e-05,  8.0807437621346330e-06}  /*     24 */
,   {      664, -2.7772391338347379e-05, -1.2933056019872690e-05}  /*     25 */
,   {    17476,  2.6116611382398249e-05, -5.2512707918252350e-06}  /*     26 */
,   {      283, -2.5901894192097461e-05,  7.6613799535167943e-07}  /*     27 */
,   {    28407, -1.3727149833568929e-06, -2.5571899427495732e-05}  /*     28 */
,   {       12,  1.7840323650284740e-05,  1.4017491541647661e-05}  /*     29 */
,   {     2875, -1.2752702589739361e-05,  1.8674367648851481e-05}  /*     30 */
,   {       63, -2.1141833621778831e-05, -6.5305966080253417e-06}  /*     31 */
,   {     2733,  2.1699001266656308e-05, -4.2325482084098446e-06}  /*     32 */
,   {     1048,  1.1071845430997370e-05,  1.3823004095404599e-05}  /*     33 */
,   {      137, -1.3154331920868999e-05, -9.3564722901572934e-06}  /*     34 */
,   {     1044,  1.3602483536936790e-05,  6.6866532885834768e-06}  /*     35 */
,   {     1614,  3.0521160394305390e-06, -1.4295846223173739e-05}  /*     36 */
,   {      239,  1.2491201410693211e-05,  6.3469005652658392e-06}  /*     37 */
,   {      133, -1.3409179117794041e-05, -2.6113349355778870e-06}  /*     38 */
,   {     1492, -4.9926812258289137e-06, -9.3442347365693383e-06}  /*     39 */
,   {     1454,  7.1182846121023137e-06,  7.7835699509488784e-06}  /*     40 */
,   {     4418,  9.5085270250543767e-06,  2.5720788995216492e-06}  /*     41 */
,   {      354, -7.7105865781446408e-06,  3.5195735852526461e-06}  /*     42 */
,   {     2663,  8.0823599962064017e-06,  8.8289336639692517e-07}  /*     43 */
,   {       75, -6.1089854320482681e-07, -8.0284473434343767e-06}  /*     44 */
,   {    17334,  2.4666413292489791e-06, -7.6105265820617346e-06}  /*     45 */
,   {      974,  6.9130746927758178e-06,  3.5048585089334658e-06}  /*     46 */
,   {    28266, -6.8235793370187026e-06, -3.3487790585021339e-06}  /*     47 */
,   {     3016, -7.0427593206860361e-06, -2.6580682576526168e-06}  /*     48 */
,   {       59, -5.5456494164143571e-06,  4.9667283828556126e-06}  /*     49 */
,   {      204, -7.0017722378178446e-06, -1.6817751033546979e-06}  /*     50 */
,   {     1115, -6.4024937184914472e-06,  2.1192315231268140e-06}  /*     51 */
,   {      978,  2.3623585675585379e-06,  5.3260537432133881e-06}  /*     52 */
,   {      129, -3.9867297079366098e-06,  4.2341276287979132e-06}  /*     53 */
,   {      574, -4.8314043929159330e-06, -2.9544972252560498e-06}  /*     54 */
,   {      612,  4.0663883671886490e-06,  3.8452439907704313e-06}  /*     55 */
,   {      200, -3.2782990554685469e-06,  4.3746397303989681e-06}  /*     56 */
,   {      247,  4.4086088444075100e-06, -3.0262449322858620e-06}  /*     57 */
,   {     1685,  2.7665561507113492e-06, -4.4616466628081997e-06}  /*     58 */
,   {     1335, -2.3134382826680641e-06,  4.5013914052228881e-06}  /*     59 */
,   {     1327,  3.1175433103076271e-06,  3.7924860549293981e-06}  /*     60 */
,   {      903,  4.4741994730303496e-06,  1.7822080022821150e-06}  /*     61 */
,   {       16, -1.4927342901310379e-06,  4.4884086825602869e-06}  /*     62 */
,   {      169,  2.5354018835911621e-06,  3.5730000800789471e-06}  /*     63 */
,   {      271, -2.1688629242033941e-06,  3.7820832945120280e-06}  /*     64 */
,   {       79,  3.8979055894631164e-06, -1.8923407502835230e-06}  /*     65 */
,   {      416,  3.9376438761867171e-06, -5.2538703904847715e-07}  /*     66 */
,   {    17405,  1.0904344065136349e-06,  3.5409195078133210e-06}  /*     67 */
,   {    28337,  3.4736718287656798e-06, -5.2407682497120616e-07}  /*     68 */
,   {     1351, -3.4436468577006500e-06, -7.4140986208801818e-08}  /*     69 */
,   {     1312,  3.3433376977429141e-06, -7.9547674738811048e-07}  /*     70 */
,   {      832,  3.0647904819735712e-06,  1.2969722855848811e-06}  /*     71 */
,   {      275,  1.5394331893610239e-06, -2.7751033895625021e-06}  /*     72 */
,   {     2592,  2.8273974520595848e-06,  1.3333677990888271e-06}  /*     73 */
,   {    17264,  1.6047868553221060e-06, -2.5731990436490451e-06}  /*     74 */
};

static const astro_top_term_t topterms_8_3_1[] =
{
    {     1331, -6.0044174932694893e-05, -5.6841814115066798e-04}  /*      0 */
,   {     1473,  2.8299001111328808e-04, -1.9899683117934169e-04}  /*      1 */
,   {       71,  1.3109847892346901e-04, -1.9645547029607899e-04}  /*      2 */
,   {     1402,  9.8753503746441094e-05,  1.1624141830784459e-04}  /*      3 */
,   {     1261,  2.0857393493693379e-05, -1.4333407898431410e-04}  /*      4 */
,   {     2945,  5.7600265965337831e-05, -7.1668848804472786e-05}  /*      5 */
,   {     1190,  3.1158128562030883e-05, -7.4221508541133025e-05}  /*      6 */
,   {      593,  4.9568824746840162e-05, -5.9443745710499677e-05}  /*      7 */
,   {        0, -6.3916210233836622e-05,  0.0000000000000000e+00}  /*      8 */
,   {        8,  3.2326459763508557e-05,  4.1097067153751401e-05}  /*      9 */
,   {        4,  3.8505374914800768e-05, -2.4690750688454580e-05}  /*     10 */
,   {     1543, -2.6127774810849349e-05,  2.9873942333388641e-05}  /*     11 */
,   {      381, -4.4772545431527864e-06, -3.7396119251182150e-05}  /*     12 */
,   {     1119,  2.2430310753365549e-05, -2.9557761386193919e-05}  /*     13 */
,   {      522,  2.9886224978272480e-05,  2.0259791365064171e-05}  /*     14 */
,   {      452, -1.3106498259692980e-05, -3.2637289766253343e-05}  /*     15 */
,   {     2804, -1.1508149132666279e-05, -2.7316790371599908e-05}  /*     16 */
,   {      310,  3.4738350608151140e-06, -2.1017437982180379e-05}  /*     17 */
,   {      212,  2.1402716283134239e-06, -1.9178262093378069e-05}  /*     18 */
,   {     2733, -2.3283762042050202e-06, -1.5882216399278110e-05}  /*     19 */
,   {     1048,  1.2246829111767810e-05, -9.7512980427413643e-06}  /*     20 */
,   {       12, -7.6135340345472092e-06,  1.2831304341856990e-05}  /*     21 */
,   {       35, -5.7732549173295510e-06,  1.0016273100035810e-05}  /*     22 */
,   {      283, -2.2363814820487809e-06, -9.6218108123675468e-06}  /*     23 */
,   {      239,  4.6765179736328208e-06, -8.5111533123445213e-06}  /*     24 */
,   {      106, -1.9836445213635969e-06,  7.8016311549488939e-06}  /*     25 */
,   {     2875,  5.8334802150818700e-06,  5.1997278670085084e-06}  /*     26 */
,   {     1044,  2.4978011849958161e-06, -6.6276491176842909e-06}  /*     27 */
};

static const astro_top_term_t topterms_8_3_2[] =
{
    {     1331, -1.3993900805601181e-04, -2.9006879879181548e-05}  /*      0 */
,   {        0,  6.7829783611580237e-05,  0.0000000000000000e+00}  /*      1 */
,   {     1261, -1.8612445466711130e-05, -2.8908877766737650e-05}  /*      2 */
,   {       71,  2.9853088727300880e-05, -7.1012915175069879e-06}  /*      3 */
,   {        4,  2.6321188906063949e-05,  1.3464234343799510e-05}  /*      4 */
,   {     1402,  2.3284325058403039e-05, -6.8144038337101704e-06}  /*      5 */
,   {     2945, -1.7707706192590389e-05, -1.5833741011362669e-05}  /*      6 */
,   {     1190, -1.4044910483808361e-05, -1.8816964596903641e-05}  /*      7 */
,   {        8, -3.2757789019448882e-06,  2.1268129579065079e-05}  /*      8 */
,   {     1473, -9.8933125150553904e-06, -1.3083286541079480e-05}  /*      9 */
,   {     1119, -7.2140288335202669e-06, -1.1729524860644351e-05}  /*     10 */
};

static const astro_top_term_t topterms_8_3_3[] =
{
    {     1331, -1.7737816039767501e-05,  2.8055420335966871e-05}  /*      0 */
};

static const astro_top_series_t topseries_8_3[] =
{
    {     75, topterms_8_3_0 }
,   {     28, topterms_8_3_1 }
,   {     11, topterms_8_3_2 }
,   {      1, topterms_8_3_3 }
};

static const astro_top_term_t topterms_8_4_0[] =
{
    {        0, -5.1702307822780000e-02,  0.0000000000000000e+00}  /*      0 */
,   {     1402, -1.3296787157385281e-04,  1.3523784099823399e-04}  /*      1 */
,   {     1543,  1.6300899640352639e-04,  5.5927755421839437e-05}  /*      2 */
,   {     1473, -2.3179975765177140e-05, -9.8184622098344091e-05}  /*      3 */
,   {      522, -1.9097390976957099e-05,  3.6607666661423369e-05}  /*      4 */
,   {      664,  3.2543564829276270e-05,  1.0672339397019160e-06}  /*      5 */
,   {     1331, -2.0841818514170151e-05,  1.2410223144740100e-05}  /*      6 */
,   {        4, -8.0356073429830634e-06,  1.9857752253856971e-05}  /*      7 */
,   {      593, -1.0459691293338759e-05, -1.7857031017673810e-05}  /*      8 */
,   {     1614,  2.0002336741252269e-05,  1.4700260033157791e-06}  /*      9 */
,   {     2875, -3.5880881863000781e-06,  8.0256130260543578e-06}  /*     10 */
,   {     3016,  8.5254320868063325e-06, -1.8479247218425171e-07}  /*     11 */
,   {        8, -5.4541260810660720e-06,  3.9805913226076404e-06}  /*     12 */
,   {      452, -3.7355080445397231e-06,  4.0025013123176398e-06}  /*     13 */
,   {       35, -2.1277230177098351e-06, -4.9607027782605714e-06}  /*     14 */
,   {      137, -3.8233051639183096e-06, -3.1849844866680141e-06}  /*     15 */
};

static const astro_top_term_t topterms_8_4_1[] =
{
    {        0,  1.9166844047183211e-04,  0.0000000000000000e+00}  /*      0 */
,   {     1402,  3.9019701387765488e-05,  3.8671378812923788e-05}  /*      1 */
,   {     1543,  1.5119207079965939e-05, -4.3332491992581899e-05}  /*      2 */
,   {     1331,  5.8922122899993431e-06,  1.0029385340576861e-05}  /*      3 */
,   {      522,  1.0280806202111021e-05,  5.2916455231996232e-06}  /*      4 */
};

static const astro_top_series_t topseries_8_4[] =
{
    {     16, topterms_8_4_0 }
,   {      5, topterms_8_4_1 }
};

static const astro_top_term_t topterms_8_5_0[] =
{
    {        0,  1.3977992515640000e-01,  0.0000000000000000e+00}  /*      0 */
,   {     1402,  1.2883200572372361e-04,  1.3471156433694759e-04}  /*      1 */
,   {     1543, -4.9979310098879121e-05,  1.6180896042292819e-04}  /*      2 */
,   {     1473, -2.2060691195269469e-05, -9.3668360843677731e-05}  /*      3 */
,   {      522,  3.5537595493916502e-05,  1.9797390613745271e-05}  /*      4 */
,   {      664, -7.4709042862279015e-08,  3.1983280221550278e-05}  /*      5 */
,   {     1331,  1.1695376261855580e-05,  2.0798562232570440e-05}  /*      6 */
,   {      593, -9.9523808317349038e-06, -1.7074224315439731e-05}  /*      7 */
,   {     1614, -9.4118113483089596e-07,  1.9684674172624949e-05}  /*      8 */
,   {        4, -1.5932093431444279e-05,  9.4977011164352140e-06}  /*      9 */
,   {     2875,  7.8417800319575014e-06,  3.7228551547920001e-06}  /*     10 */
,   {     3016,  4.3444197029476981e-07,  8.3694483911743903e-06}  /*     11 */
,   {        8, -5.9902893836260643e-06, -4.3498982251067908e-06}  /*     12 */
,   {      137, -3.0046395968887028e-06,  4.7462972937701063e-06}  /*     13 */
,   {      452,  3.9819715036126689e-06,  3.7540995120279268e-06}  /*     14 */
,   {       35, -4.0365965056395762e-08,  5.3159364396026402e-06}  /*     15 */
,   {      141,  1.2649492629014499e-07,  4.8536627302669076e-06}  /*     16 */
,   {     1261,  1.1901325673451860e-06,  4.3994722705730970e-06}  /*     17 */
,   {     2945, -2.7264749792923272e-06, -3.6454310686215340e-06}  /*     18 */
,   {     1685,  6.9410338646847810e-07,  3.3566050467992739e-06}  /*     19 */
,   {      734,  8.2952840171332113e-07,  3.2356209378486682e-06}  /*     20 */
,   {      208,  2.4256846662067760e-06, -9.5181484992359002e-07}  /*     21 */
,   {      279, -2.5271926303381422e-06, -3.7922031983877281e-07}  /*     22 */
,   {     1115,  5.6515925473070844e-07,  2.1078436053637770e-06}  /*     23 */
,   {     1257, -1.3558265037946960e-06,  1.4215822110290420e-06}  /*     24 */
,   {       12,  1.3204389683748120e-07, -1.8835278008457151e-06}  /*     25 */
,   {       71, -1.9042145379339269e-07, -1.8630215116559800e-06}  /*     26 */
,   {      212,  1.6053042001407380e-06,  8.4671630560881122e-07}  /*     27 */
,   {    17405,  1.3113174080732559e-06, -4.9311123653462177e-07}  /*     28 */
,   {    17547,  1.0602949344425471e-06,  8.9833952832208403e-07}  /*     29 */
,   {     2804,  9.2561678972169579e-07,  9.7236819416524451e-07}  /*     30 */
,   {    28337, -9.4657252626238933e-08, -1.1160255453211769e-06}  /*     31 */
,   {    28478,  9.2430781502463316e-07, -6.2265716348343233e-07}  /*     32 */
,   {       75, -9.2130387409465058e-07,  5.9403283989994324e-07}  /*     33 */
,   {      106, -2.0304770919600219e-07, -1.0388792453459521e-06}  /*     34 */
,   {      381,  5.9915368756519242e-07,  8.3180916441107633e-07}  /*     35 */
,   {     1190,  2.7201162062442490e-08,  1.0190741100380801e-06}  /*     36 */
,   {     3087,  2.9873098731161120e-07,  9.6483981362737291e-07}  /*     37 */
,   {       63,  5.8956418822020129e-07, -7.7745585971202028e-07}  /*     38 */
,   {     1186,  2.4117256641652101e-07, -9.2801537744071498e-07}  /*     39 */
,   {       67, -3.2941946279216570e-07,  7.9963658829793748e-07}  /*     40 */
,   {    17476, -7.4361213537121858e-07, -1.2891468267067801e-07}  /*     41 */
,   {      283,  2.8076344536095831e-08,  7.3853891759583075e-07}  /*     42 */
,   {     1756,  3.1161152193772379e-07,  6.3104215205486464e-07}  /*     43 */
,   {      247, -3.0658465478480348e-07, -6.2524043128321752e-07}  /*     44 */
,   {    28407, -2.6107491970896922e-07,  5.4499817140536954e-07}  /*     45 */
,   {      354,  4.9460726403829815e-07,  2.9447867230664199e-07}  /*     46 */
,   {     1421, -3.1957781759109811e-07, -4.3033184225872391e-07}  /*     47 */
,   {     1383,  4.0796989384481828e-07,  3.2647659230997358e-07}  /*     48 */
,   {      133,  1.1815490710962760e-07, -5.0524398677283852e-07}  /*     49 */
,   {      177, -4.6755777925983522e-07,  2.0641392675487111e-07}  /*     50 */
,   {      805,  2.4510909183821792e-07,  4.3308098936875489e-07}  /*     51 */
,   {     1563,  2.0041064699998999e-07, -4.4445832345851582e-07}  /*     52 */
,   {       59, -1.2658294965005670e-07, -4.6671870382330612e-07}  /*     53 */
,   {     1524, -7.8614887212231827e-08,  4.6760439016103602e-07}  /*     54 */
,   {     4348,  4.6767825285878958e-07,  3.5290700804681861e-08}  /*     55 */
,   {       16,  4.2127049178471667e-07, -1.8684252468076170e-07}  /*     56 */
,   {     4489,  1.8055735116211371e-07,  4.2107641051871891e-07}  /*     57 */
,   {       79, -4.0548063730912672e-07, -2.1065225409216059e-07}  /*     58 */
,   {     1398, -1.2435454796683390e-07,  3.9919922981210330e-07}  /*     59 */
,   {      145, -3.0487589061275071e-07,  2.8316787420643652e-07}  /*     60 */
,   {      345,  3.4373900983948560e-07, -2.2193498580565619e-07}  /*     61 */
,   {      275, -7.6720036516022037e-09, -4.0650745658300768e-07}  /*     62 */
,   {     1406, -3.7816142635595027e-07,  9.1996049592574050e-08}  /*     63 */
,   {     1044,  7.4671401050977241e-08,  3.7193907363226299e-07}  /*     64 */
,   {      318,  7.1927751245935557e-08,  3.6928457890682560e-07}  /*     65 */
,   {      416,  4.1798174498739052e-08, -3.6796074151811958e-07}  /*     66 */
,   {      204, -1.5624747021437499e-08, -3.6071195814331651e-07}  /*     67 */
,   {     1539,  3.3199296973509958e-07, -1.0438655407426570e-07}  /*     68 */
,   {     1547,  2.1465474690281490e-07,  2.6396972342688539e-07}  /*     69 */
,   {     1327, -1.5719732197123969e-07,  2.8693854391332578e-07}  /*     70 */
,   {      541,  2.6198971229951878e-07,  1.7511268865594321e-07}  /*     71 */
};

static const astro_top_term_t topterms_8_5_1[] =
{
    {        0,  7.7993297015907634e-05,  0.0000000000000000e+00}  /*      0 */
,   {     1402,  3.9152687150310170e-05, -3.7172486736475712e-05}  /*      1 */
,   {     1543,  4.3031136375758243e-05,  1.3517906141260201e-05}  /*      2 */
,   {     1331,  1.0005173914454110e-05, -5.5510992766576859e-06}  /*      3 */
,   {      522,  5.4799934835759310e-06, -9.9803123122233115e-06}  /*      4 */
,   {     1473, -9.2379228749259814e-06,  2.0049455605464851e-06}  /*      5 */
,   {      664, -3.2895419619690672e-06, -4.5974891215164072e-08}  /*      6 */
,   {        8,  1.3908106265782320e-06, -2.9292856582721550e-06}  /*      7 */
,   {     2875,  1.5562630279974180e-06, -2.5756858183805818e-06}  /*      8 */
,   {        4, -2.4245388783587939e-06, -1.5068482062721909e-06}  /*      9 */
,   {     3016,  2.7144173403334680e-06,  1.5168858882407361e-07}  /*     10 */
,   {     2945, -1.9556025171954529e-06,  1.2707620314125440e-06}  /*     11 */
,   {       35,  1.3344840689621260e-06,  1.4001396639490341e-06}  /*     12 */
,   {      593, -1.5045227426363480e-06,  9.0293863845346919e-07}  /*     13 */
,   {     1614,  1.4931074066475380e-06,  1.3901882598918031e-07}  /*     14 */
,   {     1261,  1.3653087138211129e-06, -3.4669288959066708e-07}  /*     15 */
,   {       12,  1.2315835838874919e-06, -5.9561839689487617e-08}  /*     16 */
,   {      137,  9.0428610682743405e-07,  5.7865365381345505e-07}  /*     17 */
,   {     2804,  5.6172637460901586e-07, -4.7331838000806132e-07}  /*     18 */
,   {      212, -6.0813020695312143e-07,  1.1354550736721950e-07}  /*     19 */
};

static const astro_top_term_t topterms_8_5_2[] =
{
    {     1402, -2.9412500484509639e-06, -8.0598049136548207e-06}  /*      0 */
,   {     1543, -1.1369975198335059e-06, -6.6800431913703534e-06}  /*      1 */
,   {     1331, -5.6564967655275266e-07, -2.8456118566588939e-06}  /*      2 */
,   {        0,  2.4849999999999999e-06,  0.0000000000000000e+00}  /*      3 */
,   {      522, -1.0188268325841671e-06, -1.4049711115470319e-06}  /*      4 */
,   {        4, -1.2789041308060199e-06, -7.7265243668823065e-07}  /*      5 */
,   {        8,  1.2656119872224440e-06, -4.3448350220923409e-07}  /*      6 */
};

static const astro_top_series_t topseries_8_5[] =
{
    {     72, topterms_8_5_0 }
,   {     20, topterms_8_5_1 }
,   {      7, topterms_8_5_2 }
};

static const astro_top_model_t topmodel_8[] =
{
    {  3, topseries_8_0 }
,   {  5, topseries_8_1 }
,   {  3, topseries_8_2 }
,   {  4, topseries_8_3 }
,   {  2, topseries_8_4 }
,   {  3, topseries_8_5 }
};




static top_elliptical_t TopCalcElliptical(int planet, const astro_top_model_t *model, double tt)
{
    /* Translated from: TOP2013.f */
    /* See: https://github.com/cosinekitty/ephemeris/tree/master/top2013 */
    /* Copied from: ftp://ftp.imcce.fr/pub/ephem/planets/top2013 */
    static const double freq[] =
    {
        0.5296909622785881e+03,
        0.2132990811942489e+03,
        0.7478166163181234e+02,
        0.3813297236217556e+02,
        0.2533566020437000e+02
    };
    int i, f, s, t;
    double time[TOP_NSERIES];
    double el[6];
    double arg, dmu, xl;
    top_elliptical_t ellip;

    /* Time */
    time[0] = 1.0;
    time[1] = tt / 365250.0;
    for (i=1; i < TOP_NSERIES; ++i)
        time[i] = time[i-1] * time[1];

    dmu = (freq[0] - freq[1]) / 880.0;

    for (f=0; f < TOP_NCOORDS; ++f)
    {
        el[f] = 0.0;
        for (s=0; s < model[f].nseries; ++s)
        {
            const astro_top_series_t *series = &model[f].series[s];
            for (t=0; t < series->nterms; ++t)
            {
                const astro_top_term_t *term = &series->terms[t];
                if (f==1 && s==1 && term->k==0)
                    continue;
                arg = term->k * dmu * time[1];
                el[f] += time[s] * (term->c*cos(arg) + term->s*sin(arg));
            }
        }
    }

    xl = el[1] + freq[planet - 5] * time[1];
    xl = fmod(xl, PI2);
    if (xl < 0.0)
        xl += PI2;
    el[1] = xl;

    /* Convert elliptical elements from array 'el' to friendly struct layout. */
    ellip.a      = el[0];
    ellip.lambda = el[1];
    ellip.k      = el[2];
    ellip.h      = el[3];
    ellip.q      = el[4];
    ellip.p      = el[5];
    return ellip;
}


static void TopEcliptic(const top_elliptical_t *ellip, astro_vector_t *ecl)
{
    double xa, xl, xk, xh, xq, xp;
    double xfi, xki, u, ex, ex2, ex3;
    double zr, zi;
    double z1r, z1i;
    double z2r, z2i;
    double z3r, z3i;
    double zteta_r, zteta_i;
    double zto_r, zto_i;
    double gl, gm, e, dl, rsa;
    double xcw, xsw, xm, xr;

    xa = ellip->a;
    xl = ellip->lambda;
    xk = ellip->k;
    xh = ellip->h;
    xq = ellip->q;
    xp = ellip->p;

    xfi = sqrt(1.0 - xk*xk - xh*xh);
    xki = sqrt(1.0 - xq*xq - xp*xp);
    zr = xk; zi = xh;       /* z = dcmplx(xk,xh) */
    u = 1.0 / (1.0 + xfi);
    ex2 = zr*zr + zi*zi;
    ex = sqrt(ex2);         /* ex = cdabs(z) */
    ex3 = ex * ex2;
    z1r = zr; z1i = -zi;

    gl = fmod(xl, PI2);
    gm = gl - atan2(xh, xk);
    e = gl + (ex - 0.125*ex3)*sin(gm) + 0.5*ex2*sin(2.0*gm) + 0.375*ex3*sin(3.0*gm);

    do
    {
        z2r = 0.0; z2i = e;
        zteta_r = cos(z2i);
        zteta_i = sin(z2i);
        z3r = z1r*zteta_r - z1i*zteta_i;
        z3i = z1r*zteta_i + z1i*zteta_r;
        dl = gl - e + z3i;
        rsa = 1.0 - z3r;
        e += dl/rsa;
    } while (fabs(dl) >= 1.0e-15);

    z1r = z3i * u * zr;
    z1i = z3i * u * zi;
    z2r = +z1i;
    z2i = -z1r;
    zto_r = (-zr + zteta_r + z2r) / rsa;
    zto_i = (-zi + zteta_i + z2i) / rsa;
    xcw = zto_r;
    xsw = zto_i;
    xm = xp*xcw - xq*xsw;
    xr = xa*rsa;

    ecl->x = xr*(xcw - 2.0*xp*xm);
    ecl->y = xr*(xsw + 2.0*xq*xm);
    ecl->z = -2.0*xr*xki*xm;
}



static void TopEquatorial(const astro_vector_t *ecl, astro_vector_t *equ)
{
    static int initialized;
    static double rot[3][3];

    if (!initialized)
    {
        const double sdrad = DEG2RAD / 3600.0;
        const double eps = (23.0 + 26.0/60.0 + 21.41136/3600.0)*DEG2RAD;
        const double phi = -0.05188 * sdrad;
        const double ceps = cos(eps);
        const double seps = sin(eps);
        const double cphi = cos(phi);
        const double sphi = sin(phi);

        rot[0][0] =  cphi;
        rot[0][1] = -sphi*ceps;
        rot[0][2] =  sphi*seps;
        rot[1][0] =  sphi;
        rot[1][1] =  cphi*ceps;
        rot[1][2] = -cphi*seps;
        rot[2][0] =  0.0;
        rot[2][1] =  seps;
        rot[2][2] =  ceps;

        initialized = 1;
    }

    equ->status = ecl->status;
    equ->t = ecl->t;
    equ->x = (rot[0][0] * ecl->x) + (rot[0][1] * ecl->y) + (rot[0][2] * ecl->z);
    equ->y = (rot[1][0] * ecl->x) + (rot[1][1] * ecl->y) + (rot[1][2] * ecl->z);
    equ->z = (rot[2][0] * ecl->x) + (rot[2][1] * ecl->y) + (rot[2][2] * ecl->z);
}


static astro_vector_t TopPosition(const astro_top_model_t *model, int planet, astro_time_t time)
{
    top_elliptical_t ellip;
    astro_vector_t ecl, equ;

    ellip = TopCalcElliptical(planet, model, time.tt);
    TopEcliptic(&ellip, &ecl);
    TopEquatorial(&ecl, &equ);
    equ.status = ASTRO_SUCCESS;
    equ.t = time;

    return equ;
}


/** @cond DOXYGEN_SKIP */
#define CalcPluto(time)    (TopPosition(topmodel_8, 9, (time)))
/** @endcond */

/*------------------ end of generated code ------------------*/

static void AdjustBarycenter(astro_vector_t *ssb, astro_time_t time, astro_body_t body, double pmass)
{
    astro_vector_t planet;
    double shift;

    shift = pmass / (pmass + SUN_MASS);
    planet = CalcVsop(&vsop[body], time);
    ssb->x += shift * planet.x;
    ssb->y += shift * planet.y;
    ssb->z += shift * planet.z;
}

static astro_vector_t CalcSolarSystemBarycenter(astro_time_t time)
{
    astro_vector_t ssb;

    ssb.status = ASTRO_SUCCESS;
    ssb.t = time;
    ssb.x = ssb.y = ssb.z = 0.0;

    AdjustBarycenter(&ssb, time, BODY_JUPITER, JUPITER_MASS);
    AdjustBarycenter(&ssb, time, BODY_SATURN,  SATURN_MASS);
    AdjustBarycenter(&ssb, time, BODY_URANUS,  URANUS_MASS);
    AdjustBarycenter(&ssb, time, BODY_NEPTUNE, NEPTUNE_MASS);

    return ssb;
}

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
 * If given an invalid value for `body`, this function will fail. The caller should always check
 * the `status` field inside the returned #astro_vector_t for `ASTRO_SUCCESS` (success)
 * or any other value (failure) before trusting the resulting vector.
 *
 * @param body
 *      A body for which to calculate a heliocentric position: the Sun, Moon, any of the planets,
 *      the Solar System Barycenter (SSB), or the Earth Moon Barycenter (EMB).
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

    case BODY_EMB:
        vector = Astronomy_GeoMoon(time);
        earth = CalcEarth(time);
        vector.x = earth.x + (vector.x / (1.0 + EARTH_MOON_MASS_RATIO));
        vector.y = earth.y + (vector.y / (1.0 + EARTH_MOON_MASS_RATIO));
        vector.z = earth.z + (vector.z / (1.0 + EARTH_MOON_MASS_RATIO));
        return vector;

    case BODY_SSB:
        return CalcSolarSystemBarycenter(time);

    default:
        return VecError(ASTRO_INVALID_BODY, time);
    }
}

/**
 * @brief Calculates the distance from a body to the Sun at a given time.
 *
 * Given a date and time, this function calculates the distance between
 * the center of `body` and the center of the Sun.
 * For the planets Mercury through Neptune, this function is significantly
 * more efficient than calling #Astronomy_HelioVector followed by #Astronomy_VectorLength.
 *
 * @param body
 *      A body for which to calculate a heliocentric distance: the Sun, Moon, or any of the planets.
 *
 * @param time
 *      The date and time for which to calculate the heliocentric distance.
 *
 * @return
 *      If successful, an #astro_func_result_t structure whose `status` is `ASTRO_SUCCESS`
 *      and whose `value` holds the heliocentric distance in AU.
 *      Otherwise, `status` reports an error condition.
 */
astro_func_result_t Astronomy_HelioDistance(astro_body_t body, astro_time_t time)
{
    astro_vector_t vector;
    astro_func_result_t result;

    switch (body)
    {
    case BODY_SUN:
        result.status = ASTRO_SUCCESS;
        result.value = 0.0;
        return result;

    case BODY_MERCURY:
    case BODY_VENUS:
    case BODY_EARTH:
    case BODY_MARS:
    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
        result.status = ASTRO_SUCCESS;
        result.value = VsopHelioDistance(&vsop[body], time);
        return result;

    default:
        /* For non-VSOP objects, fall back to taking the length of the heliocentric vector. */
        vector = Astronomy_HelioVector(body, time);
        if (vector.status != ASTRO_SUCCESS)
            return FuncError(vector.status);
        result.status = ASTRO_SUCCESS;
        result.value = Astronomy_VectorLength(vector);
        return result;
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
 * If given an invalid value for `body`, this function will fail. The caller should always check
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
 * @return              A geocentric position vector of the center of the given body.
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
    astro_time_t *time,
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
    gc = Astronomy_GeoVector(body, *time, aberration);
    if (gc.status != ASTRO_SUCCESS)
        return EquError(gc.status);

    j2000[0] = gc.x - gc_observer[0];
    j2000[1] = gc.y - gc_observer[1];
    j2000[2] = gc.z - gc_observer[2];

    switch (equdate)
    {
    case EQUATOR_OF_DATE:
        precession(0.0, j2000, time->tt, temp);
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
 *      See function remarks for more details.
 *
 * @param dec
 *      The declination of the body in degrees. See function remarks for more details.
 *
 * @param refraction
 *      Selects whether to correct for atmospheric refraction, and if so, which model to use.
 *      The recommended value for most uses is `REFRACTION_NORMAL`.
 *      See function remarks for more details.
 *
 * @return
 *      The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.
 */
astro_horizon_t Astronomy_Horizon(
    astro_time_t *time, astro_observer_t observer, double ra, double dec, astro_refraction_t refraction)
{
    astro_horizon_t hor;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double p[3], pz, pn, pw, proj;
    double az, zd;
    double spin_angle;

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

    spin_angle = -15.0 * sidereal_time(time);
    spin(spin_angle, uze, uz);
    spin(spin_angle, une, un);
    spin(spin_angle, uwe, uw);

    p[0] = cosdc * cosra;
    p[1] = cosdc * sinra;
    p[2] = sindc;

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];

    proj = sqrt(pn*pn + pw*pw);
    az = 0.0;
    if (proj > 0.0)
    {
        az = -atan2(pw, pn) * RAD2DEG;
        if (az < 0)
            az += 360;
        else if (az >= 360)
            az -= 360;
    }
    zd = atan2(proj, pz) * RAD2DEG;
    hor.ra = ra;
    hor.dec = dec;

    if (refraction == REFRACTION_NORMAL || refraction == REFRACTION_JPLHOR)
    {
        double zd0, refr;

        zd0 = zd;
        refr = Astronomy_Refraction(refraction, 90.0 - zd);
        zd -= refr;

        if (refr > 0.0 && zd > 3.0e-4)
        {
            int j;
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
                hor.ra = atan2(pr[1], pr[0]) * (RAD2DEG / 15.0);
                if (hor.ra < 0.0)
                    hor.ra += 24.0;
                else if (hor.ra >= 24.0)
                    hor.ra -= 24.0;
            }
            else
            {
                hor.ra = 0.0;
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
    nutation(&adjusted_time, 0, stemp, sun_ofdate);

    /* Convert equatorial coordinates to ecliptic coordinates. */
    true_obliq = DEG2RAD * e_tilt(&adjusted_time).tobl;
    return RotateEquatorialToEcliptic(sun_ofdate, true_obliq);
}

/**
 * @brief Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.
 *
 * Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
 * on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
 * which are relative to the plane of the Earth's orbit around the Sun.
 *
 * @param equ
 *      Equatorial coordinates in the J2000 frame of reference.
 *      You can call #Astronomy_GeoVector to obtain suitable equatorial coordinates.
 *
 * @return
 *      Ecliptic coordinates in the J2000 frame of reference.
 */
astro_ecliptic_t Astronomy_Ecliptic(astro_vector_t equ)
{
    /* Based on NOVAS functions equ2ecl() and equ2ecl_vec(). */
    static const double ob2000 = 0.40909260059599012;   /* mean obliquity of the J2000 ecliptic in radians */
    double pos[3];

    if (equ.status != ASTRO_SUCCESS)
        return EclError(equ.status);

    pos[0] = equ.x;
    pos[1] = equ.y;
    pos[2] = equ.z;

    return RotateEquatorialToEcliptic(pos, ob2000);
}

/**
 * @brief   Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.
 *
 * This function calculates the angle around the plane of the Earth's orbit
 * of a celestial body, as seen from the center of the Sun.
 * The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
 * in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).
 *
 * @param body
 *      A body other than the Sun.
 *
 * @param time
 *      The date and time at which the body's ecliptic longitude is to be calculated.
 *
 * @return
 *      On success, returns a structure whose `status` is `ASTRO_SUCCESS` and whose
 *      `angle` holds the ecliptic longitude in degrees.
 *      On failure, `status` holds a value other than `ASTRO_SUCCESS`.
 */
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

/**
 * @brief
 *      Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.
 *
 * This function finds the moment in time, if any exists in the given time window,
 * that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.
 *
 * This function can be used to determine equinoxes and solstices.
 * However, it is usually more convenient and efficient to call #Astronomy_Seasons
 * to calculate all equinoxes and solstices for a given calendar year.
 *
 * The function searches the window of time specified by `startTime` and `startTime+limitDays`.
 * The search will return an error if the Sun never reaches the longitude `targetLon` or
 * if the window is so large that the longitude ranges more than 180 degrees within it.
 * It is recommended to keep the window smaller than 10 days when possible.
 *
 * @param targetLon
 *      The desired ecliptic longitude in degrees, relative to the true equinox of date.
 *      This may be any value in the range [0, 360), although certain values have
 *      conventional meanings:
 *      0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice.
 *
 * @param startTime
 *      The date and time for starting the search for the desired longitude event.
 *
 * @param limitDays
 *      The real-valued number of days, which when added to `startTime`, limits the
 *      range of time over which the search looks.
 *      It is recommended to keep this value between 1 and 10 days.
 *      See function remarks for more details.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the `time` field will contain the date and time the Sun reaches the target longitude.
 *      Any other value indicates an error.
 *      See remarks in #Astronomy_Search (which this function calls) for more information about possible error codes.
 */
astro_search_result_t Astronomy_SearchSunLongitude(
    double targetLon,
    astro_time_t startTime,
    double limitDays)
{
    astro_time_t t2 = Astronomy_AddDays(startTime, limitDays);
    return Astronomy_Search(sun_offset, &targetLon, startTime, t2, 1.0);
}

/** @cond DOXYGEN_SKIP */
#define CALLFUNC(f,t)  \
    do { \
        funcres = func(context, (t)); \
        if (funcres.status != ASTRO_SUCCESS) return SearchError(funcres.status); \
        (f) = funcres.value; \
    } while(0)
/** @endcond */

/**
 * @brief Searches for a time at which a function's value increases through zero.
 *
 * Certain astronomy calculations involve finding a time when an event occurs.
 * Often such events can be defined as the root of a function:
 * the time at which the function's value becomes zero.
 *
 * `Astronomy_Search` finds the *ascending root* of a function: the time at which
 * the function's value becomes zero while having a positive slope. That is, as time increases,
 * the function transitions from a negative value, through zero at a specific moment,
 * to a positive value later. The goal of the search is to find that specific moment.
 *
 * The search function is specified by two parameters: `func` and `context`.
 * The `func` parameter is a pointer to the function itself, which accepts a time
 * and a context containing any other arguments needed to evaluate the function.
 * The `context` parameter supplies that context for the given search.
 * As an example, a caller may wish to find the moment a celestial body reaches a certain
 * ecliptic longitude. In that case, the caller might create a structure that contains
 * an #astro_body_t member to specify the body and a `double` to hold the target longitude.
 * The function would cast the pointer `context` passed in as a pointer to that structure type.
 * It could subtract the target longitude from the actual longitude at a given time;
 * thus the difference would equal zero at the moment in time the planet reaches the
 * desired longitude.
 *
 * The `func` returns an #astro_func_result_t structure every time it is called.
 * If the returned structure has a value of `status` other than `ASTRO_SUCCESS`,
 * the search immediately fails and reports that same error code in the `status`
 * returned by `Astronomy_Search`. Otherwise, `status` is `ASTRO_SUCCESS` and
 * `value` is the value of the function, and the search proceeds until it either
 * finds the ascending root or fails for some reason.
 *
 * The search calls `func` repeatedly to rapidly narrow in on any ascending
 * root within the time window specified by `t1` and `t2`. The search never
 * reports a solution outside this time window.
 *
 * `Astronomy_Search` uses a combination of bisection and quadratic interpolation
 * to minimize the number of function calls. However, it is critical that the
 * supplied time window be small enough that there cannot be more than one root
 * (ascedning or descending) within it; otherwise the search can fail.
 * Beyond that, it helps to make the time window as small as possible, ideally
 * such that the function itself resembles a smooth parabolic curve within that window.
 *
 * If an ascending root is not found, or more than one root
 * (ascending and/or descending) exists within the window `t1`..`t2`,
 * the search will fail with status code `ASTRO_SEARCH_FAILURE`.
 *
 * If the search does not converge within 20 iterations, it will fail
 * with status code `ASTRO_NO_CONVERGE`.
 *
 * @param func
 *      The function for which to find the time of an ascending root.
 *      See function remarks for more details.
 *
 * @param context
 *      Any ancillary data needed by the function `func` to calculate a value.
 *      The data type varies depending on the function passed in.
 *      For example, the function may involve a specific celestial body that
 *      must be specified somehow.
 *
 * @param t1
 *      The lower time bound of the search window.
 *      See function remarks for more details.
 *
 * @param t2
 *      The upper time bound of the search window.
 *      See function remarks for more details.
 *
 * @param dt_tolerance_seconds
 *      Specifies an amount of time in seconds within which a bounded ascending root
 *      is considered accurate enough to stop. A typical value is 1 second.
 *
 * @return
 *      If successful, the returned structure has `status` equal to `ASTRO_SUCCESS`
 *      and `time` set to a value within `dt_tolerance_seconds` of an ascending root.
 *      On success, the `time` value will always be in the inclusive range [`t1`, `t2`].
 *      If the search fails, `status` will be set to a value other than `ASTRO_SUCCESS`.
 *      See function remarks for more details.
 */
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
    double f1, f2, fmid=0.0, fq, dt_days, dt, dt_guess;
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
            tq = Astronomy_TimeFromDays(q_ut);
            CALLFUNC(fq, tq);
            if (q_df_dt != 0.0)
            {
                dt_guess = fabs(fq / q_df_dt);
                if (dt_guess < dt_days)
                {
                    /* The estimated time error is small enough that we can quit now. */
                    result.time = tq;
                    result.status = ASTRO_SUCCESS;
                    return result;
                }

                /* Try guessing a tighter boundary with the interpolated root at the center. */
                dt_guess *= 1.2;
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
    astro_time_t startTime = Astronomy_MakeTime(year, month, day, 0, 0, 0.0);
    astro_search_result_t result = Astronomy_SearchSunLongitude(targetLon, startTime, 4.0);
    *time = result.time;
    return result.status;
}

/**
 * @brief Finds both equinoxes and both solstices for a given calendar year.
 *
 * The changes of seasons are defined by solstices and equinoxes.
 * Given a calendar year number, this function calculates the
 * March and September equinoxes and the June and December solstices.
 *
 * The equinoxes are the moments twice each year when the plane of the
 * Earth's equator passes through the center of the Sun. In other words,
 * the Sun's declination is zero at both equinoxes.
 * The March equinox defines the beginning of spring in the northern hemisphere
 * and the beginning of autumn in the southern hemisphere.
 * The September equinox defines the beginning of autumn in the northern hemisphere
 * and the beginning of spring in the southern hemisphere.
 *
 * The solstices are the moments twice each year when one of the Earth's poles
 * is most tilted toward the Sun. More precisely, the Sun's declination reaches
 * its minimum value at the December solstice, which defines the beginning of
 * winter in the northern hemisphere and the beginning of summer in the southern
 * hemisphere. The Sun's declination reaches its maximum value at the June solstice,
 * which defines the beginning of summer in the northern hemisphere and the beginning
 * of winter in the southern hemisphere.
 *
 * @param year
 *      The calendar year number for which to calculate equinoxes and solstices.
 *      The value may be any integer, but only the years 1800 through 2100 have been
 *      validated for accuracy: unit testing against data from the
 *      United States Naval Observatory confirms that all equinoxes and solstices
 *      for that range of years are within 2 minutes of the correct time.
 *
 * @return
 *      The times of the four seasonal changes in the given calendar year.
 *      This function should always succeed. However, to be safe, callers
 *      should check the `status` field of the returned structure to make sure
 *      it contains `ASTRO_SUCCESS`. Any failures indicate a bug in the algorithm
 *      and should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues).
 */
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

/**
 * @brief   Returns the angle between the given body and the Sun, as seen from the Earth.
 *
 * This function calculates the angular separation between the given body and the Sun,
 * as seen from the center of the Earth. This angle is helpful for determining how
 * easy it is to see the body away from the glare of the Sun.
 *
 * @param body
 *      The celestial body whose angle from the Sun is to be measured.
 *      Not allowed to be `BODY_EARTH`.
 *
 * @param time
 *      The time at which the observation is made.
 *
 * @return
 *      If successful, the returned structure contains `ASTRO_SUCCESS` in the `status` field
 *      and `angle` holds the angle in degrees between the Sun and the specified body as
 *      seen from the center of the Earth.
 *      If an error occurs, the `status` field contains a value other than `ASTRO_SUCCESS`
 *      that indicates the error condition.
 */
astro_angle_result_t Astronomy_AngleFromSun(astro_body_t body, astro_time_t time)
{
    astro_vector_t sv, bv;

    if (body == BODY_EARTH)
        return AngleError(ASTRO_EARTH_NOT_ALLOWED);

    sv = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (sv.status != ASTRO_SUCCESS)
        return AngleError(sv.status);

    bv = Astronomy_GeoVector(body, time, ABERRATION);
    if (bv.status != ASTRO_SUCCESS)
        return AngleError(bv.status);

    return AngleBetween(sv, bv);
}

/**
 * @brief
 *      Determines visibility of a celestial body relative to the Sun, as seen from the Earth.
 *
 * This function returns an #astro_elongation_t structure, which provides the following
 * information about the given celestial body at the given time:
 *
 * - `visibility` is an enumerated type that specifies whether the body is more easily seen
 *    in the morning before sunrise, or in the evening after sunset.
 *
 * - `elongation` is the angle in degrees between two vectors: one from the center of the Earth to the
 *    center of the Sun, the other from the center of the Earth to the center of the specified body.
 *    This angle indicates how far away the body is from the glare of the Sun.
 *    The elongation angle is always in the range [0, 180].
 *
 * - `ecliptic_separation` is the absolute value of the difference between the body's ecliptic longitude
 *   and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures
 *   around the plane of the Earth's orbit, and ignores how far above or below that plane the body is.
 *   The ecliptic separation is measured in degrees and is always in the range [0, 180].
 *
 * @param body
 *      The celestial body whose visibility is to be calculated.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @return
 *      If successful, the `status` field in the returned structure contains `ASTRO_SUCCESS`
 *      and all the other fields in the structure are valid. On failure, `status` contains
 *      some other value as an error code and the other fields contain invalid values.
 */
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
        result.ecliptic_separation = 360.0 - angres.angle;
    }
    else
    {
        result.visibility = VISIBLE_EVENING;
        result.ecliptic_separation = angres.angle;
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

/**
 * @brief
 *      Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.
 *
 * Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
 * Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
 * The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
 * a telescope without atmospheric interference, are when these planets reach maximum elongation.
 * These are events where the planets reach the maximum angle from the Sun as seen from the Earth.
 *
 * This function solves for those times, reporting the next maximum elongation event's date and time,
 * the elongation value itself, the relative longitude with the Sun, and whether the planet is best
 * observed in the morning or evening. See #Astronomy_Elongation for more details about the returned structure.
 *
 * @param body
 *      Either `BODY_MERCURY` or `BODY_VENUS`. Any other value will fail with the error `ASTRO_INVALID_BODY`.
 *      To find the best viewing opportunites for planets farther from the Sun than the Earth is (Mars through Pluto)
 *      use #Astronomy_SearchRelativeLongitude to find the next opposition event.
 *
 * @param startTime
 *      The date and time at which to begin the search. The maximum elongation event found will always
 *      be the first one that occurs after this date and time.
 *
 * @return
 *      If successful, the `status` field of the returned structure will be `ASTRO_SUCCESS`
 *      and the other structure fields will be valid. Otherwise, `status` will contain
 *      some other value indicating an error.
 */
astro_elongation_t Astronomy_SearchMaxElongation(astro_body_t body, astro_time_t startTime)
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
        plon = Astronomy_EclipticLongitude(body, startTime);
        if (plon.status != ASTRO_SUCCESS)
            return ElongError(plon.status);

        elon = Astronomy_EclipticLongitude(BODY_EARTH, startTime);
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

        t_start = Astronomy_AddDays(startTime, adjust_days);

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

        if (searchx.time.tt >= startTime.tt)
            return Astronomy_Elongation(body, searchx.time);

        /* This event is in the past (earlier than startTime). */
        /* We need to search forward from t2 to find the next possible window. */
        /* We never need to search more than twice. */
        startTime = Astronomy_AddDays(t2, 1.0);
    }

    return ElongError(ASTRO_SEARCH_FAILURE);
}

/**
 * @brief
 *      Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.
 *
 * This function can be used to determine where a planet appears around the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) as seen from the Earth,
 * relative to the Sun's apparent position.
 *
 * The angle starts at 0 when the body and the Sun are at the same ecliptic longitude
 * as seen from the Earth. The angle increases in the prograde direction
 * (the direction that the planets orbit the Sun and the Moon orbits the Earth).
 *
 * When the angle is 180 degrees, it means the Sun and the body appear on opposite sides
 * of the sky for an Earthly observer. When `body` is a planet whose orbit around the
 * Sun is farther than the Earth's, 180 degrees indicates opposition. For the Moon,
 * it indicates a full moon.
 *
 * The angle keeps increasing up to 360 degrees as the body's apparent prograde
 * motion continues relative to the Sun. When the angle reaches 360 degrees, it starts
 * over at 0 degrees.
 *
 * Values between 0 and 180 degrees indicate that the body is visible in the evening sky
 * after sunset.  Values between 180 degrees and 360 degrees indicate that the body
 * is visible in the morning sky before sunrise.
 *
 * @param body
 *      The celestial body for which to find longitude from the Sun.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @return
 *      On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and
 *      the `angle` field holds a value in the range [0, 360).
 *      On failure, the `status` field contains some other value indicating an error condition.
 */
astro_angle_result_t Astronomy_LongitudeFromSun(astro_body_t body, astro_time_t time)
{
    astro_vector_t sv, bv;
    astro_ecliptic_t se, be;
    astro_angle_result_t result;

    if (body == BODY_EARTH)
        return AngleError(ASTRO_EARTH_NOT_ALLOWED);

    sv = Astronomy_GeoVector(BODY_SUN, time, NO_ABERRATION);
    se = Astronomy_Ecliptic(sv);        /* checks for errors in sv */
    if (se.status != ASTRO_SUCCESS)
        return AngleError(se.status);

    bv = Astronomy_GeoVector(body, time, NO_ABERRATION);
    be = Astronomy_Ecliptic(bv);        /* checks for errors in bv */
    if (be.status != ASTRO_SUCCESS)
        return AngleError(be.status);

    result.status = ASTRO_SUCCESS;
    result.angle = NormalizeLongitude(be.elon - se.elon);
    return result;
}

/**
 * @brief
 *      Returns the Moon's phase as an angle from 0 to 360 degrees.
 *
 * This function determines the phase of the Moon using its apparent
 * ecliptic longitude relative to the Sun, as seen from the center of the Earth.
 * Certain values of the angle have conventional definitions:
 *
 * - 0 = new moon
 * - 90 = first quarter
 * - 180 = full moon
 * - 270 = third quarter
 *
 * @param time
 *      The date and time of the observation.
 *
 * @return
 *      On success, the function returns the angle as described in the function remarks
 *      in the `angle` field and `ASTRO_SUCCESS` in the `status` field.
 *      The function should always succeed, but it is a good idea for callers to check
 *      the `status` field in the returned structure.
 *      Any other value in `status` indicates a failure that should be
 *      [reported as an issue](https://github.com/cosinekitty/astronomy/issues).
 */
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

/**
 * @brief
 *      Searches for the time that the Moon reaches a specified phase.
 *
 * Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
 * longitude with respect to the Sun's geocentric ecliptic longitude.
 * When the Moon and the Sun have the same longitude, that is defined as a new moon.
 * When their longitudes are 180 degrees apart, that is defined as a full moon.
 *
 * This function searches for any value of the lunar phase expressed as an
 * angle in degrees in the range [0, 360).
 *
 * If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
 * it is much easier to call the functions #Astronomy_SearchMoonQuarter and #Astronomy_NextMoonQuarter.
 * This function is useful for finding general phase angles outside those four quarters.
 *
 * @param targetLon
 *      The difference in geocentric longitude between the Sun and Moon
 *      that specifies the lunar phase being sought. This can be any value
 *      in the range [0, 360).  Certain values have conventional names:
 *      0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter.
 *
 * @param startTime
 *      The beginning of the time window in which to search for the Moon reaching the specified phase.
 *
 * @param limitDays
 *      The number of days after `startTime` that limits the time window for the search.
 *
 * @return
 *      On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and
 *      the `time` field holds the date and time when the Moon reaches the target longitude.
 *      On failure, `status` holds some other value as an error code.
 *      One possible error code is `ASTRO_NO_MOON_QUARTER` if `startTime` and `limitDays`
 *      do not enclose the desired event. See remarks in #Astronomy_Search for other possible
 *      error codes.
 */
astro_search_result_t Astronomy_SearchMoonPhase(double targetLon, astro_time_t startTime, double limitDays)
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
        Return ASTRO_NO_MOON_QUARTER if the final result goes beyond limitDays after startTime.
    */
    const double uncertainty = 0.9;
    astro_func_result_t funcres;
    double ya, est_dt, dt1, dt2;
    astro_time_t t1, t2;

    funcres = moon_offset(&targetLon, startTime);
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
    t1 = Astronomy_AddDays(startTime, dt1);
    t2 = Astronomy_AddDays(startTime, dt2);
    return Astronomy_Search(moon_offset, &targetLon, t1, t2, 1.0);
}

/**
 * @brief
 *      Finds the first lunar quarter after the specified date and time.
 *
 * A lunar quarter is one of the following four lunar phase events:
 * new moon, first quarter, full moon, third quarter.
 * This function finds the lunar quarter that happens soonest
 * after the specified date and time.
 *
 * To continue iterating through consecutive lunar quarters, call this function once,
 * followed by calls to #Astronomy_NextMoonQuarter as many times as desired.
 *
 * @param startTime
 *      The date and time at which to start the search.
 *
 * @return
 *      This function should always succeed, indicated by the `status` field
 *      in the returned structure holding `ASTRO_SUCCESS`. Any other value indicates
 *      an internal error, which should be [reported as an issue](https://github.com/cosinekitty/astronomy/issues).
 *      To be safe, calling code should always check the `status` field for errors.
 */
astro_moon_quarter_t Astronomy_SearchMoonQuarter(astro_time_t startTime)
{
    astro_moon_quarter_t mq;
    astro_angle_result_t angres;
    astro_search_result_t srchres;

    /* Determine what the next quarter phase will be. */
    angres = Astronomy_MoonPhase(startTime);
    if (angres.status != ASTRO_SUCCESS)
        return MoonQuarterError(angres.status);

    mq.quarter = (1 + (int)floor(angres.angle / 90.0)) % 4;
    srchres = Astronomy_SearchMoonPhase(90.0 * mq.quarter, startTime, 10.0);
    if (srchres.status != ASTRO_SUCCESS)
        return MoonQuarterError(srchres.status);

    mq.status = ASTRO_SUCCESS;
    mq.time = srchres.time;
    return mq;
}

/**
 * @brief
 *      Continues searching for lunar quarters from a previous search.
 *
 * After calling #Astronomy_SearchMoonQuarter, this function can be called
 * one or more times to continue finding consecutive lunar quarters.
 * This function finds the next consecutive moon quarter event after the one passed in as the parameter `mq`.
 *
 * @param mq
 *      A value returned by a prior call to #Astronomy_SearchMoonQuarter or #Astronomy_NextMoonQuarter.
 *
 * @return
 *      If `mq` is valid, this function should always succeed, indicated by the `status` field
 *      in the returned structure holding `ASTRO_SUCCESS`. Any other value indicates
 *      an internal error, which (after confirming that `mq` is valid) should be
 *      [reported as an issue](https://github.com/cosinekitty/astronomy/issues).
 *      To be safe, calling code should always check the `status` field for errors.
 */
astro_moon_quarter_t Astronomy_NextMoonQuarter(astro_moon_quarter_t mq)
{
    astro_time_t time;
    astro_moon_quarter_t next_mq;

    if (mq.status != ASTRO_SUCCESS)
        return MoonQuarterError(ASTRO_INVALID_PARAMETER);

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

/**
 * @brief
 *      Searches for the time when the Earth and another planet are separated by a specified angle
 *      in ecliptic longitude, as seen from the Sun.
 *
 * A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit
 * (the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored.
 * If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around
 * that plane from one body to the other in the direction the planets orbit the Sun, you will get an
 * angle somewhere between 0 and 360 degrees. This is the relative longitude.
 *
 * Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
 * this function searches for the next time that the relative longitude measured from the planet
 * to the Earth is `targetRelLon`.
 *
 * Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:
 *
 * - When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun.
 *   For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*,
 *   a time when the other planet becomes very difficult to see because of being lost in the Sun's glare.
 *   (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing
 *   between the Earth and the Sun.)
 *
 * - When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
 *   this is known as *opposition*.  Opposition is when the planet is closest to the Earth, and
 *   also when it is visible for most of the night, so it is considered the best time to observe the planet.
 *
 * - When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun
 *   from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is
 *   very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.
 *
 * @param body
 *      A planet other than the Earth. If `body` is not a planet other than the Earth, an error occurs.
 *
 * @param targetRelLon
 *      The desired relative longitude, expressed in degrees. Must be in the range [0, 360).
 *
 * @param startTime
 *      The date and time at which to begin the search.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and `time` will hold the date and time of the relative longitude event.
 *      Otherwise `status` will hold some other value that indicates an error condition.
 */
astro_search_result_t Astronomy_SearchRelativeLongitude(astro_body_t body, double targetRelLon, astro_time_t startTime)
{
    astro_search_result_t result;
    astro_func_result_t syn;
    astro_func_result_t error_angle;
    double prev_angle;
    astro_time_t time;
    int iter, direction;

    if (body == BODY_EARTH)
        return SearchError(ASTRO_EARTH_NOT_ALLOWED);

    if (body == BODY_MOON || body == BODY_SUN)
        return SearchError(ASTRO_INVALID_BODY);

    syn = SynodicPeriod(body);
    if (syn.status != ASTRO_SUCCESS)
        return SearchError(syn.status);

    direction = IsSuperiorPlanet(body) ? +1 : -1;

    /* Iterate until we converge on the desired event. */
    /* Calculate the error angle, which will be a negative number of degrees, */
    /* meaning we are "behind" the target relative longitude. */

    error_angle = rlon_offset(body, startTime, direction, targetRelLon);
    if (error_angle.status != ASTRO_SUCCESS)
        return SearchError(error_angle.status);

    if (error_angle.value > 0)
        error_angle.value -= 360;    /* force searching forward in time */

    time = startTime;
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

/**
 * @brief
 *      Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.
 *
 * The *hour angle* of a celestial body indicates its position in the sky with respect
 * to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
 * The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
 * The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
 * to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
 * the number of hours that have passed since the most recent time that the body has culminated,
 * or reached its highest point.
 *
 * This function searches for the next time a celestial body reaches the given hour angle
 * after the date and time specified by `startTime`.
 * To find when a body culminates, pass 0 for `hourAngle`.
 * To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.
 *
 * Note that, especially close to the Earth's poles, a body as seen on a given day
 * may always be above the horizon or always below the horizon, so the caller cannot
 * assume that a culminating object is visible nor that an object is below the horizon
 * at its minimum altitude.
 *
 * On success, the function reports the date and time, along with the horizontal coordinates
 * of the body at that time, as seen by the given observer.
 *
 * @param body
 *      The celestial body, which can the Sun, the Moon, or any planet other than the Earth.
 *
 * @param observer
 *      Indicates a location on or near the surface of the Earth where the observer is located.
 *      Call #Astronomy_MakeObserver to create an observer structure.
 *
 * @param hourAngle
 *      An hour angle value in the range [0, 24) indicating the number of sidereal hours after the
 *      body's most recent culmination.
 *
 * @param startTime
 *      The date and time at which to start the search.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and the other structure fields are valid. Otherwise, `status` holds some other value
 *      that indicates an error condition.
 */
astro_hour_angle_t Astronomy_SearchHourAngle(
    astro_body_t body,
    astro_observer_t observer,
    double hourAngle,
    astro_time_t startTime)
{
    int iter = 0;
    astro_time_t time;
    astro_equatorial_t ofdate;
    astro_hour_angle_t result;
    double delta_sidereal_hours, delta_days, gast;

    if (body < MIN_BODY || body > MAX_BODY)
        return HourAngleError(ASTRO_INVALID_BODY);

    if (body == BODY_EARTH)
        return HourAngleError(ASTRO_EARTH_NOT_ALLOWED);

    if (hourAngle < 0.0 || hourAngle >= 24.0)
        return HourAngleError(ASTRO_INVALID_PARAMETER);

    time = startTime;
    for(;;)
    {
        ++iter;

        /* Calculate Greenwich Apparent Sidereal Time (GAST) at the given time. */
        gast = sidereal_time(&time);

        /* Obtain equatorial coordinates of date for the body. */
        ofdate = Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION);
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
            result.hor = Astronomy_Horizon(&time, observer, ofdate.ra, ofdate.dec, REFRACTION_NORMAL);
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

    ofdate = Astronomy_Equator(p->body, &time, p->observer, EQUATOR_OF_DATE, ABERRATION);
    if (ofdate.status != ASTRO_SUCCESS)
        return FuncError(ofdate.status);

    /* We calculate altitude without refraction, then add fixed refraction near the horizon. */
    /* This gives us the time of rise/set without the extra work. */
    hor = Astronomy_Horizon(&time, p->observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
    result.value = p->direction * (hor.altitude + RAD2DEG*(p->body_radius_au / ofdate.dist) + REFRACTION_NEAR_HORIZON);
    result.status = ASTRO_SUCCESS;
    return result;
}

/**
 * @brief
 *      Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
 *
 * This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
 * Rise time is when the body first starts to be visible above the horizon.
 * For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
 * Set time is the moment when the body appears to vanish below the horizon.
 *
 * This function corrects for typical atmospheric refraction, which causes celestial
 * bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
 * It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).
 *
 * Note that rise or set may not occur in every 24 hour period.
 * For example, near the Earth's poles, there are long periods of time where
 * the Sun stays below the horizon, never rising.
 * Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
 * This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
 * significant amount during each rotation of the Earth.
 * Therefore callers must not assume that the function will always succeed.
 *
 * @param body
 *      The Sun, Moon, or any planet other than the Earth.
 *
 * @param observer
 *      The location where observation takes place.
 *      You can create an observer structure by calling #Astronomy_MakeObserver.
 *
 * @param direction
 *      Either `DIRECTION_RISE` to find a rise time or `DIRECTION_SET` to find a set time.
 *
 * @param startTime
 *      The date and time at which to start the search.
 *
 * @param limitDays
 *      Limits how many days to search for a rise or set time.
 *      To limit a rise or set time to the same day, you can use a value of 1 day.
 *      In cases where you want to find the next rise or set time no matter how far
 *      in the future (for example, for an observer near the south pole), you can
 *      pass in a larger value like 365.
 *
 * @return
 *      On success, the `status` field in the returned structure contains `ASTRO_SUCCESS`
 *      and the `time` field contains the date and time of the rise or set time as requested.
 *      If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the rise or set
 *      event does not occur within `limitDays` days of `startTime`. This is a normal condition,
 *      not an error. Any other value of `status` indicates an error of some kind.
 */
astro_search_result_t Astronomy_SearchRiseSet(
    astro_body_t body,
    astro_observer_t observer,
    astro_direction_t direction,
    astro_time_t startTime,
    double limitDays)
{
    context_peak_altitude_t context;
    double ha_before, ha_after;
    astro_time_t time_start, time_before;
    astro_func_result_t alt_before, alt_after;
    astro_hour_angle_t evt_before, evt_after;

    if (body == BODY_EARTH)
        return SearchError(ASTRO_EARTH_NOT_ALLOWED);

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
    case BODY_SUN:  context.body_radius_au = SUN_RADIUS_AU;                 break;
    case BODY_MOON: context.body_radius_au = MOON_EQUATORIAL_RADIUS_AU;     break;
    default:        context.body_radius_au = 0.0;                           break;
    }

    /*
        See if the body is currently above/below the horizon.
        If we are looking for next rise time and the body is below the horizon,
        we use the current time as the lower time bound and the next culmination
        as the upper bound.
        If the body is above the horizon, we search for the next bottom and use it
        as the lower bound and the next culmination after that bottom as the upper bound.
        The same logic applies for finding set times, only we swap the hour angles.
    */

    time_start = startTime;
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

        /* If we didn't find the desired event, use evt_after.time to find the next before-event. */
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

/**
 * @brief
 *      Finds visual magnitude, phase angle, and other illumination information about a celestial body.
 *
 * This function calculates information about how bright a celestial body appears from the Earth,
 * reported as visual magnitude, which is a smaller (or even negative) number for brighter objects
 * and a larger number for dimmer objects.
 *
 * For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
 * the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
 * of the body appears illuminated as seen from the Earth. For example, when the phase angle is
 * near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
 * 180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
 * of 90 degrees means the body appears "half full".
 * For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
 * so it doesn't have a phase angle.
 *
 * When the body is Saturn, the returned structure contains a field `ring_tilt` that holds
 * the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
 * the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
 * 0 for all bodies other than Saturn.
 *
 * @param body
 *      The Sun, Moon, or any planet other than the Earth.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @return
 *      On success, the `status` field of the return structure holds `ASTRO_SUCCESS`
 *      and the other structure fields are valid.
 *      Any other value indicates an error, in which case the remaining structure fields are not valid.
 */
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

/**
 * @brief
 *      Searches for the date and time Venus will next appear brightest as seen from the Earth.
 *
 * This function searches for the date and time Venus appears brightest as seen from the Earth.
 * Currently only Venus is supported for the `body` parameter, though this could change in the future.
 * Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth,
 * so peak magnitude events have little practical value for that planet.
 * Planets other than Venus and Mercury reach peak magnitude at opposition, which can
 * be found using #Astronomy_SearchRelativeLongitude.
 * The Moon reaches peak magnitude at full moon, which can be found using
 * #Astronomy_SearchMoonQuarter or #Astronomy_SearchMoonPhase.
 * The Sun reaches peak magnitude at perihelion, which occurs each year in January.
 * However, the difference is minor and has little practical value.
 *
 * @param body
 *      Currently only `BODY_VENUS` is allowed. Any other value results in the error `ASTRO_INVALID_BODY`.
 *      See function remarks for more details.
 *
 * @param startTime
 *      The date and time to start searching for the next peak magnitude event.
 *
 * @return
 *      See documentation about the return value from #Astronomy_Illumination.
 */
astro_illum_t Astronomy_SearchPeakMagnitude(astro_body_t body, astro_time_t startTime)
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
        plon = Astronomy_EclipticLongitude(body, startTime);
        if (plon.status != ASTRO_SUCCESS)
            return IllumError(plon.status);

        elon = Astronomy_EclipticLongitude(BODY_EARTH, startTime);
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
        else if (rlon >= +s2 || rlon < -s2)
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
        t_start = Astronomy_AddDays(startTime, adjust_days);
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

        if (tx.time.tt >= startTime.tt)
            return Astronomy_Illumination(body, tx.time);

        /* This event is in the past (earlier than startTime). */
        /* We need to search forward from t2 to find the next possible window. */
        /* We never need to search more than twice. */
        startTime = Astronomy_AddDays(t2.time, 1.0);
    }

    return IllumError(ASTRO_SEARCH_FAILURE);
}

static double MoonDistance(astro_time_t t)
{
    double lon, lat, dist;
    CalcMoon(t.tt / 36525.0, &lon, &lat, &dist);
    return dist;
}

static astro_func_result_t moon_distance_slope(void *context, astro_time_t time)
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

/**
 * @brief
 *      Finds the date and time of the Moon's closest distance (perigee)
 *      or farthest distance (apogee) with respect to the Earth.
 *
 * Given a date and time to start the search in `startTime`, this function finds the
 * next date and time that the center of the Moon reaches the closest or farthest point
 * in its orbit with respect to the center of the Earth, whichever comes first
 * after `startTime`.
 *
 * The closest point is called *perigee* and the farthest point is called *apogee*.
 * The word *apsis* refers to either event.
 *
 * To iterate through consecutive alternating perigee and apogee events, call `Astronomy_SearchLunarApsis`
 * once, then use the return value to call #Astronomy_NextLunarApsis. After that,
 * keep feeding the previous return value from `Astronomy_NextLunarApsis` into another
 * call of `Astronomy_NextLunarApsis` as many times as desired.
 *
 * @param startTime
 *      The date and time at which to start searching for the next perigee or apogee.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`,
 *      `time` holds the date and time of the next lunar apsis, `kind` holds either
 *      `APSIS_PERICENTER` for perigee or `APSIS_APOCENTER` for apogee, and the distance
 *      values `dist_au` (astronomical units) and `dist_km` (kilometers) are valid.
 *      If the function fails, `status` holds some value other than `ASTRO_SUCCESS` that
 *      indicates what went wrong, and the other structure fields are invalid.
 */
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
    m1 = moon_distance_slope(&positive_direction, t1);
    if (m1.status != ASTRO_SUCCESS)
        return ApsisError(m1.status);

    for (iter=0; iter * increment < 2.0 * MEAN_SYNODIC_MONTH; ++iter)
    {
        t2 = Astronomy_AddDays(t1, increment);
        m2 = moon_distance_slope(&positive_direction, t2);
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
                search = Astronomy_Search(moon_distance_slope, &positive_direction, t1, t2, 1.0);
                result.kind = APSIS_PERICENTER;
            }
            else if (m1.value > 0.0 || m2.value < 0.0)
            {
                /* We found a maximum-distance event: apogee. */
                /* Search the time range for the time when the slope goes from positive to negative. */
                search = Astronomy_Search(moon_distance_slope, &negative_direction, t1, t2, 1.0);
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

/**
 * @brief
 *      Finds the next lunar perigee or apogee event in a series.
 *
 * This function requires an #astro_apsis_t value obtained from a call
 * to #Astronomy_SearchLunarApsis or `Astronomy_NextLunarApsis`. Given
 * an apogee event, this function finds the next perigee event, and vice versa.
 *
 * See #Astronomy_SearchLunarApsis for more details.
 *
 * @param apsis
 *      An apsis event obtained from a call to #Astronomy_SearchLunarApsis or `Astronomy_NextLunarApsis`.
 *      See #Astronomy_SearchLunarApsis for more details.
 *
 * @return
 *      Same as the return value for #Astronomy_SearchLunarApsis.
 */
astro_apsis_t Astronomy_NextLunarApsis(astro_apsis_t apsis)
{
    static const double skip = 11.0;    /* number of days to skip to start looking for next apsis event */
    astro_apsis_t next;
    astro_time_t time;

    if (apsis.status != ASTRO_SUCCESS)
        return ApsisError(ASTRO_INVALID_PARAMETER);

    if (apsis.kind != APSIS_APOCENTER && apsis.kind != APSIS_PERICENTER)
        return ApsisError(ASTRO_INVALID_PARAMETER);

    time = Astronomy_AddDays(apsis.time, skip);
    next = Astronomy_SearchLunarApsis(time);
    if (next.status == ASTRO_SUCCESS)
    {
        /* Verify that we found the opposite apsis from the previous one. */
        if (next.kind + apsis.kind != 1)
            return ApsisError(ASTRO_INTERNAL_ERROR);
    }
    return next;
}


/** @cond DOXYGEN_SKIP */
typedef struct
{
    int direction;
    astro_body_t body;
}
planet_distance_context_t;
/** @endcond */


static astro_func_result_t planet_distance_slope(void *context, astro_time_t time)
{
    static const double dt = 0.001;
    const planet_distance_context_t *pc = context;
    astro_time_t t1 = Astronomy_AddDays(time, -dt/2.0);
    astro_time_t t2 = Astronomy_AddDays(time, +dt/2.0);
    astro_func_result_t dist1, dist2, result;

    dist1 = Astronomy_HelioDistance(pc->body, t1);
    if (dist1.status != ASTRO_SUCCESS)
        return dist1;

    dist2 = Astronomy_HelioDistance(pc->body, t2);
    if (dist2.status != ASTRO_SUCCESS)
        return dist2;

    result.value = pc->direction * (dist2.value - dist1.value) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}

static astro_apsis_t PlanetExtreme(
    astro_body_t body,
    astro_apsis_kind_t kind,
    astro_time_t start_time,
    double dayspan)
{
    astro_apsis_t apsis;
    const double direction = (kind == APSIS_APOCENTER) ? +1.0 : -1.0;
    const int npoints = 10;
    int i, best_i;
    double interval;
    double dist, best_dist;
    astro_time_t time;
    astro_func_result_t result;

    for(;;)
    {
        interval = dayspan / (npoints - 1);

        if (interval < 1.0 / 1440.0)    /* iterate until uncertainty is less than one minute */
        {
            apsis.status = ASTRO_SUCCESS;
            apsis.kind = kind;
            apsis.time = Astronomy_AddDays(start_time, interval / 2.0);
            result = Astronomy_HelioDistance(body, apsis.time);
            if (result.status != ASTRO_SUCCESS)
                return ApsisError(result.status);
            apsis.dist_au = result.value;
            apsis.dist_km = apsis.dist_au * KM_PER_AU;
            return apsis;
        }

        best_i = -1;
        best_dist = 0.0;
        for (i=0; i < npoints; ++i)
        {
            time = Astronomy_AddDays(start_time, i * interval);
            result = Astronomy_HelioDistance(body, time);
            if (result.status != ASTRO_SUCCESS)
                return ApsisError(result.status);
            dist = direction * result.value;
            if (i==0 || dist > best_dist)
            {
                best_i = i;
                best_dist = dist;
            }
        }

        /* Narrow in on the extreme point. */
        start_time = Astronomy_AddDays(start_time, (best_i - 1) * interval);
        dayspan = 2.0 * interval;
    }
}


static astro_apsis_t BruteSearchPlanetApsis(astro_body_t body, astro_time_t startTime)
{
    const int npoints = 100;
    int i;
    astro_time_t t1, t2, time, t_min, t_max;
    double dist, max_dist, min_dist;
    astro_apsis_t perihelion, aphelion;
    double interval;
    double period;
    astro_func_result_t result;

    /*
        Neptune is a special case for two reasons:
        1. Its orbit is nearly circular (low orbital eccentricity).
        2. It is so distant from the Sun that the orbital period is very long.
        Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
        to be so significant that there are 3 local minima in the distance-vs-time curve
        near each apsis. Therefore, unlike for other planets, we can't use an optimized
        algorithm for finding dr/dt = 0.
        Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
        heliocentric distance.

        There is a similar problem in the TOP2013 model for Pluto:
        Its position vector has high-frequency oscillations that confuse the
        slope-based determination of apsides.
    */

    /*
        Rewind approximately 30 degrees in the orbit,
        then search forward for 270 degrees.
        This is a very cautious way to prevent missing an apsis.
        Typically we will find two apsides, and we pick whichever
        apsis is ealier, but after startTime.
        Sample points around this orbital arc and find when the distance
        is greatest and smallest.
    */
    period = PlanetOrbitalPeriod(body);
    t1 = Astronomy_AddDays(startTime, period * ( -30.0 / 360.0));
    t2 = Astronomy_AddDays(startTime, period * (+270.0 / 360.0));
    t_min = t_max = t1;
    min_dist = max_dist = -1.0;     /* prevent warning about uninitialized variables */
    interval = (t2.ut - t1.ut) / (npoints - 1.0);

    for (i=0; i < npoints; ++i)
    {
        double ut = t1.ut + (i * interval);
        time = Astronomy_TimeFromDays(ut);
        result = Astronomy_HelioDistance(body, time);
        if (result.status != ASTRO_SUCCESS)
            return ApsisError(result.status);
        dist = result.value;
        if (i == 0)
        {
            max_dist = min_dist = dist;
        }
        else
        {
            if (dist > max_dist)
            {
                max_dist = dist;
                t_max = time;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                t_min = time;
            }
        }
    }

    t1 = Astronomy_AddDays(t_min, -2 * interval);
    perihelion = PlanetExtreme(body, APSIS_PERICENTER, t1, 4 * interval);

    t1 = Astronomy_AddDays(t_max, -2 * interval);
    aphelion = PlanetExtreme(body, APSIS_APOCENTER, t1, 4 * interval);

    if (perihelion.status == ASTRO_SUCCESS && perihelion.time.tt >= startTime.tt)
    {
        if (aphelion.status == ASTRO_SUCCESS && aphelion.time.tt >= startTime.tt)
        {
            /* Perihelion and aphelion are both valid. Pick the one that comes first. */
            if (aphelion.time.tt < perihelion.time.tt)
                return aphelion;
        }
        return perihelion;
    }

    if (aphelion.status == ASTRO_SUCCESS && aphelion.time.tt >= startTime.tt)
        return aphelion;

    return ApsisError(ASTRO_FAIL_APSIS);
}


/**
 * @brief
 *      Finds the date and time of a planet's perihelion (closest approach to the Sun)
 *      or aphelion (farthest distance from the Sun) after a given time.
 *
 * Given a date and time to start the search in `startTime`, this function finds the
 * next date and time that the center of the specified planet reaches the closest or farthest point
 * in its orbit with respect to the center of the Sun, whichever comes first
 * after `startTime`.
 *
 * The closest point is called *perihelion* and the farthest point is called *aphelion*.
 * The word *apsis* refers to either event.
 *
 * To iterate through consecutive alternating perihelion and aphelion events,
 * call `Astronomy_SearchPlanetApsis` once, then use the return value to call
 * #Astronomy_NextPlanetApsis. After that, keep feeding the previous return value
 * from `Astronomy_NextPlanetApsis` into another call of `Astronomy_NextPlanetApsis`
 * as many times as desired.
 *
 * @param body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `BODY_SUN` or `BODY_MOON`.
 *
 * @param startTime
 *      The date and time at which to start searching for the next perihelion or aphelion.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`,
 *      `time` holds the date and time of the next planetary apsis, `kind` holds either
 *      `APSIS_PERICENTER` for perihelion or `APSIS_APOCENTER` for aphelion, and the distance
 *      values `dist_au` (astronomical units) and `dist_km` (kilometers) are valid.
 *      If the function fails, `status` holds some value other than `ASTRO_SUCCESS` that
 *      indicates what went wrong, and the other structure fields are invalid.
 */
astro_apsis_t Astronomy_SearchPlanetApsis(astro_body_t body, astro_time_t startTime)
{
    astro_time_t t1, t2;
    astro_search_result_t search;
    astro_func_result_t m1, m2;
    planet_distance_context_t context;
    astro_apsis_t result;
    int iter;
    double orbit_period_days;
    double increment;   /* number of days to skip in each iteration */
    astro_func_result_t dist;

    if (body == BODY_NEPTUNE || body == BODY_PLUTO)
        return BruteSearchPlanetApsis(body, startTime);

    orbit_period_days = PlanetOrbitalPeriod(body);
    if (orbit_period_days == 0.0)
        return ApsisError(ASTRO_INVALID_BODY);      /* The body must be a planet. */

    increment = orbit_period_days / 6.0;

    context.body = body;

    t1 = startTime;
    context.direction = +1;
    m1 = planet_distance_slope(&context, t1);
    if (m1.status != ASTRO_SUCCESS)
        return ApsisError(m1.status);

    for (iter=0; iter * increment < 2.0 * orbit_period_days; ++iter)
    {
        t2 = Astronomy_AddDays(t1, increment);
        context.direction = +1;
        m2 = planet_distance_slope(&context, t2);
        if (m2.status != ASTRO_SUCCESS)
            return ApsisError(m2.status);

        if (m1.value * m2.value <= 0.0)
        {
            /* There is a change of slope polarity within the time range [t1, t2]. */
            /* Therefore this time range contains an apsis. */
            /* Figure out whether it is perihelion or aphelion. */

            if (m1.value < 0.0 || m2.value > 0.0)
            {
                /* We found a minimum-distance event: perihelion. */
                /* Search the time range for the time when the slope goes from negative to positive. */
                context.direction = +1;
                result.kind = APSIS_PERICENTER;
            }
            else if (m1.value > 0.0 || m2.value < 0.0)
            {
                /* We found a maximum-distance event: aphelion. */
                /* Search the time range for the time when the slope goes from positive to negative. */
                context.direction = -1;
                result.kind = APSIS_APOCENTER;
            }
            else
            {
                /* This should never happen. It should not be possible for both slopes to be zero. */
                return ApsisError(ASTRO_INTERNAL_ERROR);
            }

            search = Astronomy_Search(planet_distance_slope, &context, t1, t2, 1.0);
            if (search.status != ASTRO_SUCCESS)
                return ApsisError(search.status);

            dist = Astronomy_HelioDistance(body, search.time);
            if (dist.status != ASTRO_SUCCESS)
                return ApsisError(dist.status);

            result.status = ASTRO_SUCCESS;
            result.time = search.time;
            result.dist_au = dist.value;
            result.dist_km = dist.value * KM_PER_AU;
            return result;
        }

        /* We have not yet found a slope polarity change. Keep searching. */
        t1 = t2;
        m1 = m2;
    }

    /* It should not be possible to fail to find an apsis within 2 orbits. */
    return ApsisError(ASTRO_INTERNAL_ERROR);
}

/**
 * @brief
 *      Finds the next planetary perihelion or aphelion event in a series.
 *
 * This function requires an #astro_apsis_t value obtained from a call
 * to #Astronomy_SearchPlanetApsis or `Astronomy_NextPlanetApsis`.
 * Given an aphelion event, this function finds the next perihelion event, and vice versa.
 *
 * See #Astronomy_SearchPlanetApsis for more details.
 *
 * @param body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `BODY_SUN` or `BODY_MOON`.
 *      Must match the body passed into the call that produced the `apsis` parameter.
 *
 * @param apsis
 *      An apsis event obtained from a call to #Astronomy_SearchPlanetApsis or `Astronomy_NextPlanetApsis`.
 *
 * @return
 *      Same as the return value for #Astronomy_SearchPlanetApsis.
 */
astro_apsis_t Astronomy_NextPlanetApsis(astro_body_t body, astro_apsis_t apsis)
{
    double skip;    /* number of days to skip to start looking for next apsis event */
    astro_apsis_t next;
    astro_time_t time;

    if (apsis.status != ASTRO_SUCCESS)
        return ApsisError(ASTRO_INVALID_PARAMETER);

    if (apsis.kind != APSIS_APOCENTER && apsis.kind != APSIS_PERICENTER)
        return ApsisError(ASTRO_INVALID_PARAMETER);

    skip = 0.25 * PlanetOrbitalPeriod(body);        /* skip 1/4 of an orbit before starting search again */
    if (skip <= 0.0)
        return ApsisError(ASTRO_INVALID_BODY);      /* body must be a planet */

    time = Astronomy_AddDays(apsis.time, skip);
    next = Astronomy_SearchPlanetApsis(body, time);
    if (next.status == ASTRO_SUCCESS)
    {
        /* Verify that we found the opposite apsis from the previous one. */
        if (next.kind + apsis.kind != 1)
            return ApsisError(ASTRO_INTERNAL_ERROR);
    }
    return next;
}


/**
 * @brief Calculates the inverse of a rotation matrix.
 *
 * Given a rotation matrix that performs some coordinate transform,
 * this function returns the matrix that reverses that trasnform.
 *
 * @param rotation
 *      The rotation matrix to be inverted.
 *
 * @return
 *      A rotation matrix that performs the opposite transformation.
 */
astro_rotation_t Astronomy_InverseRotation(astro_rotation_t rotation)
{
    astro_rotation_t inverse;

    if (rotation.status != ASTRO_SUCCESS)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    inverse.status = ASTRO_SUCCESS;
    inverse.rot[0][0] = rotation.rot[0][0];
    inverse.rot[0][1] = rotation.rot[1][0];
    inverse.rot[0][2] = rotation.rot[2][0];
    inverse.rot[1][0] = rotation.rot[0][1];
    inverse.rot[1][1] = rotation.rot[1][1];
    inverse.rot[1][2] = rotation.rot[2][1];
    inverse.rot[2][0] = rotation.rot[0][2];
    inverse.rot[2][1] = rotation.rot[1][2];
    inverse.rot[2][2] = rotation.rot[2][2];

    return inverse;
}

/**
 * @brief Creates a rotation based on applying one rotation followed by another.
 *
 * Given two rotation matrices, returns a combined rotation matrix that is
 * equivalent to rotating based on the first matrix, followed by the second.
 *
 * @param a
 *      The first rotation to apply.
 *
 * @param b
 *      The second rotation to apply.
 *
 * @return
 *      The combined rotation matrix.
 */
astro_rotation_t Astronomy_CombineRotation(astro_rotation_t a, astro_rotation_t b)
{
    astro_rotation_t c;

    if (a.status != ASTRO_SUCCESS || b.status != ASTRO_SUCCESS)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    /*
        Use matrix multiplication: c = b*a.
        We put 'b' on the left and 'a' on the right because,
        just like when you use a matrix M to rotate a vector V,
        you put the M on the left in the product M*V.
        We can think of this as 'b' rotating all the 3 column vectors in 'a'.
    */
    c.rot[0][0] = b.rot[0][0]*a.rot[0][0] + b.rot[1][0]*a.rot[0][1] + b.rot[2][0]*a.rot[0][2];
    c.rot[1][0] = b.rot[0][0]*a.rot[1][0] + b.rot[1][0]*a.rot[1][1] + b.rot[2][0]*a.rot[1][2];
    c.rot[2][0] = b.rot[0][0]*a.rot[2][0] + b.rot[1][0]*a.rot[2][1] + b.rot[2][0]*a.rot[2][2];
    c.rot[0][1] = b.rot[0][1]*a.rot[0][0] + b.rot[1][1]*a.rot[0][1] + b.rot[2][1]*a.rot[0][2];
    c.rot[1][1] = b.rot[0][1]*a.rot[1][0] + b.rot[1][1]*a.rot[1][1] + b.rot[2][1]*a.rot[1][2];
    c.rot[2][1] = b.rot[0][1]*a.rot[2][0] + b.rot[1][1]*a.rot[2][1] + b.rot[2][1]*a.rot[2][2];
    c.rot[0][2] = b.rot[0][2]*a.rot[0][0] + b.rot[1][2]*a.rot[0][1] + b.rot[2][2]*a.rot[0][2];
    c.rot[1][2] = b.rot[0][2]*a.rot[1][0] + b.rot[1][2]*a.rot[1][1] + b.rot[2][2]*a.rot[1][2];
    c.rot[2][2] = b.rot[0][2]*a.rot[2][0] + b.rot[1][2]*a.rot[2][1] + b.rot[2][2]*a.rot[2][2];

    c.status = ASTRO_SUCCESS;
    return c;
}

/**
 * @brief Converts spherical coordinates to Cartesian coordinates.
 *
 * Given spherical coordinates and a time at which they are valid,
 * returns a vector of Cartesian coordinates. The returned value
 * includes the time, as required by the type #astro_vector_t.
 *
 * @param sphere
 *      Spherical coordinates to be converted.
 *
 * @param time
 *      The time that should be included in the return value.
 *
 * @return
 *      The vector form of the supplied spherical coordinates.
 */
astro_vector_t Astronomy_VectorFromSphere(astro_spherical_t sphere, astro_time_t time)
{
    astro_vector_t vector;
    double radlat, radlon, rcoslat;

    if (sphere.status != ASTRO_SUCCESS)
        return VecError(ASTRO_INVALID_PARAMETER, time);

    radlat = sphere.lat * DEG2RAD;
    radlon = sphere.lon * DEG2RAD;
    rcoslat = sphere.dist * cos(radlat);

    vector.status = ASTRO_SUCCESS;
    vector.t = time;
    vector.x = rcoslat * cos(radlon);
    vector.y = rcoslat * sin(radlon);
    vector.z = sphere.dist * sin(radlat);

    return vector;
}


/**
 * @brief Converts Cartesian coordinates to spherical coordinates.
 *
 * Given a Cartesian vector, returns latitude, longitude, and distance.
 *
 * @param vector
 *      Cartesian vector to be converted to spherical coordinates.
 *
 * @return
 *      Spherical coordinates that are equivalent to the given vector.
 */
astro_spherical_t Astronomy_SphereFromVector(astro_vector_t vector)
{
    double xyproj;
    astro_spherical_t sphere;

    xyproj = vector.x*vector.x + vector.y*vector.y;
    sphere.dist = sqrt(xyproj + vector.z*vector.z);
    if (xyproj == 0.0)
    {
        if (vector.z == 0.0)
        {
            /* Indeterminate coordinates; pos vector has zero length. */
            return SphereError(ASTRO_INVALID_PARAMETER);
        }

        sphere.lon = 0.0;
        sphere.lat = (vector.z < 0.0) ? -90.0 : +90.0;
    }
    else
    {
        sphere.lon = RAD2DEG * atan2(vector.y, vector.x);
        if (sphere.lon < 0.0)
            sphere.lon += 360.0;

        sphere.lat = RAD2DEG * atan2(vector.z, sqrt(xyproj));
    }

    sphere.status = ASTRO_SUCCESS;
    return sphere;
}


/**
 * @brief
 *      Given angular equatorial coordinates in `equ`, calculates equatorial vector.
 *
 * @param equ
 *      Angular equatorial coordinates to be converted to a vector.
 *
 * @param time
 *      The date and time of the observation. This is needed because the returned
 *      vector requires a valid time value when passed to certain other functions.
 *
 * @return
 *      A vector in the equatorial system.
 */
astro_vector_t Astronomy_VectorFromEquator(astro_equatorial_t equ, astro_time_t time)
{
    astro_spherical_t sphere;

    if (equ.status != ASTRO_SUCCESS)
        return VecError(ASTRO_INVALID_PARAMETER, time);

    sphere.status = ASTRO_SUCCESS;
    sphere.lat = equ.dec;
    sphere.lon = 15.0 * equ.ra;     /* convert sidereal hours to degrees */
    sphere.dist = equ.dist;

    return Astronomy_VectorFromSphere(sphere, time);
}


/**
 * @brief
 *      Given an equatorial vector, calculates equatorial angular coordinates.
 *
 * @param vector
 *      A vector in an equatorial coordinate system.
 *
 * @return
 *      Angular coordinates expressed in the same equatorial system as `vector`.
 */
astro_equatorial_t Astronomy_EquatorFromVector(astro_vector_t vector)
{
    astro_equatorial_t equ;
    astro_spherical_t sphere;

    sphere = Astronomy_SphereFromVector(vector);
    if (sphere.status != ASTRO_SUCCESS)
        return EquError(sphere.status);

    equ.status = ASTRO_SUCCESS;
    equ.dec = sphere.lat;
    equ.ra = sphere.lon / 15.0;     /* convert degrees to sidereal hours */
    equ.dist = sphere.dist;

    return equ;
}


static double ToggleAzimuthDirection(double az)
{
    az = 360.0 - az;
    if (az >= 360.0)
        az -= 360.0;
    else if (az < 0.0)
        az += 360.0;
    return az;
}

/**
 * @brief Converts Cartesian coordinates to horizontal coordinates.
 *
 * Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
 *
 * *IMPORTANT:* This function differs from #Astronomy_SphereFromVector in two ways:
 * - `Astronomy_SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
 *   from north (e.g., west = +90), but this function represents a clockwise rotation
 *   (e.g., east = +90). The difference is because `Astronomy_SphereFromVector` is intended
 *   to preserve the vector "right-hand rule", while this function defines azimuth in a more
 *   traditional way as used in navigation and cartography.
 * - This function optionally corrects for atmospheric refraction, while `Astronomy_SphereFromVector`
 *   does not.
 *
 * The returned structure contains the azimuth in `lon`.
 * It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
 *
 * The altitude is stored in `lat`.
 *
 * The distance to the observed object is stored in `dist`,
 * and is expressed in astronomical units (AU).
 *
 * @param vector
 *      Cartesian vector to be converted to horizontal coordinates.
 *
 * @param refraction
 *      `REFRACTION_NORMAL`: correct altitude for atmospheric refraction (recommended).
 *      `REFRACTION_NONE`: no atmospheric refraction correction is performed.
 *      `REFRACTION_JPLHOR`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *
 * @return
 *      If successful, `status` holds `ASTRO_SUCCESS` and the other fields are valid as described
 *      in the function remarks.
 *      Otherwise `status` holds an error code and the other fields are undefined.
 */
astro_spherical_t Astronomy_HorizonFromVector(astro_vector_t vector, astro_refraction_t refraction)
{
    astro_spherical_t sphere;

    sphere = Astronomy_SphereFromVector(vector);
    if (sphere.status == ASTRO_SUCCESS)
    {
        /* Convert azimuth from counterclockwise-from-north to clockwise-from-north. */
        sphere.lon = ToggleAzimuthDirection(sphere.lon);
        sphere.lat += Astronomy_Refraction(refraction, sphere.lat);
    }

    return sphere;
}


/**
 * @brief
 *      Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.
 *
 * @param sphere
 *      A structure that contains apparent horizontal coordinates:
 *      `lat` holds the refracted azimuth angle,
 *      `lon` holds the azimuth in degrees clockwise from north,
 *      and `dist` holds the distance from the observer to the object in AU.
 *
 * @param time
 *      The date and time of the observation. This is needed because the returned
 *      #astro_vector_t structure requires a valid time value when passed to certain other functions.
 *
 * @param refraction
 *      The refraction option used to model atmospheric lensing. See #Astronomy_Refraction.
 *      This specifies how refraction is to be removed from the altitude stored in `sphere.lat`.
 *
 * @return
 *      A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
 */
astro_vector_t Astronomy_VectorFromHorizon(astro_spherical_t sphere, astro_time_t time, astro_refraction_t refraction)
{
    if (sphere.status != ASTRO_SUCCESS)
        return VecError(ASTRO_INVALID_PARAMETER, time);

    /* Convert azimuth from clockwise-from-north to counterclockwise-from-north. */
    sphere.lon = ToggleAzimuthDirection(sphere.lon);

    /* Reverse any applied refraction. */
    sphere.lat += Astronomy_InverseRefraction(refraction, sphere.lat);

    return Astronomy_VectorFromSphere(sphere, time);
}


/**
 * @brief
 *      Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
 *
 * Given an altitude angle and a refraction option, calculates
 * the amount of "lift" caused by atmospheric refraction.
 * This is the number of degrees higher in the sky an object appears
 * due to the lensing of the Earth's atmosphere.
 *
 * @param refraction
 *      The option selecting which refraction correction to use.
 *      If `REFRACTION_NORMAL`, uses a well-behaved refraction model that works well for
 *      all valid values (-90 to +90) of `altitude`.
 *      If `REFRACTION_JPLHOR`, this function returns a compatible value with the JPL Horizons tool.
 *      If any other value (including `REFRACTION_NONE`), this function returns 0.
 *
 * @param altitude
 *      An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90.
 *
 * @return
 *      The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.
 */
double Astronomy_Refraction(astro_refraction_t refraction, double altitude)
{
    double refr, hd;

    if (altitude < -90.0 || altitude > +90.0)
        return 0.0;     /* no attempt to correct an invalid altitude */

    if (refraction == REFRACTION_NORMAL || refraction == REFRACTION_JPLHOR)
    {
        // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        // JPL Horizons says it uses refraction algorithm from
        // Meeus "Astronomical Algorithms", 1991, p. 101-102.
        // I found the following Go implementation:
        // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        // This is a translation from the function "Saemundsson" there.
        // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        // This is important because the 'refr' formula below goes crazy near hd = -5.11.
        hd = altitude;
        if (hd < -1.0)
            hd = -1.0;

        refr = (1.02 / tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

        if (refraction == REFRACTION_NORMAL && altitude < -1.0)
        {
            // In "normal" mode we gradually reduce refraction toward the nadir
            // so that we never get an altitude angle less than -90 degrees.
            // When horizon angle is -1 degrees, the factor is exactly 1.
            // As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            refr *= (altitude + 90.0) / 89.0;
        }
    }
    else
    {
        /* No refraction, or the refraction option is invalid. */
        refr = 0.0;
    }

    return refr;
}


/**
 * @brief
 *      Calculates the inverse of an atmospheric refraction angle.
 *
 * Given an observed altitude angle that includes atmospheric refraction,
 * calculate the negative angular correction to obtain the unrefracted
 * altitude. This is useful for cases where observed horizontal
 * coordinates are to be converted to another orientation system,
 * but refraction first must be removed from the observed position.
 *
 * @param refraction
 *      The option selecting which refraction correction to use.
 *      See #Astronomy_Refraction.
 *
 * @param bent_altitude
 *      The apparent altitude that includes atmospheric refraction.
 *
 * @return
 *      The angular adjustment in degrees to be added to the
 *      altitude angle to correct for atmospheric lensing.
 *      This will be less than or equal to zero.
 */
double Astronomy_InverseRefraction(astro_refraction_t refraction, double bent_altitude)
{
    double altitude, diff;

    if (bent_altitude < -90.0 || bent_altitude > +90.0)
        return 0.0;     /* no attempt to correct an invalid altitude */

    /* Find the pre-adjusted altitude whose refraction correction leads to 'altitude'. */
    altitude = bent_altitude - Astronomy_Refraction(refraction, bent_altitude);
    for(;;)
    {
        /* See how close we got. */
        diff = (altitude + Astronomy_Refraction(refraction, altitude)) - bent_altitude;
        if (fabs(diff) < 1.0e-14)
            return altitude - bent_altitude;

        altitude -= diff;
    }
}

/**
 * @brief
 *      Applies a rotation to a vector, yielding a rotated vector.
 *
 * This function transforms a vector in one orientation to a vector
 * in another orientation.
 *
 * @param rotation
 *      A rotation matrix that specifies how the orientation of the vector is to be changed.
 *
 * @param vector
 *      The vector whose orientation is to be changed.
 *
 * @return
 *      A vector in the orientation specified by `rotation`.
 */
astro_vector_t Astronomy_RotateVector(astro_rotation_t rotation, astro_vector_t vector)
{
    astro_vector_t target;

    if (rotation.status != ASTRO_SUCCESS || vector.status != ASTRO_SUCCESS)
        return VecError(ASTRO_INVALID_PARAMETER, vector.t);

    target.status = ASTRO_SUCCESS;
    target.t = vector.t;
    target.x = rotation.rot[0][0]*vector.x + rotation.rot[1][0]*vector.y + rotation.rot[2][0]*vector.z;
    target.y = rotation.rot[0][1]*vector.x + rotation.rot[1][1]*vector.y + rotation.rot[2][1]*vector.z;
    target.z = rotation.rot[0][2]*vector.x + rotation.rot[1][2]*vector.y + rotation.rot[2][2]*vector.z;

    return target;
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @return
 *      A rotation matrix that converts EQJ to ECL.
 */
astro_rotation_t Astronomy_Rotation_EQJ_ECL(void)
{
    /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
    static const double c = 0.9174821430670688;    /* cos(ob) */
    static const double s = 0.3977769691083922;    /* sin(ob) */
    astro_rotation_t r;

    r.status = ASTRO_SUCCESS;
    r.rot[0][0] = 1.0;  r.rot[1][0] = 0.0;  r.rot[2][0] = 0.0;
    r.rot[0][1] = 0.0;  r.rot[1][1] = +c;   r.rot[2][1] = +s;
    r.rot[0][2] = 0.0;  r.rot[1][2] = -s;   r.rot[2][2] = +c;
    return r;
}

/**
 * @brief
 *      Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @return
 *      A rotation matrix that converts ECL to EQJ.
 */
astro_rotation_t Astronomy_Rotation_ECL_EQJ(void)
{
    /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
    static const double c = 0.9174821430670688;    /* cos(ob) */
    static const double s = 0.3977769691083922;    /* sin(ob) */
    astro_rotation_t r;

    r.status = ASTRO_SUCCESS;
    r.rot[0][0] = 1.0;  r.rot[1][0] = 0.0;  r.rot[2][0] = 0.0;
    r.rot[0][1] = 0.0;  r.rot[1][1] = +c;   r.rot[2][1] = -s;
    r.rot[0][2] = 0.0;  r.rot[1][2] = +s;   r.rot[2][2] = +c;
    return r;
}

/**
 * @brief
 *      Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the target orientation.
 *
 * @return
 *      A rotation matrix that converts EQJ to EQD at `time`.
 */
astro_rotation_t Astronomy_Rotation_EQJ_EQD(astro_time_t time)
{
    astro_rotation_t prec, nut;

    prec = precession_rot(0.0, time.tt);
    nut = nutation_rot(&time, 0);
    return Astronomy_CombineRotation(prec, nut);
}

/**
 * @brief
 *      Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the source orientation.
 *
 * @return
 *      A rotation matrix that converts EQD at `time` to EQJ.
 */
astro_rotation_t Astronomy_Rotation_EQD_EQJ(astro_time_t time)
{
    astro_rotation_t prec, nut;

    nut = nutation_rot(&time, 1);
    prec = precession_rot(time.tt, 0.0);
    return Astronomy_CombineRotation(nut, prec);
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time at which the Earth's equator applies.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts EQD to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
astro_rotation_t Astronomy_Rotation_EQD_HOR(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t rot;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double spin_angle;

    double sinlat = sin(observer.latitude * DEG2RAD);
    double coslat = cos(observer.latitude * DEG2RAD);
    double sinlon = sin(observer.longitude * DEG2RAD);
    double coslon = cos(observer.longitude * DEG2RAD);

    uze[0] = coslat * coslon;
    uze[1] = coslat * sinlon;
    uze[2] = sinlat;

    une[0] = -sinlat * coslon;
    une[1] = -sinlat * sinlon;
    une[2] = coslat;

    uwe[0] = sinlon;
    uwe[1] = -coslon;
    uwe[2] = 0.0;

    spin_angle = -15.0 * sidereal_time(&time);
    spin(spin_angle, uze, uz);
    spin(spin_angle, une, un);
    spin(spin_angle, uwe, uw);

    rot.rot[0][0] = un[0]; rot.rot[1][0] = un[1]; rot.rot[2][0] = un[2];
    rot.rot[0][1] = uw[0]; rot.rot[1][1] = uw[1]; rot.rot[2][1] = uw[2];
    rot.rot[0][2] = uz[0]; rot.rot[1][2] = uz[1]; rot.rot[2][2] = uz[2];

    rot.status = ASTRO_SUCCESS;
    return rot;
}


/**
 * @brief
 *      Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param time
 *      The date and time at which the Earth's equator applies.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
astro_rotation_t Astronomy_Rotation_HOR_EQD(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t rot = Astronomy_Rotation_EQD_HOR(time, observer);
    return Astronomy_InverseRotation(rot);
}


/**
 * @brief
 *      Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQJ = equatorial system, using equator at the J2000 epoch.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
astro_rotation_t Astronomy_Rotation_HOR_EQJ(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t hor_eqd, eqd_eqj;

    hor_eqd = Astronomy_Rotation_HOR_EQD(time, observer);
    eqd_eqj = Astronomy_Rotation_EQD_EQJ(time);
    return Astronomy_CombineRotation(hor_eqd, eqd_eqj);
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time of the desired horizontal orientation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
astro_rotation_t Astronomy_Rotation_EQJ_HOR(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t rot = Astronomy_Rotation_HOR_EQJ(time, observer);
    return Astronomy_InverseRotation(rot);
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of date.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time of the source equator.
 *
 * @return
 *      A rotation matrix that converts EQD to ECL.
 */
astro_rotation_t Astronomy_Rotation_EQD_ECL(astro_time_t time)
{
    astro_rotation_t eqd_eqj;
    astro_rotation_t eqj_ecl;

    eqd_eqj = Astronomy_Rotation_EQD_EQJ(time);
    eqj_ecl = Astronomy_Rotation_EQJ_ECL();
    return Astronomy_CombineRotation(eqd_eqj, eqj_ecl);
}


/**
 * @brief
 *      Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of date.
 *
 * @param time
 *      The date and time of the desired equator.
 *
 * @return
 *      A rotation matrix that converts ECL to EQD.
 */
astro_rotation_t Astronomy_Rotation_ECL_EQD(astro_time_t time)
{
    astro_rotation_t rot = Astronomy_Rotation_EQD_ECL(time);
    return Astronomy_InverseRotation(rot);
}

/**
 * @brief
 *      Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time of the desired horizontal orientation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts ECL to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
astro_rotation_t Astronomy_Rotation_ECL_HOR(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t ecl_eqd = Astronomy_Rotation_ECL_EQD(time);
    astro_rotation_t eqd_hor = Astronomy_Rotation_EQD_HOR(time, observer);
    return Astronomy_CombineRotation(ecl_eqd, eqd_hor);
}

/**
 * @brief
 *      Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time of the horizontal observation.
 *
 * @param observer
 *      The location of the horizontal observer.
 *
 * @return
 *      A rotation matrix that converts HOR to ECL.
 */
astro_rotation_t Astronomy_Rotation_HOR_ECL(astro_time_t time, astro_observer_t observer)
{
    astro_rotation_t rot = Astronomy_Rotation_ECL_HOR(time, observer);
    return Astronomy_InverseRotation(rot);
}


/** @cond DOXYGEN_SKIP */
typedef struct
{
    const char *symbol;
    const char *name;
}
constel_info_t;


typedef struct
{
    int    index;
    double ra_lo;
    double ra_hi;
    double dec_lo;
}
constel_boundary_t;
/** @endcond */

#define NUM_CONSTELLATIONS   88

static const constel_info_t ConstelInfo[] = {
    /*  0 */ { "And", "Andromeda"            }
,   /*  1 */ { "Ant", "Antila"               }
,   /*  2 */ { "Aps", "Apus"                 }
,   /*  3 */ { "Aql", "Aquila"               }
,   /*  4 */ { "Aqr", "Aquarius"             }
,   /*  5 */ { "Ara", "Ara"                  }
,   /*  6 */ { "Ari", "Aries"                }
,   /*  7 */ { "Aur", "Auriga"               }
,   /*  8 */ { "Boo", "Bootes"               }
,   /*  9 */ { "Cae", "Caelum"               }
,   /* 10 */ { "Cam", "Camelopardis"         }
,   /* 11 */ { "Cap", "Capricornus"          }
,   /* 12 */ { "Car", "Carina"               }
,   /* 13 */ { "Cas", "Cassiopeia"           }
,   /* 14 */ { "Cen", "Centaurus"            }
,   /* 15 */ { "Cep", "Cepheus"              }
,   /* 16 */ { "Cet", "Cetus"                }
,   /* 17 */ { "Cha", "Chamaeleon"           }
,   /* 18 */ { "Cir", "Circinus"             }
,   /* 19 */ { "CMa", "Canis Major"          }
,   /* 20 */ { "CMi", "Canis Minor"          }
,   /* 21 */ { "Cnc", "Cancer"               }
,   /* 22 */ { "Col", "Columba"              }
,   /* 23 */ { "Com", "Coma Berenices"       }
,   /* 24 */ { "CrA", "Corona Australis"     }
,   /* 25 */ { "CrB", "Corona Borealis"      }
,   /* 26 */ { "Crt", "Crater"               }
,   /* 27 */ { "Cru", "Crux"                 }
,   /* 28 */ { "Crv", "Corvus"               }
,   /* 29 */ { "CVn", "Canes Venatici"       }
,   /* 30 */ { "Cyg", "Cygnus"               }
,   /* 31 */ { "Del", "Delphinus"            }
,   /* 32 */ { "Dor", "Dorado"               }
,   /* 33 */ { "Dra", "Draco"                }
,   /* 34 */ { "Equ", "Equuleus"             }
,   /* 35 */ { "Eri", "Eridanus"             }
,   /* 36 */ { "For", "Fornax"               }
,   /* 37 */ { "Gem", "Gemini"               }
,   /* 38 */ { "Gru", "Grus"                 }
,   /* 39 */ { "Her", "Hercules"             }
,   /* 40 */ { "Hor", "Horologium"           }
,   /* 41 */ { "Hya", "Hydra"                }
,   /* 42 */ { "Hyi", "Hydrus"               }
,   /* 43 */ { "Ind", "Indus"                }
,   /* 44 */ { "Lac", "Lacerta"              }
,   /* 45 */ { "Leo", "Leo"                  }
,   /* 46 */ { "Lep", "Lepus"                }
,   /* 47 */ { "Lib", "Libra"                }
,   /* 48 */ { "LMi", "Leo Minor"            }
,   /* 49 */ { "Lup", "Lupus"                }
,   /* 50 */ { "Lyn", "Lynx"                 }
,   /* 51 */ { "Lyr", "Lyra"                 }
,   /* 52 */ { "Men", "Mensa"                }
,   /* 53 */ { "Mic", "Microscopium"         }
,   /* 54 */ { "Mon", "Monoceros"            }
,   /* 55 */ { "Mus", "Musca"                }
,   /* 56 */ { "Nor", "Norma"                }
,   /* 57 */ { "Oct", "Octans"               }
,   /* 58 */ { "Oph", "Ophiuchus"            }
,   /* 59 */ { "Ori", "Orion"                }
,   /* 60 */ { "Pav", "Pavo"                 }
,   /* 61 */ { "Peg", "Pegasus"              }
,   /* 62 */ { "Per", "Perseus"              }
,   /* 63 */ { "Phe", "Phoenix"              }
,   /* 64 */ { "Pic", "Pictor"               }
,   /* 65 */ { "PsA", "Pisces Austrinus"     }
,   /* 66 */ { "Psc", "Pisces"               }
,   /* 67 */ { "Pup", "Puppis"               }
,   /* 68 */ { "Pyx", "Pyxis"                }
,   /* 69 */ { "Ret", "Reticulum"            }
,   /* 70 */ { "Scl", "Sculptor"             }
,   /* 71 */ { "Sco", "Scorpius"             }
,   /* 72 */ { "Sct", "Scutum"               }
,   /* 73 */ { "Ser", "Serpens"              }
,   /* 74 */ { "Sex", "Sextans"              }
,   /* 75 */ { "Sge", "Sagitta"              }
,   /* 76 */ { "Sgr", "Sagittarius"          }
,   /* 77 */ { "Tau", "Taurus"               }
,   /* 78 */ { "Tel", "Telescopium"          }
,   /* 79 */ { "TrA", "Triangulum Australe"  }
,   /* 80 */ { "Tri", "Triangulum"           }
,   /* 81 */ { "Tuc", "Tucana"               }
,   /* 82 */ { "UMa", "Ursa Major"           }
,   /* 83 */ { "UMi", "Ursa Minor"           }
,   /* 84 */ { "Vel", "Vela"                 }
,   /* 85 */ { "Vir", "Virgo"                }
,   /* 86 */ { "Vol", "Volans"               }
,   /* 87 */ { "Vul", "Vulpecula"            }
};

static const constel_boundary_t ConstelBounds[] = {
    { 83,  0.00000000000000, 24.00000000000000, 88.00000000000000 }    /* UMi */
,   { 83,  8.00000000000000, 14.50000000000000, 86.50000000000000 }    /* UMi */
,   { 83, 21.00000000000000, 23.00000000000000, 86.16666666666667 }    /* UMi */
,   { 83, 18.00000000000000, 21.00000000000000, 86.00000000000000 }    /* UMi */
,   { 15,  0.00000000000000,  8.00000000000000, 85.00000000000000 }    /* Cep */
,   { 10,  9.16666666666667, 10.66666666666667, 82.00000000000000 }    /* Cam */
,   { 15,  0.00000000000000,  5.00000000000000, 80.00000000000000 }    /* Cep */
,   { 10, 10.66666666666667, 14.50000000000000, 80.00000000000000 }    /* Cam */
,   { 83, 17.50000000000000, 18.00000000000000, 80.00000000000000 }    /* UMi */
,   { 33, 20.16666666666667, 21.00000000000000, 80.00000000000000 }    /* Dra */
,   { 15,  0.00000000000000,  3.50833333333333, 77.00000000000000 }    /* Cep */
,   { 10, 11.50000000000000, 13.58333333333333, 77.00000000000000 }    /* Cam */
,   { 83, 16.53333333333333, 17.50000000000000, 75.00000000000000 }    /* UMi */
,   { 15, 20.16666666666667, 20.66666666666667, 75.00000000000000 }    /* Cep */
,   { 10,  7.96666666666667,  9.16666666666667, 73.50000000000000 }    /* Cam */
,   { 33,  9.16666666666667, 11.33333333333333, 73.50000000000000 }    /* Dra */
,   { 83, 13.00000000000000, 16.53333333333333, 70.00000000000000 }    /* UMi */
,   { 13,  3.10000000000000,  3.41666666666667, 68.00000000000000 }    /* Cas */
,   { 33, 20.41666666666667, 20.66666666666667, 67.00000000000000 }    /* Dra */
,   { 33, 11.33333333333333, 12.00000000000000, 66.50000000000000 }    /* Dra */
,   { 15,  0.00000000000000,  0.33333333333333, 66.00000000000000 }    /* Cep */
,   { 83, 14.00000000000000, 15.66666666666667, 66.00000000000000 }    /* UMi */
,   { 15, 23.58333333333333, 24.00000000000000, 66.00000000000000 }    /* Cep */
,   { 33, 12.00000000000000, 13.50000000000000, 64.00000000000000 }    /* Dra */
,   { 33, 13.50000000000000, 14.41666666666667, 63.00000000000000 }    /* Dra */
,   { 15, 23.16666666666667, 23.58333333333333, 63.00000000000000 }    /* Cep */
,   { 10,  6.10000000000000,  7.00000000000000, 62.00000000000000 }    /* Cam */
,   { 33, 20.00000000000000, 20.41666666666667, 61.50000000000000 }    /* Dra */
,   { 15, 20.53666666666667, 20.60000000000000, 60.91666666666666 }    /* Cep */
,   { 10,  7.00000000000000,  7.96666666666667, 60.00000000000000 }    /* Cam */
,   { 82,  7.96666666666667,  8.41666666666667, 60.00000000000000 }    /* UMa */
,   { 33, 19.76666666666667, 20.00000000000000, 59.50000000000000 }    /* Dra */
,   { 15, 20.00000000000000, 20.53666666666667, 59.50000000000000 }    /* Cep */
,   { 15, 22.86666666666667, 23.16666666666667, 59.08333333333334 }    /* Cep */
,   { 13,  0.00000000000000,  2.43333333333333, 58.50000000000000 }    /* Cas */
,   { 33, 19.41666666666667, 19.76666666666667, 58.00000000000000 }    /* Dra */
,   { 13,  1.70000000000000,  1.90833333333333, 57.50000000000000 }    /* Cas */
,   { 13,  2.43333333333333,  3.10000000000000, 57.00000000000000 }    /* Cas */
,   { 10,  3.10000000000000,  3.16666666666667, 57.00000000000000 }    /* Cam */
,   { 15, 22.31666666666667, 22.86666666666667, 56.25000000000000 }    /* Cep */
,   { 10,  5.00000000000000,  6.10000000000000, 56.00000000000000 }    /* Cam */
,   { 82, 14.03333333333333, 14.41666666666667, 55.50000000000000 }    /* UMa */
,   { 33, 14.41666666666667, 19.41666666666667, 55.50000000000000 }    /* Dra */
,   { 10,  3.16666666666667,  3.33333333333333, 55.00000000000000 }    /* Cam */
,   { 15, 22.13333333333333, 22.31666666666667, 55.00000000000000 }    /* Cep */
,   { 15, 20.60000000000000, 21.96666666666667, 54.83333333333334 }    /* Cep */
,   { 13,  0.00000000000000,  1.70000000000000, 54.00000000000000 }    /* Cas */
,   { 50,  6.10000000000000,  6.50000000000000, 54.00000000000000 }    /* Lyn */
,   { 82, 12.08333333333333, 13.50000000000000, 53.00000000000000 }    /* UMa */
,   { 33, 15.25000000000000, 15.75000000000000, 53.00000000000000 }    /* Dra */
,   { 15, 21.96666666666667, 22.13333333333333, 52.75000000000000 }    /* Cep */
,   { 10,  3.33333333333333,  5.00000000000000, 52.50000000000000 }    /* Cam */
,   { 13, 22.86666666666667, 23.33333333333333, 52.50000000000000 }    /* Cas */
,   { 33, 15.75000000000000, 17.00000000000000, 51.50000000000000 }    /* Dra */
,   { 62,  2.04166666666667,  2.51666666666667, 50.50000000000000 }    /* Per */
,   { 33, 17.00000000000000, 18.23333333333333, 50.50000000000000 }    /* Dra */
,   { 13,  0.00000000000000,  1.36666666666667, 50.00000000000000 }    /* Cas */
,   { 62,  1.36666666666667,  1.66666666666667, 50.00000000000000 }    /* Per */
,   { 50,  6.50000000000000,  6.80000000000000, 50.00000000000000 }    /* Lyn */
,   { 13, 23.33333333333333, 24.00000000000000, 50.00000000000000 }    /* Cas */
,   { 82, 13.50000000000000, 14.03333333333333, 48.50000000000000 }    /* UMa */
,   { 13,  0.00000000000000,  1.11666666666667, 48.00000000000000 }    /* Cas */
,   { 13, 23.58333333333333, 24.00000000000000, 48.00000000000000 }    /* Cas */
,   { 39, 18.17500000000000, 18.23333333333333, 47.50000000000000 }    /* Her */
,   { 33, 18.23333333333333, 19.08333333333333, 47.50000000000000 }    /* Dra */
,   { 30, 19.08333333333333, 19.16666666666667, 47.50000000000000 }    /* Cyg */
,   { 62,  1.66666666666667,  2.04166666666667, 47.00000000000000 }    /* Per */
,   { 82,  8.41666666666667,  9.16666666666667, 47.00000000000000 }    /* UMa */
,   { 13,  0.16666666666667,  0.86666666666667, 46.00000000000000 }    /* Cas */
,   { 82, 12.00000000000000, 12.08333333333333, 45.00000000000000 }    /* UMa */
,   { 50,  6.80000000000000,  7.36666666666667, 44.50000000000000 }    /* Lyn */
,   { 30, 21.90833333333333, 21.96666666666667, 44.00000000000000 }    /* Cyg */
,   { 30, 21.87500000000000, 21.90833333333333, 43.75000000000000 }    /* Cyg */
,   { 30, 19.16666666666667, 19.40000000000000, 43.50000000000000 }    /* Cyg */
,   { 82,  9.16666666666667, 10.16666666666667, 42.00000000000000 }    /* UMa */
,   { 82, 10.16666666666667, 10.78333333333333, 40.00000000000000 }    /* UMa */
,   {  8, 15.43333333333333, 15.75000000000000, 40.00000000000000 }    /* Boo */
,   { 39, 15.75000000000000, 16.33333333333333, 40.00000000000000 }    /* Her */
,   { 50,  9.25000000000000,  9.58333333333333, 39.75000000000000 }    /* Lyn */
,   {  0,  0.00000000000000,  2.51666666666667, 36.75000000000000 }    /* And */
,   { 62,  2.51666666666667,  2.56666666666667, 36.75000000000000 }    /* Per */
,   { 51, 19.35833333333333, 19.40000000000000, 36.50000000000000 }    /* Lyr */
,   { 62,  4.50000000000000,  4.69166666666667, 36.00000000000000 }    /* Per */
,   { 30, 21.73333333333333, 21.87500000000000, 36.00000000000000 }    /* Cyg */
,   { 44, 21.87500000000000, 22.00000000000000, 36.00000000000000 }    /* Lac */
,   {  7,  6.53333333333333,  7.36666666666667, 35.50000000000000 }    /* Aur */
,   { 50,  7.36666666666667,  7.75000000000000, 35.50000000000000 }    /* Lyn */
,   {  0,  0.00000000000000,  2.00000000000000, 35.00000000000000 }    /* And */
,   { 44, 22.00000000000000, 22.81666666666667, 35.00000000000000 }    /* Lac */
,   { 44, 22.81666666666667, 22.86666666666667, 34.50000000000000 }    /* Lac */
,   {  0, 22.86666666666667, 23.50000000000000, 34.50000000000000 }    /* And */
,   { 62,  2.56666666666667,  2.71666666666667, 34.00000000000000 }    /* Per */
,   { 82, 10.78333333333333, 11.00000000000000, 34.00000000000000 }    /* UMa */
,   { 29, 12.00000000000000, 12.33333333333333, 34.00000000000000 }    /* CVn */
,   { 50,  7.75000000000000,  9.25000000000000, 33.50000000000000 }    /* Lyn */
,   { 48,  9.25000000000000,  9.88333333333333, 33.50000000000000 }    /* LMi */
,   {  0,  0.71666666666667,  1.40833333333333, 33.00000000000000 }    /* And */
,   {  8, 15.18333333333333, 15.43333333333333, 33.00000000000000 }    /* Boo */
,   {  0, 23.50000000000000, 23.75000000000000, 32.08333333333334 }    /* And */
,   { 29, 12.33333333333333, 13.25000000000000, 32.00000000000000 }    /* CVn */
,   {  0, 23.75000000000000, 24.00000000000000, 31.33333333333333 }    /* And */
,   { 29, 13.95833333333333, 14.03333333333333, 30.75000000000000 }    /* CVn */
,   { 80,  2.41666666666667,  2.71666666666667, 30.66666666666667 }    /* Tri */
,   { 62,  2.71666666666667,  4.50000000000000, 30.66666666666667 }    /* Per */
,   {  7,  4.50000000000000,  4.75000000000000, 30.00000000000000 }    /* Aur */
,   { 51, 18.17500000000000, 19.35833333333333, 30.00000000000000 }    /* Lyr */
,   { 82, 11.00000000000000, 12.00000000000000, 29.00000000000000 }    /* UMa */
,   { 30, 19.66666666666667, 20.91666666666667, 29.00000000000000 }    /* Cyg */
,   {  7,  4.75000000000000,  5.88333333333333, 28.50000000000000 }    /* Aur */
,   { 48,  9.88333333333333, 10.50000000000000, 28.50000000000000 }    /* LMi */
,   { 29, 13.25000000000000, 13.95833333333333, 28.50000000000000 }    /* CVn */
,   {  0,  0.00000000000000,  0.06666666666667, 28.00000000000000 }    /* And */
,   { 80,  1.40833333333333,  1.66666666666667, 28.00000000000000 }    /* Tri */
,   {  7,  5.88333333333333,  6.53333333333333, 28.00000000000000 }    /* Aur */
,   { 37,  7.88333333333333,  8.00000000000000, 28.00000000000000 }    /* Gem */
,   { 30, 20.91666666666667, 21.73333333333333, 28.00000000000000 }    /* Cyg */
,   { 30, 19.25833333333333, 19.66666666666667, 27.50000000000000 }    /* Cyg */
,   { 80,  1.91666666666667,  2.41666666666667, 27.25000000000000 }    /* Tri */
,   { 25, 16.16666666666667, 16.33333333333333, 27.00000000000000 }    /* CrB */
,   {  8, 15.08333333333333, 15.18333333333333, 26.00000000000000 }    /* Boo */
,   { 25, 15.18333333333333, 16.16666666666667, 26.00000000000000 }    /* CrB */
,   { 51, 18.36666666666667, 18.86666666666667, 26.00000000000000 }    /* Lyr */
,   { 48, 10.75000000000000, 11.00000000000000, 25.50000000000000 }    /* LMi */
,   { 51, 18.86666666666667, 19.25833333333333, 25.50000000000000 }    /* Lyr */
,   { 80,  1.66666666666667,  1.91666666666667, 25.00000000000000 }    /* Tri */
,   { 66,  0.71666666666667,  0.85000000000000, 23.75000000000000 }    /* Psc */
,   { 48, 10.50000000000000, 10.75000000000000, 23.50000000000000 }    /* LMi */
,   { 87, 21.25000000000000, 21.41666666666667, 23.50000000000000 }    /* Vul */
,   { 77,  5.70000000000000,  5.88333333333333, 22.83333333333333 }    /* Tau */
,   {  0,  0.06666666666667,  0.14166666666667, 22.00000000000000 }    /* And */
,   { 73, 15.91666666666667, 16.03333333333333, 22.00000000000000 }    /* Ser */
,   { 37,  5.88333333333333,  6.21666666666667, 21.50000000000000 }    /* Gem */
,   { 87, 19.83333333333333, 20.25000000000000, 21.25000000000000 }    /* Vul */
,   { 87, 18.86666666666667, 19.25000000000000, 21.08333333333333 }    /* Vul */
,   {  0,  0.14166666666667,  0.85000000000000, 21.00000000000000 }    /* And */
,   { 87, 20.25000000000000, 20.56666666666667, 20.50000000000000 }    /* Vul */
,   { 37,  7.80833333333333,  7.88333333333333, 20.00000000000000 }    /* Gem */
,   { 87, 20.56666666666667, 21.25000000000000, 19.50000000000000 }    /* Vul */
,   { 87, 19.25000000000000, 19.83333333333333, 19.16666666666667 }    /* Vul */
,   {  6,  3.28333333333333,  3.36666666666667, 19.00000000000000 }    /* Ari */
,   { 75, 18.86666666666667, 19.00000000000000, 18.50000000000000 }    /* Sge */
,   { 59,  5.70000000000000,  5.76666666666667, 18.00000000000000 }    /* Ori */
,   { 37,  6.21666666666667,  6.30833333333333, 17.50000000000000 }    /* Gem */
,   { 75, 19.00000000000000, 19.83333333333333, 16.16666666666667 }    /* Sge */
,   { 77,  4.96666666666667,  5.33333333333333, 16.00000000000000 }    /* Tau */
,   { 39, 15.91666666666667, 16.08333333333333, 16.00000000000000 }    /* Her */
,   { 75, 19.83333333333333, 20.25000000000000, 15.75000000000000 }    /* Sge */
,   { 77,  4.61666666666667,  4.96666666666667, 15.50000000000000 }    /* Tau */
,   { 77,  5.33333333333333,  5.60000000000000, 15.50000000000000 }    /* Tau */
,   { 23, 12.83333333333333, 13.50000000000000, 15.00000000000000 }    /* Com */
,   { 39, 17.25000000000000, 18.25000000000000, 14.33333333333333 }    /* Her */
,   { 23, 11.86666666666667, 12.83333333333333, 14.00000000000000 }    /* Com */
,   { 37,  7.50000000000000,  7.80833333333333, 13.50000000000000 }    /* Gem */
,   { 39, 16.75000000000000, 17.25000000000000, 12.83333333333333 }    /* Her */
,   { 61,  0.00000000000000,  0.14166666666667, 12.50000000000000 }    /* Peg */
,   { 77,  5.60000000000000,  5.76666666666667, 12.50000000000000 }    /* Tau */
,   { 37,  7.00000000000000,  7.50000000000000, 12.50000000000000 }    /* Gem */
,   { 61, 21.11666666666667, 21.33333333333333, 12.50000000000000 }    /* Peg */
,   { 37,  6.30833333333333,  6.93333333333333, 12.00000000000000 }    /* Gem */
,   { 39, 18.25000000000000, 18.86666666666667, 12.00000000000000 }    /* Her */
,   { 31, 20.87500000000000, 21.05000000000000, 11.83333333333333 }    /* Del */
,   { 61, 21.05000000000000, 21.11666666666667, 11.83333333333333 }    /* Peg */
,   { 45, 11.51666666666667, 11.86666666666667, 11.00000000000000 }    /* Leo */
,   { 59,  6.24166666666667,  6.30833333333333, 10.00000000000000 }    /* Ori */
,   { 37,  6.93333333333333,  7.00000000000000, 10.00000000000000 }    /* Gem */
,   { 21,  7.80833333333333,  7.92500000000000, 10.00000000000000 }    /* Cnc */
,   { 61, 23.83333333333333, 24.00000000000000, 10.00000000000000 }    /* Peg */
,   {  6,  1.66666666666667,  3.28333333333333,  9.91666666666667 }    /* Ari */
,   { 31, 20.14166666666667, 20.30000000000000,  8.50000000000000 }    /* Del */
,   {  8, 13.50000000000000, 15.08333333333333,  8.00000000000000 }    /* Boo */
,   { 61, 22.75000000000000, 23.83333333333333,  7.50000000000000 }    /* Peg */
,   { 21,  7.92500000000000,  9.25000000000000,  7.00000000000000 }    /* Cnc */
,   { 45,  9.25000000000000, 10.75000000000000,  7.00000000000000 }    /* Leo */
,   { 58, 18.25000000000000, 18.66222222222222,  6.25000000000000 }    /* Oph */
,   {  3, 18.66222222222222, 18.86666666666667,  6.25000000000000 }    /* Aql */
,   { 31, 20.83333333333333, 20.87500000000000,  6.00000000000000 }    /* Del */
,   { 20,  7.00000000000000,  7.01666666666667,  5.50000000000000 }    /* CMi */
,   { 73, 18.25000000000000, 18.42500000000000,  4.50000000000000 }    /* Ser */
,   { 39, 16.08333333333333, 16.75000000000000,  4.00000000000000 }    /* Her */
,   { 58, 18.25000000000000, 18.42500000000000,  3.00000000000000 }    /* Oph */
,   { 61, 21.46666666666667, 21.66666666666667,  2.75000000000000 }    /* Peg */
,   { 66,  0.00000000000000,  2.00000000000000,  2.00000000000000 }    /* Psc */
,   { 73, 18.58333333333333, 18.86666666666667,  2.00000000000000 }    /* Ser */
,   { 31, 20.30000000000000, 20.83333333333333,  2.00000000000000 }    /* Del */
,   { 34, 20.83333333333333, 21.33333333333333,  2.00000000000000 }    /* Equ */
,   { 61, 21.33333333333333, 21.46666666666667,  2.00000000000000 }    /* Peg */
,   { 61, 22.00000000000000, 22.75000000000000,  2.00000000000000 }    /* Peg */
,   { 61, 21.66666666666667, 22.00000000000000,  1.75000000000000 }    /* Peg */
,   { 20,  7.01666666666667,  7.20000000000000,  1.50000000000000 }    /* CMi */
,   { 77,  3.58333333333333,  4.61666666666667,  0.00000000000000 }    /* Tau */
,   { 59,  4.61666666666667,  4.66666666666667,  0.00000000000000 }    /* Ori */
,   { 20,  7.20000000000000,  8.08333333333333,  0.00000000000000 }    /* CMi */
,   { 85, 14.66666666666667, 15.08333333333333,  0.00000000000000 }    /* Vir */
,   { 58, 17.83333333333333, 18.25000000000000,  0.00000000000000 }    /* Oph */
,   { 16,  2.65000000000000,  3.28333333333333, -1.75000000000000 }    /* Cet */
,   { 77,  3.28333333333333,  3.58333333333333, -1.75000000000000 }    /* Tau */
,   { 73, 15.08333333333333, 16.26666666666667, -3.25000000000000 }    /* Ser */
,   { 59,  4.66666666666667,  5.08333333333333, -4.00000000000000 }    /* Ori */
,   { 59,  5.83333333333333,  6.24166666666667, -4.00000000000000 }    /* Ori */
,   { 73, 17.83333333333333, 17.96666666666667, -4.00000000000000 }    /* Ser */
,   { 73, 18.25000000000000, 18.58333333333333, -4.00000000000000 }    /* Ser */
,   {  3, 18.58333333333333, 18.86666666666667, -4.00000000000000 }    /* Aql */
,   { 66, 22.75000000000000, 23.83333333333333, -4.00000000000000 }    /* Psc */
,   { 45, 10.75000000000000, 11.51666666666667, -6.00000000000000 }    /* Leo */
,   { 85, 11.51666666666667, 11.83333333333333, -6.00000000000000 }    /* Vir */
,   { 66,  0.00000000000000,  0.33333333333333, -7.00000000000000 }    /* Psc */
,   { 66, 23.83333333333333, 24.00000000000000, -7.00000000000000 }    /* Psc */
,   { 85, 14.25000000000000, 14.66666666666667, -8.00000000000000 }    /* Vir */
,   { 58, 15.91666666666667, 16.26666666666667, -8.00000000000000 }    /* Oph */
,   {  3, 20.00000000000000, 20.53333333333333, -9.00000000000000 }    /* Aql */
,   {  4, 21.33333333333333, 21.86666666666667, -9.00000000000000 }    /* Aqr */
,   { 58, 17.16666666666667, 17.96666666666667, -10.00000000000000 }    /* Oph */
,   { 54,  5.83333333333333,  8.08333333333333, -11.00000000000000 }    /* Mon */
,   { 35,  4.91666666666667,  5.08333333333333, -11.00000000000000 }    /* Eri */
,   { 59,  5.08333333333333,  5.83333333333333, -11.00000000000000 }    /* Ori */
,   { 41,  8.08333333333333,  8.36666666666667, -11.00000000000000 }    /* Hya */
,   { 74,  9.58333333333333, 10.75000000000000, -11.00000000000000 }    /* Sex */
,   { 85, 11.83333333333333, 12.83333333333333, -11.00000000000000 }    /* Vir */
,   { 58, 17.58333333333333, 17.66666666666667, -11.66666666666667 }    /* Oph */
,   {  3, 18.86666666666667, 20.00000000000000, -12.03333333333333 }    /* Aql */
,   { 35,  4.83333333333333,  4.91666666666667, -14.50000000000000 }    /* Eri */
,   {  4, 20.53333333333333, 21.33333333333333, -15.00000000000000 }    /* Aqr */
,   { 73, 17.16666666666667, 18.25000000000000, -16.00000000000000 }    /* Ser */
,   { 72, 18.25000000000000, 18.86666666666667, -16.00000000000000 }    /* Sct */
,   { 41,  8.36666666666667,  8.58333333333333, -17.00000000000000 }    /* Hya */
,   { 58, 16.26666666666667, 16.37500000000000, -18.25000000000000 }    /* Oph */
,   { 41,  8.58333333333333,  9.08333333333333, -19.00000000000000 }    /* Hya */
,   { 26, 10.75000000000000, 10.83333333333333, -19.00000000000000 }    /* Crt */
,   { 71, 16.26666666666667, 16.37500000000000, -19.25000000000000 }    /* Sco */
,   { 47, 15.66666666666667, 15.91666666666667, -20.00000000000000 }    /* Lib */
,   { 28, 12.58333333333333, 12.83333333333333, -22.00000000000000 }    /* Crv */
,   { 85, 12.83333333333333, 14.25000000000000, -22.00000000000000 }    /* Vir */
,   { 41,  9.08333333333333,  9.75000000000000, -24.00000000000000 }    /* Hya */
,   { 16,  1.66666666666667,  2.65000000000000, -24.38333333333333 }    /* Cet */
,   { 35,  2.65000000000000,  3.75000000000000, -24.38333333333333 }    /* Eri */
,   { 26, 10.83333333333333, 11.83333333333333, -24.50000000000000 }    /* Crt */
,   { 28, 11.83333333333333, 12.58333333333333, -24.50000000000000 }    /* Crv */
,   { 47, 14.25000000000000, 14.91666666666667, -24.50000000000000 }    /* Lib */
,   { 58, 16.26666666666667, 16.75000000000000, -24.58333333333333 }    /* Oph */
,   { 16,  0.00000000000000,  1.66666666666667, -25.50000000000000 }    /* Cet */
,   { 11, 21.33333333333333, 21.86666666666667, -25.50000000000000 }    /* Cap */
,   {  4, 21.86666666666667, 23.83333333333333, -25.50000000000000 }    /* Aqr */
,   { 16, 23.83333333333333, 24.00000000000000, -25.50000000000000 }    /* Cet */
,   { 41,  9.75000000000000, 10.25000000000000, -26.50000000000000 }    /* Hya */
,   { 35,  4.70000000000000,  4.83333333333333, -27.25000000000000 }    /* Eri */
,   { 46,  4.83333333333333,  6.11666666666667, -27.25000000000000 }    /* Lep */
,   { 11, 20.00000000000000, 21.33333333333333, -28.00000000000000 }    /* Cap */
,   { 41, 10.25000000000000, 10.58333333333333, -29.16666666666667 }    /* Hya */
,   { 41, 12.58333333333333, 14.91666666666667, -29.50000000000000 }    /* Hya */
,   { 47, 14.91666666666667, 15.66666666666667, -29.50000000000000 }    /* Lib */
,   { 71, 15.66666666666667, 16.00000000000000, -29.50000000000000 }    /* Sco */
,   { 35,  4.58333333333333,  4.70000000000000, -30.00000000000000 }    /* Eri */
,   { 58, 16.75000000000000, 17.60000000000000, -30.00000000000000 }    /* Oph */
,   { 76, 17.60000000000000, 17.83333333333333, -30.00000000000000 }    /* Sgr */
,   { 41, 10.58333333333333, 10.83333333333333, -31.16666666666667 }    /* Hya */
,   { 19,  6.11666666666667,  7.36666666666667, -33.00000000000000 }    /* CMa */
,   { 41, 12.25000000000000, 12.58333333333333, -33.00000000000000 }    /* Hya */
,   { 41, 10.83333333333333, 12.25000000000000, -35.00000000000000 }    /* Hya */
,   { 36,  3.50000000000000,  3.75000000000000, -36.00000000000000 }    /* For */
,   { 68,  8.36666666666667,  9.36666666666667, -36.75000000000000 }    /* Pyx */
,   { 35,  4.26666666666667,  4.58333333333333, -37.00000000000000 }    /* Eri */
,   { 76, 17.83333333333333, 19.16666666666667, -37.00000000000000 }    /* Sgr */
,   { 65, 21.33333333333333, 23.00000000000000, -37.00000000000000 }    /* PsA */
,   { 70, 23.00000000000000, 23.33333333333333, -37.00000000000000 }    /* Scl */
,   { 36,  3.00000000000000,  3.50000000000000, -39.58333333333334 }    /* For */
,   {  1,  9.36666666666667, 11.00000000000000, -39.75000000000000 }    /* Ant */
,   { 70,  0.00000000000000,  1.66666666666667, -40.00000000000000 }    /* Scl */
,   { 36,  1.66666666666667,  3.00000000000000, -40.00000000000000 }    /* For */
,   { 35,  3.86666666666667,  4.26666666666667, -40.00000000000000 }    /* Eri */
,   { 70, 23.33333333333333, 24.00000000000000, -40.00000000000000 }    /* Scl */
,   { 14, 14.16666666666667, 14.91666666666667, -42.00000000000000 }    /* Cen */
,   { 49, 15.66666666666667, 16.00000000000000, -42.00000000000000 }    /* Lup */
,   { 71, 16.00000000000000, 16.42083333333333, -42.00000000000000 }    /* Sco */
,   {  9,  4.83333333333333,  5.00000000000000, -43.00000000000000 }    /* Cae */
,   { 22,  5.00000000000000,  6.58333333333333, -43.00000000000000 }    /* Col */
,   { 67,  8.00000000000000,  8.36666666666667, -43.00000000000000 }    /* Pup */
,   { 35,  3.41666666666667,  3.86666666666667, -44.00000000000000 }    /* Eri */
,   { 71, 16.42083333333333, 17.83333333333333, -45.50000000000000 }    /* Sco */
,   { 24, 17.83333333333333, 19.16666666666667, -45.50000000000000 }    /* CrA */
,   { 76, 19.16666666666667, 20.33333333333333, -45.50000000000000 }    /* Sgr */
,   { 53, 20.33333333333333, 21.33333333333333, -45.50000000000000 }    /* Mic */
,   { 35,  3.00000000000000,  3.41666666666667, -46.00000000000000 }    /* Eri */
,   {  9,  4.50000000000000,  4.83333333333333, -46.50000000000000 }    /* Cae */
,   { 49, 15.33333333333333, 15.66666666666667, -48.00000000000000 }    /* Lup */
,   { 63,  0.00000000000000,  2.33333333333333, -48.16666666666666 }    /* Phe */
,   { 35,  2.66666666666667,  3.00000000000000, -49.00000000000000 }    /* Eri */
,   { 40,  4.08333333333333,  4.26666666666667, -49.00000000000000 }    /* Hor */
,   {  9,  4.26666666666667,  4.50000000000000, -49.00000000000000 }    /* Cae */
,   { 38, 21.33333333333333, 22.00000000000000, -50.00000000000000 }    /* Gru */
,   { 67,  6.00000000000000,  8.00000000000000, -50.75000000000000 }    /* Pup */
,   { 84,  8.00000000000000,  8.16666666666667, -50.75000000000000 }    /* Vel */
,   { 35,  2.41666666666667,  2.66666666666667, -51.00000000000000 }    /* Eri */
,   { 40,  3.83333333333333,  4.08333333333333, -51.00000000000000 }    /* Hor */
,   { 63,  0.00000000000000,  1.83333333333333, -51.50000000000000 }    /* Phe */
,   { 12,  6.00000000000000,  6.16666666666667, -52.50000000000000 }    /* Car */
,   { 84,  8.16666666666667,  8.45000000000000, -53.00000000000000 }    /* Vel */
,   { 40,  3.50000000000000,  3.83333333333333, -53.16666666666666 }    /* Hor */
,   { 32,  3.83333333333333,  4.00000000000000, -53.16666666666666 }    /* Dor */
,   { 63,  0.00000000000000,  1.58333333333333, -53.50000000000000 }    /* Phe */
,   { 35,  2.16666666666667,  2.41666666666667, -54.00000000000000 }    /* Eri */
,   { 64,  4.50000000000000,  5.00000000000000, -54.00000000000000 }    /* Pic */
,   { 49, 15.05000000000000, 15.33333333333333, -54.00000000000000 }    /* Lup */
,   { 84,  8.45000000000000,  8.83333333333333, -54.50000000000000 }    /* Vel */
,   { 12,  6.16666666666667,  6.50000000000000, -55.00000000000000 }    /* Car */
,   { 14, 11.83333333333333, 12.83333333333333, -55.00000000000000 }    /* Cen */
,   { 49, 14.16666666666667, 15.05000000000000, -55.00000000000000 }    /* Lup */
,   { 56, 15.05000000000000, 15.33333333333333, -55.00000000000000 }    /* Nor */
,   { 32,  4.00000000000000,  4.33333333333333, -56.50000000000000 }    /* Dor */
,   { 84,  8.83333333333333, 11.00000000000000, -56.50000000000000 }    /* Vel */
,   { 14, 11.00000000000000, 11.25000000000000, -56.50000000000000 }    /* Cen */
,   {  5, 17.50000000000000, 18.00000000000000, -57.00000000000000 }    /* Ara */
,   { 78, 18.00000000000000, 20.33333333333333, -57.00000000000000 }    /* Tel */
,   { 38, 22.00000000000000, 23.33333333333333, -57.00000000000000 }    /* Gru */
,   { 40,  3.20000000000000,  3.50000000000000, -57.50000000000000 }    /* Hor */
,   { 64,  5.00000000000000,  5.50000000000000, -57.50000000000000 }    /* Pic */
,   { 12,  6.50000000000000,  6.83333333333333, -58.00000000000000 }    /* Car */
,   { 63,  0.00000000000000,  1.33333333333333, -58.50000000000000 }    /* Phe */
,   { 35,  1.33333333333333,  2.16666666666667, -58.50000000000000 }    /* Eri */
,   { 63, 23.33333333333333, 24.00000000000000, -58.50000000000000 }    /* Phe */
,   { 32,  4.33333333333333,  4.58333333333333, -59.00000000000000 }    /* Dor */
,   { 56, 15.33333333333333, 16.42083333333333, -60.00000000000000 }    /* Nor */
,   { 43, 20.33333333333333, 21.33333333333333, -60.00000000000000 }    /* Ind */
,   { 64,  5.50000000000000,  6.00000000000000, -61.00000000000000 }    /* Pic */
,   { 18, 15.16666666666667, 15.33333333333333, -61.00000000000000 }    /* Cir */
,   {  5, 16.42083333333333, 16.58333333333333, -61.00000000000000 }    /* Ara */
,   { 18, 14.91666666666667, 15.16666666666667, -63.58333333333334 }    /* Cir */
,   {  5, 16.58333333333333, 16.75000000000000, -63.58333333333334 }    /* Ara */
,   { 64,  6.00000000000000,  6.83333333333333, -64.00000000000000 }    /* Pic */
,   { 12,  6.83333333333333,  9.03333333333333, -64.00000000000000 }    /* Car */
,   { 14, 11.25000000000000, 11.83333333333333, -64.00000000000000 }    /* Cen */
,   { 27, 11.83333333333333, 12.83333333333333, -64.00000000000000 }    /* Cru */
,   { 14, 12.83333333333333, 14.53333333333333, -64.00000000000000 }    /* Cen */
,   { 18, 13.50000000000000, 13.66666666666667, -65.00000000000000 }    /* Cir */
,   {  5, 16.75000000000000, 16.83333333333333, -65.00000000000000 }    /* Ara */
,   { 40,  2.16666666666667,  3.20000000000000, -67.50000000000000 }    /* Hor */
,   { 69,  3.20000000000000,  4.58333333333333, -67.50000000000000 }    /* Ret */
,   { 18, 14.75000000000000, 14.91666666666667, -67.50000000000000 }    /* Cir */
,   {  5, 16.83333333333333, 17.50000000000000, -67.50000000000000 }    /* Ara */
,   { 60, 17.50000000000000, 18.00000000000000, -67.50000000000000 }    /* Pav */
,   { 81, 22.00000000000000, 23.33333333333333, -67.50000000000000 }    /* Tuc */
,   { 32,  4.58333333333333,  6.58333333333333, -70.00000000000000 }    /* Dor */
,   { 18, 13.66666666666667, 14.75000000000000, -70.00000000000000 }    /* Cir */
,   { 79, 14.75000000000000, 17.00000000000000, -70.00000000000000 }    /* TrA */
,   { 81,  0.00000000000000,  1.33333333333333, -75.00000000000000 }    /* Tuc */
,   { 42,  3.50000000000000,  4.58333333333333, -75.00000000000000 }    /* Hyi */
,   { 86,  6.58333333333333,  9.03333333333333, -75.00000000000000 }    /* Vol */
,   { 12,  9.03333333333333, 11.25000000000000, -75.00000000000000 }    /* Car */
,   { 55, 11.25000000000000, 13.66666666666667, -75.00000000000000 }    /* Mus */
,   { 60, 18.00000000000000, 21.33333333333333, -75.00000000000000 }    /* Pav */
,   { 43, 21.33333333333333, 23.33333333333333, -75.00000000000000 }    /* Ind */
,   { 81, 23.33333333333333, 24.00000000000000, -75.00000000000000 }    /* Tuc */
,   { 81,  0.75000000000000,  1.33333333333333, -76.00000000000000 }    /* Tuc */
,   { 42,  0.00000000000000,  3.50000000000000, -82.50000000000000 }    /* Hyi */
,   { 17,  7.66666666666667, 13.66666666666667, -82.50000000000000 }    /* Cha */
,   {  2, 13.66666666666667, 18.00000000000000, -82.50000000000000 }    /* Aps */
,   { 52,  3.50000000000000,  7.66666666666667, -85.00000000000000 }    /* Men */
,   { 57,  0.00000000000000, 24.00000000000000, -90.00000000000000 }    /* Oct */
};

#define NUM_CONSTEL_BOUNDARIES  357



/**
 * @brief
 *      Determines the constellation that contains the given point in the sky.
 *
 * Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
 * constellation that contains that point.
 *
 * @param ra
 *      The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
 *
 * @param dec
 *      The declination (DEC) of a point in the sky, using the J2000 equatorial system.
 *
 * @return
 *      If successful, `status` holds `ASTRO_SUCCESS`,
 *      `symbol` holds a pointer to a 3-character string like "Ori", and
 *      `name` holds a pointer to the full constellation name like "Orion".
 */
astro_constellation_t Astronomy_Constellation(double ra, double dec)
{
    static astro_time_t epoch2000;
    static astro_rotation_t rot = { ASTRO_NOT_INITIALIZED };
    astro_constellation_t constel;
    astro_equatorial_t j2000, b1875;
    astro_vector_t vec2000, vec1875;
    int i, c;

    if (dec < -90.0 || dec > +90.0)
        return ConstelErr(ASTRO_INVALID_PARAMETER);

    /* Allow right ascension to "wrap around". Clamp to [0, 24) sidereal hours. */
    ra = fmod(ra, 24.0);
    if (ra < 0.0)
        ra += 24.0;

    /* Lazy-initialize the rotation matrix for converting J2000 to B1875. */
    if (rot.status != ASTRO_SUCCESS)
    {
        /*
            Need to calculate the B1875 epoch. Based on this:
            https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years
            B = 1900 + (JD - 2415020.31352) / 365.242198781
            I'm interested in using TT instead of JD, giving:
            B = 1900 + ((TT+2451545) - 2415020.31352) / 365.242198781
            B = 1900 + (TT + 36524.68648) / 365.242198781
            TT = 365.242198781*(B - 1900) - 36524.68648 = -45655.741449525
            But Astronomy_TimeFromDays() wants UT, not TT.
            Near that date, I get a historical correction of ut-tt = 3.2 seconds.
            That gives UT = -45655.74141261017 for the B1875 epoch,
            or 1874-12-31T18:12:21.950Z.
        */
        astro_time_t time = Astronomy_TimeFromDays(-45655.74141261017);
        rot = Astronomy_Rotation_EQJ_EQD(time);
        if (rot.status != ASTRO_SUCCESS)
            return ConstelErr(rot.status);

        epoch2000 = Astronomy_TimeFromDays(0.0);
    }

    /* Convert coordinates from J2000 to year 1875. */
    j2000.status = ASTRO_SUCCESS;
    j2000.ra = ra;
    j2000.dec = dec;
    j2000.dist = 1.0;
    vec2000 = Astronomy_VectorFromEquator(j2000, epoch2000);
    if (vec2000.status != ASTRO_SUCCESS)
        return ConstelErr(vec2000.status);

    vec1875 = Astronomy_RotateVector(rot, vec2000);
    if (vec1875.status != ASTRO_SUCCESS)
        return ConstelErr(vec1875.status);

    b1875 = Astronomy_EquatorFromVector(vec1875);
    if (b1875.status != ASTRO_SUCCESS)
        return ConstelErr(b1875.status);

    /* Search for the constellation using the B1875 coordinates. */
    c = -1;     /* constellation not (yet) found */
    for (i=0; i < NUM_CONSTEL_BOUNDARIES; ++i)
    {
        const constel_boundary_t *b = &ConstelBounds[i];
        if ((b->dec_lo <= b1875.dec) && (b->ra_hi > b1875.ra) && (b->ra_lo <= b1875.ra))
        {
            c = b->index;
            break;
        }
    }

    if (c < 0 || c >= NUM_CONSTELLATIONS)
        return ConstelErr(ASTRO_INTERNAL_ERROR);    /* should have been able to find the constellation */

    constel.status = ASTRO_SUCCESS;
    constel.symbol = ConstelInfo[c].symbol;
    constel.name = ConstelInfo[c].name;
    constel.ra_1875 = b1875.ra;
    constel.dec_1875 = b1875.dec;
    return constel;
}


static astro_lunar_eclipse_t LunarEclipseError(astro_status_t status)
{
    astro_lunar_eclipse_t eclipse;
    eclipse.status = status;
    eclipse.kind = ECLIPSE_NONE;
    eclipse.peak = TimeError();
    eclipse.sd_penum = eclipse.sd_partial = eclipse.sd_total = NAN;
    return eclipse;
}


/** @cond DOXYGEN_SKIP */
typedef struct
{
    astro_status_t status;
    astro_time_t time;
    double  u;              /* dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is */
    double  r;              /* km distance between center of Moon/Earth (shaded body) and the line passing through the centers of the Sun and Earth/Moon (casting body). */
    double  k;              /* umbra radius in km, at the shadow plane */
    double  p;              /* penumbra radius in km, at the shadow plane */
    astro_vector_t target;  /* coordinates of target body relative to shadow-casting body at 'time' */
    astro_vector_t dir;     /* heliocentric coordinates of shadow-casting body at 'time' */
}
shadow_t;               /* Represents alignment of the Moon/Earth with the Earth's/Moon's shadow, for finding eclipses. */

typedef struct
{
    double radius_limit;
    double direction;
}
shadow_context_t;
/** @endcond */


static shadow_t ShadowError(astro_status_t status)
{
    shadow_t shadow;
    memset(&shadow, 0, sizeof(shadow));
    shadow.status = status;
    return shadow;
}


static shadow_t CalcShadow(
    double body_radius_km,
    astro_time_t time,
    astro_vector_t target,
    astro_vector_t dir)
{
    double dx, dy, dz;
    shadow_t shadow;

    shadow.target = target;
    shadow.dir = dir;

    shadow.u = (dir.x*target.x + dir.y*target.y + dir.z*target.z) / (dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);

    dx = (shadow.u * dir.x) - target.x;
    dy = (shadow.u * dir.y) - target.y;
    dz = (shadow.u * dir.z) - target.z;
    shadow.r = KM_PER_AU * sqrt(dx*dx + dy*dy + dz*dz);

    shadow.k = +SUN_RADIUS_KM - (1.0 + shadow.u)*(SUN_RADIUS_KM - body_radius_km);
    shadow.p = -SUN_RADIUS_KM + (1.0 + shadow.u)*(SUN_RADIUS_KM + body_radius_km);
    shadow.status = ASTRO_SUCCESS;
    shadow.time = time;

    return shadow;
}


static shadow_t PlanetShadow(astro_body_t body, double planet_radius_km, astro_time_t time)
{
    astro_vector_t e, p, g;

    /* Calculate light-travel-corrected vector from Earth to planet. */
    g = Astronomy_GeoVector(body, time, NO_ABERRATION);
    if (g.status != ASTRO_SUCCESS)
        return ShadowError(g.status);

    /* Calculate light-travel-corrected vector from Earth to Sun. */
    e = Astronomy_GeoVector(BODY_SUN, time, NO_ABERRATION);
    if (e.status != ASTRO_SUCCESS)
        return ShadowError(e.status);

    /* Deduce light-travel-corrected vector from Sun to planet. */
    p.t = time;
    p.x = g.x - e.x;
    p.y = g.y - e.y;
    p.z = g.z - e.z;

    /* Calcluate Earth's position from the planet's point of view. */
    e.x = -g.x;
    e.y = -g.y;
    e.z = -g.z;

    return CalcShadow(planet_radius_km, time, e, p);
}


static shadow_t EarthShadow(astro_time_t time)
{
    /* This function helps find when the Earth's shadow falls upon the Moon. */
    astro_vector_t e, m;

    e = CalcEarth(time);            /* This function never fails; no need to check return value */
    m = Astronomy_GeoMoon(time);    /* This function never fails; no need to check return value */

    return CalcShadow(EARTH_ECLIPSE_RADIUS_KM, time, m, e);
}


static shadow_t MoonShadow(astro_time_t time)
{
    /* This function helps find when the Moon's shadow falls upon the Earth. */

    astro_vector_t h, e, m;

    /*
        This is a variation on the logic in EarthShadow().
        Instead of a heliocentric Earth and a geocentric Moon,
        we want a heliocentric Moon and a lunacentric Earth.
    */

    h = CalcEarth(time);            /* heliocentric Earth */
    m = Astronomy_GeoMoon(time);    /* geocentric Moon */

    /* Calculate lunacentric Earth. */
    e.status = m.status;
    e.x = -m.x;
    e.y = -m.y;
    e.z = -m.z;
    e.t = m.t;

    /* Convert geocentric moon to heliocentric Moon. */
    m.x += h.x;
    m.y += h.y;
    m.z += h.z;

    return CalcShadow(MOON_MEAN_RADIUS_KM, time, e, m);
}


/** @cond DOXYGEN_SKIP */
typedef shadow_t (* shadow_func_t) (astro_time_t time);
/** @endcond */


static astro_func_result_t shadow_distance_slope(void *context, astro_time_t time)
{
    const double dt = 1.0 / 86400.0;
    astro_time_t t1, t2;
    astro_func_result_t result;
    shadow_t shadow1, shadow2;
    shadow_func_t shadowfunc = context;

    t1 = Astronomy_AddDays(time, -dt);
    t2 = Astronomy_AddDays(time, +dt);

    shadow1 = shadowfunc(t1);
    if (shadow1.status != ASTRO_SUCCESS)
        return FuncError(shadow1.status);

    shadow2 = shadowfunc(t2);
    if (shadow2.status != ASTRO_SUCCESS)
        return FuncError(shadow2.status);

    result.value = (shadow2.r - shadow1.r) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}


static shadow_t PeakEarthShadow(astro_time_t search_center_time)
{
    /* Search for when the Earth's shadow axis is closest to the center of the Moon. */

    astro_time_t t1, t2;
    astro_search_result_t result;
    const double window = 0.03;        /* days before/after full moon to search for minimum shadow distance */

    t1 = Astronomy_AddDays(search_center_time, -window);
    t2 = Astronomy_AddDays(search_center_time, +window);

    result = Astronomy_Search(shadow_distance_slope, EarthShadow, t1, t2, 1.0);
    if (result.status != ASTRO_SUCCESS)
        return ShadowError(result.status);

    return EarthShadow(result.time);
}


static shadow_t PeakMoonShadow(astro_time_t search_center_time)
{
    /* Search for when the Moon's shadow axis is closest to the center of the Earth. */

    astro_time_t t1, t2;
    astro_search_result_t result;
    const double window = 0.03;     /* days before/after new moon to search for minimum shadow distance */

    t1 = Astronomy_AddDays(search_center_time, -window);
    t2 = Astronomy_AddDays(search_center_time, +window);

    result = Astronomy_Search(shadow_distance_slope, MoonShadow, t1, t2, 1.0);
    if (result.status != ASTRO_SUCCESS)
        return ShadowError(result.status);

    return MoonShadow(result.time);
}


/** @cond DOXYGEN_SKIP */
typedef struct
{
    astro_body_t    body;
    double          planet_radius_km;
    double          direction;          /* used for transit start/finish search only */
}
planet_shadow_context_t;
/** @endcond */


static astro_func_result_t planet_shadow_distance_slope(void *context, astro_time_t time)
{
    const double dt = 1.0 / 86400.0;
    astro_time_t t1, t2;
    astro_func_result_t result;
    shadow_t shadow1, shadow2;
    const planet_shadow_context_t *p = context;

    t1 = Astronomy_AddDays(time, -dt);
    t2 = Astronomy_AddDays(time, +dt);

    shadow1 = PlanetShadow(p->body, p->planet_radius_km, t1);
    if (shadow1.status != ASTRO_SUCCESS)
        return FuncError(shadow1.status);

    shadow2 = PlanetShadow(p->body, p->planet_radius_km, t2);
    if (shadow2.status != ASTRO_SUCCESS)
        return FuncError(shadow2.status);

    result.value = (shadow2.r - shadow1.r) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}


static shadow_t PeakPlanetShadow(astro_body_t body, double planet_radius_km, astro_time_t search_center_time)
{
    /* Search for when the body's shadow is closest to the center of the Earth. */

    astro_time_t t1, t2;
    astro_search_result_t result;
    planet_shadow_context_t context;
    const double window = 1.0;     /* days before/after inferior conjunction to search for minimum shadow distance */

    t1 = Astronomy_AddDays(search_center_time, -window);
    t2 = Astronomy_AddDays(search_center_time, +window);

    context.body = body;
    context.planet_radius_km = planet_radius_km;
    context.direction = 0.0;    /* not used in this search */

    result = Astronomy_Search(planet_shadow_distance_slope, &context, t1, t2, 1.0);
    if (result.status != ASTRO_SUCCESS)
        return ShadowError(result.status);

    return PlanetShadow(body, planet_radius_km, result.time);
}


static astro_func_result_t shadow_distance(void *context, astro_time_t time)
{
    astro_func_result_t result;
    const shadow_context_t *p = context;
    shadow_t shadow = EarthShadow(time);
    if (shadow.status != ASTRO_SUCCESS)
        return FuncError(shadow.status);

    result.value = p->direction * (shadow.r - p->radius_limit);
    result.status = ASTRO_SUCCESS;
    return result;
}


static double ShadowSemiDurationMinutes(astro_time_t center_time, double radius_limit, double window_minutes)
{
    /* Search backwards and forwards from the center time until shadow axis distance crosses radius limit. */
    double window = window_minutes / (24.0 * 60.0);
    shadow_context_t context;
    astro_search_result_t s1, s2;
    astro_time_t before, after;

    before = Astronomy_AddDays(center_time, -window);
    after  = Astronomy_AddDays(center_time, +window);

    context.radius_limit = radius_limit;
    context.direction = -1.0;
    s1 = Astronomy_Search(shadow_distance, &context, before, center_time, 1.0);

    context.direction = +1.0;
    s2 = Astronomy_Search(shadow_distance, &context, center_time, after, 1.0);

    if (s1.status != ASTRO_SUCCESS || s2.status != ASTRO_SUCCESS)
        return -1.0;    /* something went wrong! */

    return (s2.time.ut - s1.time.ut) * ((24.0 * 60.0) / 2.0);       /* convert days to minutes and average the semi-durations. */
}


/**
 * @brief Searches for a lunar eclipse.
 *
 * This function finds the first lunar eclipse that occurs after `startTime`.
 * A lunar eclipse may be penumbral, partial, or total.
 * See #astro_lunar_eclipse_t for more information.
 * To find a series of lunar eclipses, call this function once,
 * then keep calling #Astronomy_NextLunarEclipse as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * @param startTime
 *      The date and time for starting the search for a lunar eclipse.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields will be valid.
 *      Any other value indicates an error.
 */
astro_lunar_eclipse_t Astronomy_SearchLunarEclipse(astro_time_t startTime)
{
    const double PruneLatitude = 1.8;   /* full Moon's ecliptic latitude above which eclipse is impossible */
    astro_time_t fmtime;
    astro_lunar_eclipse_t eclipse;
    astro_search_result_t fullmoon;
    shadow_t shadow;
    int fmcount;
    double eclip_lat, eclip_lon, distance;

    /* Iterate through consecutive full moons until we find any kind of lunar eclipse. */
    fmtime = startTime;
    for (fmcount=0; fmcount < 12; ++fmcount)
    {
        /* Search for the next full moon. Any eclipse will be near it. */
        fullmoon = Astronomy_SearchMoonPhase(180.0, fmtime, 40.0);
        if (fullmoon.status != ASTRO_SUCCESS)
            return LunarEclipseError(fullmoon.status);

        /* Pruning: if the full Moon's ecliptic latitude is too large, a lunar eclipse is not possible. */
        CalcMoon(fullmoon.time.tt / 36525.0, &eclip_lon, &eclip_lat, &distance);
        if (RAD2DEG * fabs(eclip_lat) < PruneLatitude)
        {
            /* Search near the full moon for the time when the center of the Moon */
            /* is closest to the line passing through the centers of the Sun and Earth. */
            shadow = PeakEarthShadow(fullmoon.time);
            if (shadow.status != ASTRO_SUCCESS)
                return LunarEclipseError(shadow.status);

            if (shadow.r < shadow.p + MOON_MEAN_RADIUS_KM)
            {
                /* This is at least a penumbral eclipse. We will return a result. */
                eclipse.status = ASTRO_SUCCESS;
                eclipse.kind = ECLIPSE_PENUMBRAL;
                eclipse.peak = shadow.time;
                eclipse.sd_total = 0.0;
                eclipse.sd_partial = 0.0;
                eclipse.sd_penum = ShadowSemiDurationMinutes(shadow.time, shadow.p + MOON_MEAN_RADIUS_KM, 200.0);
                if (eclipse.sd_penum <= 0.0)
                    return LunarEclipseError(ASTRO_SEARCH_FAILURE);

                if (shadow.r < shadow.k + MOON_MEAN_RADIUS_KM)
                {
                    /* This is at least a partial eclipse. */
                    eclipse.kind = ECLIPSE_PARTIAL;
                    eclipse.sd_partial = ShadowSemiDurationMinutes(shadow.time, shadow.k + MOON_MEAN_RADIUS_KM, eclipse.sd_penum);
                    if (eclipse.sd_partial <= 0.0)
                        return LunarEclipseError(ASTRO_SEARCH_FAILURE);

                    if (shadow.r + MOON_MEAN_RADIUS_KM < shadow.k)
                    {
                        /* This is a total eclipse. */
                        eclipse.kind = ECLIPSE_TOTAL;
                        eclipse.sd_total = ShadowSemiDurationMinutes(shadow.time, shadow.k - MOON_MEAN_RADIUS_KM, eclipse.sd_partial);
                        if (eclipse.sd_total <= 0.0)
                            return LunarEclipseError(ASTRO_SEARCH_FAILURE);
                    }
                }
                return eclipse;
            }
        }

        /* We didn't find an eclipse on this full moon, so search for the next one. */
        fmtime = Astronomy_AddDays(fullmoon.time, 10.0);
    }

    /* Safety valve to prevent infinite loop. */
    /* This should never happen, because at least 2 lunar eclipses happen per year. */
    return LunarEclipseError(ASTRO_INTERNAL_ERROR);
}

/**
 * @brief Searches for the next lunar eclipse in a series.
 *
 * After using #Astronomy_SearchLunarEclipse to find the first lunar eclipse
 * in a series, you can call this function to find the next consecutive lunar eclipse.
 * Pass in the `peak` value from the #astro_lunar_eclipse_t returned by the
 * previous call to `Astronomy_SearchLunarEclipse` or `Astronomy_NextLunarEclipse`
 * to find the next lunar eclipse.
 *
 * @param prevEclipseTime
 *      A date and time near a full moon. Lunar eclipse search will start at the next full moon.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields will be valid.
 *      Any other value indicates an error.
 */
astro_lunar_eclipse_t Astronomy_NextLunarEclipse(astro_time_t prevEclipseTime)
{
    astro_time_t startTime = Astronomy_AddDays(prevEclipseTime, 10.0);
    return Astronomy_SearchLunarEclipse(startTime);
}


static astro_global_solar_eclipse_t GlobalSolarEclipseError(astro_status_t status)
{
    astro_global_solar_eclipse_t eclipse;

    eclipse.status = status;
    eclipse.kind = ECLIPSE_NONE;
    eclipse.peak = TimeError();
    eclipse.distance = eclipse.latitude = eclipse.longitude = NAN;

    return eclipse;
}

/* The umbra radius tells us what kind of eclipse the observer sees. */
/* If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular. */
/* HACK: I added a tiny bias (14 meters) to match Espenak test data. */
#define EclipseKindFromUmbra(k)     (((k) > 0.014) ? ECLIPSE_TOTAL : ECLIPSE_ANNULAR)

static astro_global_solar_eclipse_t GeoidIntersect(shadow_t shadow)
{
    astro_global_solar_eclipse_t eclipse;
    astro_rotation_t rot, inv;
    astro_vector_t v, e, o;
    shadow_t surface;
    double A, B, C, radic, u, R;
    double px, py, pz, proj;
    double gast;

    eclipse.status = ASTRO_SUCCESS;
    eclipse.kind = ECLIPSE_PARTIAL;
    eclipse.peak = shadow.time;
    eclipse.distance = shadow.r;
    eclipse.latitude = eclipse.longitude = NAN;

    /*
        We want to calculate the intersection of the shadow axis with the Earth's geoid.
        First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
        coordinates that are perfectly aligned with the Earth's equator at this
        moment in time.
    */
    rot = Astronomy_Rotation_EQJ_EQD(shadow.time);
    if (rot.status != ASTRO_SUCCESS)
        return GlobalSolarEclipseError(rot.status);

    v = Astronomy_RotateVector(rot, shadow.dir);        /* shadow-axis vector in equator-of-date coordinates */
    if (v.status != ASTRO_SUCCESS)
        return GlobalSolarEclipseError(v.status);

    e = Astronomy_RotateVector(rot, shadow.target);     /* lunacentric Earth in equator-of-date coordinates */
    if (e.status != ASTRO_SUCCESS)
        return GlobalSolarEclipseError(e.status);

    /*
        Convert all distances from AU to km.
        But dilate the z-coordinates so that the Earth becomes a perfect sphere.
        Then find the intersection of the vector with the sphere.
        See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
    */
    v.x *= KM_PER_AU;
    v.y *= KM_PER_AU;
    v.z *= KM_PER_AU / EARTH_FLATTENING;

    e.x *= KM_PER_AU;
    e.y *= KM_PER_AU;
    e.z *= KM_PER_AU / EARTH_FLATTENING;

    /*
        Solve the quadratic equation that finds whether and where
        the shadow axis intersects with the Earth in the dilated coordinate system.
    */
    R = EARTH_EQUATORIAL_RADIUS_KM;
    A = v.x*v.x + v.y*v.y + v.z*v.z;
    B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z);
    C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R;
    radic = B*B - 4*A*C;

    if (radic > 0.0)
    {
        /* Calculate the closer of the two intersection points. */
        /* This will be on the day side of the Earth. */
        u = (-B - sqrt(radic)) / (2 * A);

        /* Convert lunacentric dilated coordinates to geocentric coordinates. */
        px = u*v.x - e.x;
        py = u*v.y - e.y;
        pz = (u*v.z - e.z) * EARTH_FLATTENING;

        /* Convert cartesian coordinates into geodetic latitude/longitude. */
        proj = sqrt(px*px + py*py) * (EARTH_FLATTENING * EARTH_FLATTENING);
        if (proj == 0.0)
            eclipse.latitude = (pz > 0.0) ? +90.0 : -90.0;
        else
            eclipse.latitude = RAD2DEG * atan(pz / proj);

        /* Adjust longitude for Earth's rotation at the given UT. */
        gast = sidereal_time(&eclipse.peak);
        eclipse.longitude = fmod((RAD2DEG*atan2(py, px)) - (15*gast), 360.0);
        if (eclipse.longitude <= -180.0)
            eclipse.longitude += 360.0;
        else if (eclipse.longitude > +180.0)
            eclipse.longitude -= 360.0;

        /* We want to determine whether the observer sees a total eclipse or an annular eclipse. */
        /* We need to perform a series of vector calculations... */
        /* Calculate the inverse rotation matrix, so we can convert EQD to EQJ. */
        inv = Astronomy_InverseRotation(rot);
        if (inv.status != ASTRO_SUCCESS)
            return GlobalSolarEclipseError(inv.status);

        /* Put the EQD geocentric coordinates of the observer into the vector 'o'. */
        /* Also convert back from kilometers to astronomical units. */
        o.status = ASTRO_SUCCESS;
        o.t = shadow.time;
        o.x = px / KM_PER_AU;
        o.y = py / KM_PER_AU;
        o.z = pz / KM_PER_AU;

        /* Rotate the observer's geocentric EQD back to the EQJ system. */
        o = Astronomy_RotateVector(inv, o);

        /* Convert geocentric vector to lunacentric vector. */
        o.x += shadow.target.x;
        o.y += shadow.target.y;
        o.z += shadow.target.z;

        /* Recalculate the shadow using a vector from the Moon's center toward the observer. */
        surface = CalcShadow(MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir);

        /* If we did everything right, the shadow distance should be very close to zero. */
        /* That's because we already determined the observer 'o' is on the shadow axis! */
        if (surface.r > 1.0e-9 || surface.r < 0.0)
            return GlobalSolarEclipseError(ASTRO_INTERNAL_ERROR);

        eclipse.kind = EclipseKindFromUmbra(surface.k);
    }

    return eclipse;
}


/**
 * @brief Searches for a solar eclipse visible anywhere on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See #astro_global_solar_eclipse_t for more information.
 * To find a series of solar eclipses, call this function once,
 * then keep calling #Astronomy_NextGlobalSolarEclipse as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * @param startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields are as described in #astro_global_solar_eclipse_t.
 *      Any other value indicates an error.
 */
astro_global_solar_eclipse_t Astronomy_SearchGlobalSolarEclipse(astro_time_t startTime)
{
    const double PruneLatitude = 1.8;   /* Moon's ecliptic latitude beyond which eclipse is impossible */
    astro_time_t nmtime;
    astro_search_result_t newmoon;
    shadow_t shadow;
    int nmcount;
    double eclip_lat, eclip_lon, distance;

    /* Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth. */
    nmtime = startTime;
    for (nmcount=0; nmcount < 12; ++nmcount)
    {
        /* Search for the next new moon. Any eclipse will be near it. */
        newmoon = Astronomy_SearchMoonPhase(0.0, nmtime, 40.0);
        if (newmoon.status != ASTRO_SUCCESS)
            return GlobalSolarEclipseError(newmoon.status);

        /* Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible. */
        CalcMoon(newmoon.time.tt / 36525.0, &eclip_lon, &eclip_lat, &distance);
        if (RAD2DEG * fabs(eclip_lat) < PruneLatitude)
        {
            /* Search near the new moon for the time when the center of the Earth */
            /* is closest to the line passing through the centers of the Sun and Moon. */
            shadow = PeakMoonShadow(newmoon.time);
            if (shadow.status != ASTRO_SUCCESS)
                return GlobalSolarEclipseError(shadow.status);

            if (shadow.r < shadow.p + EARTH_MEAN_RADIUS_KM)
            {
                /* This is at least a partial solar eclipse visible somewhere on Earth. */
                /* Try to find an intersection between the shadow axis and the Earth's oblate geoid. */
                return GeoidIntersect(shadow);
            }
        }

        /* We didn't find an eclipse on this new moon, so search for the next one. */
        nmtime = Astronomy_AddDays(newmoon.time, 10.0);
    }

    /* Safety valve to prevent infinite loop. */
    /* This should never happen, because at least 2 solar eclipses happen per year. */
    return GlobalSolarEclipseError(ASTRO_INTERNAL_ERROR);
}


/**
 * @brief Searches for the next global solar eclipse in a series.
 *
 * After using #Astronomy_SearchGlobalSolarEclipse to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the #astro_global_solar_eclipse_t returned by the
 * previous call to `Astronomy_SearchGlobalSolarEclipse` or `Astronomy_NextGlobalSolarEclipse`
 * to find the next solar eclipse.
 *
 * @param prevEclipseTime
 *      A date and time near a new moon. Solar eclipse search will start at the next new moon.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields are as described in #astro_global_solar_eclipse_t.
 *      Any other value indicates an error.
 */
astro_global_solar_eclipse_t Astronomy_NextGlobalSolarEclipse(astro_time_t prevEclipseTime)
{
    astro_time_t startTime = Astronomy_AddDays(prevEclipseTime, 10.0);
    return Astronomy_SearchGlobalSolarEclipse(startTime);
}


static astro_eclipse_event_t EclipseEventError(void)
{
    astro_eclipse_event_t evt;
    evt.time = TimeError();
    evt.altitude = NAN;
    return evt;
}


static astro_local_solar_eclipse_t LocalSolarEclipseError(astro_status_t status)
{
    astro_local_solar_eclipse_t eclipse;

    eclipse.status = status;
    eclipse.kind = ECLIPSE_NONE;

    eclipse.partial_begin = EclipseEventError();
    eclipse.total_begin   = EclipseEventError();
    eclipse.peak          = EclipseEventError();
    eclipse.total_end     = EclipseEventError();
    eclipse.partial_end   = EclipseEventError();

    return eclipse;
}


static shadow_t LocalMoonShadow(astro_time_t time, astro_observer_t observer)
{
    astro_vector_t h, o, m;
    double pos[3];

    /* Calculate observer's geocentric position. */
    /* For efficiency, do this first, to populate the earth rotation parameters in 'time'. */
    /* That way they can be recycled instead of recalculated. */
    geo_pos(&time, observer, pos);

    h = CalcEarth(time);            /* heliocentric Earth */
    m = Astronomy_GeoMoon(time);    /* geocentric Moon */

    /* Calculate lunacentric location of an observer on the Earth's surface. */
    o.status = m.status;
    o.x = pos[0] - m.x;
    o.y = pos[1] - m.y;
    o.z = pos[2] - m.z;
    o.t = m.t;

    /* Convert geocentric moon to heliocentric Moon. */
    m.x += h.x;
    m.y += h.y;
    m.z += h.z;

    return CalcShadow(MOON_MEAN_RADIUS_KM, time, o, m);
}


static astro_func_result_t local_shadow_distance_slope(void *context, astro_time_t time)
{
    const double dt = 1.0 / 86400.0;
    astro_time_t t1, t2;
    astro_func_result_t result;
    shadow_t shadow1, shadow2;
    const astro_observer_t *observer = context;

    t1 = Astronomy_AddDays(time, -dt);
    t2 = Astronomy_AddDays(time, +dt);

    shadow1 = LocalMoonShadow(t1, *observer);
    if (shadow1.status != ASTRO_SUCCESS)
        return FuncError(shadow1.status);

    shadow2 = LocalMoonShadow(t2, *observer);
    if (shadow2.status != ASTRO_SUCCESS)
        return FuncError(shadow2.status);

    result.value = (shadow2.r - shadow1.r) / dt;
    result.status = ASTRO_SUCCESS;
    return result;
}


static shadow_t PeakLocalMoonShadow(astro_time_t search_center_time, astro_observer_t observer)
{
    astro_time_t t1, t2;
    astro_search_result_t result;
    const double window = 0.2;

    /*
        Search for the time near search_center_time that the Moon's shadow comes
        closest to the given observer.
    */

    t1 = Astronomy_AddDays(search_center_time, -window);
    t2 = Astronomy_AddDays(search_center_time, +window);

    result = Astronomy_Search(local_shadow_distance_slope, &observer, t1, t2, 1.0);
    if (result.status != ASTRO_SUCCESS)
        return ShadowError(result.status);

    return LocalMoonShadow(result.time, observer);
}


static double local_partial_distance(const shadow_t *shadow)
{
    return shadow->p - shadow->r;
}

static double local_total_distance(const shadow_t *shadow)
{
    /* Must take the absolute value of the umbra radius 'k' */
    /* because it can be negative for an annular eclipse. */
    return fabs(shadow->k) - shadow->r;
}

/** @cond DOXYGEN_SKIP */
typedef double (* local_distance_func) (const shadow_t *shadow);

typedef struct
{
    local_distance_func     func;
    double                  direction;
    astro_observer_t        observer;
}
eclipse_transition_t;
/* @endcond */


static astro_func_result_t local_eclipse_func(void *context, astro_time_t time)
{
    const eclipse_transition_t *trans = context;
    shadow_t shadow;
    astro_func_result_t result;

    shadow = LocalMoonShadow(time, trans->observer);
    if (shadow.status != ASTRO_SUCCESS)
        return FuncError(shadow.status);

    result.value = trans->direction * trans->func(&shadow);
    result.status = ASTRO_SUCCESS;
    return result;
}


astro_func_result_t SunAltitude(
    astro_time_t time,
    astro_observer_t observer)
{
    astro_equatorial_t equ;
    astro_horizon_t hor;
    astro_func_result_t result;

    equ = Astronomy_Equator(BODY_SUN, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (equ.status != ASTRO_SUCCESS)
        return FuncError(equ.status);

    hor = Astronomy_Horizon(&time, observer, equ.ra, equ.dec, REFRACTION_NORMAL);
    result.value = hor.altitude;
    result.status = ASTRO_SUCCESS;
    return result;
}


static astro_status_t CalcEvent(
    astro_observer_t observer,
    astro_time_t time,
    astro_eclipse_event_t *evt)
{
    astro_func_result_t result;

    result = SunAltitude(time, observer);
    if (result.status != ASTRO_SUCCESS)
    {
        evt->time = TimeError();
        evt->altitude = NAN;
        return result.status;
    }

    evt->time = time;
    evt->altitude = result.value;
    return ASTRO_SUCCESS;
}


static astro_status_t LocalEclipseTransition(
    astro_observer_t observer,
    double direction,
    local_distance_func func,
    astro_time_t t1,
    astro_time_t t2,
    astro_eclipse_event_t *evt)
{
    eclipse_transition_t trans;
    astro_search_result_t search;

    trans.func = func;
    trans.direction = direction;
    trans.observer = observer;

    search = Astronomy_Search(local_eclipse_func, &trans, t1, t2, 1.0);
    if (search.status != ASTRO_SUCCESS)
    {
        evt->time = TimeError();
        evt->altitude = NAN;
        return search.status;
    }

    return CalcEvent(observer, search.time, evt);
}


static astro_local_solar_eclipse_t LocalEclipse(
    shadow_t shadow,
    astro_observer_t observer)
{
    const double PARTIAL_WINDOW = 0.2;
    const double TOTAL_WINDOW = 0.01;
    astro_local_solar_eclipse_t eclipse;
    astro_time_t t1, t2;
    astro_status_t status;

    status = CalcEvent(observer, shadow.time, &eclipse.peak);
    if (status != ASTRO_SUCCESS)
        return LocalSolarEclipseError(status);

    t1 = Astronomy_AddDays(shadow.time, -PARTIAL_WINDOW);
    t2 = Astronomy_AddDays(shadow.time, +PARTIAL_WINDOW);

    status = LocalEclipseTransition(observer, +1.0, local_partial_distance, t1, shadow.time, &eclipse.partial_begin);
    if (status != ASTRO_SUCCESS)
        return LocalSolarEclipseError(status);

    status = LocalEclipseTransition(observer, -1.0, local_partial_distance, shadow.time, t2, &eclipse.partial_end);
    if (status != ASTRO_SUCCESS)
        return LocalSolarEclipseError(status);

    if (shadow.r < fabs(shadow.k))      /* take absolute value of 'k' to handle annular eclipses too. */
    {
        t1 = Astronomy_AddDays(shadow.time, -TOTAL_WINDOW);
        t2 = Astronomy_AddDays(shadow.time, +TOTAL_WINDOW);

        status = LocalEclipseTransition(observer, +1.0, local_total_distance, t1, shadow.time, &eclipse.total_begin);
        if (status != ASTRO_SUCCESS)
            return LocalSolarEclipseError(status);

        status = LocalEclipseTransition(observer, -1.0, local_total_distance, shadow.time, t2, &eclipse.total_end);
        if (status != ASTRO_SUCCESS)
            return LocalSolarEclipseError(status);

        eclipse.kind = EclipseKindFromUmbra(shadow.k);
    }
    else
    {
        eclipse.total_begin = eclipse.total_end = EclipseEventError();
        eclipse.kind = ECLIPSE_PARTIAL;
    }

    eclipse.status = ASTRO_SUCCESS;
    return eclipse;
}


/**
 * @brief Searches for a solar eclipse visible at a specific location on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See #astro_local_solar_eclipse_t for more information.
 * To find a series of solar eclipses, call this function once,
 * then keep calling #Astronomy_NextLocalSolarEclipse as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * IMPORTANT: An eclipse reported by this function might be partly or
 * completely invisible to the observer due to the time of day.
 *
 * @param startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @param observer
 *      The geographic location of the observer.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields are as described in #astro_local_solar_eclipse_t.
 *      Any other value indicates an error.
 */
astro_local_solar_eclipse_t Astronomy_SearchLocalSolarEclipse(
    astro_time_t startTime,
    astro_observer_t observer)
{
    const double PruneLatitude = 1.8;   /* Moon's ecliptic latitude beyond which eclipse is impossible */
    astro_time_t nmtime;
    astro_search_result_t newmoon;
    shadow_t shadow;
    double eclip_lat, eclip_lon, distance;
    astro_local_solar_eclipse_t eclipse;

    /* Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth. */
    nmtime = startTime;
    for(;;)
    {
        /* Search for the next new moon. Any eclipse will be near it. */
        newmoon = Astronomy_SearchMoonPhase(0.0, nmtime, 40.0);
        if (newmoon.status != ASTRO_SUCCESS)
            return LocalSolarEclipseError(newmoon.status);

        /* Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible. */
        CalcMoon(newmoon.time.tt / 36525.0, &eclip_lon, &eclip_lat, &distance);
        if (RAD2DEG * fabs(eclip_lat) < PruneLatitude)
        {
            /* Search near the new moon for the time when the observer */
            /* is closest to the line passing through the centers of the Sun and Moon. */
            shadow = PeakLocalMoonShadow(newmoon.time, observer);
            if (shadow.status != ASTRO_SUCCESS)
                return LocalSolarEclipseError(shadow.status);

            if (shadow.r < shadow.p)
            {
                /* This is at least a partial solar eclipse for the observer. */
                eclipse = LocalEclipse(shadow, observer);

                /* If any error occurs, something is really wrong and we should bail out. */
                if (eclipse.status != ASTRO_SUCCESS)
                    return eclipse;

                /* Ignore any eclipse that happens completely at night. */
                /* More precisely, the center of the Sun must be above the horizon */
                /* at the beginning or the end of the eclipse, or we skip the event. */
                if (eclipse.partial_begin.altitude > 0.0 || eclipse.partial_end.altitude > 0.0)
                    return eclipse;
            }
        }

        /* We didn't find an eclipse on this new moon, so search for the next one. */
        nmtime = Astronomy_AddDays(newmoon.time, 10.0);
    }
}


/**
 * @brief Searches for the next local solar eclipse in a series.
 *
 * After using #Astronomy_SearchLocalSolarEclipse to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the #astro_local_solar_eclipse_t returned by the
 * previous call to `Astronomy_SearchLocalSolarEclipse` or `Astronomy_NextLocalSolarEclipse`
 * to find the next solar eclipse.
 *
 * @param prevEclipseTime
 *      A date and time near a new moon. Solar eclipse search will start at the next new moon.
 *
 * @param observer
 *      The geographic location of the observer.
 *
 * @return
 *      If successful, the `status` field in the returned structure will contain `ASTRO_SUCCESS`
 *      and the remaining structure fields are as described in #astro_local_solar_eclipse_t.
 *      Any other value indicates an error.
 */
astro_local_solar_eclipse_t Astronomy_NextLocalSolarEclipse(
    astro_time_t prevEclipseTime,
    astro_observer_t observer)
{
    astro_time_t startTime = Astronomy_AddDays(prevEclipseTime, 10.0);
    return Astronomy_SearchLocalSolarEclipse(startTime, observer);
}


static astro_func_result_t planet_transit_bound(void *context, astro_time_t time)
{
    shadow_t shadow;
    astro_func_result_t result;
    const planet_shadow_context_t *p = context;

    shadow = PlanetShadow(p->body, p->planet_radius_km, time);
    if (shadow.status != ASTRO_SUCCESS)
        return FuncError(shadow.status);

    result.status = ASTRO_SUCCESS;
    result.value = p->direction * (shadow.r - shadow.p);
    return result;
}


static astro_search_result_t PlanetTransitBoundary(
    astro_body_t body,
    double planet_radius_km,
    astro_time_t t1,
    astro_time_t t2,
    double direction)
{
    /* Search for the time the planet's penumbra begins/ends making contact with the center of the Earth. */
    planet_shadow_context_t context;

    context.body = body;
    context.planet_radius_km = planet_radius_km;
    context.direction = direction;

    return Astronomy_Search(planet_transit_bound, &context, t1, t2, 1.0);
}


/**
 * @brief Searches for the first transit of Mercury or Venus after a given date.
 *
 * Finds the first transit of Mercury or Venus after a specified date.
 * A transit is when an inferior planet passes between the Sun and the Earth
 * so that the silhouette of the planet is visible against the Sun in the background.
 * To continue the search, pass the `finish` time in the returned structure to
 * #Astronomy_NextTransit.
 *
 * @param body
 *      The planet whose transit is to be found. Must be `BODY_MERCURY` or `BODY_VENUS`.
 *
 * @param startTime
 *      The date and time for starting the search for a transit.
 *
 * @return
 *      If successful, the `status` field in the returned structure hold `ASTRO_SUCCESS`
 *      and the other fields are as documented in #astro_transit_t.
 *      Otherwise, `status` holds an error code and the other structure members are undefined.
 */
astro_transit_t Astronomy_SearchTransit(astro_body_t body, astro_time_t startTime)
{
    astro_time_t search_time;
    astro_transit_t transit;
    astro_search_result_t conj, search;
    astro_angle_result_t conj_separation, min_separation;
    shadow_t shadow;
    double planet_radius_km;
    astro_time_t tx;
    const double threshold_angle = 0.4;     /* maximum angular separation to attempt transit calculation */
    const double dt_days = 1.0;

    /* Validate the planet and find its mean radius. */
    switch (body)
    {
    case BODY_MERCURY:  planet_radius_km = 2439.7;  break;
    case BODY_VENUS:    planet_radius_km = 6051.8;  break;
    default:
        return TransitErr(ASTRO_INVALID_BODY);
    }

    search_time = startTime;
    for(;;)
    {
        /*
            Search for the next inferior conjunction of the given planet.
            This is the next time the Earth and the other planet have the same
            ecliptic longitude as seen from the Sun.
        */
        conj = Astronomy_SearchRelativeLongitude(body, 0.0, search_time);
        if (conj.status != ASTRO_SUCCESS)
            return TransitErr(conj.status);

        /* Calculate the angular separation between the body and the Sun at this time. */
        conj_separation = Astronomy_AngleFromSun(body, conj.time);
        if (conj_separation.status != ASTRO_SUCCESS)
            return TransitErr(conj_separation.status);

        if (conj_separation.angle < threshold_angle)
        {
            /*
                The planet's angular separation from the Sun is small enough
                to consider it a transit candidate.
                Search for the moment when the line passing through the Sun
                and planet are closest to the Earth's center.
            */
            shadow = PeakPlanetShadow(body, planet_radius_km, conj.time);
            if (shadow.status != ASTRO_SUCCESS)
                return TransitErr(shadow.status);

            if (shadow.r < shadow.p)        /* does the planet's penumbra touch the Earth's center? */
            {
                /* Find the beginning and end of the penumbral contact. */
                tx = Astronomy_AddDays(shadow.time, -dt_days);
                search = PlanetTransitBoundary(body, planet_radius_km, tx, shadow.time, -1.0);
                if (search.status != ASTRO_SUCCESS)
                    return TransitErr(search.status);
                transit.start = search.time;

                tx = Astronomy_AddDays(shadow.time, +dt_days);
                search = PlanetTransitBoundary(body, planet_radius_km, shadow.time, tx, +1.0);
                if (search.status != ASTRO_SUCCESS)
                    return TransitErr(search.status);
                transit.finish = search.time;
                transit.status = ASTRO_SUCCESS;
                transit.peak = shadow.time;

                min_separation = Astronomy_AngleFromSun(body, shadow.time);
                if (min_separation.status != ASTRO_SUCCESS)
                    return TransitErr(min_separation.status);

                transit.separation = 60.0 * min_separation.angle;   /* convert degrees to arcminutes */
                return transit;
            }
        }

        /* This inferior conjunction was not a transit. Try the next inferior conjunction. */
        search_time = Astronomy_AddDays(conj.time, 10.0);
    }
}


/**
 * @brief Searches for another transit of Mercury or Venus.
 *
 * After calling #Astronomy_SearchTransit to find a transit of Mercury or Venus,
 * this function finds the next transit after that.
 * Keep calling this function as many times as you want to keep finding more transits.
 *
 * @param body
 *      The planet whose transit is to be found. Must be `BODY_MERCURY` or `BODY_VENUS`.
 *
 * @param prevTransitTime
 *      A date and time near the previous transit.
 *
 * @return
 *      If successful, the `status` field in the returned structure hold `ASTRO_SUCCESS`
 *      and the other fields are as documented in #astro_transit_t.
 *      Otherwise, `status` holds an error code and the other structure members are undefined.
 */
astro_transit_t Astronomy_NextTransit(astro_body_t body, astro_time_t prevTransitTime)
{
    astro_time_t startTime;

    startTime = Astronomy_AddDays(prevTransitTime, 100.0);
    return Astronomy_SearchTransit(body, startTime);
}


#ifdef __cplusplus
}
#endif
