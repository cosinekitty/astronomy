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
 * @return If `name` is one of the strings (case-sensitive) listed above, the returned value is the corresponding #astro_body_t value, otherwise it is `BODY_INVALID`.
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
    transit.start = transit.finish = TimeError();
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
$ASTRO_IAU_DATA()
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
            case 1: ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2: ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3: ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            case 4: ARG=D;  MAX=6; FAC=1.0;                       break;
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
$ASTRO_ADDSOL()
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
 * If given an invalid value for `body`, or the body is `BODY_PLUTO` and the `time` is outside
 * the year range 1700..2200, this function will fail. The caller should always check
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
 *      See remarks above for more details.
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
 *      See remarks above for more details.
 *
 * @param context
 *      Any ancillary data needed by the function `func` to calculate a value.
 *      The data type varies depending on the function passed in.
 *      For example, the function may involve a specific celestial body that
 *      must be specified somehow.
 *
 * @param t1
 *      The lower time bound of the search window.
 *      See remarks above for more details.
 *
 * @param t2
 *      The upper time bound of the search window.
 *      See remarks above for more details.
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
 *      See the remarks above for more details.
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

    sv = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    se = Astronomy_Ecliptic(sv);        /* checks for errors in sv */
    if (se.status != ASTRO_SUCCESS)
        return AngleError(se.status);

    bv = Astronomy_GeoVector(body, time, ABERRATION);
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
 *      On success, the function returns the angle as described above in the `angle`
 *      field and `ASTRO_SUCCESS` in the `status` field.
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
 *      See remarks above for more details.
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

#define NeptuneHelioDistance(time)      VsopHelioDistance(&vsop[BODY_NEPTUNE], (time))

static astro_apsis_t NeptuneExtreme(astro_apsis_kind_t kind, astro_time_t start_time, double dayspan)
{
    astro_apsis_t apsis;
    const double direction = (kind == APSIS_APOCENTER) ? +1.0 : -1.0;
    const int npoints = 10;
    int i, best_i;
    double interval;
    double dist, best_dist;
    astro_time_t time;

    for(;;)
    {
        interval = dayspan / (npoints - 1);

        if (interval < 1.0 / 1440.0)    /* iterate until uncertainty is less than one minute */
        {
            apsis.status = ASTRO_SUCCESS;
            apsis.kind = kind;
            apsis.time = Astronomy_AddDays(start_time, interval / 2.0);
            apsis.dist_au = NeptuneHelioDistance(apsis.time);
            apsis.dist_km = apsis.dist_au * KM_PER_AU;
            return apsis;
        }

        best_i = -1;
        best_dist = 0.0;
        for (i=0; i < npoints; ++i)
        {
            time = Astronomy_AddDays(start_time, i * interval);
            dist = direction * NeptuneHelioDistance(time);
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


static astro_apsis_t SearchNeptuneApsis(astro_time_t startTime)
{
    const int npoints = 100;
    int i;
    astro_time_t t1, t2, time, t_min, t_max;
    double dist, max_dist, min_dist;
    astro_apsis_t perihelion, aphelion;
    double interval;

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
    t1 = Astronomy_AddDays(startTime, NEPTUNE_ORBITAL_PERIOD * ( -30.0 / 360.0));
    t2 = Astronomy_AddDays(startTime, NEPTUNE_ORBITAL_PERIOD * (+270.0 / 360.0));
    t_min = t_max = t1;
    min_dist = max_dist = -1.0;     /* prevent warning about uninitialized variables */
    interval = (t2.ut - t1.ut) / (npoints - 1.0);

    for (i=0; i < npoints; ++i)
    {
        double ut = t1.ut + (i * interval);
        time = Astronomy_TimeFromDays(ut);
        dist = NeptuneHelioDistance(time);
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
    perihelion = NeptuneExtreme(APSIS_PERICENTER, t1, 4 * interval);

    t1 = Astronomy_AddDays(t_max, -2 * interval);
    aphelion = NeptuneExtreme(APSIS_APOCENTER, t1, 4 * interval);

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

    return ApsisError(ASTRO_FAIL_NEPTUNE_APSIS);
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

    if (body == BODY_NEPTUNE)
        return SearchNeptuneApsis(startTime);

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
 *      If successful, `status` holds `ASTRO_SUCCESS` and the other fields are valid as described above.
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

$ASTRO_CONSTEL()

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

    /* To include light travel time compensation, we use GeoVector instead of HelioVector. */

    g = Astronomy_GeoVector(body, time, ABERRATION);
    if (g.status != ASTRO_SUCCESS)
        return ShadowError(g.status);

    /* Calculate heliocentric Earth. */
    e = CalcEarth(time);

    /* Convert light-travel corrected geocentric planet to heliocentric planet. */
    p.t = time;
    p.x = e.x + g.x;
    p.y = e.y + g.y;
    p.z = e.z + g.z;

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
    astro_angle_result_t separation;
    shadow_t shadow;
    double planet_radius_km;
    astro_time_t tx;
    const double threshold_angle = 10.0;     /* maximum angular separation to attempt transit calculation */
    const double dt_days = 10.0;

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
        separation = Astronomy_AngleFromSun(body, conj.time);
        if (separation.status != ASTRO_SUCCESS)
            return TransitErr(separation.status);

        if (separation.angle < threshold_angle)
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
                return transit;
            }
        }

        /* This inferior conjunction was not a transit. Try the next inferior conjunction. */
        search_time = Astronomy_AddDays(conj.time, 10.0);
    }
}


/**
 * @brief Searches for the another transit of Mercury and Venus.
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
