/*
    Astronomy Engine for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2023 Don Cross <cosinekitty@gmail.com>

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

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined(ASTRONOMY_ENGINE_NO_CURRENT_TIME)
#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/time.h>
#elif defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <time.h>
#endif
#endif

#include "astronomy.h"

#ifdef __FAST_MATH__
#error Astronomy Engine does not support "fast math" optimization because it causes incorrect behavior. See: https://github.com/cosinekitty/astronomy/issues/245
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @cond DOXYGEN_SKIP */
#define PI      3.14159265358979323846

//$ASTRO_C_PLUTO_CONST()

typedef enum
{
    FROM_2000,
    INTO_2000
}
precess_dir_t;

typedef struct
{
    double x;
    double y;
    double z;
}
terse_vector_t;

typedef struct
{
    double          tt;   /* J2000 terrestrial time [days] */
    terse_vector_t  r;    /* position [au] */
    terse_vector_t  v;    /* velocity [au/day] */
    terse_vector_t  a;    /* acceleration [au/day^2] */
}
body_grav_calc_t;

typedef struct
{
    body_grav_calc_t   step[PLUTO_NSTEPS];
}
body_segment_t;

typedef struct
{
    double          tt;  /* Terrestrial Time in J2000 days */
    terse_vector_t  r;   /* position [au] */
    terse_vector_t  v;   /* velocity [au/day] */
}
body_state_t;

typedef struct
{
    body_state_t Sun;
    body_state_t Jupiter;
    body_state_t Saturn;
    body_state_t Uranus;
    body_state_t Neptune;
}
major_bodies_t;

typedef struct
{
    astro_time_t      time;
    body_state_t      gravitators[1 + BODY_SUN];
    body_grav_calc_t *bodies;
}
gravsim_endpoint_t;

struct astro_grav_sim_s
{
    astro_body_t        originBody;
    int                 numBodies;
    gravsim_endpoint_t  endpoint[2];
    gravsim_endpoint_t *prev;
    gravsim_endpoint_t *curr;
};

typedef struct
{
    double ra;
    double dec;
    double dist;
}
stardef_t;

/* Mean obliquity of the J2000 ecliptic in radians. */
#define OBLIQ_2000       0.40909260059599012
#define COS_OBLIQ_2000   0.9174821430670688
#define SIN_OBLIQ_2000   0.3977769691083922

/** @endcond */

#define NSTARS 8
static stardef_t StarTable[NSTARS];

#define GetStarPointer(body)    (((body) >= BODY_STAR1) && ((body) <= BODY_STAR8) ? &StarTable[(body) - BODY_STAR1] : NULL)

static stardef_t *UserDefinedStar(astro_body_t body)
{
    stardef_t *star = GetStarPointer(body);
    if (star != NULL && star->dist > 0.0)
        return star;
    return NULL;
}


/**
 * @brief Assign equatorial coordinates to a user-defined star.
 *
 * Some Astronomy Engine functions allow their `body` parameter to
 * be a user-defined fixed point in the sky, loosely called a "star".
 * This function assigns a right ascension, declination, and distance
 * to one of the eight user-defined stars `BODY_STAR1` .. `BODY_STAR8`.
 *
 * Stars are not valid until defined. Once defined, they retain their
 * definition until re-defined by another call to `Astronomy_DefineStar`.
 *
 * @param body
 *      One of the eight user-defined star identifiers: `BODY_STAR1` .. `BODY_STAR8`.
 *
 * @param ra
 *      The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ).
 *      The value is in units of sidereal hours, and must be within the half-open range [0, 24).
 *
 * @param dec
 *      The declination to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ).
 *      The value is in units of degrees north (positive) or south (negative) of the J2000 equator,
 *      and must be within the closed range [-90, +90].
 *
 * @param distanceLightYears
 *      The distance between the star and the Sun, expressed in light-years.
 *      This value is used to calculate the tiny parallax shift as seen by an observer on Earth.
 *      If you don't know the distance to the star, using a large value like 1000 will generally work well.
 *      The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations.
 *
 * @return
 *      `ASTRO_SUCCESS` indicates the star has been defined. Any other value indicates an error,
 *      in which case no change has taken place to any of the star definitions.
 */
astro_status_t Astronomy_DefineStar(astro_body_t body, double ra, double dec, double distanceLightYears)
{
    stardef_t *star = GetStarPointer(body);

    if (star == NULL)
        return ASTRO_INVALID_BODY;

    if (!isfinite(ra) || ra < 0.0 || ra >= 24.0)
        return ASTRO_INVALID_PARAMETER;

    if (!isfinite(dec) || dec < -90.0 || dec > +90.0)
        return ASTRO_INVALID_PARAMETER;

    if (!isfinite(distanceLightYears) || distanceLightYears < 1.0)
        return ASTRO_INVALID_PARAMETER;

    star->ra   = ra;
    star->dec  = dec;
    star->dist = distanceLightYears * AU_PER_LY;

    return ASTRO_SUCCESS;
}

static const terse_vector_t VecZero = { 0.0, 0.0, 0.0 };

static terse_vector_t VecAdd(terse_vector_t a, terse_vector_t b)
{
    terse_vector_t c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}

static void VecIncr(terse_vector_t *target, terse_vector_t source)
{
    target->x += source.x;
    target->y += source.y;
    target->z += source.z;
}

static void VecDecr(terse_vector_t *target, terse_vector_t source)
{
    target->x -= source.x;
    target->y -= source.y;
    target->z -= source.z;
}

static terse_vector_t VecMul(double s, terse_vector_t v)
{
    terse_vector_t p;
    p.x = s * v.x;
    p.y = s * v.y;
    p.z = s * v.z;
    return p;
}

static void VecScale(terse_vector_t *target, double scalar)
{
    target->x *= scalar;
    target->y *= scalar;
    target->z *= scalar;
}

static terse_vector_t VecRamp(terse_vector_t a, terse_vector_t b, double ramp)
{
    terse_vector_t c;
    c.x = (1-ramp)*a.x + ramp*b.x;
    c.y = (1-ramp)*a.y + ramp*b.y;
    c.z = (1-ramp)*a.z + ramp*b.z;
    return c;
}

static terse_vector_t VecMean(terse_vector_t a, terse_vector_t b)
{
    terse_vector_t c;
    c.x = (a.x + b.x) / 2;
    c.y = (a.y + b.y) / 2;
    c.z = (a.z + b.z) / 2;
    return c;
}

static astro_state_vector_t ExportState(body_state_t terse, astro_time_t time)
{
    astro_state_vector_t state;

    state.status = ASTRO_SUCCESS;
    state.x  = terse.r.x;
    state.y  = terse.r.y;
    state.z  = terse.r.z;
    state.vx = terse.v.x;
    state.vy = terse.v.y;
    state.vz = terse.v.z;
    state.t  = time;

    return state;
}

static astro_state_vector_t ExportGravCalc(body_grav_calc_t calc, astro_time_t time)
{
    astro_state_vector_t state;

    state.status = ASTRO_SUCCESS;
    state.x  = calc.r.x;
    state.y  = calc.r.y;
    state.z  = calc.r.z;
    state.vx = calc.v.x;
    state.vy = calc.v.y;
    state.vz = calc.v.z;
    state.t  = time;

    return state;
}

static const double DAYS_PER_TROPICAL_YEAR = 365.24217;
static const double ASEC360 = 1296000.0;
static const double ASEC2RAD = 4.848136811095359935899141e-6;
static const double PI2 = 2.0 * PI;
static const double ARC = 3600.0 * 180.0 / PI;          /* arcseconds per radian */
static const double SECONDS_PER_DAY = 24.0 * 3600.0;
static const double SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;
static const double MEAN_SYNODIC_MONTH = 29.530588;     /* average number of days for Moon to return to the same phase */
static const double EARTH_ORBITAL_PERIOD = 365.256;
static const double NEPTUNE_ORBITAL_PERIOD = 60189.0;

/*
    Degrees of refractive "lift" seen for objects near horizon.
    More precisely, the angle below the horizon a point has to be, at sea level,
    to appear to be exactly on the horizon.
    If the ground plane is higher than sea level, this angle
    needs to be corrected for decreased atmospheric density.
*/
static const double REFRACTION_NEAR_HORIZON = 34.0 / 60.0;


#define             SUN_RADIUS_AU  (SUN_RADIUS_KM / KM_PER_AU)

#define EARTH_MEAN_RADIUS_KM        6371.0            /* mean radius of the Earth's geoid, without atmosphere */
#define EARTH_ATMOSPHERE_KM           88.0            /* effective atmosphere thickness for lunar eclipses. see: https://eclipse.gsfc.nasa.gov/LEcat5/shadow.html */
#define EARTH_ECLIPSE_RADIUS_KM     (EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM)
#define EARTH_EQUATORIAL_RADIUS_AU  (EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU)

#define MOON_MEAN_RADIUS_KM         1737.4
#define MOON_EQUATORIAL_RADIUS_AU   (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU)
#define MOON_POLAR_RADIUS_AU        (MOON_POLAR_RADIUS_KM / KM_PER_AU)

/* The inclination of the moon's rotation axis to the ecliptic plane, in radians. */
#define MOON_AXIS_INCLINATION_RADIANS    (DEG2RAD * 1.543)

static const double ASEC180 = 180.0 * 60.0 * 60.0;      /* arcseconds per 180 degrees (or pi radians) */
static const double EARTH_MOON_MASS_RATIO = 81.30056;

/*
    Masses of the Sun and planets, used for:
    (1) Calculating the Solar System Barycenter
    (2) Integrating the movement of Pluto

    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf

    Page 10 in the above document describes the constants used in the DE405 ephemeris.
    The following are GM values (gravity constant * mass) in [au^3 / day^2].
    This side-steps issues of not knowing the exact values of G and masses M[i];
    the products GM[i] are known extremely accurately.
*/
static const double SUN_GM     = 0.2959122082855911e-03;
static const double MERCURY_GM = 0.4912547451450812e-10;
static const double VENUS_GM   = 0.7243452486162703e-09;
static const double EARTH_GM   = 0.8887692390113509e-09;
static const double MARS_GM    = 0.9549535105779258e-10;
static const double JUPITER_GM = 0.2825345909524226e-06;
static const double SATURN_GM  = 0.8459715185680659e-07;
static const double URANUS_GM  = 0.1292024916781969e-07;
static const double NEPTUNE_GM = 0.1524358900784276e-07;
static const double PLUTO_GM   = 0.2188699765425970e-11;

#define MOON_GM   (EARTH_GM / EARTH_MOON_MASS_RATIO)

/** @cond DOXYGEN_SKIP */
#define ASTRO_ARRAYSIZE(x)    (sizeof(x) / sizeof(x[0]))
#define AU_PER_PARSEC   (ASEC180 / PI)             /* exact definition of how many AU = one parsec */
#define Y2000_IN_MJD    (T0 - MJD_BASIS)
/** @endcond */

static astro_ecliptic_t RotateEquatorialToEcliptic(const double pos[3], double obliq_radians, astro_time_t time);
static int QuadInterp(
    double tm, double dt, double fa, double fm, double fb,
    double *t, double *df_dt);

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
 * @brief Returns the average number of days it takes for a planet to orbit the Sun.
 * @param body One of the planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, or Pluto.
 * @return The mean orbital period of the body, or 0.0 if the `body` parameter is not valid.
 */
double Astronomy_PlanetOrbitalPeriod(astro_body_t body)
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

static astro_vector_t VecError(astro_status_t status, astro_time_t time)
{
    astro_vector_t vec;
    vec.x = vec.y = vec.z = NAN;
    vec.t = time;
    vec.status = status;
    return vec;
}

static astro_state_vector_t StateVecError(astro_status_t status, astro_time_t time)
{
    astro_state_vector_t vec;
    vec.x = vec.y = vec.z = NAN;
    vec.vx = vec.vy = vec.vz = NAN;
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

static astro_time_t TimeError(void)
{
    astro_time_t time;
    time.tt = time.ut = time.eps = time.psi = time.st = NAN;
    return time;
}

static astro_equatorial_t EquError(astro_status_t status)
{
    astro_equatorial_t equ;
    equ.vec = VecError(status, TimeError());
    equ.ra = equ.dec = equ.dist = NAN;
    equ.status = status;
    return equ;
}

static astro_ecliptic_t EclError(astro_status_t status)
{
    astro_ecliptic_t ecl;
    ecl.status = status;
    ecl.elon = ecl.elat = NAN;
    ecl.vec = VecError(status, TimeError());
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
    result.phase_fraction = NAN;
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

static astro_axis_t AxisErr(astro_status_t status, astro_time_t time)
{
    astro_axis_t axis;
    axis.status = status;
    axis.ra = axis.dec = axis.spin = NAN;
    axis.north = VecError(status, time);
    return axis;
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

    Tp = Astronomy_PlanetOrbitalPeriod(body);
    if (Tp <= 0.0)
        return FuncError(ASTRO_INVALID_BODY);

    result.status = ASTRO_SUCCESS;
    result.value = fabs(EARTH_ORBITAL_PERIOD / (EARTH_ORBITAL_PERIOD/Tp - 1.0));
    return result;
}


/**
 * @brief Calculates the angle between two vectors.
 *
 * Given a pair of vectors, this function returns the angle in degrees
 * between the two vectors in 3D space.
 * The angle is measured in the plane that contains both vectors.
 *
 * @param a
 *      The first vector.
 *
 * @param b
 *      The second vector.
 *
 * @returns
 *      On success, the `status` field holds `ASTRO_SUCCESS` and `angle` holds
 *      a number of degrees in the range [0, 180].
 *      If either vector has a zero magnitude or contains NAN (not a number)
 *      components, the `status` will hold the error code `ASTRO_BAD_VECTOR`.
 */
astro_angle_result_t Astronomy_AngleBetween(astro_vector_t a, astro_vector_t b)
{
    double r, dot;
    astro_angle_result_t result;

    r = Astronomy_VectorLength(a) * Astronomy_VectorLength(b);
    if (r < 1.0e-8 || !isfinite(r))
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
 * @brief Converts a J2000 day value to an #astro_time_t value.
 *
 * This function can be useful for reproducing an #astro_time_t structure
 * from its `ut` field only.
 *
 * @param ut
 *      The floating point number of days since noon UTC on January 1, 2000.
 *      This time is based on UTC/UT1 civil time.
 *      See #Astronomy_TerrestrialTime if you instead want to create
 *      a time value based on atomic Terrestrial Time (TT).
 *
 * @returns
 *      An #astro_time_t value for the given `ut` value.
 */
astro_time_t Astronomy_TimeFromDays(double ut)
{
    astro_time_t  time;
    time.ut = ut;
    time.tt = TerrestrialTime(ut);
    time.psi = time.eps = time.st = NAN;
    return time;
}


/**
 * @brief Converts a terrestrial time value into an #astro_time_t value.
 *
 * This function can be used in rare cases where a time must be based
 * on Terrestrial Time (TT) rather than Universal Time (UT).
 * Most developers will want to call #Astronomy_TimeFromDays instead of
 * this function, because usually time is based on civil time adjusted
 * by leap seconds to match the Earth's rotation, rather than the uniformly
 * flowing TT used to calculate solar system dynamics. In rare cases
 * where the caller already knows TT, this function is provided to create
 * an #astro_time_t value that can be passed to Astronomy Engine functions.
 *
 * @param tt
 *      The floating point number of days of uniformly flowing
 *      Terrestrial Time since the J2000 epoch.
 *
 * @returns
 *      An #astro_time_t value for the given `tt` value.
 */
astro_time_t Astronomy_TerrestrialTime(double tt)
{
    /* Iterate to solve to find the correct ut for a given tt, and create an astro_time_t for that time. */
    astro_time_t time = Astronomy_TimeFromDays(tt);
    for(;;)
    {
        double err = tt - time.tt;
        if (fabs(err) < 1.0e-12)
            return time;
        time = Astronomy_AddDays(time, err);
    }
}

#if !defined(ASTRONOMY_ENGINE_NO_CURRENT_TIME)
/**
 * @brief Returns the computer's current date and time in the form of an #astro_time_t.
 *
 * Uses the computer's system clock to find the current UTC date and time.
 * Converts that date and time to an #astro_time_t value and returns the result.
 * Callers can pass this value to other Astronomy Engine functions to calculate
 * current observational conditions.
 *
 * On supported platforms (Linux/Unix, Mac, Windows), the time is measured with
 * microsecond resolution.
 *
 * On unsupported platforms, a compiler error will occur due to lack of
 * microsecond resolution support. However, if whole second resolution is good
 * enough for your application, you can define the preprocessor symbol
 * `ASTRONOMY_ENGINE_WHOLE_SECOND` to use the portable function `time(NULL)`.
 * Alternatively, if you do not need to use `Astronomy_CurrentTime`, you can
 * define the preprocessor symbol `ASTRONOMY_ENGINE_NO_CURRENT_TIME` to
 * exclude this function from your code.
 */
astro_time_t Astronomy_CurrentTime(void)
{
    astro_time_t t;
    double sec;         /* Seconds since midnight January 1, 1970. */

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = (double)tv.tv_sec + tv.tv_usec/1.0e+6;
#elif defined(_WIN32)
    FILETIME ft;
    ULARGE_INTEGER large;
    /* Get time in 100-nanosecond units from January 1, 1601. */
    GetSystemTimePreciseAsFileTime(&ft);
    large.u.LowPart  = ft.dwLowDateTime;
    large.u.HighPart = ft.dwHighDateTime;
    sec = (large.QuadPart - 116444736000000000ULL) / 1.0e+7;
#elif defined(ASTRONOMY_ENGINE_WHOLE_SECOND)
    sec = time(NULL);
#else
    #error Microsecond time resolution is not supported on this platform. Define ASTRONOMY_ENGINE_WHOLE_SECOND to use second resolution instead.
#endif

    /* Convert seconds to days, then subtract to get days since noon on January 1, 2000. */
    t.ut = (sec / SECONDS_PER_DAY) - 10957.5;
    t.tt = TerrestrialTime(t.ut);
    t.psi = t.eps = t.st = NAN;
    return t;
}
#endif

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
    int64_t y = (int64_t)year;
    int64_t m = (int64_t)month;
    int64_t d = (int64_t)day;
    int64_t f = (14 - m) / 12;

    /*
        This formula is adapted from NOVAS C 3.1 function julian_date(),
        which in turn comes from Henry F. Fliegel & Thomas C. Van Flendern:
        Communications of the ACM, Vol 11, No 10, October 1968, p. 657.
        See: https://dl.acm.org/doi/pdf/10.1145/364096.364097

        [Don Cross - 2023-02-25] I modified the formula so that it will
        work correctly with years as far back as -999999.
    */
    int64_t y2000 = (
        (d - 365972956)
        + (1461*(y + 1000000 - f))/4
        + (367*(m - 2 + 12*f))/12
        - (3*((y + 1000100 - f) / 100))/4
    );

    time.ut = (y2000 - 0.5) + (hour / 24.0) + (minute / 1440.0) + (second / 86400.0);
    time.tt = TerrestrialTime(time.ut);
    time.psi = time.eps = time.st = NAN;

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
    sum.eps = sum.psi = sum.st = NAN;

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
    int64_t jd, k, m, n;
    double djd, x;
    const int64_t c = 2500;

    djd = time.ut + 2451545.5;
    jd = (int64_t)floor(djd);

    x = 24.0 * fmod(djd, 1.0);
    if (x < 0.0)
        x += 24.0;
    utc.hour = (int)x;
    x = 60.0 * fmod(x, 1.0);
    utc.minute = (int)x;
    utc.second = 60.0 * fmod(x, 1.0);

    /*
        This is my own adjustment to the NOVAS cal_date logic
        so that it can handle dates much farther back in the past.
        I add c*400 years worth of days at the front,
        then subtract c*400 years at the back,
        which avoids negative values in the formulas that mess up
        the calendar date calculations.
        Any multiple of 400 years has the same number of days,
        because it eliminates all the special cases for leap years.
    */
    k = jd + (68569 + c*146097);
    n = (4 * k) / 146097;
    k = k - (146097 * n + 3)/4;
    m = (4000 * (k+1)) / 1461001;
    k = k - (1461 * m)/4 + 31;

    utc.month = (int) ((80 * k) / 2447);
    utc.day = (int) (k - (2447*utc.month)/80);
    k = utc.month / 11;

    utc.month = (int) (utc.month + 2 - 12*k);
    utc.year = (int) (100 * (n - 49) + m + k - 400*c);

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
 *      If the civil time is outside the year range -999999 to +999999, the function fails
 *      and returns `ASTRO_BAD_TIME`. Years prior to 1583 are treated as if they are
 *      using the modern Gregorian calendar, even when the Julian calendar was actually in effect.
 *      The year before 1 AD, commonly known as 1 BC, is represented by the value 0.
 *      The year 2 BC is represented by -1, etc.
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
 *      A buffer that is `TIME_TEXT_BYTES` (28) bytes or larger is always large enough for this function.
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
    char ytext[20];   /* worst case: "+999999" = 8 characters including terminal '\0'. But gcc 12.2 still complains! */

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

    /* Perform rounding. */
    time.ut += rounding;

    /* Convert linear J2000 days to Gregorian UTC date/time. */
    utc = Astronomy_UtcFromTime(time);
    if (utc.year < -999999 || utc.year > +999999)
        return ASTRO_BAD_TIME;

    if (utc.year < 0)
    {
        snprintf(ytext, sizeof(ytext), "-%06d", -utc.year);
        min_size += 3;  /* '-' prefix and two extra year digits. */
    }
    else if (utc.year <= 9999)
    {
        snprintf(ytext, sizeof(ytext), "%04d", utc.year);
    }
    else
    {
        snprintf(ytext, sizeof(ytext), "+%06d", utc.year);
        min_size += 3;  /* '+' prefix and two extra year digits. */
    }

    /* Check for insufficient buffer size. */
    if (size < min_size)
        return ASTRO_BUFFER_TOO_SMALL;

    /* Format the string. */
    switch (format)
    {
    case TIME_FORMAT_DAY:
        nprinted = snprintf(text, size, "%s-%02d-%02d",
            ytext, utc.month, utc.day);
        break;

    case TIME_FORMAT_MINUTE:
        nprinted = snprintf(text, size, "%s-%02d-%02dT%02d:%02dZ",
            ytext, utc.month, utc.day,
            utc.hour, utc.minute);
        break;

    case TIME_FORMAT_SECOND:
        nprinted = snprintf(text, size, "%s-%02d-%02dT%02d:%02d:%02.0lfZ",
            ytext, utc.month, utc.day,
            utc.hour, utc.minute, floor(utc.second));
        break;

    case TIME_FORMAT_MILLI:
        nprinted = snprintf(text, size, "%s-%02d-%02dT%02d:%02d:%06.3lfZ",
            ytext, utc.month, utc.day,
            utc.hour, utc.minute, floor(1000.0 * utc.second) / 1000.0);
        break;

    default:
        /* We should have already failed for any unknown 'format' value. */
        return ASTRO_INTERNAL_ERROR;
    }

    if (nprinted < 0)
        return ASTRO_INTERNAL_ERROR;    /* should not be possible for snprintf to return a negative number */

    if ((size_t)(1+nprinted) != min_size)
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
    /* Truncated and hand-optimized nutation model. */

    if ((time != NULL) && isnan(time->psi))
    {
        double t, elp, f, d, om, arg, dp, de, sarg, carg;

        t = time->tt / 36525.0;
        elp = fmod(1287104.79305 + t * 129596581.0481,  ASEC360) * ASEC2RAD;
        f   = fmod(335779.526232 + t * 1739527262.8478, ASEC360) * ASEC2RAD;
        d   = fmod(1072260.70369 + t * 1602961601.2090, ASEC360) * ASEC2RAD;
        om  = fmod(450160.398036 - t * 6962890.5431,    ASEC360) * ASEC2RAD;

        sarg = sin(om);
        carg = cos(om);
        dp = (-172064161.0 - 174666.0*t)*sarg + 33386.0*carg;
        de = (92052331.0 + 9086.0*t)*carg + 15377.0*sarg;

        arg = 2.0*(f - d + om);
        sarg = sin(arg);
        carg = cos(arg);
        dp += (-13170906.0 - 1675.0*t)*sarg - 13696.0*carg;
        de += (5730336.0 - 3015.0*t)*carg - 4587.0*sarg;

        arg = 2.0*(f + om);
        sarg = sin(arg);
        carg = cos(arg);
        dp += (-2276413.0 - 234.0*t)*sarg + 2796.0*carg;
        de += (978459.0 - 485.0*t)*carg + 1374.0*sarg;

        arg = 2.0*om;
        sarg = sin(arg);
        carg = cos(arg);
        dp += (2074554.0 + 207.0*t)*sarg - 698.0*carg;
        de += (-897492.0 + 470.0*t)*carg - 291.0*sarg;

        sarg = sin(elp);
        carg = cos(elp);
        dp += (1475877.0 - 3633.0*t)*sarg + 11817.0*carg;
        de += (73871.0 - 184.0*t)*carg - 1924.0*sarg;

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

    /* There is no good answer for what to do if time==NULL. Callers must prevent! */

    iau2000b(time);
    et.dpsi = time->psi;
    et.deps = time->eps;
    et.mobl = mean_obliq(time->tt);
    et.tobl = et.mobl + (et.deps / 3600.0);
    et.tt = time->tt;
    et.ee = et.dpsi * cos(et.mobl * DEG2RAD) / 15.0;

    return et;
}

static void obl_ecl2equ_vec(double obl, astro_time_t time, const double ecl[3], double equ[3])
{
    double obl_rad = obl * DEG2RAD;
    double cos_obl = cos(obl_rad);
    double sin_obl = sin(obl_rad);

    equ[0] = ecl[0];
    equ[1] = ecl[1]*cos_obl - ecl[2]*sin_obl;
    equ[2] = ecl[1]*sin_obl + ecl[2]*cos_obl;
}

static void ecl2equ_vec(astro_time_t time, const double ecl[3], double equ[3])
{
    double obl = mean_obliq(time.tt);
    obl_ecl2equ_vec(obl, time, ecl, equ);
}

static astro_rotation_t precession_rot(astro_time_t time, precess_dir_t dir)
{
    /*
        dir==INTO_2000: converts mean equator of date (EQM) to J2000 mean equator (EQJ).
        dir==FROM_2000: converts J2000 mean equator (EQJ) to mean equator of date (EQM).
    */

    astro_rotation_t rotation;
    double xx, yx, zx, xy, yy, zy, xz, yz, zz;
    double t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;
    double eps0 = 84381.406;

    t = time.tt / 36525;

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

    if (dir == INTO_2000)
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


static void rotate(const double invec[3], const double rot[3][3], double outvec[3])
{
    outvec[0] = rot[0][0]*invec[0] + rot[1][0]*invec[1] + rot[2][0]*invec[2];
    outvec[1] = rot[0][1]*invec[0] + rot[1][1]*invec[1] + rot[2][1]*invec[2];
    outvec[2] = rot[0][2]*invec[0] + rot[1][2]*invec[1] + rot[2][2]*invec[2];
}

static void precession(
    const double pos1[3],
    astro_time_t time,
    precess_dir_t dir,
    double pos2[3])
{
    astro_rotation_t r = precession_rot(time, dir);
    rotate(pos1, r.rot, pos2);
}

static void precession_posvel(
    const double pos1[3],
    const double vel1[3],
    astro_time_t time,
    precess_dir_t dir,
    double pos2[3],
    double vel2[3])
{
    astro_rotation_t r = precession_rot(time, dir);
    rotate(pos1, r.rot, pos2);
    rotate(vel1, r.rot, vel2);
}


static astro_equatorial_t vector2radec(const double pos[3], astro_time_t time)
{
    astro_equatorial_t equ;
    double xyproj;

    /* Copy the cartesian coordinates from the input into the returned structure. */
    equ.vec.status = ASTRO_SUCCESS;
    equ.vec.t = time;
    equ.vec.x = pos[0];
    equ.vec.y = pos[1];
    equ.vec.z = pos[2];

    /* Calculate spherical coordinates: RA, DEC, distance. */
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
        equ.ra = RAD2HOUR * atan2(pos[1], pos[0]);
        if (equ.ra < 0)
            equ.ra += 24.0;

        equ.dec = RAD2DEG * atan2(pos[2], sqrt(xyproj));
    }

    return equ;
}


static astro_rotation_t nutation_rot(astro_time_t *time, precess_dir_t dir)
{
    /*
        Creates a rotation matrix that adds/removes nutation from
        an equatorial vector of date:

        The `dir` parameter is a little misleading, but reflects the
        common task of working together with precession_rot() to
        convert mean equator of 2000 to/from true equator of date.
        Here is the actual result of `dir`:

        dir==INTO_2000: Subtract nutation from true equator of date (EQD)
        to produce mean equator of date (EQM).

        dir==FROM_2000: Add nutation from mean equator of date (EQM) to
        produce true equator of date (EQD).
    */

    astro_rotation_t rotation;
    earth_tilt_t tilt;
    double oblm, oblt, psi, cobm, sobm, cobt, sobt, cpsi, spsi;
    double xx, yx, zx, xy, yy, zy, xz, yz, zz;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    tilt = e_tilt(time);
    oblm = tilt.mobl * DEG2RAD;
    oblt = tilt.tobl * DEG2RAD;
    psi = tilt.dpsi * ASEC2RAD;
    cobm = cos(oblm);
    sobm = sin(oblm);
    cobt = cos(oblt);
    sobt = sin(oblt);
    cpsi = cos(psi);
    spsi = sin(psi);

    xx = cpsi;
    yx = -spsi * cobm;
    zx = -spsi * sobm;
    xy = spsi * cobt;
    yy = cpsi * cobm * cobt + sobm * sobt;
    zy = cpsi * sobm * cobt - cobm * sobt;
    xz = spsi * sobt;
    yz = cpsi * cobm * sobt - sobm * cobt;
    zz = cpsi * sobm * sobt + cobm * cobt;

    if (dir == FROM_2000)
    {
        /* convert J2000 to of-date */
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
        /* convert of-date to J2000 */
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

static void nutation(
    const double inpos[3],
    astro_time_t *time,
    precess_dir_t dir,
    double outpos[3])
{
    astro_rotation_t r = nutation_rot(time, dir);
    rotate(inpos, r.rot, outpos);
}

static void nutation_posvel(
    const double inpos[3],
    const double invel[3],
    astro_time_t *time,
    precess_dir_t dir,
    double outpos[3],
    double outvel[3])
{
    astro_rotation_t r = nutation_rot(time, dir);
    rotate(inpos, r.rot, outpos);
    rotate(invel, r.rot, outvel);
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

/**
 * @brief Calculates Greenwich Apparent Sidereal Time (GAST).
 *
 * Given a date and time, this function calculates the rotation of the
 * Earth, represented by the equatorial angle of the Greenwich prime meridian
 * with respect to distant stars (not the Sun, which moves relative to background
 * stars by almost one degree per day).
 * This angle is called Greenwich Apparent Sidereal Time (GAST).
 * GAST is measured in sidereal hours in the half-open range [0, 24).
 * When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
 * corrected at that time for precession and nutation of the Earth's axis.
 * In this context, the "equinox" is the direction in space where the Earth's
 * orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
 * at the location on the Earth's orbit of the (seasonal) March equinox.
 * As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
 * then starts over at 0.
 * To convert to degrees, multiply the return value by 15.
 *
 * @param time
 *      The date and time for which to find GAST.
 *      The parameter is passed by address because it can be modified by the call:
 *      As an optimization, this function caches the sidereal time value in `time`,
 *      unless it has already been cached, in which case the cached value is reused.
 *      If the `time` pointer is NULL, this function returns a NAN value.
 *
 * @returns {number}
 */
double Astronomy_SiderealTime(astro_time_t *time)
{
    if (time == NULL)
        return NAN;

    if (isnan(time->st))
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

        time->st = gst;
    }

    return time->st;     /* return sidereal hours in the half-open range [0, 24). */
}

static astro_observer_t inverse_terra(const double ovec[3], double st)
{
    double x, y, z, p, F, W, D, c, s, c2, s2;
    double lon_deg, lat_deg, lat, radicand, factor, denom, adjust;
    double height_km, stlocl;
    astro_observer_t observer;
    int count;

    /* Convert from AU to kilometers. */
    x = ovec[0] * KM_PER_AU;
    y = ovec[1] * KM_PER_AU;
    z = ovec[2] * KM_PER_AU;
    p = hypot(x, y);
    if (p < 1.0e-6)
    {
        /* Special case: within 1 millimeter of a pole! */
        /* Use arbitrary longitude, and latitude determined by polarity of z. */
        lon_deg = 0.0;
        lat_deg = (z > 0.0) ? +90.0 : -90.0;
        /* Elevation is calculated directly from z */
        height_km = fabs(z) - EARTH_POLAR_RADIUS_KM;
    }
    else
    {
        stlocl = atan2(y, x);
        /* Calculate exact longitude. */
        lon_deg = RAD2DEG*stlocl - (15.0 * st);
        /* Normalize longitude to the range (-180, +180]. */
        while (lon_deg <= -180.0)
            lon_deg += 360.0;
        while (lon_deg > +180.0)
            lon_deg -= 360.0;
        /* Numerically solve for exact latitude, using Newton's Method. */
        F = EARTH_FLATTENING * EARTH_FLATTENING;
        /* Start with initial latitude estimate, based on a spherical Earth. */
        lat = atan2(z, p);
        for (count = 0; ; ++count)
        {
            if (count > 10)
            {
                fprintf(stderr, "\nFATAL(inverse_terra): did not converge!\n");
                exit(1);
            }

            /* Calculate the error function W(lat). */
            /* We try to find the root of W, meaning where the error is 0. */
            c = cos(lat);
            s = sin(lat);
            factor = (F-1)*EARTH_EQUATORIAL_RADIUS_KM;
            c2 = c*c;
            s2 = s*s;
            radicand = c2 + F*s2;
            denom = sqrt(radicand);
            W = (factor*s*c)/denom - z*c + p*s;
            if (fabs(W) < 1.0e-8)
                break;  /* The error is now negligible. */
            /* Error is still too large. Find the next estimate. */
            /* Calculate D = the derivative of W with respect to lat. */
            D = factor*((c2 - s2)/denom - s2*c2*(F-1)/(factor*radicand)) + z*s + p*c;
            lat -= W/D;
        }
        /* We now have a solution for the latitude in radians. */
        lat_deg = lat * RAD2DEG;
        /* Solve for exact height in meters. */
        /* There are two formulas I can use. Use whichever has the less risky denominator. */
        adjust = EARTH_EQUATORIAL_RADIUS_KM / denom;
        if (fabs(s) > fabs(c))
            height_km = z/s - F*adjust;
        else
            height_km = p/c - adjust;
    }

    observer.latitude = lat_deg;
    observer.longitude = lon_deg;
    observer.height = 1000.0 * height_km;
    return observer;
}

static void terra(astro_observer_t observer, double st, double pos[3], double vel[3])
{
    static const double ANGVEL = 7.2921150e-5;

    double phi = observer.latitude * DEG2RAD;
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double c = 1.0 / hypot(cosphi, sinphi*EARTH_FLATTENING);
    double s = c * (EARTH_FLATTENING * EARTH_FLATTENING);
    double ht_km = observer.height / 1000.0;
    double ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
    double ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
    double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
    double sinst = sin(stlocl);
    double cosst = cos(stlocl);

    if (pos != NULL)
    {
        pos[0] = ach * cosphi * cosst / KM_PER_AU;
        pos[1] = ach * cosphi * sinst / KM_PER_AU;
        pos[2] = ash * sinphi / KM_PER_AU;
    }

    if (vel != NULL)
    {
        vel[0] = -(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * sinst;
        vel[1] = +(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * cosst;
        vel[2] = 0.0;
    }
}

static void geo_pos(astro_time_t *time, astro_observer_t observer, double pos[3])
{
    double gast;
    double pos1[3], pos2[3];

    if (time == NULL)
    {
        pos[0] = pos[1] = pos[2] = NAN;
    }
    else
    {
        gast = Astronomy_SiderealTime(time);
        terra(observer, gast, pos1, NULL);
        nutation(pos1, time, INTO_2000, pos2);
        precession(pos2, *time, INTO_2000, pos);
    }
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

static void Term(const MoonContext *ctx, int p, int q, int r, int s, double *x, double *y)
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
//$ASTRO_ADDSOL()
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
 * @brief Calculates equatorial geocentric position of the Moon at a given time.
 *
 * Given a time of observation, calculates the Moon's position as a vector.
 * The vector gives the location of the Moon's center relative to the Earth's center
 * with x-, y-, and z-components measured in astronomical units.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 *
 * This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
 * which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
 * It is adapted from Turbo Pascal code from the book
 * [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
 * by Montenbruck and Pfleger.
 *
 * To calculate ecliptic spherical coordinates instead, see #Astronomy_EclipticGeoMoon.
 *
 * @param time  The date and time for which to calculate the Moon's position.
 * @return The Moon's position as a vector in J2000 Cartesian equatorial (EQJ) coordinates.
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

    /* Convert equatorial coordinates from mean equinox of date to J2000 mean equinox. */
    precession(mpos1, time, INTO_2000, mpos2);

    vector.status = ASTRO_SUCCESS;
    vector.x = mpos2[0];
    vector.y = mpos2[1];
    vector.z = mpos2[2];
    vector.t = time;
    return vector;
}


/**
 * @brief Calculates spherical ecliptic geocentric position of the Moon.
 *
 * Given a time of observation, calculates the Moon's geocentric position
 * in ecliptic spherical coordinates. Provides the ecliptic latitude and
 * longitude in degrees, and the geocentric distance in astronomical units (AU).
 *
 * The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and
 * equatorial plane at the specified time. This means the Earth's equator
 * is corrected for precession and nutation, and the plane of the Earth's
 * orbit is corrected for gradual obliquity drift.
 *
 * This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
 * which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
 * It is adapted from Turbo Pascal code from the book
 * [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
 * by Montenbruck and Pfleger.
 *
 * To calculate a J2000 mean equator vector instead, use #Astronomy_GeoMoon.
 *
 * @param time  The date and time for which to calculate the Moon's position.
 * @return The Moon's position expressed in ecliptic coordinates using the true equinox of date (ECT).
 */
astro_spherical_t Astronomy_EclipticGeoMoon(astro_time_t time)
{
    astro_spherical_t sphere;
    astro_ecliptic_t eclip;
    earth_tilt_t et;
    double dist_cos_lat, ecm[3], eqm[3], eqd[3];

    /* CalcMoon produces ecliptic coordinates in mean equinox of date (ECM). */
    CalcMoon(time.tt / 36525.0, &sphere.lon, &sphere.lat, &sphere.dist);

    /* Calculate vector in ecliptic coordinates (ECM). */
    dist_cos_lat = sphere.dist * cos(sphere.lat);
    ecm[0] = dist_cos_lat * cos(sphere.lon);
    ecm[1] = dist_cos_lat * sin(sphere.lon);
    ecm[2] = sphere.dist * sin(sphere.lat);

    /* Obtain true and mean obliquity angles for the given time. */
    /* This serves to pre-calculate the nutation also, and cache it in `time`. */
    et = e_tilt(&time);

    /* Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date. */
    /* In other words, convert ECM to EQM. */
    obl_ecl2equ_vec(et.mobl, time, ecm, eqm);

    /* Add nutation to convert ECM to true equatorial coordinates of date (EQD). */
    nutation(eqm, &time, FROM_2000, eqd);

    /* Convert back to ecliptic, this time in true equinox of date (ECT). */
    eclip = RotateEquatorialToEcliptic(eqd, DEG2RAD * et.tobl, time);

    /* Package the return value. */
    /* CalcMoon() already set sphere.dist to the correct value. */
    sphere.status = eclip.status;
    sphere.lat = eclip.elat;
    sphere.lon = eclip.elon;
    return sphere;
}


/**
 * @brief Calculates equatorial geocentric position and velocity of the Moon at a given time.
 *
 * Given a time of observation, calculates the Moon's position and velocity vectors.
 * The position and velocity are of the Moon's center relative to the Earth's center.
 * The position (x, y, z) components are expressed in AU (astronomical units).
 * The velocity (vx, vy, vz) components are expressed in AU/day.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 *
 * If you need the Moon's position only, and not its velocity,
 * it is much more efficient to use #Astronomy_GeoMoon instead.
 *
 * @param time  The date and time for which to calculate the Moon's position and velocity.
 * @return The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).
 */
astro_state_vector_t Astronomy_GeoMoonState(astro_time_t time)
{
    /*
        This is a hack, because trying to figure out how to derive a time
        derivative for CalcMoon() would be extremely painful!
        Calculate just before and just after the given time.
        Average to find position, subtract to find velocity.
    */
    const double dt = 1.0e-5;   /* 0.864 seconds */
    astro_vector_t r1, r2;
    astro_time_t t1, t2;
    astro_state_vector_t s;

    t1 = Astronomy_AddDays(time, -dt);
    t2 = Astronomy_AddDays(time, +dt);

    r1 = Astronomy_GeoMoon(t1);
    r2 = Astronomy_GeoMoon(t2);

    /* The desired position is the average of the two calculated positions. */
    s.x = (r1.x + r2.x) / 2;
    s.y = (r1.y + r2.y) / 2;
    s.z = (r1.z + r2.z) / 2;

    /* The difference of the position vectors divided by the time span gives the velocity vector. */
    s.vx = (r2.x - r1.x) / (2 * dt);
    s.vy = (r2.y - r1.y) / (2 * dt);
    s.vz = (r2.z - r1.z) / (2 * dt);
    s.t = time;
    s.status = ASTRO_SUCCESS;

    return s;
}


/**
 * @brief Calculates the geocentric position and velocity of the Earth/Moon barycenter.
 *
 * Given a time of observation, calculates the geocentric position and velocity vectors
 * of the Earth/Moon barycenter (EMB).
 * The position (x, y, z) components are expressed in AU (astronomical units).
 * The velocity (vx, vy, vz) components are expressed in AU/day.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 *
 * @param time  The date and time for which to calculate the EMB vectors.
 * @return The EMB's position and velocity vectors in geocentric J2000 equatorial coordinates.
 */
astro_state_vector_t Astronomy_GeoEmbState(astro_time_t time)
{
    astro_state_vector_t s = Astronomy_GeoMoonState(time);
    const double d = 1.0 + EARTH_MOON_MASS_RATIO;
    s.x /= d;
    s.y /= d;
    s.z /= d;
    s.vx /= d;
    s.vy /= d;
    s.vz /= d;
    return s;
}


/**
 * @brief Calculates the Moon's libration angles at a given moment in time.
 *
 * Libration is an observed back-and-forth wobble of the portion of the
 * Moon visible from the Earth. It is caused by the imperfect tidal locking
 * of the Moon's fixed rotation rate, compared to its variable angular speed
 * of orbit around the Earth.
 *
 * This function calculates a pair of perpendicular libration angles,
 * one representing rotation of the Moon in ecliptic longitude `elon`, the other
 * in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.
 *
 * This function also returns the geocentric position of the Moon
 * expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
 * distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
 * and the apparent angular diameter of the Moon `diam_deg`.
 *
 * @param time  The date and time for which to calculate libration angles.
 * @return The Moon's ecliptic position and libration angles as seen from the Earth.
 */
astro_libration_t Astronomy_Libration(astro_time_t time)
{
    astro_libration_t lib;
    double t, t2, t3, t4;
    double f, omega, w, a, ldash, ldash2, bdash, bdash2;
    double k1, k2, m, mdash, d, e, rho, sigma, tau;
    const double sin_I = sin(MOON_AXIS_INCLINATION_RADIANS);
    const double cos_I = cos(MOON_AXIS_INCLINATION_RADIANS);

    t = time.tt / 36525.0;
    t2 = t * t;
    t3 = t2 * t;
    t4 = t2 * t2;

    double mlon;    /* Moon's ecliptic longitude in radians. */
    double mlat;    /* Moon's ecliptic latitude in radians. */
    CalcMoon(t, &mlon, &mlat, &lib.dist_km);
    lib.mlon = RAD2DEG * mlon;
    lib.mlat = RAD2DEG * mlat;
    lib.dist_km *= KM_PER_AU;
    lib.diam_deg = (2.0 * RAD2DEG) * atan(MOON_MEAN_RADIUS_KM / sqrt(lib.dist_km*lib.dist_km - MOON_MEAN_RADIUS_KM*MOON_MEAN_RADIUS_KM));

    /* Moon's argument of latitude in radians. */
    f = DEG2RAD * NormalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000);

    /* Moon's ascending node's mean longitude in radians. */
    omega = DEG2RAD * NormalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000);

    /* Sun's mean anomaly. */
    m = DEG2RAD * NormalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000);

    /* Moon's mean anomaly. */
    mdash = DEG2RAD * NormalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000);

    /* Moon's mean elongation. */
    d = DEG2RAD * NormalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000);

    /* Eccentricity of the Earth's orbit. */
    e = 1.0 - 0.002516*t - 0.0000074*t2;

    /* Optical librations */
    w = mlon - omega;
    a = atan2(sin(w)*cos(mlat)*cos_I - sin(mlat)*sin_I, cos(w)*cos(mlat));
    ldash = LongitudeOffset(RAD2DEG * (a - f));
    bdash = asin(-sin(w)*cos(mlat)*sin_I - sin(mlat)*cos_I);

    /* Physical librations */
    k1 = DEG2RAD*(119.75 + 131.849*t);
    k2 = DEG2RAD*(72.56 + 20.186*t);

    rho = (
        -0.02752*cos(mdash) +
        -0.02245*sin(f) +
        +0.00684*cos(mdash - 2*f) +
        -0.00293*cos(2*f) +
        -0.00085*cos(2*f - 2*d) +
        -0.00054*cos(mdash - 2*d) +
        -0.00020*sin(mdash + f) +
        -0.00020*cos(mdash + 2*f) +
        -0.00020*cos(mdash - f) +
        +0.00014*cos(mdash + 2*f - 2*d)
    );

    sigma = (
        -0.02816*sin(mdash) +
        +0.02244*cos(f) +
        -0.00682*sin(mdash - 2*f) +
        -0.00279*sin(2*f) +
        -0.00083*sin(2*f - 2*d) +
        +0.00069*sin(mdash - 2*d) +
        +0.00040*cos(mdash + f) +
        -0.00025*sin(2*mdash) +
        -0.00023*sin(mdash + 2*f) +
        +0.00020*cos(mdash - f) +
        +0.00019*sin(mdash - f) +
        +0.00013*sin(mdash + 2*f - 2*d) +
        -0.00010*cos(mdash - 3*f)
    );

    tau = (
        +0.02520*e*sin(m) +
        +0.00473*sin(2*mdash - 2*f) +
        -0.00467*sin(mdash) +
        +0.00396*sin(k1) +
        +0.00276*sin(2*mdash - 2*d) +
        +0.00196*sin(omega) +
        -0.00183*cos(mdash - f) +
        +0.00115*sin(mdash - 2*d) +
        -0.00096*sin(mdash - d) +
        +0.00046*sin(2*f - 2*d) +
        -0.00039*sin(mdash - f) +
        -0.00032*sin(mdash - m - d) +
        +0.00027*sin(2*mdash - m - 2*d) +
        +0.00023*sin(k2) +
        -0.00014*sin(2*d) +
        +0.00014*cos(2*mdash - 2*f) +
        -0.00012*sin(mdash - 2*f) +
        -0.00012*sin(2*mdash) +
        +0.00011*sin(2*mdash - 2*m - 2*d)
    );

    ldash2 = -tau + (rho*cos(a) + sigma*sin(a))*tan(bdash);
    bdash *= RAD2DEG;
    bdash2 = sigma*cos(a) - rho*sin(a);

    lib.elon = ldash + ldash2;
    lib.elat = bdash + bdash2;

    return lib;
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

typedef struct
{
    double mu;
    double al[2];
    vsop_series_t a;
    vsop_series_t l;
    vsop_series_t z;
    vsop_series_t zeta;
}
jupiter_moon_t;
/** @endcond */

//$ASTRO_C_VSOP(Mercury);
//$ASTRO_C_VSOP(Venus);
//$ASTRO_C_VSOP(Earth);
//$ASTRO_C_VSOP(Mars);
//$ASTRO_C_VSOP(Jupiter);
//$ASTRO_C_VSOP(Saturn);
//$ASTRO_C_VSOP(Uranus);
//$ASTRO_C_VSOP(Neptune);

/** @cond DOXYGEN_SKIP */
#define VSOPFORMULA(x)    { ASTRO_ARRAYSIZE(x), x }
/** @endcond */

static const vsop_model_t vsop[] =
{
    { { VSOPFORMULA(vsop_lon_Mercury),  VSOPFORMULA(vsop_lat_Mercury),  VSOPFORMULA(vsop_rad_Mercury) } },
    { { VSOPFORMULA(vsop_lon_Venus),    VSOPFORMULA(vsop_lat_Venus),    VSOPFORMULA(vsop_rad_Venus)   } },
    { { VSOPFORMULA(vsop_lon_Earth),    VSOPFORMULA(vsop_lat_Earth),    VSOPFORMULA(vsop_rad_Earth)   } },
    { { VSOPFORMULA(vsop_lon_Mars),     VSOPFORMULA(vsop_lat_Mars),     VSOPFORMULA(vsop_rad_Mars)    } },
    { { VSOPFORMULA(vsop_lon_Jupiter),  VSOPFORMULA(vsop_lat_Jupiter),  VSOPFORMULA(vsop_rad_Jupiter) } },
    { { VSOPFORMULA(vsop_lon_Saturn),   VSOPFORMULA(vsop_lat_Saturn),   VSOPFORMULA(vsop_rad_Saturn)  } },
    { { VSOPFORMULA(vsop_lon_Uranus),   VSOPFORMULA(vsop_lat_Uranus),   VSOPFORMULA(vsop_rad_Uranus)  } },
    { { VSOPFORMULA(vsop_lon_Neptune),  VSOPFORMULA(vsop_lat_Neptune),  VSOPFORMULA(vsop_rad_Neptune) } }
};

/** @cond DOXYGEN_SKIP */
#define CalcEarth(time)     CalcVsop(&vsop[BODY_EARTH], (time))
#define LON_INDEX 0
#define LAT_INDEX 1
#define RAD_INDEX 2
/** @endcond */

static void VsopCoords(const vsop_model_t *model, double t, double sphere[3])
{
    int k, s, i;
    double incr;

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
                sum  += term->amplitude * cos(term->phase + (t * term->frequency));
            }
            incr = tpower * sum;
            if (k == LON_INDEX)
                incr = fmod(incr, PI2);     /* improve precision for longitudes, which can be hundreds of radians */
            sphere[k] += incr;
            tpower *= t;
        }
    }
}


static terse_vector_t VsopRotate(const double ecl[3])
{
    terse_vector_t equ;

    /*
        X        +1.000000000000  +0.000000440360  -0.000000190919   X
        Y     =  -0.000000479966  +0.917482137087  -0.397776982902   Y
        Z FK5     0.000000000000  +0.397776982902  +0.917482137087   Z VSOP87A
    */

    equ.x = ecl[0] + 0.000000440360*ecl[1] - 0.000000190919*ecl[2];
    equ.y = -0.000000479966*ecl[0] + 0.917482137087*ecl[1] - 0.397776982902*ecl[2];
    equ.z = 0.397776982902*ecl[1] + 0.917482137087*ecl[2];

    return equ;
}


static void VsopSphereToRect(double lon, double lat, double radius, double pos[3])
{
    double r_coslat = radius * cos(lat);
    double coslon = cos(lon);
    double sinlon = sin(lon);
    pos[0] = r_coslat * coslon;
    pos[1] = r_coslat * sinlon;
    pos[2] = radius * sin(lat);
}

static const double DAYS_PER_MILLENNIUM = 365250.0;


static astro_vector_t CalcVsop(const vsop_model_t *model, astro_time_t time)
{
    double t = time.tt / DAYS_PER_MILLENNIUM;
    double sphere[3];       /* lon, lat, rad */
    double eclip[3];
    astro_vector_t vector;
    terse_vector_t pos;

    /* Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates. */
    VsopCoords(model, t, sphere);

    /* Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates. */
    VsopSphereToRect(sphere[LON_INDEX], sphere[LAT_INDEX], sphere[RAD_INDEX], eclip);

    /* Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates. */
    pos = VsopRotate(eclip);

    /* Package the position as astro_vector_t. */
    vector.status = ASTRO_SUCCESS;
    vector.t = time;
    vector.x = pos.x;
    vector.y = pos.y;
    vector.z = pos.z;

    return vector;
}


static void VsopDeriv(const vsop_model_t *model, double t, double deriv[3])
{
    int k, s, i;

    for (k=0; k < 3; ++k)
    {
        double tpower = 1.0;        /* t^s */
        double dpower = 0.0;        /* t^(s-1) */
        const vsop_formula_t *formula = &model->formula[k];
        deriv[k] = 0.0;
        for (s=0; s < formula->nseries; ++s)
        {
            double sin_sum = 0.0;
            double cos_sum = 0.0;
            const vsop_series_t *series = &formula->series[s];
            for (i=0; i < series->nterms; ++i)
            {
                const vsop_term_t *term = &series->term[i];
                double angle = term->phase + (t * term->frequency);
                sin_sum += term->amplitude * term->frequency * sin(angle);
                if (s > 0)
                    cos_sum += term->amplitude * cos(angle);
            }
            deriv[k] += (s * dpower * cos_sum) - (tpower * sin_sum);
            dpower = tpower;
            tpower *= t;
        }
    }
}


static body_state_t CalcVsopPosVel(const vsop_model_t *model, double tt)
{
    body_state_t state;
    double t = tt / DAYS_PER_MILLENNIUM;
    double sphere[3];   /* lon, lat, r */
    double deriv[3];    /* d(lon)/dt, d(lat)/dt, dr/dt */
    double eclip[3];
    double dr_dt, dlat_dt, dlon_dt;
    double r, coslat, coslon, sinlat, sinlon;

    state.tt = tt;
    VsopCoords(model, t, sphere);
    VsopSphereToRect(sphere[LON_INDEX], sphere[LAT_INDEX], sphere[RAD_INDEX], eclip);
    state.r = VsopRotate(eclip);

    VsopDeriv(model, t, deriv);

    /* Use spherical coords and spherical derivatives to calculate */
    /* the velocity vector in rectangular coordinates. */

    /* Calculate mnemonic variables to help keep the math straight. */
    coslon = cos(sphere[LON_INDEX]);
    sinlon = sin(sphere[LON_INDEX]);
    coslat = cos(sphere[LAT_INDEX]);
    sinlat = sin(sphere[LAT_INDEX]);
    r = sphere[RAD_INDEX];
    dlon_dt = deriv[LON_INDEX];
    dlat_dt = deriv[LAT_INDEX];
    dr_dt   = deriv[RAD_INDEX];

    /* vx = dx/dt */
    eclip[0] = (dr_dt * coslat * coslon) - (r * sinlat * coslon * dlat_dt) - (r * coslat * sinlon * dlon_dt);

    /* vy = dy/dt */
    eclip[1] = (dr_dt * coslat * sinlon) - (r * sinlat * sinlon * dlat_dt) + (r * coslat * coslon * dlon_dt);

    /* vz = dz/dt */
    eclip[2] = (dr_dt * sinlat) + (r * coslat * dlat_dt);

    /* Rotate the velocity vector from ecliptic to equatorial coordinates. */
    state.v = VsopRotate(eclip);

    /* Convert speed units from [AU/millennium] to [AU/day]. */
    VecScale(&state.v, 1 / DAYS_PER_MILLENNIUM);

    return state;
}


static double VsopHelioDistance(const vsop_model_t *model, astro_time_t time)
{
    int s, i;
    double t = time.tt / DAYS_PER_MILLENNIUM;
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


static void AdjustBarycenter(astro_vector_t *ssb, astro_time_t time, astro_body_t body, double planet_gm)
{
    astro_vector_t planet;
    double shift;

    shift = planet_gm / (planet_gm + SUN_GM);
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

    AdjustBarycenter(&ssb, time, BODY_JUPITER, JUPITER_GM);
    AdjustBarycenter(&ssb, time, BODY_SATURN,  SATURN_GM);
    AdjustBarycenter(&ssb, time, BODY_URANUS,  URANUS_GM);
    AdjustBarycenter(&ssb, time, BODY_NEPTUNE, NEPTUNE_GM);

    return ssb;
}

/*------------------ begin general gravity simulator ------------------*/

static terse_vector_t UpdatePosition(double dt, terse_vector_t r, terse_vector_t v, terse_vector_t a)
{
    r.x += (v.x + a.x*dt/2) * dt;
    r.y += (v.y + a.y*dt/2) * dt;
    r.z += (v.z + a.z*dt/2) * dt;
    return r;
}


static terse_vector_t UpdateVelocity(double dt, terse_vector_t v, terse_vector_t a)
{
    v.x += dt * a.x;
    v.y += dt * a.y;
    v.z += dt * a.z;
    return v;
}


static body_state_t AdjustBarycenterPosVel(body_state_t *ssb, double tt, astro_body_t body, double planet_gm)
{
    body_state_t planet;
    double shift;

    /*
        This function does 2 important things:
        1. Adjusts 'ssb' by the effect of one major body on the Solar System Barycenter.
        2, Returns the heliocentric position of that major body.
    */

    shift = planet_gm / (planet_gm + SUN_GM);
    planet = CalcVsopPosVel(&vsop[body], tt);
    VecIncr(&ssb->r, VecMul(shift, planet.r));
    VecIncr(&ssb->v, VecMul(shift, planet.v));

    return planet;
}


static void MajorBodyBary(major_bodies_t *bary, double tt)
{
    /* Initialize the Sun's position as a zero vector, then adjust it from pulls from the planets. */
    bary->Sun.tt = tt;
    bary->Sun.r = VecZero;
    bary->Sun.v = VecZero;

    /* Calculate heliocentric planet positions and SSB. */
    bary->Jupiter = AdjustBarycenterPosVel(&bary->Sun, tt, BODY_JUPITER, JUPITER_GM);
    bary->Saturn  = AdjustBarycenterPosVel(&bary->Sun, tt, BODY_SATURN,  SATURN_GM);
    bary->Uranus  = AdjustBarycenterPosVel(&bary->Sun, tt, BODY_URANUS,  URANUS_GM);
    bary->Neptune = AdjustBarycenterPosVel(&bary->Sun, tt, BODY_NEPTUNE, NEPTUNE_GM);

    /* Convert planet [pos, vel] from heliocentric to barycentric. */
    VecDecr(&bary->Jupiter.r, bary->Sun.r);  VecDecr(&bary->Jupiter.v, bary->Sun.v);
    VecDecr(&bary->Saturn.r,  bary->Sun.r);  VecDecr(&bary->Saturn.v,  bary->Sun.v);
    VecDecr(&bary->Uranus.r,  bary->Sun.r);  VecDecr(&bary->Uranus.v,  bary->Sun.v);
    VecDecr(&bary->Neptune.r, bary->Sun.r);  VecDecr(&bary->Neptune.v, bary->Sun.v);

    /* Convert heliocentric SSB to barycentric Sun. */
    VecScale(&bary->Sun.r, -1.0);
    VecScale(&bary->Sun.v, -1.0);
}


static void AddAcceleration(terse_vector_t *acc, terse_vector_t small_pos, double gm, terse_vector_t major_pos)
{
    double dx, dy, dz, r2, pull;

    dx = major_pos.x - small_pos.x;
    dy = major_pos.y - small_pos.y;
    dz = major_pos.z - small_pos.z;

    r2 = dx*dx + dy*dy + dz*dz;
    pull = gm / (r2 * sqrt(r2));

    acc->x += dx * pull;
    acc->y += dy * pull;
    acc->z += dz * pull;
}


static terse_vector_t SmallBodyAcceleration(terse_vector_t small_pos, const major_bodies_t *bary)
{
    terse_vector_t acc = VecZero;

    /* Use barycentric coordinates of the Sun and major planets to calculate gravitational accelerations. */
    AddAcceleration(&acc, small_pos, SUN_GM,     bary->Sun.r);
    AddAcceleration(&acc, small_pos, JUPITER_GM, bary->Jupiter.r);
    AddAcceleration(&acc, small_pos, SATURN_GM,  bary->Saturn.r);
    AddAcceleration(&acc, small_pos, URANUS_GM,  bary->Uranus.r);
    AddAcceleration(&acc, small_pos, NEPTUNE_GM, bary->Neptune.r);

    return acc;
}


body_grav_calc_t GravSim(           /* out: [pos, vel, acc] of the simulated body at time tt2 */
    major_bodies_t *bary2,          /* temp: work area for major body barycentric state */
    double tt2,                     /* in:  a target time to be calculated (either before or after tt1) */
    const body_grav_calc_t *calc1)  /* in:  [pos, vel, acc] of the simulated body at time tt1 */
{
    body_grav_calc_t calc2;
    terse_vector_t approx_pos;
    terse_vector_t acc;
    const double dt = tt2 - calc1->tt;

    /* Calculate where the major bodies (Sun, Jupiter...Neptune) will be at the next time step. */
    MajorBodyBary(bary2, tt2);

    /* Estimate position of small body as if current acceleration applies across the whole time interval. */
    /* approx_pos = pos1 + vel1*dt + (1/2)acc*dt^2 */
    approx_pos = UpdatePosition(dt, calc1->r, calc1->v, calc1->a);

    /* Calculate acceleration experienced by small body at approximate next location. */
    acc = SmallBodyAcceleration(approx_pos, bary2);

    /* Calculate the average acceleration of the endpoints. */
    /* This becomes our estimate of the mean effective acceleration over the whole interval. */
    acc = VecMean(acc, calc1->a);

    /* Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration. */
    calc2.r = UpdatePosition(dt, calc1->r, calc1->v, acc);
    calc2.v = UpdateVelocity(dt, calc1->v, acc);
    calc2.a = SmallBodyAcceleration(calc2.r, bary2);
    calc2.tt = tt2;
    return calc2;
}


static body_grav_calc_t GravFromState(major_bodies_t *bary, const body_state_t *state)
{
    body_grav_calc_t calc;

    MajorBodyBary(bary, state->tt);

    calc.tt = state->tt;
    calc.r  = VecAdd(state->r, bary->Sun.r);      /* convert heliocentric to barycentric */
    calc.v  = VecAdd(state->v, bary->Sun.v);      /* convert heliocentric to barycentric */
    calc.a  = SmallBodyAcceleration(calc.r, bary);

    return calc;
}


static void CalcSolarSystem(astro_grav_sim_t *sim)
{
    int body;
    double tt = sim->curr->time.tt;
    body_state_t *grav = sim->curr->gravitators;
    body_state_t *sun = &grav[BODY_SUN];

    /* Initialize the Sun's position/velocity as zero vectors, then adjust from pulls from the planets. */
    sun->tt = tt;
    sun->r  = VecZero;
    sun->v  = VecZero;

    /* Calculate the position of each planet, and adjust the SSB position accordingly. */
    grav[BODY_MERCURY] = AdjustBarycenterPosVel(sun, tt, BODY_MERCURY, MERCURY_GM);
    grav[BODY_VENUS  ] = AdjustBarycenterPosVel(sun, tt, BODY_VENUS,   VENUS_GM);
    grav[BODY_EARTH  ] = AdjustBarycenterPosVel(sun, tt, BODY_EARTH,   EARTH_GM + MOON_GM);
    grav[BODY_MARS   ] = AdjustBarycenterPosVel(sun, tt, BODY_MARS,    MARS_GM);
    grav[BODY_JUPITER] = AdjustBarycenterPosVel(sun, tt, BODY_JUPITER, JUPITER_GM);
    grav[BODY_SATURN ] = AdjustBarycenterPosVel(sun, tt, BODY_SATURN,  SATURN_GM);
    grav[BODY_URANUS ] = AdjustBarycenterPosVel(sun, tt, BODY_URANUS,  URANUS_GM);
    grav[BODY_NEPTUNE] = AdjustBarycenterPosVel(sun, tt, BODY_NEPTUNE, NEPTUNE_GM);

    /* Convert planet state vectors from heliocentric to barycentric. */
    for (body = BODY_MERCURY; body <= BODY_NEPTUNE; ++body)
    {
        VecDecr(&grav[body].r, sun->r);
        VecDecr(&grav[body].v, sun->v);
    }

    /* Convert heliocentric SSB to barycentric Sun. */
    VecScale(&sun->r, -1.0);
    VecScale(&sun->v, -1.0);
}


static void CalcBodyAccelerations(astro_grav_sim_t *sim)
{
    int i;
    const body_state_t *grav = sim->curr->gravitators;

    /* Calculate the gravitational acceleration experienced by the simulated bodies. */
    for (i = 0; i < sim->numBodies; ++i)
    {
        body_grav_calc_t *calc = &sim->curr->bodies[i];

        calc->a = VecZero;

        AddAcceleration(&calc->a, calc->r, SUN_GM,              grav[BODY_SUN    ].r);
        AddAcceleration(&calc->a, calc->r, MERCURY_GM,          grav[BODY_MERCURY].r);
        AddAcceleration(&calc->a, calc->r, VENUS_GM,            grav[BODY_VENUS  ].r);
        AddAcceleration(&calc->a, calc->r, EARTH_GM + MOON_GM,  grav[BODY_EARTH  ].r);
        AddAcceleration(&calc->a, calc->r, MARS_GM,             grav[BODY_MARS   ].r);
        AddAcceleration(&calc->a, calc->r, JUPITER_GM,          grav[BODY_JUPITER].r);
        AddAcceleration(&calc->a, calc->r, SATURN_GM,           grav[BODY_SATURN ].r);
        AddAcceleration(&calc->a, calc->r, URANUS_GM,           grav[BODY_URANUS ].r);
        AddAcceleration(&calc->a, calc->r, NEPTUNE_GM,          grav[BODY_NEPTUNE].r);
    }
}


static body_state_t *GravSimBodyStatePtr(astro_grav_sim_t *sim, astro_body_t body)
{
    /*
        Return a pointer to the place where we cache the barycentric
        state of the given body, or NULL if this body is not one
        that we use for calculating gravitational interactions.
    */

    if ((body == BODY_SUN) || (body >= BODY_MERCURY && body <= BODY_NEPTUNE))
        return &sim->curr->gravitators[body];

    return NULL;
}


static astro_state_vector_t GravSimOriginState(astro_grav_sim_t *sim)
{
    body_state_t *optr;
    astro_time_t time = sim->curr->time;

    if (sim->originBody == BODY_SSB)
    {
        /* The barycentric state of the SSB is zero, by definition. */
        astro_state_vector_t state;
        state.status = ASTRO_SUCCESS;
        state.t  = time;
        state.x  = 0.0;
        state.y  = 0.0;
        state.z  = 0.0;
        state.vx = 0.0;
        state.vy = 0.0;
        state.vz = 0.0;
        return state;
    }

    optr = GravSimBodyStatePtr(sim, sim->originBody);
    if (optr != NULL)
        return ExportState(*optr, time);

    /* We only support the VSOP bodies, for efficiency. */
    return StateVecError(ASTRO_INVALID_BODY, time);
}


static void GravSimDuplicate(astro_grav_sim_t *sim)
{
    /* Copy the current state into the previous state, so that both become the same moment in time. */
    sim->prev->time = sim->curr->time;
    memcpy(sim->prev->gravitators, sim->curr->gravitators, sizeof(sim->prev->gravitators));
    memcpy(sim->prev->bodies, sim->curr->bodies, ((size_t)sim->numBodies) * sizeof(sim->prev->bodies[0]));
}


/**
 * @brief Allocate and initialize a gravity step simulator.
 *
 * Prepares to simulate a series of incremental time steps,
 * simulating the movement of zero or more small bodies through the Solar System
 * acting under gravitational attraction from the Sun and planets.
 *
 * After calling this function, you can call #Astronomy_GravSimUpdate
 * as many times as desired to advance the simulation by small time steps.
 *
 * If this function succeeds (returns `ASTRO_SUCCESS`), `sim`
 * will be set to a dynamically allocated object. The caller is
 * then responsible for eventually calling #Astronomy_GravSimFree
 * to release the memory.
 *
 * @param simOut
 *      The address of a pointer to store the newly allocated simulation object.
 *      The type #astro_grav_sim_t is an opaque type, so its internal structure is not documented.
 *
 * @param originBody
 *      Specifies the origin of the reference frame.
 *      All position vectors and velocity vectors will use `originBody`
 *      as the origin of the coordinate system.
 *      This origin applies to all the input vectors provided in the
 *      `bodyStateArray` parameter of this function, along with all
 *      output vectors returned by #Astronomy_GravSimUpdate.
 *      Most callers will want to provide one of the following:
 *      `BODY_SUN` for heliocentric coordinates,
 *      `BODY_SSB` for solar system barycentric coordinates,
 *      or `BODY_EARTH` for geocentric coordinates. Note that the
 *      gravity simulator does not correct for light travel time;
 *      all state vectors are tied to a Newtonian "instantaneous" time.
 *
 * @param time
 *      The initial time at which to start the simulation.
 *
 * @param numBodies
 *      The number of small bodies to be simulated. This may be any non-negative integer.
 *
 * @param bodyStateArray
 *      An array of initial state vectors (positions and velocities) of the small bodies to be simulated.
 *      The caller must know the positions and velocities of the small bodies at an initial moment in time.
 *      Their positions and velocities are expressed with respect to `originBody`, using J2000 mean equator orientation (EQJ).
 *      Positions are expressed in astronomical units (AU). Velocities are expressed in AU/day.
 *      All the times embedded within the state vectors must be exactly equal to `time`,
 *      or this function will fail with the error `ASTRO_INCONSISTENT_TIMES`.
 *
 * @return
 *      `ASTRO_SUCCESS` on success, with `*sim` set to a non-NULL value. Otherwise an error code with `*sim` set to NULL.
 */
astro_status_t Astronomy_GravSimInit(
    astro_grav_sim_t **simOut,
    astro_body_t originBody,
    astro_time_t time,
    int numBodies,
    const astro_state_vector_t *bodyStateArray)
{
    astro_grav_sim_t *sim;
    astro_status_t status;
    body_grav_calc_t *array;
    int i;

    /* Validate parameters before attempting to allocate memory. */

    if (simOut == NULL)
        return ASTRO_INVALID_PARAMETER;

    *simOut = NULL;

    if (numBodies < 0)
        return ASTRO_INVALID_PARAMETER;

    if (numBodies > 0 && bodyStateArray == NULL)
        return ASTRO_INVALID_PARAMETER;

    if (originBody < BODY_MERCURY || originBody > BODY_SSB)
        return ASTRO_INVALID_BODY;

    /* Verify that all the state vectors are valid and have matching times. */
    for (i = 0; i < numBodies; ++i)
    {
        if (bodyStateArray[i].status != ASTRO_SUCCESS)
            return ASTRO_INVALID_PARAMETER;

        if (bodyStateArray[i].t.tt != time.tt)
            return ASTRO_INCONSISTENT_TIMES;
    }

    *simOut = sim = (astro_grav_sim_t *) calloc(1, sizeof(astro_grav_sim_t));
    if (sim == NULL)
        return ASTRO_OUT_OF_MEMORY;

    sim->originBody = originBody;
    sim->numBodies = numBodies;
    sim->prev = &(sim->endpoint[0]);
    sim->curr = &(sim->endpoint[1]);
    sim->curr->time = time;

    if (numBodies > 0)
    {
        sim->prev->bodies = (body_grav_calc_t *) calloc(numBodies, sizeof(body_grav_calc_t));
        sim->curr->bodies = (body_grav_calc_t *) calloc(numBodies, sizeof(body_grav_calc_t));
        if (sim->prev->bodies == NULL || sim->curr->bodies == NULL)
        {
            status = ASTRO_OUT_OF_MEMORY;
            goto fail;
        }
    }

    /* Remember the initial states of all the bodies as "current". */
    /* Convert from the public type astro_state_t to our internal type body_grav_calc_t. */
    array = sim->curr->bodies;
    for (i = 0; i < numBodies; ++i)
    {
        array[i].tt  = bodyStateArray[i].t.tt;
        array[i].r.x = bodyStateArray[i].x;
        array[i].r.y = bodyStateArray[i].y;
        array[i].r.z = bodyStateArray[i].z;
        array[i].v.x = bodyStateArray[i].vx;
        array[i].v.y = bodyStateArray[i].vy;
        array[i].v.z = bodyStateArray[i].vz;
    }

    /* Calculate the state of the Sun and planets. */
    CalcSolarSystem(sim);

    /*
        We need to do all the physics calculations in barycentric coordinates.
        But the caller can provide the input vectors with respect to any body.
        Correct the input body state vectors for the specified coordinate origin.
    */
    if (originBody != BODY_SSB)
    {
        /* Determine the barycentric state of the origin body. */
        astro_state_vector_t originState = GravSimOriginState(sim);
        if (originState.status != ASTRO_SUCCESS)
        {
            status = originState.status;
            goto fail;
        }

        /* Add barycentric origin to origin-centric body to obtain barycentric body. */
        for (i = 0; i < numBodies; ++i)
        {
            array[i].r.x += originState.x;
            array[i].r.y += originState.y;
            array[i].r.z += originState.z;
            array[i].v.x += originState.vx;
            array[i].v.y += originState.vy;
            array[i].v.z += originState.vz;
        }
    }

    /* Calculate the net acceleration experienced by the small bodies. */
    CalcBodyAccelerations(sim);

    /* To prepare for a possible swap operation, duplicate the current state into the previous state. */
    GravSimDuplicate(sim);

    return ASTRO_SUCCESS;

fail:
    Astronomy_GravSimFree(sim);
    *simOut = NULL;
    return status;
}


/**
 * @brief Advances a gravity simulation by a small time step.
 *
 * @param sim
 *      A simulation object that was created by a prior call to #Astronomy_GravSimInit.
 *
 * @param time
 *      A time that is a small increment away from the current simulation time.
 *      It is up to the developer to figure out an appropriate time increment.
 *      Depending on the trajectories, a smaller or larger increment
 *      may be needed for the desired accuracy. Some experimentation may be needed.
 *      Generally, bodies that stay in the outer Solar System and move slowly can
 *      use larger time steps.  Bodies that pass into the inner Solar System and
 *      move faster will need a smaller time step to maintain accuracy.
 *      The `time` value may be after or before the current simulation time
 *      to move forward or backward in time.
 *
 * @param numBodies
 *      The number of bodies whose state vectors are to be updated.
 *      This is the number of elements in the `bodyStateArray`.
 *      This parameter is passed as a sanity check, and must be equal
 *      to the value passed to #Astronomy_GravSimInit when `sim` was created.
 *
 * @param bodyStateArray
 *      An array big enough to hold `numBodies` state vectors, to receive
 *      the updated positions and velocities of the simulated small bodies.
 *      Alternatively, `bodyStateArray` may be NULL if the output of this
 *      simulation step is not needed. This makes the call slightly faster.
 *
 * @return
 *      `ASTRO_SUCCESS` if the calculation was successful.
 *      Otherwise, an error code if something went wrong, in which case
 *      the simulation should be considered "broken". This means there
 *      is no reliable output in `bodyStateArray` and that no more calculations
 *      can be performed with `sim`.
 */
astro_status_t Astronomy_GravSimUpdate(
    astro_grav_sim_t *sim,
    astro_time_t time,
    int numBodies,
    astro_state_vector_t *bodyStateArray)
{
    terse_vector_t acc;
    double dt;      /* terrestrial time increment */
    int i;

    /*
        The caller's understanding of the number of bodies must match the actual
        array size in `sim`, or we risk corrupting/accessing invalid memory.
    */
    if (numBodies != sim->numBodies)
        return ASTRO_INVALID_PARAMETER;

    dt = time.tt - sim->curr->time.tt;

    if (dt == 0.0)
    {
        /*
            Special case: the time has not changed, so skip the usual physics calculations.
            This allows a way for the caller to query the current state if desired.
            It is also necessary to avoid dividing by `dt` if `dt` is zero.
            To prepare for a possible swap operation, duplicate the current state into the previous state.
        */
        GravSimDuplicate(sim);
    }
    else
    {
        /* Swap the current state and the previous state. Then calculate the new current state. */
        Astronomy_GravSimSwap(sim);

        /* Update the current time. This is the only place we have a full (tt,ut) pair. */
        /* All of the Newtonian dynamics are calculated using tt only. */
        sim->curr->time = time;

        /* Now that sim->time is set, it is safe to call `CalcSolarSystem`. */
        CalcSolarSystem(sim);

        for (i = 0; i < numBodies; ++i)
        {
            /*
                Estimate the positions of the small bodies as if their
                current accelerations apply across the whole time interval.
                approx_pos = pos1 + vel1*dt + (1/2)acc*dt^2
            */
            const body_grav_calc_t *prev = &sim->prev->bodies[i];
            sim->curr->bodies[i].r = UpdatePosition(dt, prev->r, prev->v, prev->a);
        }

        /*
            Calculate the acceleration experienced by the small bodies
            at their respective approximate next locations.
        */
        CalcBodyAccelerations(sim);

        for (i = 0; i < numBodies; ++i)
        {
            const body_grav_calc_t *prev = &sim->prev->bodies[i];
            body_grav_calc_t *curr = &sim->curr->bodies[i];

            /*
                Calculate the average of the acceleration vectors
                experienced by the previous body positions and
                their estimated next positions.
                These become estimates of the mean effective accelerations over the whole interval.
            */
            acc = VecMean(prev->a, curr->a);

            /*
                Refine the estimates of position and velocity at the next time step,
                using the mean acceleration as a better approximation of the
                continuously changing acceleration acting on each body.
            */
            curr->tt = time.tt;
            curr->r = UpdatePosition(dt, prev->r, prev->v, acc);
            curr->v = UpdateVelocity(dt, prev->v, acc);
        }

        /*
            Re-calculate accelerations experienced by each body.
            These will be needed for the next simulation step (if any).
            Also, they will be potentially useful if some day we add
            a function to query the acceleration vectors for the bodies.
        */
        CalcBodyAccelerations(sim);
    }

    /*
        Translate our internal calculations of body positions
        and velocities into state vectors that the caller can understand.
        But if the output buffer `bodyStateArray` is NULL, it means
        the caller wanted us to update the simulation state without
        returning any output.
    */
    if (bodyStateArray != NULL)
    {
        for (i = 0; i < numBodies; ++i)
            bodyStateArray[i] = ExportGravCalc(sim->curr->bodies[i], time);

        if (sim->originBody != BODY_SSB)
        {
            /* Determine the barycentric state of the origin body. */
            astro_state_vector_t originState = GravSimOriginState(sim);
            if (originState.status != ASTRO_SUCCESS)
                return originState.status;

            /* Subtract vectors to convert barycentric states to origin-centric states. */
            for (i = 0; i < numBodies; ++i)
            {
                bodyStateArray[i].x  -= originState.x;
                bodyStateArray[i].y  -= originState.y;
                bodyStateArray[i].z  -= originState.z;
                bodyStateArray[i].vx -= originState.vx;
                bodyStateArray[i].vy -= originState.vy;
                bodyStateArray[i].vz -= originState.vz;
            }
        }
    }

    return ASTRO_SUCCESS;
}


/**
 * @brief Get the position and velocity of a Solar System body included in the simulation.
 *
 * In order to simulate the movement of small bodies through the Solar System,
 * the simulator needs to calculate the state vectors for the Sun and planets.
 *
 * If an application wants to know the positions of one or more of the planets
 * in addition to the small bodies, this function provides a way to obtain
 * their state vectors. This is provided for the sake of efficiency, to avoid
 * redundant calculations.
 *
 * @param sim
 *      A gravity simulator object created by a successful call to #Astronomy_GravSimInit.
 *
 * @param body
 *      The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune.
 *
 * @return
 *      If the given body is part of the set of calculated bodies
 *      (Sun and planets), returns the current time step's state
 *      vector for that body, expressed in the coordinate system
 *      that was specified by the `originBody` parameter to #Astronomy_GravSimInit.
 *      Success is indicated by the returned structure's `status` field holding
 *      `ASTRO_SUCCESS`. Any other `status` value indicates an error, meaning
 *      the returned state vector is invalid.
 */
astro_state_vector_t Astronomy_GravSimBodyState(
    astro_grav_sim_t *sim,
    astro_body_t body)
{
    astro_state_vector_t state;
    body_state_t *bptr;

    bptr = GravSimBodyStatePtr(sim, body);
    if (bptr == NULL)
        return StateVecError(ASTRO_INVALID_BODY, sim->curr->time);

    state = GravSimOriginState(sim);
    if (state.status != ASTRO_SUCCESS)
        return state;

    /*
        Subtract the origin's state from the body's barycentric
        state to get the body in the desired reference frame.
        state.t and state.status are already initialized correctly.
    */

    state.x  = bptr->r.x - state.x;
    state.y  = bptr->r.y - state.y;
    state.z  = bptr->r.z - state.z;
    state.vx = bptr->v.x - state.vx;
    state.vy = bptr->v.y - state.vy;
    state.vz = bptr->v.z - state.vz;

    return state;
}


/**
 * @brief Exchange the current time step with the previous time step.
 *
 * Sometimes it is helpful to "explore" various times near a given
 * simulation time step, while repeatedly returning to the original
 * time step. For example, when backdating a position for light travel
 * time, the caller may wish to repeatedly try different amounts of
 * backdating. When the backdating solver has converged, the caller
 * wants to leave the simulation in its original state.
 *
 * This function allows a single "undo" of a simulation, and does so
 * very efficiently.
 *
 * Usually this function will be called immediately after a matching
 * call to #Astronomy_GravSimUpdate. It has the effect of rolling
 * back the most recent update. If called twice in a row, it reverts
 * the swap and thus has no net effect.
 *
 * #Astronomy_GravSimInit initializes the current state and previous
 * state to be identical. Both states represent the `time` parameter that was
 * passed into the initializer. Therefore, `Astronomy_GravSimSwap` will
 * have no effect from the caller's point of view when passed a simulator
 * that has not yet been updated by a call to #Astronomy_GravSimUpdate.
 *
 * @param sim
 *      A gravity simulator object that was created by a prior call to #Astronomy_GravSimInit.
 */
void Astronomy_GravSimSwap(astro_grav_sim_t *sim)
{
    gravsim_endpoint_t *swap = sim->prev;
    sim->prev = sim->curr;
    sim->curr = swap;
}


/**
 * @brief Returns the time of the current simulation step.
 *
 * @param sim
 *      A gravity simulator object that was created by a prior call to #Astronomy_GravSimInit.
 */
astro_time_t Astronomy_GravSimTime(const astro_grav_sim_t *sim)
{
    return sim->curr->time;
}


/**
 * @brief Returns the number of small bodies represented in this simulation.
 *
 * When a simulation is created by a call to #Astronomy_GravSimInit,
 * the caller specifies the number of small bodies.
 * This function returns that same number, which may be convenient for a caller,
 * so that it does not need to track the body count separately.
 *
 * @param sim
 *      A gravity simulator object that was created by a prior call to #Astronomy_GravSimInit.
 */
int Astronomy_GravSimNumBodies(const astro_grav_sim_t *sim)
{
    return sim->numBodies;
}


/**
 * @brief Returns the body whose center is the coordinate origin that small bodies are referenced to.
 *
 * When a simulation is created by a call to #Astronomy_GravSimInit,
 * the caller specifies an `originBody` to indicate the coordinate origin
 * used to represent the small bodies being simulated.
 * This function returns that same #astro_body_t value.
 *
 * @param sim
 *      A gravity simulator object that was created by a prior call to #Astronomy_GravSimInit.
 */
astro_body_t Astronomy_GravSimOrigin(const astro_grav_sim_t *sim)
{
    return sim->originBody;
}


/**
 * @brief Releases memory allocated to a gravity simulator object.
 *
 * To avoid memory leaks, any successful call to #Astronomy_GravSimInit
 * must be paired with a matching call to `Astronomy_GravSimFree`.
 *
 * @param sim
 *      A gravity simulator object that was created by a prior call to #Astronomy_GravSimInit.
 */
void Astronomy_GravSimFree(astro_grav_sim_t *sim)
{
    if (sim != NULL)
    {
        free(sim->endpoint[0].bodies);
        free(sim->endpoint[1].bodies);
        free(sim);
    }
}


/*------------------ begin Pluto integrator ------------------*/

//$ASTRO_PLUTO_TABLE();

/* FIXFIXFIX - Using a global is not thread-safe. Either add thread-locks or change API to accept a cache pointer. */
static body_segment_t *pluto_cache[PLUTO_NUM_STATES-1];


static int ClampIndex(double frac, int nsteps)
{
    int index = (int) floor(frac);
    if (index < 0)
        return 0;
    if (index >= nsteps)
        return nsteps-1;
    return index;
}


static astro_status_t GetSegment(int *seg_index, body_segment_t *cache[], double tt)
{
    int i;
    body_segment_t reverse;
    body_segment_t *seg;
    major_bodies_t bary;
    double step_tt, ramp;

    if (tt < PlutoStateTable[0].tt || tt > PlutoStateTable[PLUTO_NUM_STATES-1].tt)
    {
        /* We don't bother calculating a segment. Let the caller crawl backward/forward to this time. */
        *seg_index = -1;
        return ASTRO_SUCCESS;
    }

    /* See if we have a segment that straddles the requested time. */
    /* If so, return it. Otherwise, calculate it and return it. */

    *seg_index = ClampIndex((tt - PlutoStateTable[0].tt) / PLUTO_TIME_STEP, PLUTO_NUM_STATES-1);
    if (cache[*seg_index] == NULL)
    {
        /* Allocate memory for the segment (about 11K each). */
        seg = cache[*seg_index] = (body_segment_t *) calloc(1, sizeof(body_segment_t));
        if (seg == NULL)
            return ASTRO_OUT_OF_MEMORY;

        /* Calculate the segment. */
        /* Pick the pair of bracketing body states to fill the segment. */

        /* Each endpoint is exact. */
        seg->step[0] = GravFromState(&bary, &PlutoStateTable[*seg_index]);
        seg->step[PLUTO_NSTEPS-1] = GravFromState(&bary, &PlutoStateTable[*seg_index + 1]);

        /* Simulate forwards from the lower time bound. */
        step_tt = seg->step[0].tt;
        for (i=1; i < PLUTO_NSTEPS-1; ++i)
            seg->step[i] = GravSim(&bary, step_tt += PLUTO_DT, &seg->step[i-1]);

        /* Simulate backwards from the upper time bound. */
        step_tt = seg->step[PLUTO_NSTEPS-1].tt;
        reverse.step[PLUTO_NSTEPS-1] = seg->step[PLUTO_NSTEPS-1];
        for (i=PLUTO_NSTEPS-2; i > 0; --i)
            reverse.step[i] = GravSim(&bary, step_tt -= PLUTO_DT, &reverse.step[i+1]);

        /* Fade-mix the two series so that there are no discontinuities. */
        for (i=PLUTO_NSTEPS-2; i > 0; --i)
        {
            ramp = (double)i / (PLUTO_NSTEPS-1);
            seg->step[i].r = VecRamp(seg->step[i].r, reverse.step[i].r, ramp);
            seg->step[i].v = VecRamp(seg->step[i].v, reverse.step[i].v, ramp);
            seg->step[i].a = VecRamp(seg->step[i].a, reverse.step[i].a, ramp);
        }
    }

    return ASTRO_SUCCESS;
}


static body_grav_calc_t CalcPlutoOneWay(major_bodies_t *bary, const body_state_t *init_state, double target_tt, double dt)
{
    body_grav_calc_t calc;
    int i, n;

    calc = GravFromState(bary, init_state);
    n = (int) ceil((target_tt - calc.tt) / dt);
    for (i=0; i < n; ++i)
        calc = GravSim(bary, (i+1 == n) ? target_tt : (calc.tt + dt), &calc);

    return calc;
}


static astro_status_t CalcPluto(body_state_t *bstate, astro_time_t time, int helio)
{
    terse_vector_t acc, ra, rb, va, vb;
    major_bodies_t bary;
    const body_segment_t *seg;
    int seg_index, left;
    const body_grav_calc_t *s1;
    const body_grav_calc_t *s2;
    body_grav_calc_t calc;
    astro_status_t status;
    double ramp;

    memset(bstate, 0, sizeof(body_state_t));
    bstate->tt = time.tt;

    status = GetSegment(&seg_index, pluto_cache, time.tt);
    if (status != ASTRO_SUCCESS)
        return status;

    if (seg_index < 0)
    {
        /* The target time is outside the year range 0000..4000. */
        /* Calculate it by crawling backward from 0000 or forward from 4000. */
        /* FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed. */
        if (time.tt < PlutoStateTable[0].tt)
            calc = CalcPlutoOneWay(&bary, &PlutoStateTable[0], time.tt, -PLUTO_DT);
        else
            calc = CalcPlutoOneWay(&bary, &PlutoStateTable[PLUTO_NUM_STATES-1], time.tt, +PLUTO_DT);

        bstate->r  = calc.r;
        bstate->v  = calc.v;
    }
    else
    {
        seg = pluto_cache[seg_index];
        left = ClampIndex((time.tt - seg->step[0].tt) / PLUTO_DT, PLUTO_NSTEPS-1);
        s1 = &seg->step[left];
        s2 = &seg->step[left+1];

        /* Find mean acceleration vector over the interval. */
        acc = VecMean(s1->a, s2->a);

        /* Use Newtonian mechanics to extrapolate away from t1 in the positive time direction. */
        ra = UpdatePosition(time.tt - s1->tt, s1->r, s1->v, acc);
        va = UpdateVelocity(time.tt - s1->tt, s1->v, acc);

        /* Use Newtonian mechanics to extrapolate away from t2 in the negative time direction. */
        rb = UpdatePosition(time.tt - s2->tt, s2->r, s2->v, acc);
        vb = UpdateVelocity(time.tt - s2->tt, s2->v, acc);

        /* Use fade in/out idea to blend the two position estimates. */
        ramp = (time.tt - s1->tt)/PLUTO_DT;
        bstate->r = VecRamp(ra, rb, ramp);
        bstate->v = VecRamp(va, vb, ramp);

        if (helio)
            MajorBodyBary(&bary, time.tt);
    }

    if (helio)
    {
        /* Convert barycentric coordinates back to heliocentric coordinates. */
        VecDecr(&bstate->r, bary.Sun.r);
        VecDecr(&bstate->v, bary.Sun.v);
    }

    return ASTRO_SUCCESS;
}

/*------------------ end Pluto integrator ------------------*/


/*---------------------- begin Jupiter moons ----------------------*/

//$ASTRO_JUPITER_MOONS();

static astro_state_vector_t JupiterMoon_elem2pv(astro_time_t time, double mu, const double elem[6])
{
    /* Translation of FORTRAN subroutine ELEM2PV from: */
    /* https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/ */
    astro_state_vector_t state;
    double EE, DE, CE, SE, DLE, RSAM1, ASR, PHI, PSI, X1, Y1, VX1, VY1, F2, P2, Q2, PQ;

    const double A  = elem[0];
    const double AL = elem[1];
    const double K  = elem[2];
    const double H  = elem[3];
    const double Q  = elem[4];
    const double P  = elem[5];

    const double AN = sqrt(mu / (A*A*A));

    EE = AL + K*sin(AL) - H*cos(AL);
    do
    {
        CE = cos(EE);
        SE = sin(EE);
        DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE);
        EE += DE;
    }
    while (fabs(DE) >= 1.0e-12);

    CE = cos(EE);
    SE = sin(EE);
    DLE = H*CE - K*SE;
    RSAM1 = -K*CE - H*SE;
    ASR = 1.0/(1.0 + RSAM1);
    PHI = sqrt(1.0 - K*K - H*H);
    PSI = 1.0/(1.0 + PHI);
    X1 = A*(CE - K - PSI*H*DLE);
    Y1 = A*(SE - H + PSI*K*DLE);
    VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1);
    VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1);
    F2 = 2.0*sqrt(1.0 - Q*Q - P*P);
    P2 = 1.0 - 2.0*P*P;
    Q2 = 1.0 - 2.0*Q*Q;
    PQ = 2.0*P*Q;

    state.x = X1*P2 + Y1*PQ;
    state.y = X1*PQ + Y1*Q2;
    state.z = (Q*Y1 - X1*P)*F2;

    state.vx = VX1*P2 + VY1*PQ;
    state.vy = VX1*PQ + VY1*Q2;
    state.vz = (Q*VY1 - VX1*P)*F2;

    state.t = time;
    state.status = ASTRO_SUCCESS;
    return state;
}

static astro_state_vector_t CalcJupiterMoon(astro_time_t time, int mindex)
{
    /* This is a translation of FORTRAN code by Duriez, Lainey, and Vienne: */
    /* https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/ */

    astro_state_vector_t state;
    int k;
    double arg;
    double elem[6];
    const jupiter_moon_t *m = &JupiterMoonModel[mindex];
    const double t = time.tt + 18262.5;     /* t = time since 1950-01-01T00:00:00Z */

    /* Calculate 6 orbital elements at the given time t. */

    elem[0] = 0.0;
    for (k = 0; k < m->a.nterms; ++k)
    {
        arg = m->a.term[k].phase + (t * m->a.term[k].frequency);
        elem[0] += m->a.term[k].amplitude * cos(arg);
    }

    elem[1] = m->al[0] + (t * m->al[1]);
    for (k = 0; k < m->l.nterms; ++k)
    {
        arg = m->l.term[k].phase + (t * m->l.term[k].frequency);
        elem[1] += m->l.term[k].amplitude * sin(arg);
    }
    elem[1] = fmod(elem[1], PI2);
    if (elem[1] < 0.0)
        elem[1] += PI2;

    elem[2] = elem[3] = 0.0;
    for (k = 0; k < m->z.nterms; ++k)
    {
        arg = m->z.term[k].phase + (t * m->z.term[k].frequency);
        elem[2] += m->z.term[k].amplitude * cos(arg);
        elem[3] += m->z.term[k].amplitude * sin(arg);
    }

    elem[4] = elem[5] = 0.0;
    for (k = 0; k < m->zeta.nterms; ++k)
    {
        arg = m->zeta.term[k].phase + (t * m->zeta.term[k].frequency);
        elem[4] += m->zeta.term[k].amplitude * cos(arg);
        elem[5] += m->zeta.term[k].amplitude * sin(arg);
    }

    /* Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP). */
    state = JupiterMoon_elem2pv(time, m->mu, elem);

    /* Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ). */
    return Astronomy_RotateState(Rotation_JUP_EQJ, state);
}


/**
 * @brief Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.
 *
 * Calculates position and velocity vectors for Jupiter's moons
 * Io, Europa, Ganymede, and Callisto, at the given date and time.
 * The vectors are jovicentric (relative to the center of Jupiter).
 * Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
 * The position components are expressed in astronomical units (AU), and the
 * velocity components are in AU/day.
 *
 * To convert to heliocentric position vectors, call #Astronomy_HelioVector
 * with `BODY_JUPITER` to get Jupiter's heliocentric position, then
 * add the jovicentric positions.
 *
 * Likewise, you can call #Astronomy_GeoVector
 * with `BODY_JUPITER` to convert to geocentric positions.
 * However, you will have to manually correct for light travel time
 * from the Jupiter system to Earth to figure out what time to pass
 * to `Astronomy_JupiterMoons` to get an accurate picture
 * of how Jupiter and its moons look from Earth.
 *
 * @param time  The date and time for which to calculate the position vectors.
 * @return Position vectors of Jupiter's largest 4 moons, as described above.
 */
astro_jupiter_moons_t Astronomy_JupiterMoons(astro_time_t time)
{
    astro_jupiter_moons_t jm;

    jm.io       = CalcJupiterMoon(time, 0);
    jm.europa   = CalcJupiterMoon(time, 1);
    jm.ganymede = CalcJupiterMoon(time, 2);
    jm.callisto = CalcJupiterMoon(time, 3);

    return jm;
}

/*---------------------- end Jupiter moons ----------------------*/


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
 *      Can also be a star defined by #Astronomy_DefineStar.
 * @param time  The date and time for which to calculate the position.
 * @return      A heliocentric position vector of the center of the given body.
 */
astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time)
{
    astro_vector_t vector, earth;
    body_state_t bstate;
    stardef_t *star;

    star = UserDefinedStar(body);
    if (star != NULL)
    {
        astro_spherical_t sphere;
        sphere.lat = star->dec;
        sphere.lon = 15.0 * star->ra;
        sphere.dist = star->dist;
        sphere.status = ASTRO_SUCCESS;
        return Astronomy_VectorFromSphere(sphere, time);
    }

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
        vector.t = time;
        vector.status = CalcPluto(&bstate, time, 1);
        if (vector.status != ASTRO_SUCCESS)
        {
            vector.x = vector.y = vector.z = NAN;
        }
        else
        {
            vector.x = bstate.r.x;
            vector.y = bstate.r.y;
            vector.z = bstate.r.z;
        }
        return vector;

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
 *      A body for which to calculate a heliocentric distance:
 *      the Sun, Moon, any of the planets, or a user-defined star.
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
    stardef_t *star;

    star = UserDefinedStar(body);
    if (star != NULL)
    {
        result.status = ASTRO_SUCCESS;
        result.value = star->dist;
        return result;
    }

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
 * @brief Solve for light travel time of a vector function.
 *
 * When observing a distant object, for example Jupiter as seen from Earth,
 * the amount of time it takes for light to travel from the object to the
 * observer can significantly affect the object's apparent position.
 * This function is a generic solver that figures out how long in the
 * past light must have left the observed object to reach the observer
 * at the specified observation time. It uses a context/function pair
 * as a generic interface that expresses an arbitrary position vector
 * as a function of time.
 *
 * This function repeatedly calls `func`, passing `context` and a series of time
 * estimates in the past. Then `func` must return a relative position vector between
 * the observer and the target. `Astronomy_CorrectLightTravel` keeps calling
 * `func` with more and more refined estimates of the time light must have
 * left the target to arrive at the observer.
 *
 * For common use cases, it is simpler to use #Astronomy_BackdatePosition
 * for calculating the light travel time correction of one body observing another body.
 *
 * For geocentric calculations, #Astronomy_GeoVector also backdates the returned
 * position vector for light travel time, only it returns the observation time in
 * the returned vector's `t` field rather than the backdated time.
 *
 * @param context   Holds any parameters needed by `func`.
 * @param func      Pointer to a function that returns a relative position vector as a function of time.
 * @param time      The observation time for which to solve for light travel delay.
 * @return
 *      The position vector returned by `func` at the solved backdated time.
 *      On success, the vector will hold `ASTRO_SUCCESS` in its `status` field,
 *      the backdated time in its `t` field, along with the apparent relative position.
 *      If an error occurs, `status` will hold an error code and the remaining fields
 *      should be ignored.
 */
astro_vector_t Astronomy_CorrectLightTravel(
    void *context,
    astro_position_func_t func,
    astro_time_t time)
{
    int iter;
    astro_time_t ltime, ltime2;
    astro_vector_t pos;
    double distance, dt;

    ltime = time;
    for (iter = 0; iter < 10; ++iter)
    {
        pos = func(context, ltime);
        if (pos.status != ASTRO_SUCCESS)
            return pos;

        distance = Astronomy_VectorLength(pos);

        /*
            This solver does not support more than one light-day of distance,
            because that would cause convergence problems and inaccurate
            values for stellar aberration angles.
        */
        if (distance > C_AUDAY)
            return VecError(ASTRO_INVALID_PARAMETER, time);

        ltime2 = Astronomy_AddDays(time, -distance/C_AUDAY);
        dt = fabs(ltime2.tt - ltime.tt);
        if (dt < 1.0e-9)        /* 86.4 microseconds */
            return pos;

        ltime = ltime2;
    }
    return VecError(ASTRO_NO_CONVERGE, time);   /* light travel time solver did not converge */
}


/** @cond DOXYGEN_SKIP */
typedef struct
{
    astro_body_t        observerBody;
    astro_body_t        targetBody;
    astro_aberration_t  aberration;
    astro_vector_t      observerPos;          /* used only when aberration == NO_ABERRATION */
}
backdate_context_t;
/** @endcond */


static astro_vector_t BodyPosition(void *context, astro_time_t time)
{
    const backdate_context_t *b = (const backdate_context_t *)context;
    astro_vector_t observerPos, pos;

    if (b->aberration == NO_ABERRATION)
    {
        /* No aberration, so use the pre-calculated initial position of the observer body. */
        observerPos = b->observerPos;
    }
    else
    {
        /*
            The following discussion is worded with the observer body being the Earth,
            which is often the case. However, the same reasoning applies to any observer body
            without loss of generality.

            To include aberration, make a good first-order approximation
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
        observerPos = Astronomy_HelioVector(b->observerBody, time);
    }

    if (observerPos.status != ASTRO_SUCCESS)
        return observerPos;

    pos = Astronomy_HelioVector(b->targetBody, time);
    if (pos.status == ASTRO_SUCCESS)
    {
        /* Convert heliocentric body position to observer-centric position. */
        pos.x -= observerPos.x;
        pos.y -= observerPos.y;
        pos.z -= observerPos.z;
    }
    return pos;
}


/**
 * @brief Solve for light travel time correction of apparent position.
 *
 * When observing a distant object, for example Jupiter as seen from Earth,
 * the amount of time it takes for light to travel from the object to the
 * observer can significantly affect the object's apparent position.
 *
 * This function solves the light travel time correction for the apparent
 * relative position vector of a target body as seen by an observer body
 * at a given observation time.
 *
 * For geocentric calculations, #Astronomy_GeoVector also includes light
 * travel time correction, but the time `t` embedded in its returned vector
 * refers to the observation time, not the backdated time that light left
 * the observed body. Thus `Astronomy_BackdatePosition` provides direct
 * access to the light departure time for callers that need it.
 *
 * For a more generalized light travel correction solver, see #Astronomy_CorrectLightTravel.
 *
 * @param time          The time of observation.
 * @param observerBody  The body to be used as the observation location.
 * @param targetBody    The body to be observed.
 * @param aberration    `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected.
 *
 * @return
 *      On success, the position vector at the solved backdated time.
 *      The returned vector will hold `ASTRO_SUCCESS` in its `status` field,
 *      the backdated time in its `t` field, along with the apparent relative position.
 *      If an error occurs, `status` will hold an error code and the remaining fields should be ignored.
 */
astro_vector_t Astronomy_BackdatePosition(
    astro_time_t time,
    astro_body_t observerBody,
    astro_body_t targetBody,
    astro_aberration_t aberration)
{
    if (UserDefinedStar(targetBody))
    {
        /*
            This is a user-defined star, which must be treated as a special case.
            First, we assume its heliocentric position does not change with time.
            Second, we assume its heliocentric position has already been corrected
            for light-travel time, its coordinates given as it appears on Earth at the present.
            Therefore, no backdating is applied.
        */
        astro_vector_t ovec, tvec, vec;
        double rx, ry, rz, s;
        astro_state_vector_t ostate;

        tvec = Astronomy_HelioVector(targetBody, time);
        if (tvec.status != ASTRO_SUCCESS)
            return tvec;

        switch (aberration)
        {
        case NO_ABERRATION:
            /* Return the star's position as seen from the observer. */
            ovec = Astronomy_HelioVector(observerBody, time);
            if (ovec.status != ASTRO_SUCCESS)
                return ovec;
            vec.x = tvec.x - ovec.x;
            vec.y = tvec.y - ovec.y;
            vec.z = tvec.z - ovec.z;
            vec.t = time;
            vec.status = ASTRO_SUCCESS;
            return vec;

        case ABERRATION:
            /*
                (Observer velocity) - (light vector) = (Aberration-corrected direction to target body).
                Note that this is an approximation, because technically the light vector should
                be measured in barycentric coordinates, not heliocentric. The error is very small.
            */
            ostate = Astronomy_HelioState(observerBody, time);
            if (ostate.status != ASTRO_SUCCESS)
                return VecError(ostate.status, time);

            rx = tvec.x - ostate.x;
            ry = tvec.y - ostate.y;
            rz = tvec.z - ostate.z;
            s = C_AUDAY / sqrt(rx*rx + ry*ry + rz*rz);

            vec.x = rx + ostate.vx/s;
            vec.y = ry + ostate.vy/s;
            vec.z = rz + ostate.vz/s;
            vec.t = time;
            vec.status = ASTRO_SUCCESS;
            return vec;

        default:
            return VecError(ASTRO_INVALID_PARAMETER, time);
        }
    }
    else
    {
        backdate_context_t context;

        context.observerBody = observerBody;
        context.targetBody   = targetBody;
        context.aberration   = aberration;
        switch (aberration)
        {
        case NO_ABERRATION:
            /* Without aberration, we need the observer body position at the observation time only. */
            /* For efficiency, calculate it once and hold onto it, so `BodyPosition` can keep using it. */
            context.observerPos = Astronomy_HelioVector(observerBody, time);
            break;

        case ABERRATION:
            /* With aberration, `BackdatePosition` will calculate the observer body state at different times. */
            /* Therefore, do not waste time calculating it at the observation time. */
            /* Initialize the memory with an explicitly invalid value. */
            context.observerPos = VecError(ASTRO_NOT_INITIALIZED, time);
            break;

        default:
            return VecError(ASTRO_INVALID_PARAMETER, time);
        }

        return Astronomy_CorrectLightTravel(&context, BodyPosition, time);
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
 * Unlike #Astronomy_HelioVector, this function corrects for light travel time.
 * This means the position of the body is "back-dated" by the amount of time it takes
 * light to travel from that body to an observer on the Earth.
 *
 * Also, the position can optionally be corrected for
 * [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
 * causing the apparent direction of the body to be shifted due to transverse
 * movement of the Earth with respect to the rays of light coming from that body.
 *
 * @param body
 *      A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
 *      Can also be a star defined by #Astronomy_DefineStar.
 * @param time          The date and time for which to calculate the position.
 * @param aberration    `ABERRATION` to correct for aberration, or `NO_ABERRATION` to leave uncorrected.
 * @return              A geocentric position vector of the center of the given body.
 */
astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time, astro_aberration_t aberration)
{
    astro_vector_t vector;

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
        /* The moon is so close, aberration and light travel time don't matter. */
        vector = Astronomy_GeoMoon(time);
        break;

    default:
        /* For all other bodies, apply light travel time correction. */
        vector = Astronomy_BackdatePosition(time, BODY_EARTH, body, aberration);
        break;
    }

    vector.t = time;    /* tricky: return the observation time, not the backdated time */
    return vector;
}


/**
 * @brief  Calculates barycentric position and velocity vectors for the given body.
 *
 * Given a body and a time, calculates the barycentric position and velocity
 * vectors for the center of that body at that time.
 * The vectors are expressed in J2000 mean equator coordinates (EQJ).
 *
 * @param body
 *      The celestial body whose barycentric state vector is to be calculated.
 *      Supported values are `BODY_SUN`, `BODY_MOON`, `BODY_EMB`, `BODY_SSB`, and all planets:
 *      `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`, `BODY_JUPITER`,
 *      `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`.
 * @param time
 *      The date and time for which to calculate position and velocity.
 * @return
 *      A structure that contains barycentric position and velocity vectors.
 */
astro_state_vector_t Astronomy_BaryState(astro_body_t body, astro_time_t time)
{
    astro_state_vector_t state;
    major_bodies_t bary;
    body_state_t planet, earth;

    if (body == BODY_SSB)
    {
        /* Trivial case: the solar system barycenter itself. */
        state.status = ASTRO_SUCCESS;
        state.x = state.y = state.z = 0.0;
        state.vx = state.vy = state.vz = 0.0;
        state.t = time;
        return state;
    }

    if (body == BODY_PLUTO)
    {
        astro_status_t status = CalcPluto(&planet, time, 0);
        if (status != ASTRO_SUCCESS)
            return StateVecError(status, time);
        return ExportState(planet, time);
    }

    MajorBodyBary(&bary, time.tt);

    switch (body)
    {
    /* If the caller is asking for one of the major bodies, we can immediately return the answer. */
    case BODY_SUN:      return ExportState(bary.Sun,     time);
    case BODY_JUPITER:  return ExportState(bary.Jupiter, time);
    case BODY_SATURN:   return ExportState(bary.Saturn,  time);
    case BODY_URANUS:   return ExportState(bary.Uranus,  time);
    case BODY_NEPTUNE:  return ExportState(bary.Neptune, time);

    /* Handle the remaining VSOP bodies: Mercury, Venus, Earth, Mars. */
    case BODY_MERCURY:
    case BODY_VENUS:
    case BODY_EARTH:
    case BODY_MARS:
        planet = CalcVsopPosVel(&vsop[body], time.tt);
        /* BarySun + HelioBody = BaryBody */
        state.x  = bary.Sun.r.x + planet.r.x;
        state.y  = bary.Sun.r.y + planet.r.y;
        state.z  = bary.Sun.r.z + planet.r.z;
        state.vx = bary.Sun.v.x + planet.v.x;
        state.vy = bary.Sun.v.y + planet.v.y;
        state.vz = bary.Sun.v.z + planet.v.z;
        state.t  = time;
        state.status = ASTRO_SUCCESS;
        return state;

    case BODY_MOON:
    case BODY_EMB:
        earth = CalcVsopPosVel(&vsop[BODY_EARTH], time.tt);
        if (body == BODY_MOON)
            state = Astronomy_GeoMoonState(time);
        else
            state = Astronomy_GeoEmbState(time);
        state.x  += bary.Sun.r.x + earth.r.x;
        state.y  += bary.Sun.r.y + earth.r.y;
        state.z  += bary.Sun.r.z + earth.r.z;
        state.vx += bary.Sun.v.x + earth.v.x;
        state.vy += bary.Sun.v.y + earth.v.y;
        state.vz += bary.Sun.v.z + earth.v.z;
        return state;

    default:
        return StateVecError(ASTRO_INVALID_BODY, time);
    }
}


/**
 * @brief  Calculates heliocentric position and velocity vectors for the given body.
 *
 * Given a body and a time, calculates the position and velocity
 * vectors for the center of that body at that time, relative to the center of the Sun.
 * The vectors are expressed in J2000 mean equator coordinates (EQJ).
 * If you need the position vector only, it is more efficient to call #Astronomy_HelioVector.
 * The Sun's center is a non-inertial frame of reference. In other words, the Sun
 * experiences acceleration due to gravitational forces, mostly from the larger
 * planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
 * kinetic energy, or other quantities that require a non-accelerating frame
 * of reference, consider using #Astronomy_BaryState instead.
 *
 * @param body
 *      The celestial body whose heliocentric state vector is to be calculated.
 *      Supported values are `BODY_SUN`, `BODY_MOON`, `BODY_EMB`, `BODY_SSB`, and all planets:
 *      `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`, `BODY_JUPITER`,
 *      `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`.
 *      Also allowed to be a user-defined star created by #Astronomy_DefineStar.
 * @param time
 *      The date and time for which to calculate position and velocity.
 * @return
 *      A structure that contains heliocentric position and velocity vectors.
 */
astro_state_vector_t Astronomy_HelioState(astro_body_t body, astro_time_t time)
{
    astro_status_t status;
    astro_state_vector_t state;
    major_bodies_t bary;
    body_state_t planet, earth;

    if (UserDefinedStar(body))
    {
        astro_vector_t vec = Astronomy_HelioVector(body, time);
        state.x = vec.x;
        state.y = vec.y;
        state.z = vec.z;
        state.vx = state.vy = state.vz = 0.0;
        state.t = time;
        state.status = vec.status;
        return state;
    }

    switch (body)
    {
    case BODY_SUN:
        /* Trivial case: the Sun is the origin of the heliocentric frame. */
        state.status = ASTRO_SUCCESS;
        state.x = state.y = state.z = 0.0;
        state.vx = state.vy = state.vz = 0.0;
        state.t = time;
        return state;

    case BODY_SSB:
        /* Calculate the barycentric Sun. Then the negative of that is the heliocentric SSB. */
        MajorBodyBary(&bary, time.tt);
        state.x  = -bary.Sun.r.x;
        state.y  = -bary.Sun.r.y;
        state.z  = -bary.Sun.r.z;
        state.vx = -bary.Sun.v.x;
        state.vy = -bary.Sun.v.y;
        state.vz = -bary.Sun.v.z;
        state.t  = time;
        state.status = ASTRO_SUCCESS;
        return state;

    case BODY_MERCURY:
    case BODY_VENUS:
    case BODY_EARTH:
    case BODY_MARS:
    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
        /* Planets included in the VSOP87 model. */
        planet = CalcVsopPosVel(&vsop[body], time.tt);
        return ExportState(planet, time);

    case BODY_PLUTO:
        status = CalcPluto(&planet, time, 1);
        if (status != ASTRO_SUCCESS)
            return StateVecError(status, time);
        return ExportState(planet, time);

    case BODY_MOON:
    case BODY_EMB:
        earth = CalcVsopPosVel(&vsop[BODY_EARTH], time.tt);
        if (body == BODY_MOON)
            state = Astronomy_GeoMoonState(time);
        else
            state = Astronomy_GeoEmbState(time);
        state.x  += earth.r.x;
        state.y  += earth.r.y;
        state.z  += earth.r.z;
        state.vx += earth.v.x;
        state.vy += earth.v.y;
        state.vz += earth.v.z;
        return state;

    default:
        return StateVecError(ASTRO_INVALID_BODY, time);
    }
}


/**
 * @brief Returns the product of mass and universal gravitational constant of a Solar System body.
 *
 * For problems involving the gravitational interactions of Solar System bodies,
 * it is helpful to know the product GM, where G = the universal gravitational constant
 * and M = the mass of the body. In practice, GM is known to a higher precision than
 * either G or M alone, and thus using the product results in the most accurate results.
 * This function returns the product GM in the units au^3/day^2, or 0 for invalid bodies.
 * The values come from page 10 of a
 * [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).
 *
 * @param body      The body for which to find the GM product.
 * @return          The mass product of the given body in au^3/day^2.
 */
double Astronomy_MassProduct(astro_body_t body)
{
    switch (body)
    {
    case BODY_SUN:      return SUN_GM;
    case BODY_MERCURY:  return MERCURY_GM;
    case BODY_VENUS:    return VENUS_GM;
    case BODY_EARTH:    return EARTH_GM;
    case BODY_MOON:     return MOON_GM;
    case BODY_EMB:      return EARTH_GM + MOON_GM;
    case BODY_MARS:     return MARS_GM;
    case BODY_JUPITER:  return JUPITER_GM;
    case BODY_SATURN:   return SATURN_GM;
    case BODY_URANUS:   return URANUS_GM;
    case BODY_NEPTUNE:  return NEPTUNE_GM;
    case BODY_PLUTO:    return PLUTO_GM;
    default:            return 0.0;         /* invalid body */
    }
}


/**
 * @brief Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.
 *
 * Given a more massive "major" body and a much less massive "minor" body,
 * calculates one of the five Lagrange points in relation to the minor body's
 * orbit around the major body. The parameter `point` is an integer that
 * selects the Lagrange point as follows:
 *
 * 1 = the Lagrange point between the major body and minor body.
 * 2 = the Lagrange point on the far side of the minor body.
 * 3 = the Lagrange point on the far side of the major body.
 * 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
 * 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
 *
 * The function returns the state vector for the selected Lagrange point
 * in J2000 mean equator coordinates (EQJ), with respect to the center of the
 * major body.
 *
 * To calculate Sun/Earth Lagrange points, pass in `BODY_SUN` for `major_body`
 * and `BODY_EMB` (Earth/Moon barycenter) for `minor_body`.
 * For Lagrange points of the Sun and any other planet, pass in that planet
 * (e.g. `BODY_JUPITER`) for `minor_body`.
 * To calculate Earth/Moon Lagrange points, pass in `BODY_EARTH` and `BODY_MOON`
 * for the major and minor bodies respectively.
 *
 * In some cases, it may be more efficient to call #Astronomy_LagrangePointFast,
 * especially when the state vectors have already been calculated, or are needed
 * for some other purpose.
 *
 * @param point         A value 1..5 that selects which of the Lagrange points to calculate.
 * @param time          The time at which the Lagrange point is to be calculated.
 * @param major_body    The more massive of the co-orbiting bodies: `BODY_SUN` or `BODY_EARTH`.
 * @param minor_body    The less massive of the co-orbiting bodies. See main remarks.
 * @return              The position and velocity of the selected Lagrange point with respect to the major body's center.
 */
astro_state_vector_t Astronomy_LagrangePoint(
    int point,
    astro_time_t time,
    astro_body_t major_body,
    astro_body_t minor_body)
{
    astro_state_vector_t major_state, minor_state;
    double major_mass, minor_mass;

    major_mass = Astronomy_MassProduct(major_body);
    if (major_mass <= 0.0)
        return StateVecError(ASTRO_INVALID_BODY, time);

    minor_mass = Astronomy_MassProduct(minor_body);
    if (minor_mass <= 0.0)
        return StateVecError(ASTRO_INVALID_BODY, time);

    /* Calculate the state vectors for the major and minor bodies. */
    if (major_body == BODY_EARTH && minor_body == BODY_MOON)
    {
        /* Use geocentric calculations for more precision. */

        /* The Earth's geocentric state is trivial. */
        major_state.status = ASTRO_SUCCESS;
        major_state.t = time;
        major_state.x = major_state.y = major_state.z = 0.0;
        major_state.vx = major_state.vy = major_state.vz = 0.0;

        minor_state = Astronomy_GeoMoonState(time);
        if (minor_state.status != ASTRO_SUCCESS)
            return minor_state;
    }
    else
    {
        major_state = Astronomy_HelioState(major_body, time);
        if (major_state.status != ASTRO_SUCCESS)
            return major_state;

        minor_state = Astronomy_HelioState(minor_body, time);
        if (minor_state.status != ASTRO_SUCCESS)
            return minor_state;
    }

    return Astronomy_LagrangePointFast(
        point,
        major_state,
        major_mass,
        minor_state,
        minor_mass
    );
}


/**
 * @brief Calculates one of the 5 Lagrange points from body masses and state vectors.
 *
 * Given a more massive "major" body and a much less massive "minor" body,
 * calculates one of the five Lagrange points in relation to the minor body's
 * orbit around the major body. The parameter `point` is an integer that
 * selects the Lagrange point as follows:
 *
 * 1 = the Lagrange point between the major body and minor body.
 * 2 = the Lagrange point on the far side of the minor body.
 * 3 = the Lagrange point on the far side of the major body.
 * 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
 * 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
 *
 * The caller passes in the state vector and mass for both bodies.
 * The state vectors can be in any orientation and frame of reference.
 * The body masses are expressed as GM products, where G = the universal
 * gravitation constant and M = the body's mass. Thus the units for
 * `major_mass` and `minor_mass` must be au^3/day^2.
 * Use #Astronomy_MassProduct to obtain GM values for various solar system bodies.
 *
 * The function returns the state vector for the selected Lagrange point
 * using the same orientation as the state vector parameters `major_state` and `minor_state`,
 * and the position and velocity components are with respect to the major body's center.
 *
 * Consider calling #Astronomy_LagrangePoint, instead of this function, for simpler usage in most cases.
 *
 * @param point         A value 1..5 that selects which of the Lagrange points to calculate.
 * @param major_state   The state vector of the major (more massive) of the pair of bodies.
 * @param major_mass    The mass product GM of the major body.
 * @param minor_state   The state vector of the minor (less massive) of the pair of bodies.
 * @param minor_mass    The mass product GM of the minor body.
 * @return              The position and velocity of the selected Lagrange point with respect to the major body's center.
 */
astro_state_vector_t Astronomy_LagrangePointFast(
    int point,
    astro_state_vector_t major_state,
    double major_mass,
    astro_state_vector_t minor_state,
    double minor_mass)
{
    const double cos_60 = 0.5;
    const double sin_60 = 0.8660254037844386;   /* sqrt(3) / 2 */
    double scale, dx, dy, dz;
    double vx, vy, vz;
    double R2, R, r1, r2, x, deltax, dr1, dr2, numer1, numer2, omega2, accel, deriv;
    astro_state_vector_t  p;

    if (point < 1 || point > 5)
        return StateVecError(ASTRO_INVALID_PARAMETER, major_state.t);

    if (major_state.status != ASTRO_SUCCESS || minor_state.status != ASTRO_SUCCESS)
        return StateVecError(ASTRO_INVALID_PARAMETER, major_state.t);

    if (!isfinite(major_mass) || major_mass <= 0.0)
        return StateVecError(ASTRO_INVALID_PARAMETER, major_state.t);

    if (!isfinite(minor_mass) || minor_mass <= 0.0)
        return StateVecError(ASTRO_INVALID_PARAMETER, major_state.t);

    /* Find the relative position vector <dx, dy, dz>. */
    dx = minor_state.x - major_state.x;
    dy = minor_state.y - major_state.y;
    dz = minor_state.z - major_state.z;
    R2 = (dx*dx + dy*dy + dz*dz);

    /* R = Total distance between the bodies. */
    R = sqrt(R2);

    /* Find the velocity vector <vx, vy, vz>. */
    vx = minor_state.vx - major_state.vx;
    vy = minor_state.vy - major_state.vy;
    vz = minor_state.vz - major_state.vz;

    if (point == 4 || point == 5)
    {
        double nx, ny, nz;
        double ux, uy, uz, U;
        double Dx, Dy, Dz;
        double Ux, Uy, Uz;
        double vert, vrad, vtan;

        /*
            For L4 and L5, we need to find points 60 degrees away from the
            line connecting the two bodies and in the instantaneous orbital plane.
            Define the instantaneous orbital plane as the unique plane that contains
            both the relative position vector and the relative velocity vector.
        */

        /* Take the cross product of position and velocity to find a normal vector <nx, ny, nz>. */
        nx = dy*vz - dz*vy;
        ny = dz*vx - dx*vz;
        nz = dx*vy - dy*vx;

        /* Take the cross product normal*position to get a tangential vector <ux, uy, uz>. */
        ux = ny*dz - nz*dy;
        uy = nz*dx - nx*dz;
        uz = nx*dy - ny*dx;

        /* Convert the tangential direction vector to a unit vector. */
        U = sqrt(ux*ux + uy*uy + uz*uz);
        ux /= U;
        uy /= U;
        uz /= U;

        /* Convert the relative position vector into a unit vector. */
        dx /= R;
        dy /= R;
        dz /= R;

        /* Now we have two perpendicular unit vectors in the orbital plane: 'd' and 'u'. */

        /* Create new unit vectors rotated (+/-)60 degrees from the radius/tangent directions. */
        vert = (point == 4) ? +sin_60 : -sin_60;

        /* Rotated radial vector */
        Dx = cos_60*dx + vert*ux;
        Dy = cos_60*dy + vert*uy;
        Dz = cos_60*dz + vert*uz;

        /* Rotated tangent vector */
        Ux = cos_60*ux - vert*dx;
        Uy = cos_60*uy - vert*dy;
        Uz = cos_60*uz - vert*dz;

        /* Calculate L4/L5 positions relative to the major body. */
        p.x = R * Dx;
        p.y = R * Dy;
        p.z = R * Dz;

        /* Use dot products to find radial and tangential components of the relative velocity. */
        vrad = vx*dx + vy*dy + vz*dz;
        vtan = vx*ux + vy*uy + vz*uz;

        /* Calculate L4/L5 velocities. */
        p.vx = vrad*Dx + vtan*Ux;
        p.vy = vrad*Dy + vtan*Uy;
        p.vz = vrad*Dz + vtan*Uz;
    }
    else
    {
        /*
            Calculate the distances of each body from their mutual barycenter.
            r1 = negative distance of major mass from barycenter (e.g. Sun to the left of barycenter)
            r2 = positive distance of minor mass from barycenter (e.g. Earth to the right of barycenter)
        */
        r1 = -R * (minor_mass / (major_mass + minor_mass));
        r2 = +R * (major_mass / (major_mass + minor_mass));

        /* Calculate the square of the angular orbital speed in [rad^2 / day^2]. */
        omega2 = (major_mass + minor_mass) / (R2*R);

        /*
            Use Newton's Method to numerically solve for the location where
            outward centrifugal acceleration in the rotating frame of reference
            is equal to net inward gravitational acceleration.
            First derive a good initial guess based on approximate analysis.
        */
        if (point == 1 || point == 2)
        {
            scale = (major_mass / (major_mass + minor_mass)) * cbrt(minor_mass / (3.0 * major_mass));
            numer1 = -major_mass;    /* The major mass is to the left of L1 and L2 */
            if (point == 1)
            {
                scale = 1.0 - scale;
                numer2 = +minor_mass;    /* The minor mass is to the right of L1. */
            }
            else
            {
                scale = 1.0 + scale;
                numer2 = -minor_mass;    /* The minor mass is to the left of L2. */
            }
        }
        else /* point == 3 */
        {
            scale = ((7.0/12.0)*minor_mass - major_mass) / (minor_mass + major_mass);
            numer1 = +major_mass;    /* major mass is to the right of L3. */
            numer2 = +minor_mass;    /* minor mass is to the right of L3. */
        }

        /* Iterate Newton's Method until it converges. */
        x = R*scale - r1;
        do
        {
            dr1 = x - r1;
            dr2 = x - r2;
            accel = omega2*x + numer1/(dr1*dr1) + numer2/(dr2*dr2);
            deriv = omega2 - 2*numer1/(dr1*dr1*dr1) - 2*numer2/(dr2*dr2*dr2);
            deltax = accel/deriv;
            x -= deltax;
        }
        while (fabs(deltax/R) > 1.0e-14);
        scale = (x - r1) / R;

        p.x  = scale * dx;
        p.y  = scale * dy;
        p.z  = scale * dz;
        p.vx = scale * vx;
        p.vy = scale * vy;
        p.vz = scale * vz;
    }
    p.t = major_state.t;
    p.status = ASTRO_SUCCESS;
    return p;
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
 * @return              Topocentric equatorial coordinates of the celestial body.
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

    if (time == NULL)
        return EquError(ASTRO_INVALID_PARAMETER);

    /* Calculate the geocentric location of the observer. */
    geo_pos(time, observer, gc_observer);

    /* Calculate the geocentric location of the body. */
    gc = Astronomy_GeoVector(body, *time, aberration);
    if (gc.status != ASTRO_SUCCESS)
        return EquError(gc.status);

    /* Convert geocentric coordinates to topocentric coordinates. */
    j2000[0] = gc.x - gc_observer[0];
    j2000[1] = gc.y - gc_observer[1];
    j2000[2] = gc.z - gc_observer[2];

    switch (equdate)
    {
    case EQUATOR_OF_DATE:
        precession(j2000, *time, FROM_2000, temp);
        nutation(temp, time, FROM_2000, datevect);
        equ = vector2radec(datevect, *time);
        return equ;

    case EQUATOR_J2000:
        equ = vector2radec(j2000, *time);
        return equ;

    default:
        return EquError(ASTRO_INVALID_PARAMETER);
    }
}

/**
 * @brief Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.
 *
 * This function calculates a vector from the center of the Earth to
 * a point on or near the surface of the Earth, expressed in equatorial
 * coordinates. It takes into account the rotation of the Earth at the given
 * time, along with the given latitude, longitude, and elevation of the observer.
 *
 * The caller may pass a value in `equdate` to select either `EQUATOR_J2000`
 * for using J2000 coordinates, or `EQUATOR_OF_DATE` for using coordinates relative
 * to the Earth's equator at the specified time.
 *
 * The returned vector has components expressed in astronomical units (AU).
 * To convert to kilometers, multiply the `x`, `y`, and `z` values by
 * the constant value #KM_PER_AU.
 *
 * The inverse of this function is also available: #Astronomy_VectorObserver.
 *
 * @param time
 *      The date and time for which to calculate the observer's position vector.
 *
 * @param observer
 *      The geographic location of a point on or near the surface of the Earth.
 *
 * @param equdate
 *      Selects the date of the Earth's equator in which to express the equatorial coordinates.
 *      The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by the `time` parameter.
 *      Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `time`
 *      as the orientation.
 *
 * @return
 *      If successful, the returned vector holds `ASTRO_SUCCESS` in its `status` field,
 *      and is an equatorial vector from the center of the Earth to the specified location
 *      on (or near) the Earth's surface. Otherwise, `status` holds an error code.
 */
astro_vector_t Astronomy_ObserverVector(
    astro_time_t *time,
    astro_observer_t observer,
    astro_equator_date_t equdate)
{
    astro_vector_t vec;
    double gast, pos[3], temp[3];

    if (time == NULL)
        return VecError(ASTRO_INVALID_PARAMETER, TimeError());

    gast = Astronomy_SiderealTime(time);
    terra(observer, gast, pos, NULL);

    switch (equdate)
    {
    case EQUATOR_OF_DATE:
        /* 'pos' already contains equator-of-date coordinates. */
        break;

    case EQUATOR_J2000:
        /* Convert 'pos' from equator-of-date to J2000. */
        nutation(pos, time, INTO_2000, temp);
        precession(temp, *time, INTO_2000, pos);
        break;

    default:
        /* This is not a valid value of the 'equdate' parameter. */
        return VecError(ASTRO_INVALID_PARAMETER, *time);
    }

    vec.x = pos[0];
    vec.y = pos[1];
    vec.z = pos[2];
    vec.t = *time;
    vec.status = ASTRO_SUCCESS;
    return vec;
}


/**
 * @brief Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.
 *
 * This function calculates position and velocity vectors of an observer
 * on or near the surface of the Earth, expressed in equatorial
 * coordinates. It takes into account the rotation of the Earth at the given
 * time, along with the given latitude, longitude, and elevation of the observer.
 *
 * The caller may pass a value in `equdate` to select either `EQUATOR_J2000`
 * for using J2000 coordinates, or `EQUATOR_OF_DATE` for using coordinates relative
 * to the Earth's equator at the specified time.
 *
 * The returned position vector has components expressed in astronomical units (AU).
 * To convert to kilometers, multiply the `x`, `y`, and `z` values by
 * the constant value #KM_PER_AU.
 *
 * The returned velocity vector is measured in AU/day.
 *
 * If you need the position only, and not the velocity, #Astronomy_ObserverVector
 * is slightly more efficient than this function.
 *
 * @param time
 *      The date and time for which to calculate the observer's geocentric state vector.
 *
 * @param observer
 *      The geographic location of a point on or near the surface of the Earth.
 *
 * @param equdate
 *      Selects the date of the Earth's equator in which to express the equatorial coordinates.
 *      The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by the `time` parameter.
 *      Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `time`
 *      as the orientation.
 *
 * @return
 *      If successful, the returned state vector holds `ASTRO_SUCCESS` in its `status` field,
 *      and the position (x, y, z) and velocity (vx, vy, vz) vectors are valid.
 *      Otherwise, `status` holds an error code.
 */
astro_state_vector_t Astronomy_ObserverState(
    astro_time_t *time,
    astro_observer_t observer,
    astro_equator_date_t equdate)
{
    astro_state_vector_t state;
    double gast, pos[3], vel[3], postemp[3], veltemp[3];

    if (time == NULL)
        return StateVecError(ASTRO_INVALID_PARAMETER, TimeError());

    gast = Astronomy_SiderealTime(time);
    terra(observer, gast, pos, vel);

    switch (equdate)
    {
    case EQUATOR_OF_DATE:
        /* 'pos' and 'vel' already contain equator-of-date coordinates. */
        break;

    case EQUATOR_J2000:
        /* Convert 'pos' from equator-of-date to J2000. */
        nutation_posvel(pos, vel, time, INTO_2000, postemp, veltemp);
        precession_posvel(postemp, veltemp, *time, INTO_2000, pos, vel);
        break;

    default:
        /* This is not a valid value of the 'equdate' parameter. */
        return StateVecError(ASTRO_INVALID_PARAMETER, *time);
    }

    state.x  = pos[0];
    state.y  = pos[1];
    state.z  = pos[2];
    state.vx = vel[0];
    state.vy = vel[1];
    state.vz = vel[2];
    state.t  = *time;
    state.status = ASTRO_SUCCESS;
    return state;
}


/**
 * @brief Calculates the geographic location corresponding to an equatorial vector.
 *
 * This is the inverse function of #Astronomy_ObserverVector.
 * Given a geocentric equatorial vector, it returns the geographic
 * latitude, longitude, and elevation for that vector.
 *
 * @param vector
 *      The geocentric equatorial position vector for which to find geographic coordinates.
 *      The components are expressed in Astronomical Units (AU).
 *      You can calculate AU by dividing kilometers by the constant #KM_PER_AU.
 *      The time `vector.t` determines the Earth's rotation. The caller must set `vector.t` to a valid time.
 *      The vector is passed by reference (using a pointer) so that nutation calculations
 *      can be cached inside `vector.t` as an optimization.
 *
 * @param equdate
 *      Selects the date of the Earth's equator in which `vector` is expressed.
 *      The caller may select `EQUATOR_J2000` to use the orientation of the Earth's equator
 *      at noon UTC on January 1, 2000, in which case this function corrects for precession
 *      and nutation of the Earth as it was at the moment specified by `vector.t`.
 *      Or the caller may select `EQUATOR_OF_DATE` to use the Earth's equator at `vector.t`
 *      as the orientation.
 *
 * @return
 *      The geographic latitude, longitude, and elevation above sea level
 *      that corresponds to the given equatorial vector.
 */
astro_observer_t Astronomy_VectorObserver(
    astro_vector_t *vector,
    astro_equator_date_t equdate)
{
    double gast;
    double pos1[3];
    double pos2[3];

    gast = Astronomy_SiderealTime(&vector->t);
    pos1[0] = vector->x;
    pos1[1] = vector->y;
    pos1[2] = vector->z;
    if (equdate == EQUATOR_J2000)
    {
        precession(pos1, vector->t, FROM_2000, pos2);
        nutation(pos2, &vector->t, FROM_2000, pos1);
    }
    return inverse_terra(pos1, gast);
}


/**
 * @brief Calculates the gravitational acceleration experienced by an observer on the Earth.
 *
 * This function implements the WGS 84 Ellipsoidal Gravity Formula.
 * The result is a combination of inward gravitational acceleration
 * with outward centrifugal acceleration, as experienced by an observer
 * in the Earth's rotating frame of reference.
 * The resulting value increases toward the Earth's poles and decreases
 * toward the equator, consistent with changes of the weight measured
 * by a spring scale of a fixed mass moved to different latitudes and heights
 * on the Earth.
 *
 * @param latitude
 *      The latitude of the observer in degrees north or south of the equator.
 *      By formula symmetry, positive latitudes give the same answer as negative
 *      latitudes, so the sign does not matter.
 *
 * @param height
 *      The height above the sea level geoid in meters.
 *      No range checking is done; however, accuracy is only valid in the
 *      range 0 to 100000 meters.
 *
 * @return
 *      The effective gravitational acceleration expressed in meters per second squared [m/s^2].
 */
double Astronomy_ObserverGravity(double latitude, double height)
{
    double s = sin(latitude * DEG2RAD);
    double s2 = s*s;
    double g0 = 9.7803253359 * (1.0 + 0.00193185265241*s2) / sqrt(1.0 - 0.00669437999013*s2);
    return g0 * (1.0 - (3.15704e-07 - 2.10269e-09*s2)*height + 7.37452e-14*height*height);
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
    astro_time_t *time,
    astro_observer_t observer,
    double ra,
    double dec,
    astro_refraction_t refraction)
{
    astro_horizon_t hor;
    double latrad, lonrad, decrad, rarad;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double p[3], pz, pn, pw, proj;
    double az, zd;
    double spin_angle;

    if (time == NULL)
    {
        /* The best we can do is return an invalid state. */
        /* It would break external dependencies to expect them to check for errors. */
        hor.altitude = NAN;
        hor.azimuth = NAN;
        hor.ra = NAN;
        hor.dec = NAN;
        return hor;
    }

    latrad = observer.latitude * DEG2RAD;
    lonrad = observer.longitude * DEG2RAD;
    decrad = dec * DEG2RAD;
    rarad = ra * HOUR2RAD;

    double sinlat = sin(latrad);
    double coslat = cos(latrad);
    double sinlon = sin(lonrad);
    double coslon = cos(lonrad);
    double sindc = sin(decrad);
    double cosdc = cos(decrad);
    double sinra = sin(rarad);
    double cosra = cos(rarad);

    /*
        Calculate three mutually perpendicular unit vectors
        in equatorial coordinates: uze, une, uwe.

        uze = The direction of the observer's local zenith (straight up).
        une = The direction toward due north on the observer's horizon.
        uwe = The direction toward due west on the observer's horizon.

        HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.

        The components of these 3 vectors are as follows:
        [0] = x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
        [1] = y = direction from center of Earth toward 90 degrees west longitude on equator.
        [2] = z = direction from center of Earth toward the north pole.
    */

    uze[0] = coslat * coslon;
    uze[1] = coslat * sinlon;
    uze[2] = sinlat;

    une[0] = -sinlat * coslon;
    une[1] = -sinlat * sinlon;
    une[2] = coslat;

    uwe[0] = sinlon;
    uwe[1] = -coslon;
    uwe[2] = 0.0;

    /*
        Correct the vectors uze, une, uwe for the Earth's rotation by calculating
        sidereal time. Call spin() for each uncorrected vector to rotate about
        the Earth's axis to yield corrected unit vectors uz, un, uw.
        Multiply sidereal hours by -15 to convert to degrees and flip eastward
        rotation of the Earth to westward apparent movement of objects with time.
    */

    spin_angle = -15.0 * Astronomy_SiderealTime(time);
    spin(spin_angle, uze, uz);
    spin(spin_angle, une, un);
    spin(spin_angle, uwe, uw);

    /*
        Convert angular equatorial coordinates (RA, DEC) to
        cartesian equatorial coordinates in 'p', using the
        same orientation system as uze, une, uwe.
    */

    p[0] = cosdc * cosra;
    p[1] = cosdc * sinra;
    p[2] = sindc;

    /*
        Use dot products of p with the zenith, north, and west
        vectors to obtain the cartesian coordinates of the body in
        the observer's horizontal orientation system.

        pn = north  component [-1, +1]
        pw = west   component [-1, +1]
        pz = zenith component [-1, +1]
    */

    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];
    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];

    /* proj is the "shadow" of the body vector along the observer's flat ground. */
    proj = hypot(pn, pw);
    if (proj > 0.0)
    {
        /* If the body is not exactly straight up/down, it has an azimuth. */
        /* Invert the angle to produce degrees eastward from north. */
        az = -atan2(pw, pn) * RAD2DEG;
        if (az < 0.0)
            az += 360;
    }
    else
    {
        /* The body is straight up/down, so it does not have an azimuth. */
        /* Report an arbitrary but reasonable value. */
        az = 0.0;
    }

    /* zd = the angle of the body away from the observer's zenith, in degrees. */
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

            proj = hypot(pr[0], pr[1]);
            if (proj > 0)
            {
                hor.ra = RAD2HOUR * atan2(pr[1], pr[0]);
                if (hor.ra < 0.0)
                    hor.ra += 24.0;
            }
            else
            {
                hor.ra = 0.0;
            }
            hor.dec = RAD2DEG * atan2(pr[2], proj);
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
    precession(sun2000, adjusted_time, FROM_2000, stemp);
    nutation(stemp, &adjusted_time, FROM_2000, sun_ofdate);

    /* Convert equatorial coordinates to ecliptic coordinates. */
    true_obliq = DEG2RAD * e_tilt(&adjusted_time).tobl;
    return RotateEquatorialToEcliptic(sun_ofdate, true_obliq, time);
}

/**
 * @brief Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.
 *
 * Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
 * on 1 January 2000), this function converts those coordinates to true ecliptic coordinates
 * that are relative to the plane of the Earth's orbit around the Sun on that date.
 *
 * @param eqj
 *      Equatorial coordinates in the EQJ frame of reference.
 *      You can call #Astronomy_GeoVector to obtain suitable equatorial coordinates.
 *
 * @return
 *      Spherical and vector coordinates expressed in true ecliptic coordinates of date (ECT).
 */
astro_ecliptic_t Astronomy_Ecliptic(astro_vector_t eqj)
{
    earth_tilt_t et;
    double eqj_pos[3];
    double mean_pos[3];
    double eqd_pos[3];

    if (eqj.status != ASTRO_SUCCESS)
        return EclError(eqj.status);

    /* Calculate nutation and obliquity for this time. */
    /* As an optimization, the nutation angles are cached in `time`, */
    /* and reused below when the `nutation` function is called. */
    et = e_tilt(&eqj.t);

    /* Convert mean J2000 equator (EQJ) to true equator of date (EQD). */
    eqj_pos[0] = eqj.x;
    eqj_pos[1] = eqj.y;
    eqj_pos[2] = eqj.z;
    precession(eqj_pos, eqj.t, FROM_2000, mean_pos);
    nutation(mean_pos, &eqj.t, FROM_2000, eqd_pos);

    /* Rotate from EQD to true ecliptic of date (ECT). */
    return RotateEquatorialToEcliptic(eqd_pos, et.tobl * DEG2RAD, eqj.t);
}

/**
 * @brief   Calculates heliocentric ecliptic longitude of a body.
 *
 * This function calculates the angle around the plane of the Earth's orbit
 * of a celestial body, as seen from the center of the Sun.
 * The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
 * in degrees from the true equinox of date. The ecliptic longitude is always in the range [0, 360).
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

static astro_ecliptic_t RotateEquatorialToEcliptic(const double pos[3], double obliq_radians, astro_time_t time)
{
    astro_ecliptic_t ecl;
    double cos_ob, sin_ob;
    double xyproj;

    cos_ob = cos(obliq_radians);
    sin_ob = sin(obliq_radians);

    ecl.vec.status = ASTRO_SUCCESS;
    ecl.vec.t = time;
    ecl.vec.x = +pos[0];
    ecl.vec.y = +pos[1]*cos_ob + pos[2]*sin_ob;
    ecl.vec.z = -pos[1]*sin_ob + pos[2]*cos_ob;

    xyproj = hypot(ecl.vec.x, ecl.vec.y);
    if (xyproj > 0.0)
    {
        ecl.elon = RAD2DEG * atan2(ecl.vec.y, ecl.vec.x);
        if (ecl.elon < 0.0)
            ecl.elon += 360.0;
    }
    else
        ecl.elon = 0.0;

    ecl.elat = RAD2DEG * atan2(ecl.vec.z, xyproj);
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
    return Astronomy_Search(sun_offset, &targetLon, startTime, t2, 0.01);
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
    double q_ut, q_df_dt;
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

        if (QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2, &q_ut, &q_df_dt))
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
    double *out_t, double *out_df_dt)
{
    double Q, R, S;
    double x, u, ru, x1, x2;

    Q = (fb + fa)/2.0 - fm;
    R = (fb - fa)/2.0;
    S = fm;

    if (Q == 0.0)
    {
        /* This is a line, not a parabola. */
        if (R == 0.0)
            return 0;       /* This is a HORIZONTAL line... can't make progress! */
        x = -S / R;
        if (x < -1.0 || x > +1.0)
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
            x = x1;
        }
        else if (-1.0 <= x2 && x2 <= +1.0)
            x = x2;
        else
            return 0;   /* neither root is within bounds */
    }

    *out_t = tm + x*dt;
    *out_df_dt = (2*Q*x + R) / dt;
    return 1;   /* success */
}

static astro_status_t FindSeasonChange(double targetLon, int year, int month, int day, astro_time_t *time)
{
    astro_time_t startTime;
    astro_search_result_t result;

    startTime = Astronomy_MakeTime(year, month, day, 0, 0, 0.0);
    result = Astronomy_SearchSunLongitude(targetLon, startTime, 20.0);
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

    /*
        https://github.com/cosinekitty/astronomy/issues/187
        Solstices and equinoxes drift over long spans of time,
        due to precession of the Earth's axis.
        Therefore, we have to search a wider range of time than
        one might expect. It turns out this has very little
        effect on efficiency, thanks to the quick convergence
        of quadratic interpolation inside Astronomy_Search().
    */

    status = FindSeasonChange(  0, year,  3, 10, &seasons.mar_equinox);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange( 90, year,  6, 10, &seasons.jun_solstice);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange(180, year,  9, 10, &seasons.sep_equinox);
    if (status != ASTRO_SUCCESS) seasons.status = status;

    status = FindSeasonChange(270, year, 12, 10, &seasons.dec_solstice);
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

    return Astronomy_AngleBetween(sv, bv);
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

    angres = Astronomy_PairLongitude(body, BODY_SUN, time);
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
 * @brief Returns one body's ecliptic longitude with respect to another, as seen from the Earth.
 *
 * This function determines where one body appears around the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) as seen from the Earth,
 * relative to the another body's apparent position.
 * The function returns an angle in the half-open range [0, 360) degrees.
 * The value is the ecliptic longitude of `body1` relative to the ecliptic
 * longitude of `body2`.
 *
 * The angle is 0 when the two bodies are at the same ecliptic longitude
 * as seen from the Earth. The angle increases in the prograde direction
 * (the direction that the planets orbit the Sun and the Moon orbits the Earth).
 *
 * When the angle is 180 degrees, it means the two bodies appear on opposite sides
 * of the sky for an Earthly observer.
 *
 * Neither `body1` nor `body2` is allowed to be `BODY_EARTH`.
 * If this happens, the function fails with the error code `ASTRO_EARTH_NOT_ALLOWED`.
 *
 * @param body1
 *      The first body, whose longitude is to be found relative to the second body.
 *
 * @param body2
 *      The second body, relative to which the longitude of the first body is to be found.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @return
 *      On success, the `status` field in the returned structure holds `ASTRO_SUCCESS` and
 *      the `angle` field holds a value in the range [0, 360).
 *      On failure, the `status` field contains some other value indicating an error condition.
 */
astro_angle_result_t Astronomy_PairLongitude(
    astro_body_t body1,
    astro_body_t body2,
    astro_time_t time)
{
    astro_vector_t vector1, vector2;
    astro_ecliptic_t eclip1, eclip2;
    astro_angle_result_t result;

    if (body1 == BODY_EARTH || body2 == BODY_EARTH)
        return AngleError(ASTRO_EARTH_NOT_ALLOWED);

    vector1 = Astronomy_GeoVector(body1, time, NO_ABERRATION);
    eclip1 = Astronomy_Ecliptic(vector1);        /* checks for errors in vector1 */
    if (eclip1.status != ASTRO_SUCCESS)
        return AngleError(eclip1.status);

    vector2 = Astronomy_GeoVector(body2, time, NO_ABERRATION);
    eclip2 = Astronomy_Ecliptic(vector2);        /* checks for errors in vector2 */
    if (eclip2.status != ASTRO_SUCCESS)
        return AngleError(eclip2.status);

    result.status = ASTRO_SUCCESS;
    result.angle = NormalizeLongitude(eclip1.elon - eclip2.elon);
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
    return Astronomy_PairLongitude(BODY_MOON, BODY_SUN, time);
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
 *      The number of days away from `startTime` that limits the time window for the search.
 *      If the value is negative, the search is performed into the past from `startTime`.
 *      Otherwise, the search is performed into the future from `startTime`.
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
        I have seen more than 0.9 days away from the simple prediction.
        To be safe, we take the predicted time of the event and search
        +/-1.5 days around it (a 3-day wide window).
        Return ASTRO_NO_MOON_QUARTER if the final result goes beyond limitDays after startTime.
    */
    const double uncertainty = 1.5;
    astro_func_result_t funcres;
    double ya, est_dt, dt1, dt2;
    astro_time_t t1, t2;

    funcres = moon_offset(&targetLon, startTime);
    if (funcres.status != ASTRO_SUCCESS)
        return SearchError(funcres.status);

    ya = funcres.value;
    if (limitDays < 0.0)
    {
        /* Search backward in time. */
        if (ya < 0.0) ya += 360.0;
        est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
        dt1 = est_dt - uncertainty;
        dt2 = est_dt + uncertainty;
        if (dt2 < limitDays)
            return SearchError(ASTRO_NO_MOON_QUARTER);    /* not possible for moon phase to occur within specified window (too short) */
        if (dt1 < limitDays)
            dt1 = limitDays;
    }
    else
    {
        /* Search forward in time. */
        if (ya > 0.0) ya -= 360.0;
        est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
        dt1 = est_dt - uncertainty;
        dt2 = est_dt + uncertainty;
        if (dt1 > limitDays)
            return SearchError(ASTRO_NO_MOON_QUARTER);    /* not possible for moon phase to occur within specified window (too short) */
        if (dt2 > limitDays)
            dt2 = limitDays;
    }
    t1 = Astronomy_AddDays(startTime, dt1);
    t2 = Astronomy_AddDays(startTime, dt2);
    return Astronomy_Search(moon_offset, &targetLon, t1, t2, 0.1);
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
 * @brief Searches for the time when the center of a body reaches a specified hour angle as seen by an observer on the Earth.
 *
 * The *hour angle* of a celestial body indicates its position in the sky with respect
 * to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
 * The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
 * The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
 * to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
 * the number of hours that have passed since the most recent time that the body has culminated,
 * or reached its highest point.
 *
 * This function searches for the next or previous time a celestial body reaches the given hour angle
 * relative to the date and time specified by `startTime`.
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
 *      The Sun, Moon, any planet other than the Earth,
 *      or a user-defined star that was created by a call to #Astronomy_DefineStar.
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
 * @param direction
 *      The direction in time to perform the search: a positive value
 *      searches forward in time, a negative value searches backward in time.
 *      The function will fail with `ASTRO_INVALID_PARAMETER` if `direction` is zero.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and the other structure fields are valid. Otherwise, `status` holds some other value
 *      that indicates an error condition.
 */
astro_hour_angle_t Astronomy_SearchHourAngleEx(
    astro_body_t body,
    astro_observer_t observer,
    double hourAngle,
    astro_time_t startTime,
    int direction)
{
    int iter = 0;
    astro_time_t time;
    astro_equatorial_t ofdate;
    astro_hour_angle_t result;
    double delta_sidereal_hours, delta_days, gast;

    if (body == BODY_EARTH)
        return HourAngleError(ASTRO_EARTH_NOT_ALLOWED);

    if (hourAngle < 0.0 || hourAngle >= 24.0)
        return HourAngleError(ASTRO_INVALID_PARAMETER);

    if (direction == 0)
        return HourAngleError(ASTRO_INVALID_PARAMETER);

    time = startTime;
    for(;;)
    {
        ++iter;

        /* Calculate Greenwich Apparent Sidereal Time (GAST) at the given time. */
        gast = Astronomy_SiderealTime(&time);

        /* Obtain equatorial coordinates of date for the body. */
        ofdate = Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION);
        if (ofdate.status != ASTRO_SUCCESS)
            return HourAngleError(ofdate.status);

        /* Calculate the adjustment needed in sidereal time */
        /* to bring the hour angle to the desired value. */

        delta_sidereal_hours = fmod((hourAngle + ofdate.ra - observer.longitude/15) - gast, 24.0);
        if (iter == 1)
        {
            /* On the first iteration, always search the requested time direction. */
            if (direction > 0)
            {
                /* Search forward in time. */
                if (delta_sidereal_hours < 0.0)
                    delta_sidereal_hours += 24.0;
            }
            else
            {
                /* Search backward in time. */
                if (delta_sidereal_hours > 0.0)
                    delta_sidereal_hours -= 24.0;
            }
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


/**
 * @brief Finds the hour angle of a body for a given observer and time.
 *
 * The *hour angle* of a celestial body indicates its position in the sky with respect
 * to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
 * The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
 * The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
 * to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
 * the number of hours that have passed since the most recent time that the body has culminated,
 * or reached its highest point.
 *
 * @param body
 *      The body whose observed hour angle is to be found.
 *
 * @param time
 *      The time of the observation.
 *
 * @param observer
 *      The geographic location where the observation takes place.
 *
 * @return astro_func_result_t
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and `value` holds the hour angle in the half-open range [0, 24).
 *      Otherwise, `status` is an error code that indicates failure.
 */
astro_func_result_t Astronomy_HourAngle(astro_body_t body, astro_time_t *time, astro_observer_t observer)
{
    astro_func_result_t result;
    astro_equatorial_t ofdate;
    double gast;

    if (time == NULL)
        return FuncError(ASTRO_INVALID_PARAMETER);

    /* Calculate Greenwich Apparent Sidereal Time (GAST) at the given time. */
    gast = Astronomy_SiderealTime(time);

    /* Obtain equatorial coordinates of date for the body. */
    ofdate = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (ofdate.status != ASTRO_SUCCESS)
        return FuncError(ofdate.status);

    result.value = fmod(observer.longitude/15 + gast - ofdate.ra, 24.0);
    if (result.value < 0.0)
        result.value += 24.0;
    result.status = ASTRO_SUCCESS;

    return result;
}


/** @cond DOXYGEN_SKIP */

typedef struct
{
    astro_body_t        body;
    int                 direction;          // search option: +1 = rise, -1 = set
    astro_observer_t    observer;
    double              body_radius_au;
    double              target_altitude;
}
context_altitude_t;

static const double RISE_SET_DT = 0.42;    /* 10.08 hours: Nyquist-safe for 22-hour period. */

typedef struct
{
    astro_status_t  status;
    astro_time_t    tx;
    astro_time_t    ty;
    double          ax;
    double          ay;
}
ascent_t;

int _AltitudeDiffCallCount;
int _FindAscentMaxRecursionDepth;

/** @endcond */

static astro_func_result_t altitude_diff(void *context, astro_time_t time)
{
    astro_func_result_t result;
    astro_equatorial_t ofdate;
    astro_horizon_t hor;
    double altitude;
    const context_altitude_t *p = (const context_altitude_t *)context;

    ++_AltitudeDiffCallCount;   /* for internal performance testing */

    ofdate = Astronomy_Equator(p->body, &time, p->observer, EQUATOR_OF_DATE, ABERRATION);
    if (ofdate.status != ASTRO_SUCCESS)
        return FuncError(ofdate.status);

    hor = Astronomy_Horizon(&time, p->observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
    altitude = hor.altitude + RAD2DEG*asin(p->body_radius_au / ofdate.dist);
    result.value = p->direction*(altitude - p->target_altitude);
    result.status = ASTRO_SUCCESS;
    return result;
}


static ascent_t AscentError(astro_status_t status)
{
    ascent_t ascent;
    ascent.ax = ascent.ay = NAN;
    ascent.tx = ascent.ty = TimeError();
    ascent.status = status;
    return ascent;
}


static ascent_t FindAscent(
    int depth,
    context_altitude_t *context,
    double max_deriv_alt,
    astro_time_t t1,
    astro_time_t t2,
    double a1,
    double a2)
{
    ascent_t ascent;
    double da, dt, abs_a1, abs_a2;
    astro_time_t tm;
    astro_func_result_t alt;

    /* For internal performance testing. */
    if (depth > _FindAscentMaxRecursionDepth)
        _FindAscentMaxRecursionDepth = depth;

    /* See if we can find any time interval where the altitude-diff function */
    /* rises from non-positive to positive. */
    /* Return ASTRO_SUCCESS if we do, ASTRO_SEARCH_FAILURE if we don't, or some other status for error cases. */

    if (a1 < 0.0 && a2 >= 0.0)
    {
        /* Trivial success case: the endpoints already rise through zero. */
        ascent.status = ASTRO_SUCCESS;
        ascent.tx = t1;
        ascent.ty = t2;
        ascent.ax = a1;
        ascent.ay = a2;
        return ascent;
    }

    if (a1 >= 0.0 && a2 < 0.0)
    {
        /* Trivial failure case: Assume Nyquist condition prevents an ascent. */
        return AscentError(ASTRO_SEARCH_FAILURE);
    }

    if (depth > 17)
    {
        /*
            Safety valve: do not allow unlimited recursion.
            This should never happen if the rest of the logic is working correctly,
            so fail the whole search if it does happen. It's a bug!
        */
        return AscentError(ASTRO_NO_CONVERGE);
    }

    /*
        Both altitudes are on the same side of zero: both are negative, or both are non-negative.
        There could be a convex "hill" or a concave "valley" that passes through zero.
        In polar regions sometimes there is a rise/set or set/rise pair within minutes of each other.
        For example, the Moon can be below the horizon, then the very top of it becomes
        visible (moonrise) for a few minutes, then it moves sideways and down below
        the horizon again (moonset). We want to catch these cases.
        However, for efficiency and practicality concerns, because the rise/set search itself
        has a 0.1 second threshold, we do not worry about rise/set pairs that are less than
        one second apart. These are marginal cases that are rendered highly uncertain
        anyway, due to unpredictable atmospheric refraction conditions (air temperature and pressure).
    */
    dt = (t2.ut - t1.ut) / 2;
    if (dt * SECONDS_PER_DAY < 1.0)
        return AscentError(ASTRO_SEARCH_FAILURE);

    /* Is it possible to reach zero from the altitude that is closer to zero? */
    abs_a1 = fabs(a1);
    abs_a2 = fabs(a2);
    da = (abs_a1 < abs_a2) ? abs_a1 : abs_a2;

    /*
        Without loss of generality, assume |a1| <= |a2|.
        (Reverse the argument in the case |a2| < |a1|.)
        Imagine you have to "drive" from a1 to 0, then back to a2.
        You can't go faster than max_deriv_alt. If you can't reach 0 in half the time,
        you certainly don't have time to reach 0, turn around, and still make your way
        back up to a2 (which is at least as far from 0 than a1 is) in the time interval dt.
        Therefore, the time threshold is half the time interval, or dt/2.
    */
    if (da > max_deriv_alt*(dt / 2))
    {
        /* Prune: the altitude cannot change fast enough to reach zero. */
        return AscentError(ASTRO_SEARCH_FAILURE);
    }

    /* Bisect the time interval and evaluate the altitude at the midpoint. */
    tm = Astronomy_TimeFromDays((t1.ut + t2.ut)/2);
    alt = altitude_diff(context, tm);
    if (alt.status != ASTRO_SUCCESS)
        return AscentError(ASTRO_SEARCH_FAILURE);

    /* Recurse to the left interval. */
    ascent = FindAscent(1+depth, context, max_deriv_alt, t1, tm, a1, alt.value);
    if (ascent.status == ASTRO_SEARCH_FAILURE)
    {
        /* Recurse to the right interval. */
        ascent = FindAscent(1+depth, context, max_deriv_alt, tm, t2, alt.value, a2);
    }

    return ascent;
}


static astro_func_result_t MaxAltitudeSlope(astro_body_t body, double latitude)
{
    astro_func_result_t result;
    double deriv_ra, deriv_dec, latrad;

    if (!isfinite(latitude) || latitude < -90.0 || latitude > +90.0)
    {
        result.value = NAN;
        result.status = ASTRO_INVALID_PARAMETER;
        return result;
    }

    /*
        Calculate the maximum possible rate that this body's altitude
        could change [degrees/day] as seen by this observer.
        First use experimentally determined extreme bounds for this body
        of how much topocentric RA and DEC can ever change per rate of time.
        We need minimum possible d(RA)/dt, and maximum possible magnitude of d(DEC)/dt.
        Conservatively, we round d(RA)/dt down, d(DEC)/dt up.
        Then calculate the resulting maximum possible altitude change rate.
    */

    switch (body)
    {
    case BODY_MOON:
        deriv_ra  = +4.5;
        deriv_dec = +8.2;
        break;

    case BODY_SUN:
        deriv_ra  = +0.8;
        deriv_dec = +0.5;
        break;

    case BODY_MERCURY:
        deriv_ra  = -1.6;
        deriv_dec = +1.0;
        break;

    case BODY_VENUS:
        deriv_ra  = -0.8;
        deriv_dec = +0.6;
        break;

    case BODY_MARS:
        deriv_ra  = -0.5;
        deriv_dec = +0.4;
        break;

    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
    case BODY_PLUTO:
        deriv_ra  = -0.2;
        deriv_dec = +0.2;
        break;

    case BODY_EARTH:
        result.value = NAN;
        result.status = ASTRO_EARTH_NOT_ALLOWED;
        return result;

    default:
        if (UserDefinedStar(body))
        {
            /*
                The minimum allowed heliocentric distance of a user-defined star
                is one light-year. This can cause a tiny amount of parallax (about 0.001 degrees).
                Also, including stellar aberration (22 arcsec = 0.006 degrees), we provide a
                generous safety buffer of 0.008 degrees.
            */
            deriv_ra  = -0.008;
            deriv_dec = +0.008;
            break;
        }
        result.value = NAN;
        result.status = ASTRO_INVALID_BODY;
        return result;
    }

    latrad = DEG2RAD * latitude;
    result.value = fabs(((360.0 / SOLAR_DAYS_PER_SIDEREAL_DAY) - deriv_ra)*cos(latrad)) + fabs(deriv_dec*sin(latrad));
    result.status = isfinite(result.value) ? ASTRO_SUCCESS : ASTRO_INTERNAL_ERROR;
    return result;
}


static astro_search_result_t InternalSearchAltitude(
    astro_body_t body,
    astro_observer_t observer,
    astro_direction_t direction,
    astro_time_t startTime,
    double limitDays,
    double bodyRadiusAu,
    double targetAltitude)
{
    astro_search_result_t search_result;
    astro_func_result_t func_result;
    context_altitude_t context;
    ascent_t ascent;
    astro_time_t t1, t2;
    double a1, a2, max_deriv_alt;

    if (!isfinite(targetAltitude) || targetAltitude < -90.0 || targetAltitude > +90.0)
        return SearchError(ASTRO_INVALID_PARAMETER);

    func_result = MaxAltitudeSlope(body, observer.latitude);
    if (func_result.status != ASTRO_SUCCESS)
        return SearchError(func_result.status);
    max_deriv_alt = func_result.value;

    context.body = body;
    context.direction = (int)direction;
    context.observer = observer;
    context.body_radius_au = bodyRadiusAu;
    context.target_altitude = targetAltitude;

    /* We allow searching forward or backward in time. */
    /* But we want to keep t1 < t2, so we need a few if/else statements. */
    t1 = t2 = startTime;
    func_result = altitude_diff(&context, t2);
    if (func_result.status != ASTRO_SUCCESS)
        return SearchError(func_result.status);
    a1 = a2 = func_result.value;

    for(;;)
    {
        if (limitDays < 0.0)
        {
            t1 = Astronomy_AddDays(t2, -RISE_SET_DT);
            func_result = altitude_diff(&context, t1);
            if (func_result.status != ASTRO_SUCCESS)
                return SearchError(func_result.status);
            a1 = func_result.value;
        }
        else
        {
            t2 = Astronomy_AddDays(t1, +RISE_SET_DT);
            func_result = altitude_diff(&context, t2);
            if (func_result.status != ASTRO_SUCCESS)
                return SearchError(func_result.status);
            a2 = func_result.value;
        }

        ascent = FindAscent(0, &context, max_deriv_alt, t1, t2, a1, a2);
        if (ascent.status == ASTRO_SUCCESS)
        {
            /* We found a time interval [t1, t2] that contains an alt-diff */
            /* rising from negative a1 to non-negative a2. */
            /* Search for the time where the root occurs. */
            search_result = Astronomy_Search(altitude_diff, &context, ascent.tx, ascent.ty, 0.1);
            if (search_result.status == ASTRO_SUCCESS)
            {
                /* Now that we have a solution, we have to check whether it goes outside the time bounds. */
                if (limitDays < 0.0)
                {
                    if (search_result.time.ut < startTime.ut + limitDays)
                        return SearchError(ASTRO_SEARCH_FAILURE);
                }
                else
                {
                    if (search_result.time.ut > startTime.ut + limitDays)
                        return SearchError(ASTRO_SEARCH_FAILURE);
                }
                return search_result;  /* success! */
            }

            /* The search should have succeeded. Something is wrong with FindAscent! */
            return SearchError(ASTRO_INTERNAL_ERROR);
        }
        else if (ascent.status == ASTRO_SEARCH_FAILURE)
        {
            /* There is no ascent in this interval, so keep searching. */
        }
        else
        {
            /* An unexpected error occurred. Fail the search. */
            return SearchError(ascent.status);
        }

        if (limitDays < 0.0)
        {
            if (t1.ut < startTime.ut + limitDays)
                return SearchError(ASTRO_SEARCH_FAILURE);
            t2 = t1;
            a2 = a1;
        }
        else
        {
            if (t2.ut > startTime.ut + limitDays)
                return SearchError(ASTRO_SEARCH_FAILURE);
            t1 = t2;
            a1 = a2;
        }
    }
}


/**
 * @brief Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.
 *
 * This function calculates idealized values of pressure, temperature, and density
 * using the U.S. Standard Atmosphere (1976) model.
 *
 * See:
 * https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml
 * https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
 * https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf
 *
 * @param elevationMeters
 *      The elevation above sea level at which to calculate atmospheric variables.
 *      The value must be at least -500 to +100000, or the function will
 *      fail with status `ASTRO_INVALID_PARAMETER`.
 *
 * @return astro_atmosphere_tp0
 */
astro_atmosphere_t Astronomy_Atmosphere(double elevationMeters)
{
    astro_atmosphere_t atmos;
    const double P0 = 101325.0;     /* pressure at sea level [pascals] */
    const double T0 = 288.15;       /* temperature at sea level [kelvins] */
    const double T1 = 216.65;       /* temperature between 20 km and 32 km [kelvins] */

    /*
        Formulas for air temperature and pressure at a height of `h` meters
        were found at:
        https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml

        These in turn come from:
        1. COESA, U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, DC, 1976.
        2. Jursa, A. S., Ed., Handbook of Geophysics and the Space Environment, Air Force Geophysics Laboratory, 1985.
    */

    if (!isfinite(elevationMeters) || elevationMeters < -500.0 || elevationMeters > 100000.0)
    {
        /* Invalid elevation. */
        atmos.status = ASTRO_INVALID_PARAMETER;
        atmos.pressure = atmos.temperature = atmos.density = NAN;
    }
    else
    {
        if (elevationMeters <= 11000.0)
        {
            atmos.temperature = T0 - 0.0065*elevationMeters;
            atmos.pressure = P0 * pow(T0 / atmos.temperature, -5.25577);
        }
        else if (elevationMeters <= 20000.0)
        {
            atmos.temperature = T1;
            atmos.pressure = 22632.0 * exp(-0.00015768832 * (elevationMeters - 11000.0));
        }
        else
        {
            atmos.temperature = T1 + 0.001*(elevationMeters - 20000.0);
            atmos.pressure = 5474.87 * pow(T1 / atmos.temperature, 34.16319);
        }
        /* The density is calculated relative to the sea level value. */
        /* Using the ideal gas law PV=nRT, we deduce that density is proportional to P/T. */
        atmos.density = (atmos.pressure / atmos.temperature) / (P0 / T0);
        atmos.status = ASTRO_SUCCESS;
    }

    return atmos;
}


static double HorizonDipAngle(
    astro_observer_t observer,
    double metersAboveGround)
{
    double phi, sinphi, cosphi, c, s, ht_km, ach, ash, radius_m;
    double k, dip;

    /* Calculate the effective radius of the Earth at ground level below the observer. */
    /* Correct for the Earth's oblateness. */
    phi = observer.latitude * DEG2RAD;
    sinphi = sin(phi);
    cosphi = cos(phi);
    c = 1.0 / hypot(cosphi, sinphi*EARTH_FLATTENING);
    s = c * (EARTH_FLATTENING * EARTH_FLATTENING);
    ht_km = (observer.height - metersAboveGround) / 1000.0;     /* height of ground above sea level */
    ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
    ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
    radius_m = 1000.0 * hypot(ach*cosphi, ash*sinphi);

    /*
        Correct refraction of a ray of light traveling tangent to the Earth's surface.
        Based on: https://www.largeformatphotography.info/sunmooncalc/SMCalc.js
        which in turn derives from:
        Sweer, John. 1938.  The Path of a Ray of Light Tangent to the Surface of the Earth.
        Journal of the Optical Society of America 28 (September):327-329.
    */

    /* k = refraction index */
    k = 0.175 * pow(1.0 - (6.5e-3/283.15)*(observer.height - (2.0/3.0)*metersAboveGround), 3.256);

    /* Calculate how far below the observer's horizontal plane the observed horizon dips. */
    dip = RAD2DEG * -(sqrt(2*(1 - k)*metersAboveGround / radius_m) / (1 - k));
    return dip;
}



/**
 * @brief Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
 *
 * This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
 * Rise time is when the body first starts to be visible above the horizon.
 * For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
 * Set time is the moment when the body appears to vanish below the horizon.
 * Therefore, this function adjusts for the apparent angular radius of the observed body
 * (significant only for the Sun and Moon).
 *
 * This function corrects for a typical value of atmospheric refraction, which causes celestial
 * bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
 * Astronomy Engine uses a correction of 34 arcminutes. Real-world refraction varies based
 * on air temperature, pressure, and humidity; such weather-based conditions are outside
 * the scope of Astronomy Engine.
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
 *      The Sun, Moon, any planet other than the Earth,
 *      or a user-defined star that was created by a call to #Astronomy_DefineStar.
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
 *      Limits how many days to search for a rise or set time, and defines
 *      the direction in time to search. When `limitDays` is positive, the
 *      search is performed into the future, after `startTime`.
 *      When negative, the search is performed into the past, before `startTime`.
 *      To limit a rise or set time to the same day, you can use a value of 1 day.
 *      In cases where you want to find the next rise or set time no matter how far
 *      in the future (for example, for an observer near the south pole), you can
 *      pass in a larger value like 365.
 *
 * @param metersAboveGround
 *      Usually the observer is located at ground level. Then this parameter
 *      should be zero. But if the observer is significantly higher than ground
 *      level, for example in an airplane, this parameter should be a positive
 *      number indicating how far above the ground the observer is.
 *      An error occurs if `metersAboveGround` is negative.
 *
 * @return
 *      On success, the `status` field in the returned structure contains `ASTRO_SUCCESS`
 *      and the `time` field contains the date and time of the rise or set time as requested.
 *      If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the rise or set
 *      event does not occur within `limitDays` days of `startTime`. This is a normal condition,
 *      not an error. Any other value of `status` indicates an error of some kind.
 */
astro_search_result_t Astronomy_SearchRiseSetEx(
    astro_body_t body,
    astro_observer_t observer,
    astro_direction_t direction,
    astro_time_t startTime,
    double limitDays,
    double metersAboveGround)
{
    double altitude, dip;
    double body_radius_au;
    astro_atmosphere_t atmos;

    if (!isfinite(metersAboveGround) || (metersAboveGround < 0.0))
        return SearchError(ASTRO_INVALID_PARAMETER);

    switch (body)
    {
    case BODY_SUN:  body_radius_au = SUN_RADIUS_AU;                 break;
    case BODY_MOON: body_radius_au = MOON_EQUATORIAL_RADIUS_AU;     break;
    default:        body_radius_au = 0.0;                           break;
    }

    /* Calculate atmospheric density at ground level. */
    atmos = Astronomy_Atmosphere(observer.height - metersAboveGround);
    if (atmos.status != ASTRO_SUCCESS)
        return SearchError(atmos.status);

    /* Calculate the apparent angular dip of the horizon. */
    dip = HorizonDipAngle(observer, metersAboveGround);

    /* Correct refraction for objects near the horizon, using atmospheric density at the ground. */
    altitude = dip - (REFRACTION_NEAR_HORIZON * atmos.density);

    /* Search for the top of the body crossing the corrected altitude angle. */
    return InternalSearchAltitude(body, observer, direction, startTime, limitDays, body_radius_au, altitude);
}


/**
 * @brief Finds the next time the center of a body passes through a given altitude.
 *
 * Finds when the center of the given body ascends or descends through a given
 * altitude angle, as seen by an observer at the specified location on the Earth.
 * By using the appropriate combination of `direction` and `altitude` parameters,
 * this function can be used to find when civil, nautical, or astronomical twilight
 * begins (dawn) or ends (dusk).
 *
 * Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
 * the horizon. To find civil dawn, pass `DIRECTION_RISE` for `direction` and -6 for `altitude`.
 *
 * Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
 * To find civil dusk, pass `DIRECTION_SET` for `direction` and -6 for `altitude`.
 *
 * Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.
 *
 * Astronomical twilight uses -18 degrees as the `altitude` value.
 *
 * By convention for twilight time calculations, the altitude is not corrected for
 * atmospheric refraction. This is because the target altitudes are below the horizon,
 * and refraction is not directly observable.
 *
 * `Astronomy_SearchAltitude` is not intended to find rise/set times of a body for two reasons:
 * (1) Rise/set times of the Sun or Moon are defined by their topmost visible portion, not their centers.
 * (2) Rise/set times are affected significantly by atmospheric refraction.
 * Therefore, it is better to use #Astronomy_SearchRiseSetEx to find rise/set times, which
 * corrects for both of these considerations.
 *
 * `Astronomy_SearchAltitude` will not work reliably for altitudes at or near the body's
 * maximum or minimum altitudes. To find the time a body reaches minimum or maximum altitude
 * angles, use #Astronomy_SearchHourAngleEx.
 *
 * @param body
 *      The Sun, Moon, any planet other than the Earth,
 *      or a user-defined star that was created by a call to #Astronomy_DefineStar.
 *
 * @param observer
 *      The location where observation takes place.
 *      You can create an observer structure by calling #Astronomy_MakeObserver.
 *
 * @param direction
 *      Either `DIRECTION_RISE` to find when the body ascends through the altitude,
 *      or `DIRECTION_SET` for when the body descends through the altitude.
 *
 * @param startTime
 *      The date and time at which to start the search.
 *
 * @param limitDays
 *      Limits how many days to search for the body reaching the altitude angle,
 *      and defines the direction in time to search. When `limitDays` is positive, the
 *      search is performed into the future, after `startTime`.
 *      When negative, the search is performed into the past, before `startTime`.
 *      To limit the search to the same day, you can use a value of 1 day.
 *      In cases where you want to find the altitude event no matter how far
 *      in the future (for example, for an observer near the south pole), you can
 *      pass in a larger value like 365.
 *
 * @param altitude
 *      The desired altitude angle of the body's center above (positive)
 *      or below (negative) the observer's local horizon, expressed in degrees.
 *      Must be in the range [-90, +90].
 *
 * @return
 *      On success, the `status` field in the returned structure contains `ASTRO_SUCCESS`
 *      and the `time` field contains the date and time of the requested altitude event.
 *      If the `status` field contains `ASTRO_SEARCH_FAILURE`, it means the altitude
 *      event does not occur within `limitDays` days of `startTime`. This is a normal condition,
 *      not an error. Any other value of `status` indicates an error of some kind.
 */
astro_search_result_t Astronomy_SearchAltitude(
    astro_body_t body,
    astro_observer_t observer,
    astro_direction_t direction,
    astro_time_t startTime,
    double limitDays,
    double altitude)
{
    return InternalSearchAltitude(body, observer, direction, startTime, limitDays, 0.0, altitude);
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

        phase = Astronomy_AngleBetween(gc, hc);
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
        if (status != ASTRO_SUCCESS)
            return IllumError(status);
        break;
    }

    illum.status = ASTRO_SUCCESS;
    illum.time = time;
    illum.mag = mag;
    illum.phase_angle = phase.angle;
    illum.phase_fraction = (1.0 + cos(DEG2RAD * phase.angle)) / 2.0;
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
    const planet_distance_context_t *pc = (const planet_distance_context_t *) context;
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
    period = Astronomy_PlanetOrbitalPeriod(body);
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

    orbit_period_days = Astronomy_PlanetOrbitalPeriod(body);
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

    skip = 0.25 * Astronomy_PlanetOrbitalPeriod(body);        /* skip 1/4 of an orbit before starting search again */
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
 * this function returns the matrix that reverses that transform.
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
 * @brief Creates an identity rotation matrix.
 *
 * Returns a rotation matrix that has no effect on orientation.
 * This matrix can be the starting point for other operations,
 * such as using a series of calls to #Astronomy_Pivot to
 * create a custom rotation matrix.
 *
 * @return
 *      The identity matrix.
 */
astro_rotation_t Astronomy_IdentityMatrix(void)
{
    astro_rotation_t r;

    r.rot[0][0] = 1.0;  r.rot[1][0] = 0.0;  r.rot[2][0] = 0.0;
    r.rot[0][1] = 0.0;  r.rot[1][1] = 1.0;  r.rot[2][1] = 0.0;
    r.rot[0][2] = 0.0;  r.rot[1][2] = 0.0;  r.rot[2][2] = 1.0;

    r.status = ASTRO_SUCCESS;

    return r;
}

/**
 * @brief Re-orients a rotation matrix by pivoting it by an angle around one of its axes.
 *
 * Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
 * this function pivots the rotation matrix by that angle around that coordinate axis.
 *
 * For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
 * to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
 * of a telescope camera pointed at a given body, you can use `Astronomy_Pivot` twice:
 * (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
 * western axis by the body's altitude angle. The resulting rotation matrix will then
 * reorient ECL coordinates to the orientation of your telescope camera.
 *
 * @param rotation
 *      The input rotation matrix.
 *
 * @param axis
 *      An integer that selects which coordinate axis to rotate around:
 *      0 = x, 1 = y, 2 = z. Any other value will fail with the error code
 *      `ASTRO_INVALID_PARAMETER` in the `status` field of the return value.
 *
 * @param angle
 *      An angle in degrees indicating the amount of rotation around the specified axis.
 *      Positive angles indicate rotation counterclockwise as seen from the positive
 *      direction along that axis, looking towards the origin point of the orientation system.
 *      If `angle` is NAN or infinite, the function will fail with the error code
 *      `ASTRO_INVALID_PARAMETER`. Any finite number of degrees is allowed, but best
 *      precision will result from keeping `angle` in the range [-360, +360].
 *
 * @return
 *      If successful, the return value will have `ASTRO_SUCCESS` in the `status`
 *      field, along with a pivoted rotation matrix. Otherwise, `status` holds
 *      an appropriate error code and the rotation matrix is invalid.
 */
astro_rotation_t Astronomy_Pivot(astro_rotation_t rotation, int axis, double angle)
{
    astro_rotation_t p;
    double radians, c, s;
    int i, j, k;

    /* Check for an invalid input matrix. */
    if (rotation.status != ASTRO_SUCCESS)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    /* Check for an invalid coordinate axis. */
    if (axis < 0 || axis > 2)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    /* Check for an invalid angle value. */
    if (!isfinite(angle))
        return RotationErr(ASTRO_INVALID_PARAMETER);

    radians = angle * DEG2RAD;
    c = cos(radians);
    s = sin(radians);

    /*
        We need to maintain the "right-hand" rule, no matter which
        axis was selected. That means we pick (i, j, k) axis order
        such that the following vector cross product is satisfied:
        i x j = k
    */
    i = (axis + 1) % 3;
    j = (axis + 2) % 3;
    k = axis;

    p.rot[i][i] = c*rotation.rot[i][i] - s*rotation.rot[i][j];
    p.rot[i][j] = s*rotation.rot[i][i] + c*rotation.rot[i][j];
    p.rot[i][k] = rotation.rot[i][k];

    p.rot[j][i] = c*rotation.rot[j][i] - s*rotation.rot[j][j];
    p.rot[j][j] = s*rotation.rot[j][i] + c*rotation.rot[j][j];
    p.rot[j][k] = rotation.rot[j][k];

    p.rot[k][i] = c*rotation.rot[k][i] - s*rotation.rot[k][j];
    p.rot[k][j] = s*rotation.rot[k][i] + c*rotation.rot[k][j];
    p.rot[k][k] = rotation.rot[k][k];

    p.status = ASTRO_SUCCESS;
    return p;
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

    if (vector.status != ASTRO_SUCCESS)
        return SphereError(vector.status);

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
    equ.vec = vector;

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
 *      `lat` holds the refracted altitude angle,
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
 * This function works best near sea level.
 * To correct for higher elevations, call #Astronomy_Atmosphere for that
 * elevation and multiply the refraction angle by the resulting relative density.
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
        /*
            http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
            JPL Horizons says it uses refraction algorithm from
            Meeus "Astronomical Algorithms", 1991, p. 101-102.
            I found the following Go implementation:
            https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
            This is a translation from the function "Saemundsson" there.
            I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
            This is important because the 'refr' formula below goes crazy near hd = -5.11.
        */

        hd = altitude;
        if (hd < -1.0)
            hd = -1.0;

        refr = (1.02 / tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

        if (refraction == REFRACTION_NORMAL && altitude < -1.0)
        {
            /*
                In "normal" mode we gradually reduce refraction toward the nadir
                so that we never get an altitude angle less than -90 degrees.
                When horizon angle is -1 degrees, the factor is exactly 1.
                As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            */
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
 * calculates the negative angular correction to obtain the unrefracted
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
 * @brief Applies a rotation to a vector, yielding a rotated vector.
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
 * @brief Applies a rotation to a state vector, yielding a rotated vector.
 *
 * This function transforms a state vector in one orientation to a vector
 * in another orientation.
 *
 * @param rotation
 *      A rotation matrix that specifies how the orientation of the state vector is to be changed.
 *
 * @param state
 *      The state vector whose orientation is to be changed.
 *      Both the position and velocity components are transformed.
 *
 * @return
 *      A state vector in the orientation specified by `rotation`.
 */
astro_state_vector_t Astronomy_RotateState(astro_rotation_t rotation, astro_state_vector_t state)
{
    astro_state_vector_t target;

    if (rotation.status != ASTRO_SUCCESS || state.status != ASTRO_SUCCESS)
        return StateVecError(ASTRO_INVALID_PARAMETER, state.t);

    target.status = ASTRO_SUCCESS;
    target.t = state.t;

    target.x = rotation.rot[0][0]*state.x + rotation.rot[1][0]*state.y + rotation.rot[2][0]*state.z;
    target.y = rotation.rot[0][1]*state.x + rotation.rot[1][1]*state.y + rotation.rot[2][1]*state.z;
    target.z = rotation.rot[0][2]*state.x + rotation.rot[1][2]*state.y + rotation.rot[2][2]*state.z;

    target.vx = rotation.rot[0][0]*state.vx + rotation.rot[1][0]*state.vy + rotation.rot[2][0]*state.vz;
    target.vy = rotation.rot[0][1]*state.vx + rotation.rot[1][1]*state.vy + rotation.rot[2][1]*state.vz;
    target.vz = rotation.rot[0][2]*state.vx + rotation.rot[1][2]*state.vy + rotation.rot[2][2]*state.vz;

    return target;
}


/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).
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
    static const double c = COS_OBLIQ_2000;
    static const double s = SIN_OBLIQ_2000;
    astro_rotation_t r;

    r.status = ASTRO_SUCCESS;
    r.rot[0][0] = 1.0;  r.rot[1][0] = 0.0;  r.rot[2][0] = 0.0;
    r.rot[0][1] = 0.0;  r.rot[1][1] = +c;   r.rot[2][1] = +s;
    r.rot[0][2] = 0.0;  r.rot[1][2] = -s;   r.rot[2][2] = +c;
    return r;
}

/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).
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
    static const double c = COS_OBLIQ_2000;
    static const double s = SIN_OBLIQ_2000;
    astro_rotation_t r;

    r.status = ASTRO_SUCCESS;
    r.rot[0][0] = 1.0;  r.rot[1][0] = 0.0;  r.rot[2][0] = 0.0;
    r.rot[0][1] = 0.0;  r.rot[1][1] = +c;   r.rot[2][1] = -s;
    r.rot[0][2] = 0.0;  r.rot[1][2] = +s;   r.rot[2][2] = +c;
    return r;
}

/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean equator (EQJ) to equatorial of-date (EQD).
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
astro_rotation_t Astronomy_Rotation_EQJ_EQD(astro_time_t *time)
{
    astro_rotation_t prec, nut;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    prec = precession_rot(*time, FROM_2000);
    nut = nutation_rot(time, FROM_2000);
    return Astronomy_CombineRotation(prec, nut);
}


/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using mean equator at J2000 epoch.
 * Target: ECT = ecliptic system, using true equinox of the specified date/time.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the target orientation.
 *
 * @return
 *      A rotation matrix that converts EQJ to ECT at `time`.
 */
astro_rotation_t Astronomy_Rotation_EQJ_ECT(astro_time_t *time)
{
    astro_rotation_t rot, step;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    rot = Astronomy_Rotation_EQJ_EQD(time);
    step = Astronomy_Rotation_EQD_ECT(time);
    return Astronomy_CombineRotation(rot, step);
}


/**
 * @brief
 *      Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECT = ecliptic system, using true equinox of the specified date/time.
 * Target: EQJ = equatorial system, using mean equator at J2000 epoch.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the target orientation.
 *
 * @return
 *      A rotation matrix that converts ECT to EQJ at `time`.
 */
astro_rotation_t Astronomy_Rotation_ECT_EQJ(astro_time_t *time)
{
    astro_rotation_t rot, step;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    rot = Astronomy_Rotation_ECT_EQD(time);
    step = Astronomy_Rotation_EQD_EQJ(time);
    return Astronomy_CombineRotation(rot, step);
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean equator (EQJ).
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
astro_rotation_t Astronomy_Rotation_EQD_EQJ(astro_time_t *time)
{
    astro_rotation_t prec, nut;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    nut = nutation_rot(time, INTO_2000);
    prec = precession_rot(*time, INTO_2000);
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
astro_rotation_t Astronomy_Rotation_EQD_HOR(astro_time_t *time, astro_observer_t observer)
{
    astro_rotation_t rot;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double spin_angle;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

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

    spin_angle = -15.0 * Astronomy_SiderealTime(time);
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
astro_rotation_t Astronomy_Rotation_HOR_EQD(astro_time_t *time, astro_observer_t observer)
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
 *      A rotation matrix that converts HOR to EQJ at `time` and for `observer`.
 */
astro_rotation_t Astronomy_Rotation_HOR_EQJ(astro_time_t *time, astro_observer_t observer)
{
    astro_rotation_t hor_eqd, eqd_eqj;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    hor_eqd = Astronomy_Rotation_HOR_EQD(time, observer);
    eqd_eqj = Astronomy_Rotation_EQD_EQJ(time);
    return Astronomy_CombineRotation(hor_eqd, eqd_eqj);
}


/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR).
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
astro_rotation_t Astronomy_Rotation_EQJ_HOR(astro_time_t *time, astro_observer_t observer)
{
    astro_rotation_t rot = Astronomy_Rotation_HOR_EQJ(time, observer);
    return Astronomy_InverseRotation(rot);
}


/**
 * @brief
 *      Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).
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
astro_rotation_t Astronomy_Rotation_EQD_ECL(astro_time_t *time)
{
    astro_rotation_t eqd_eqj;
    astro_rotation_t eqj_ecl;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    eqd_eqj = Astronomy_Rotation_EQD_EQJ(time);
    eqj_ecl = Astronomy_Rotation_EQJ_ECL();
    return Astronomy_CombineRotation(eqd_eqj, eqj_ecl);
}


/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean ecliptic (ECL) to equatorial of-date (EQD).
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
astro_rotation_t Astronomy_Rotation_ECL_EQD(astro_time_t *time)
{
    astro_rotation_t rot = Astronomy_Rotation_EQD_ECL(time);
    return Astronomy_InverseRotation(rot);
}

/**
 * @brief
 *      Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR).
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
astro_rotation_t Astronomy_Rotation_ECL_HOR(astro_time_t *time, astro_observer_t observer)
{
    astro_rotation_t ecl_eqd;
    astro_rotation_t eqd_hor;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    ecl_eqd = Astronomy_Rotation_ECL_EQD(time);
    eqd_hor = Astronomy_Rotation_EQD_HOR(time, observer);
    return Astronomy_CombineRotation(ecl_eqd, eqd_hor);
}

/**
 * @brief
 *      Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL).
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
astro_rotation_t Astronomy_Rotation_HOR_ECL(astro_time_t *time, astro_observer_t observer)
{
    astro_rotation_t rot = Astronomy_Rotation_ECL_HOR(time, observer);
    return Astronomy_InverseRotation(rot);
}

/**
 * @brief
 *      Returns a rotation matrix from J2000 mean ecliptic (EQJ) to galactic (GAL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: GAL = galactic system (IAU 1958 definition).
 *
 * @return
 *      A rotation matrix that converts EQJ to GAL.
 */
astro_rotation_t Astronomy_Rotation_EQJ_GAL(void)
{
    astro_rotation_t rot;

    /*
        This rotation matrix was calculated by the following script
        in this same source code repository:
        demo/python/galeqj_matrix.py
    */

    rot.rot[0][0] = -0.0548624779711344;
    rot.rot[0][1] = +0.4941095946388765;
    rot.rot[0][2] = -0.8676668813529025;

    rot.rot[1][0] = -0.8734572784246782;
    rot.rot[1][1] = -0.4447938112296831;
    rot.rot[1][2] = -0.1980677870294097;

    rot.rot[2][0] = -0.4838000529948520;
    rot.rot[2][1] = +0.7470034631630423;
    rot.rot[2][2] = +0.4559861124470794;

    rot.status = ASTRO_SUCCESS;

    return rot;
}

/**
 * @brief
 *      Returns a rotation matrix from ecliptic galactic (GAL) to J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: GAL = galactic system (IAU 1958 definition).
 * Target: EQJ = equatorial system, using the equator at the J2000 epoch.
 *
 * @return
 *      A rotation matrix that converts GAL to EQJ.
 */
astro_rotation_t Astronomy_Rotation_GAL_EQJ(void)
{
    astro_rotation_t rot;

    /*
        This rotation matrix was calculated by the following script
        in this same source code repository:
        demo/python/galeqj_matrix.py
    */

    rot.rot[0][0] = -0.0548624779711344;
    rot.rot[0][1] = -0.8734572784246782;
    rot.rot[0][2] = -0.4838000529948520;

    rot.rot[1][0] = +0.4941095946388765;
    rot.rot[1][1] = -0.4447938112296831;
    rot.rot[1][2] = +0.7470034631630423;

    rot.rot[2][0] = -0.8676668813529025;
    rot.rot[2][1] = -0.1980677870294097;
    rot.rot[2][2] = +0.4559861124470794;

    rot.status = ASTRO_SUCCESS;

    return rot;
}


/**
 * @brief
 *      Returns a rotation matrix from true ecliptic of date (ECT) to equator of date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECT = true ecliptic of date.
 * Target: EQD = equator of date.
 *
 * @param time
 *      The date and time of the ecliptic/equator conversion.
 *
 * @return
 *      A rotation matrix that converts ECT to EQD.
 */
astro_rotation_t Astronomy_Rotation_ECT_EQD(astro_time_t *time)
{
    astro_rotation_t m;
    earth_tilt_t et;
    double tobl, cos_tobl, sin_tobl;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    /* Find true ecliptic obliquity for this time. */
    et = e_tilt(time);
    tobl = et.tobl * DEG2RAD;
    cos_tobl = cos(tobl);
    sin_tobl = sin(tobl);

    /* EQD.x = ECT.x */
    m.rot[0][0] = 1.0;
    m.rot[1][0] = 0.0;
    m.rot[2][0] = 0.0;

    /* EQD.y = +cos*ECT.y - sin*ECT.z */
    m.rot[0][1] = 0.0;
    m.rot[1][1] = +cos_tobl;
    m.rot[2][1] = -sin_tobl;

    /* EQD.z = +sin*ECT.y + cos*ECT.z */
    m.rot[0][2] = 0.0;
    m.rot[1][2] = +sin_tobl;
    m.rot[2][2] = +cos_tobl;

    m.status = ASTRO_SUCCESS;
    return m;
}


/**
 * @brief
 *      Returns a rotation matrix from equator of date (EQD) to true ecliptic of date (ECT).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equator of date.
 * Target: ECT = true ecliptic of date.
 *
 * @param time
 *      The date and time of the equator/ecliptic conversion.
 *
 * @return
 *      A rotation matrix that converts EQD to ECT.
 */
astro_rotation_t Astronomy_Rotation_EQD_ECT(astro_time_t *time)
{
    astro_rotation_t m;
    earth_tilt_t et;
    double tobl, cos_tobl, sin_tobl;

    if (time == NULL)
        return RotationErr(ASTRO_INVALID_PARAMETER);

    /* Find true ecliptic obliquity for this time. */
    et = e_tilt(time);
    tobl = et.tobl * DEG2RAD;
    cos_tobl = cos(tobl);
    sin_tobl = sin(tobl);

    /* ECT.x = EQD.x */
    m.rot[0][0] = 1.0;
    m.rot[1][0] = 0.0;
    m.rot[2][0] = 0.0;

    /* ECT.y = +cos*EQJ.y + sin*EQJ.z */
    m.rot[0][1] = 0.0;
    m.rot[1][1] = +cos_tobl;
    m.rot[2][1] = +sin_tobl;

    /* ECT.z = -sin*EQJ.y + cos*EQJ.z */
    m.rot[0][2] = 0.0;
    m.rot[1][2] = -sin_tobl;
    m.rot[2][2] = +cos_tobl;

    m.status = ASTRO_SUCCESS;
    return m;
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

//$ASTRO_CONSTEL()

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
    astro_spherical_t s2000;
    astro_equatorial_t b1875;
    astro_vector_t vec2000, vec1875;
    double x_ra, x_dec;
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
        rot = Astronomy_Rotation_EQJ_EQD(&time);
        if (rot.status != ASTRO_SUCCESS)
            return ConstelErr(rot.status);

        epoch2000 = Astronomy_TimeFromDays(0.0);
    }

    /* Convert coordinates from J2000 to year 1875. */
    s2000.status = ASTRO_SUCCESS;
    s2000.lon = ra * 15.0;
    s2000.lat = dec;
    s2000.dist = 1.0;
    vec2000 = Astronomy_VectorFromSphere(s2000, epoch2000);
    if (vec2000.status != ASTRO_SUCCESS)
        return ConstelErr(vec2000.status);

    vec1875 = Astronomy_RotateVector(rot, vec2000);
    if (vec1875.status != ASTRO_SUCCESS)
        return ConstelErr(vec1875.status);

    b1875 = Astronomy_EquatorFromVector(vec1875);
    if (b1875.status != ASTRO_SUCCESS)
        return ConstelErr(b1875.status);

    /* Convert DEC from degrees, and RA from hours, to compact angle units used in the ContelBounds table. */
    x_ra = (24.0 * 15.0) * b1875.ra;
    x_dec = 24.0 * b1875.dec;

    /* Search for the constellation using the B1875 coordinates. */
    c = -1;     /* constellation not (yet) found */
    for (i=0; i < NUM_CONSTEL_BOUNDARIES; ++i)
    {
        const constel_boundary_t *b = &ConstelBounds[i];
        if ((b->dec_lo <= x_dec) && (b->ra_hi > x_ra) && (b->ra_lo <= x_ra))
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
    eclipse.obscuration = eclipse.sd_penum = eclipse.sd_partial = eclipse.sd_total = NAN;
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


static double Obscuration(      /* returns area of intersection of the two discs, divided by area of first disc. */
    double a,   /* radius of first disc */
    double b,   /* radius of second disc */
    double c)   /* distance between the centers of the discs */
{
    double x;   /* Horizontal location of intersection point on both circumferences */
    double y;   /* Positive vertical location of intersection point on both circumferences */
    double radicand, lens1, lens2, obs;

    if (a <= 0.0 || b <= 0.0)
        return 0.0;     /* invalid radius */

    if (c < 0.0)
        return 0.0;     /* invalid distance between centers */

    if (c >= a + b)
        return 0.0;     /* the discs are too far apart to have any overlapping area */

    if (c == 0.0)
    {
        /* The discs have a common center. Therefore, one disc is inside the other. */
        return (a <= b) ? 1.0 : (b*b)/(a*a);
    }

    x = (a*a - b*b + c*c) / (2 * c);

    radicand = a*a - x*x;
    if (radicand <= 0.0)
    {
        /* The circumferences do not intersect, or are tangent. */
        /* We already ruled out the case of non-overlapping discs. */
        /* Therefore, one disc is inside the other. */
        return (a <= b) ? 1.0 : (b*b)/(a*a);
    }

    /* The discs overlap fractionally in a pair of lens-shaped areas. */

    y = sqrt(radicand);

    /* Return the overlapping fractional area. */
    /* There are two lens-shaped areas, one to the left of x, the other to the right of x. */
    /* Each part is calculated by subtracting a triangular area from a sector's area. */
    lens1 = a*a*acos(x/a) - x*y;
    lens2 = b*b*acos((c-x)/b) - (c-x)*y;

    /* Find the fractional area with respect to the first disc. */
    obs = (lens1 + lens2) / (PI*a*a);
    return obs;
}


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
    g = Astronomy_GeoVector(body, time, ABERRATION);
    if (g.status != ASTRO_SUCCESS)
        return ShadowError(g.status);

    /* Calculate light-travel-corrected vector from Earth to Sun. */
    e = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (e.status != ASTRO_SUCCESS)
        return ShadowError(e.status);

    /* Deduce light-travel-corrected vector from Sun to planet. */
    p.status = ASTRO_SUCCESS;
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

    /* Find geocentric Sun with aberration correction. */
    e = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (e.status != ASTRO_SUCCESS)
        return ShadowError(e.status);

    /* Convert geocentric Sun to heliocentric Earth. */
    /* Thus `e` points in the direction of sunlight heading toward the Earth's center. */
    e.x = -e.x;
    e.y = -e.y;
    e.z = -e.z;

    /* Find geocentric Moon. */
    m = Astronomy_GeoMoon(time);    /* This function never fails; no need to check return value */

    return CalcShadow(EARTH_ECLIPSE_RADIUS_KM, time, m, e);
}


static shadow_t MoonShadow(astro_time_t time)
{
    /* This function helps find when the Moon's shadow falls upon the Earth. */

    astro_vector_t s, e, m;

    /* Calculate geocentric Sun with aberration correction. */
    s = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (s.status != ASTRO_SUCCESS)
        return ShadowError(s.status);

    m = Astronomy_GeoMoon(time);    /* geocentric Moon */

    /* Calculate lunacentric Earth. */
    e.status = m.status;
    e.x = -m.x;
    e.y = -m.y;
    e.z = -m.z;
    e.t = m.t;

    /* Convert geocentric moon to heliocentric Moon. */
    m.x -= s.x;
    m.y -= s.y;
    m.z -= s.z;

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
    shadow_func_t shadowfunc = (shadow_func_t) context;

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

    result = Astronomy_Search(shadow_distance_slope, (void *)EarthShadow, t1, t2, 1.0);
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

    result = Astronomy_Search(shadow_distance_slope, (void *)MoonShadow, t1, t2, 1.0);
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
    const planet_shadow_context_t *p = (const planet_shadow_context_t *) context;

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
    const shadow_context_t *p = (const shadow_context_t *) context;
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
                eclipse.obscuration = 0.0;
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
                        eclipse.obscuration = 1.0;
                        eclipse.sd_total = ShadowSemiDurationMinutes(shadow.time, shadow.k - MOON_MEAN_RADIUS_KM, eclipse.sd_partial);
                        if (eclipse.sd_total <= 0.0)
                            return LunarEclipseError(ASTRO_SEARCH_FAILURE);
                    }
                    else
                    {
                        /* For lunar eclipses, we calculate the fraction of the Moon's disc covered by the Earth's umbra. */
                        eclipse.obscuration = Obscuration(MOON_MEAN_RADIUS_KM, shadow.k, shadow.r);
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
    eclipse.obscuration = eclipse.distance = eclipse.latitude = eclipse.longitude = NAN;

    return eclipse;
}

/* The umbra radius tells us what kind of eclipse the observer sees. */
/* If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular. */
/* HACK: I added a tiny bias (14 meters) to match Espenak test data. */
#define EclipseKindFromUmbra(k)     (((k) > 0.014) ? ECLIPSE_TOTAL : ECLIPSE_ANNULAR)

static double SolarEclipseObscuration(
    astro_vector_t hm,     /* heliocentric Moon */
    astro_vector_t lo)     /* lunacentric observer */
{
    astro_vector_t ho;     /* heliocentric observer */
    astro_angle_result_t sun_moon_separation;
    double obscuration, sun_au, sun_radius, moon_radius;

    /* Find heliocentric observer. */
    ho.status = ASTRO_SUCCESS;
    ho.t = lo.t;
    ho.x = hm.x + lo.x;
    ho.y = hm.y + lo.y;
    ho.z = hm.z + lo.z;

    /* Find the distance from the Sun's center to the observer. */
    sun_au = Astronomy_VectorLength(ho);

    /* Calculate the apparent angular radius of the Sun for the observer. */
    sun_radius = asin(SUN_RADIUS_AU / sun_au);

    /* Calculate the apparent angular radius of the Moon for the observer. */
    moon_radius = asin(MOON_POLAR_RADIUS_AU / Astronomy_VectorLength(lo));

    /* Calculate the apparent angular separation between the Sun's center and the Moon's center. */
    sun_moon_separation = Astronomy_AngleBetween(lo, ho);

    if (sun_moon_separation.status != ASTRO_SUCCESS)
        return -1.0;    /* internal error! should never happen. */

    obscuration = Obscuration(sun_radius, moon_radius, sun_moon_separation.angle * DEG2RAD);

    /* HACK: In marginal cases, we need to clamp obscuration to less than 1.0. */
    if (obscuration > 0.9999)
        obscuration = 0.9999;

    return obscuration;
}

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
    rot = Astronomy_Rotation_EQJ_EQD(&shadow.time);
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
        proj = hypot(px, py) * (EARTH_FLATTENING * EARTH_FLATTENING);
        if (proj == 0.0)
            eclipse.latitude = (pz > 0.0) ? +90.0 : -90.0;
        else
            eclipse.latitude = RAD2DEG * atan(pz / proj);

        /* Adjust longitude for Earth's rotation at the given UT. */
        gast = Astronomy_SiderealTime(&eclipse.peak);
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
        if (eclipse.kind == ECLIPSE_TOTAL)
            eclipse.obscuration = 1.0;
        else
            eclipse.obscuration = SolarEclipseObscuration(shadow.dir, o);
    }
    else
    {
        /* This is a partial solar eclipse. It does not make practical sense to calculate obscuration. */
        /* Anyone who wants obscuration should use Astronomy_SearchLocalSolarEclipse for a specific location on the Earth. */
        eclipse.obscuration = NAN;
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
    eclipse.obscuration = NAN;

    eclipse.partial_begin = EclipseEventError();
    eclipse.total_begin   = EclipseEventError();
    eclipse.peak          = EclipseEventError();
    eclipse.total_end     = EclipseEventError();
    eclipse.partial_end   = EclipseEventError();

    return eclipse;
}


static shadow_t LocalMoonShadow(astro_time_t time, astro_observer_t observer)
{
    astro_vector_t s, o, m;
    double pos[3];

    /* Calculate observer's geocentric position. */
    /* For efficiency, do this first, to populate the earth rotation parameters in 'time'. */
    /* That way they can be recycled instead of recalculated. */
    geo_pos(&time, observer, pos);

    /* Calculate geocentric Sun with aberration correction. */
    s = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (s.status != ASTRO_SUCCESS)
        return ShadowError(s.status);

    m = Astronomy_GeoMoon(time);    /* geocentric Moon */

    /* Calculate lunacentric location of an observer on the Earth's surface. */
    o.status = m.status;
    o.x = pos[0] - m.x;
    o.y = pos[1] - m.y;
    o.z = pos[2] - m.z;
    o.t = m.t;

    /* Convert geocentric moon to heliocentric Moon. */
    m.x -= s.x;
    m.y -= s.y;
    m.z -= s.z;

    return CalcShadow(MOON_MEAN_RADIUS_KM, time, o, m);
}


static astro_func_result_t local_shadow_distance_slope(void *context, astro_time_t time)
{
    const double dt = 1.0 / 86400.0;
    astro_time_t t1, t2;
    astro_func_result_t result;
    shadow_t shadow1, shadow2;
    const astro_observer_t *observer = (const astro_observer_t *) context;

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
    const eclipse_transition_t *trans = (const eclipse_transition_t *) context;
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
        if (eclipse.kind == ECLIPSE_TOTAL)
            eclipse.obscuration = 1.0;
        else
            eclipse.obscuration = SolarEclipseObscuration(shadow.dir, shadow.target);
    }
    else
    {
        eclipse.total_begin = eclipse.total_end = EclipseEventError();
        eclipse.kind = ECLIPSE_PARTIAL;
        eclipse.obscuration = SolarEclipseObscuration(shadow.dir, shadow.target);
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
    const planet_shadow_context_t *p = (const planet_shadow_context_t *) context;

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
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and the other fields are as documented in #astro_transit_t.
 *      Otherwise, `status` holds an error code and the other structure members are undefined.
 */
astro_transit_t Astronomy_NextTransit(astro_body_t body, astro_time_t prevTransitTime)
{
    astro_time_t startTime;

    startTime = Astronomy_AddDays(prevTransitTime, 100.0);
    return Astronomy_SearchTransit(body, startTime);
}


static astro_node_event_t NodeError(astro_status_t status)
{
    astro_node_event_t node;

    node.status = status;
    node.time = TimeError();
    node.kind = INVALID_NODE;

    return node;
}

static astro_func_result_t MoonNodeSearchFunc(void *context, astro_time_t time)
{
    astro_func_result_t result;
    astro_spherical_t eclip;
    astro_node_kind_t kind = *((astro_node_kind_t *)context);

    eclip = Astronomy_EclipticGeoMoon(time);

    result.value = eclip.lat * (double)kind;
    result.status = ASTRO_SUCCESS;
    return result;
}

static const double MOON_NODE_STEP_DAYS = +10.0;    /* a safe number of days to step without missing a Moon node */

/**
 * @brief Searches for a time when the Moon's center crosses through the ecliptic plane.
 *
 * Searches for the first ascending or descending node of the Moon after `startTime`.
 * An ascending node is when the Moon's center passes through the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) from south to north.
 * A descending node is when the Moon's center passes through the ecliptic plane
 * from north to south. Nodes indicate possible times of solar or lunar eclipses,
 * if the Moon also happens to be in the correct phase (new or full, respectively).
 *
 * Call `Astronomy_SearchMoonNode` to find the first of a series of nodes.
 * Then call #Astronomy_NextMoonNode to find as many more consecutive nodes as desired.
 *
 * @param startTime
 *      The date and time for starting the search for an ascending or descending node of the Moon.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and the other fields are as documented in #astro_node_event_t.
 *      Otherwise, `status` holds an error code and the other structure members are undefined.
 */
astro_node_event_t Astronomy_SearchMoonNode(astro_time_t startTime)
{
    astro_node_event_t node;
    astro_time_t time1, time2;
    astro_spherical_t eclip1, eclip2;
    astro_node_kind_t kind;
    astro_search_result_t result;

    /* Start at the given moment in time and sample the Moon's ecliptic latitude. */
    /* Step 10 days at a time, searching for an interval where that latitude crosses zero. */
    time1 = startTime;
    eclip1 = Astronomy_EclipticGeoMoon(time1);    /* never returns a failure code */

    for(;;)
    {
        time2 = Astronomy_AddDays(time1, MOON_NODE_STEP_DAYS);
        eclip2 = Astronomy_EclipticGeoMoon(time2);      /* never returns a failure code */
        if (eclip1.lat * eclip2.lat <= 0.0)
        {
            /* There is a node somewhere inside this closed time interval. */
            /* Figure out whether it is an ascending node or a descending node. */
            kind = (eclip2.lat > eclip1.lat) ? ASCENDING_NODE : DESCENDING_NODE;
            result = Astronomy_Search(MoonNodeSearchFunc, &kind, time1, time2, 1.0);
            if (result.status != ASTRO_SUCCESS)
                return NodeError(result.status);

            node.status = ASTRO_SUCCESS;
            node.time = result.time;
            node.kind = kind;
            return node;
        }
        time1 = time2;
        eclip1 = eclip2;
    }
}


/**
 * @brief Searches for the next time when the Moon's center crosses through the ecliptic plane.
 *
 * Call #Astronomy_SearchMoonNode to find the first of a series of nodes.
 * Then call `Astronomy_NextMoonNode` to find as many more consecutive nodes as desired.
 *
 * @param prevNode
 *      The previous node found from calling #Astronomy_SearchMoonNode or `Astronomy_NextMoonNode`.
 *
 * @return
 *      If successful, the `status` field in the returned structure holds `ASTRO_SUCCESS`
 *      and the other fields are as documented in #astro_node_event_t.
 *      Otherwise, `status` holds an error code and the other structure members are undefined.
 */
astro_node_event_t Astronomy_NextMoonNode(astro_node_event_t prevNode)
{
    astro_time_t time;
    astro_node_event_t node;

    if (prevNode.status != ASTRO_SUCCESS)
        return NodeError(ASTRO_INVALID_PARAMETER);

    if (prevNode.kind != ASCENDING_NODE && prevNode.kind != DESCENDING_NODE)
        return NodeError(ASTRO_INVALID_PARAMETER);

    time = Astronomy_AddDays(prevNode.time, MOON_NODE_STEP_DAYS);
    node = Astronomy_SearchMoonNode(time);
    if (node.status == ASTRO_SUCCESS)
    {
        /* Verify nodes are alternating as expected. */
        if (prevNode.kind == ASCENDING_NODE)
        {
            if (node.kind != DESCENDING_NODE)
                return NodeError(ASTRO_INTERNAL_ERROR);
        }
        else
        {
            if (node.kind != ASCENDING_NODE)
                return NodeError(ASTRO_INTERNAL_ERROR);
        }
    }
    return node;
}


/**
 * @brief Frees up all dynamic memory allocated by Astronomy Engine.
 *
 * Astronomy Engine uses dynamic memory allocation in only one place:
 * it makes calculation of Pluto's orbit more efficient by caching 11 KB
 * segments and recycling them. To force purging this cache and
 * freeing all the dynamic memory, you can call this function at any time.
 * It is always safe to call, although it will slow down the very next
 * calculation of Pluto's position for a nearby time value.
 * Calling this function before your program exits is optional, but
 * it will be helpful for leak-checkers like valgrind.
 */
void Astronomy_Reset(void)
{
    int i;
    for (i=0; i < PLUTO_NUM_STATES-1; ++i)
    {
        free(pluto_cache[i]);
        pluto_cache[i] = NULL;
    }
}


static astro_axis_t EarthRotationAxis(astro_time_t *time)
{
    astro_axis_t axis;
    double pos1[3];
    double pos2[3];
    astro_equatorial_t equ;

    if (time == NULL)
        return AxisErr(ASTRO_INVALID_PARAMETER, TimeError());

    /*
        Unlike the other planets, we have a model of precession and nutation
        for the Earth's axis that provides a north pole vector.
        So calculate the vector first, then derive the (RA,DEC) angles from the vector.
    */

    /* Start with a north pole vector in equator-of-date coordinates: (0,0,1). */
    pos1[0] = 0.0;
    pos1[1] = 0.0;
    pos1[2] = 1.0;

    /* Convert the vector into J2000 coordinates. */
    nutation(pos1, time, INTO_2000, pos2);
    precession(pos2, *time, INTO_2000, pos1);
    axis.north.x = pos1[0];
    axis.north.y = pos1[1];
    axis.north.z = pos1[2];
    axis.north.t = *time;
    axis.north.status = ASTRO_SUCCESS;

    /* Derive angular values: right ascension and declination. */
    equ = Astronomy_EquatorFromVector(axis.north);
    if (equ.status != ASTRO_SUCCESS)
        return AxisErr(equ.status, *time);
    axis.ra = equ.ra;
    axis.dec = equ.dec;

    /* Use a modified version of the era() function that does not trim to 0..360 degrees. */
    /* This expression is also corrected to give the correct angle at the J2000 epoch. */
    axis.spin = 190.41375788700253 + (360.9856122880876 * time->ut);

    axis.status = ASTRO_SUCCESS;

    return axis;
}


/**
 * @brief Calculates information about a body's rotation axis at a given time.
 *
 * Calculates the orientation of a body's rotation axis, along with
 * the rotation angle of its prime meridian, at a given moment in time.
 *
 * This function uses formulas standardized by the IAU Working Group
 * on Cartographics and Rotational Elements 2015 report, as described
 * in the following document:
 *
 * https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf
 *
 * See #astro_axis_t for more detailed information.
 *
 * @param body
 *      The body whose rotation axis is to be found. The supported bodies are:
 *      `BODY_SUN`, `BODY_MOON`, `BODY_MERCURY`, `BODY_VENUS`, `BODY_EARTH`, `BODY_MARS`,
 *      `BODY_JUPITER`, `BODY_SATURN`, `BODY_URANUS`, `BODY_NEPTUNE`, `BODY_PLUTO`.
 *
 * @param time
 *      The time for which the body's rotation axis is to be found.
 *
 * @return astro_axis_t
 */
astro_axis_t Astronomy_RotationAxis(astro_body_t body, astro_time_t *time)
{
    astro_axis_t axis;
    double ra, dec, w;
    double radlat, radlon, rcoslat;
    double Ja, Jb, Jc, Jd, Je, N;
    double E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13;
    double d, T;

    if (time == NULL)
        return AxisErr(ASTRO_INVALID_PARAMETER, TimeError());

    d = time->tt;
    T = d / 36525.0;

    switch (body)
    {
    case BODY_SUN:
        ra = 286.13;
        dec = 63.87;
        w = 84.176 + (14.1844 * d);
        break;

    case BODY_MERCURY:
        ra = 281.0103 - (0.0328 * T);
        dec = 61.4155 - (0.0049 * T);
        w = (
            329.5988
            + (6.1385108 * d)
            + (0.01067257 * sin(DEG2RAD*(174.7910857 + 4.092335*d)))
            - (0.00112309 * sin(DEG2RAD*(349.5821714 + 8.184670*d)))
            - (0.00011040 * sin(DEG2RAD*(164.3732571 + 12.277005*d)))
            - (0.00002539 * sin(DEG2RAD*(339.1643429 + 16.369340*d)))
            - (0.00000571 * sin(DEG2RAD*(153.9554286 + 20.461675*d)))
        );
        break;

    case BODY_VENUS:
        ra = 272.76;
        dec = 67.16;
        w = 160.20 - (1.4813688 * d);
        break;

    case BODY_EARTH:
        return EarthRotationAxis(time);

    case BODY_MOON:
        /*
            https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf
            Page 8, Table 2.
        */
        E1  = DEG2RAD * (125.045 -  0.0529921*d);
        E2  = DEG2RAD * (250.089 -  0.1059842*d);
        E3  = DEG2RAD * (260.008 + 13.0120009*d);
        E4  = DEG2RAD * (176.625 + 13.3407154*d);
        E5  = DEG2RAD * (357.529 +  0.9856003*d);
        E6  = DEG2RAD * (311.589 + 26.4057084*d);
        E7  = DEG2RAD * (134.963 + 13.0649930*d);
        E8  = DEG2RAD * (276.617 +  0.3287146*d);
        E9  = DEG2RAD * (34.226  +  1.7484877*d);
        E10 = DEG2RAD * (15.134  -  0.1589763*d);
        E11 = DEG2RAD * (119.743 +  0.0036096*d);
        E12 = DEG2RAD * (239.961 +  0.1643573*d);
        E13 = DEG2RAD * (25.053  + 12.9590088*d);

        ra = (
            269.9949 + 0.0031*T
            - 3.8787*sin(E1)
            - 0.1204*sin(E2)
            + 0.0700*sin(E3)
            - 0.0172*sin(E4)
            + 0.0072*sin(E6)
            - 0.0052*sin(E10)
            + 0.0043*sin(E13)
        );

        dec = (
            66.5392 + 0.0130*T
            + 1.5419*cos(E1)
            + 0.0239*cos(E2)
            - 0.0278*cos(E3)
            + 0.0068*cos(E4)
            - 0.0029*cos(E6)
            + 0.0009*cos(E7)
            + 0.0008*cos(E10)
            - 0.0009*cos(E13)
        );

        w = (
            38.3213 + (13.17635815 - 1.4e-12*d)*d
            + 3.5610*sin(E1)
            + 0.1208*sin(E2)
            - 0.0642*sin(E3)
            + 0.0158*sin(E4)
            + 0.0252*sin(E5)
            - 0.0066*sin(E6)
            - 0.0047*sin(E7)
            - 0.0046*sin(E8)
            + 0.0028*sin(E9)
            + 0.0052*sin(E10)
            + 0.0040*sin(E11)
            + 0.0019*sin(E12)
            - 0.0044*sin(E13)
        );
        break;

    case BODY_MARS:
        ra = (
            317.269202 - 0.10927547*T
            + 0.000068 * sin(DEG2RAD*(198.991226 + 19139.4819985*T))
            + 0.000238 * sin(DEG2RAD*(226.292679 + 38280.8511281*T))
            + 0.000052 * sin(DEG2RAD*(249.663391 + 57420.7251593*T))
            + 0.000009 * sin(DEG2RAD*(266.183510 + 76560.6367950*T))
            + 0.419057 * sin(DEG2RAD*(79.398797 + 0.5042615*T))
        );

        dec = (
            54.432516 - 0.05827105*T
            + 0.000051*cos(DEG2RAD*(122.433576 + 19139.9407476*T))
            + 0.000141*cos(DEG2RAD*(43.058401 + 38280.8753272*T))
            + 0.000031*cos(DEG2RAD*(57.663379 + 57420.7517205*T))
            + 0.000005*cos(DEG2RAD*(79.476401 + 76560.6495004*T))
            + 1.591274*cos(DEG2RAD*(166.325722 + 0.5042615*T))
        );

        w = (
            176.049863 + 350.891982443297*d
            + 0.000145*sin(DEG2RAD*(129.071773 + 19140.0328244*T))
            + 0.000157*sin(DEG2RAD*(36.352167 + 38281.0473591*T))
            + 0.000040*sin(DEG2RAD*(56.668646 + 57420.9295360*T))
            + 0.000001*sin(DEG2RAD*(67.364003 + 76560.2552215*T))
            + 0.000001*sin(DEG2RAD*(104.792680 + 95700.4387578*T))
            + 0.584542*sin(DEG2RAD*(95.391654 + 0.5042615*T))
        );
        break;

    case BODY_JUPITER:
        Ja = DEG2RAD*(99.360714 + 4850.4046*T);
        Jb = DEG2RAD*(175.895369 + 1191.9605*T);
        Jc = DEG2RAD*(300.323162 + 262.5475*T);
        Jd = DEG2RAD*(114.012305 + 6070.2476*T);
        Je = DEG2RAD*(49.511251 + 64.3000*T);

        ra = (
            268.056595 - 0.006499*T
            + 0.000117*sin(Ja)
            + 0.000938*sin(Jb)
            + 0.001432*sin(Jc)
            + 0.000030*sin(Jd)
            + 0.002150*sin(Je)
        );

        dec = (
            64.495303 + 0.002413*T
            + 0.000050*cos(Ja)
            + 0.000404*cos(Jb)
            + 0.000617*cos(Jc)
            - 0.000013*cos(Jd)
            + 0.000926*cos(Je)
        );

        w = 284.95 + 870.536*d;
        break;

    case BODY_SATURN:
        ra = 40.589 - 0.036*T;
        dec = 83.537 - 0.004*T;
        w = 38.90 + 810.7939024*d;
        break;

    case BODY_URANUS:
        ra = 257.311;
        dec = -15.175;
        w = 203.81 - 501.1600928*d;
        break;

    case BODY_NEPTUNE:
        N = DEG2RAD*(357.85 + 52.316*T);
        ra = 299.36 + 0.70*sin(N);
        dec = 43.46 - 0.51*cos(N);
        w = 249.978 + 541.1397757*d - 0.48*sin(N);
        break;

    case BODY_PLUTO:
        ra = 132.993;
        dec = -6.163;
        w = 302.695 + 56.3625225*d;
        break;

    default:
        return AxisErr(ASTRO_INVALID_BODY, *time);
    }

    axis.ra = ra / 15.0;      /* convert degrees to sidereal hours */
    axis.dec = dec;
    axis.spin = w;

    /* calculate the north pole vector using the given angles. */
    radlat = dec * DEG2RAD;
    radlon = ra * DEG2RAD;
    rcoslat = cos(radlat);
    axis.north.x = rcoslat * cos(radlon);
    axis.north.y = rcoslat * sin(radlon);
    axis.north.z = sin(radlat);
    axis.north.t = *time;
    axis.north.status = ASTRO_SUCCESS;

    axis.status = ASTRO_SUCCESS;
    return axis;
}

#ifdef __cplusplus
}
#endif
