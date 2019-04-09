/*
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
#include "eph_manager.h"
#include "novas.h"
#include "astro_vector.h"
#include "vsop.h"
#include "earth.h"
#include "novas_body.h"

static double jd_begin;
static double jd_end;
static short int de_number;

static int OpenEphem(void);
static int PrintUsage(void);
static int GenerateAllSource(void);
static int TestVsopModel(VsopModel *model, int body, double threshold, double *max_arcmin, int *trunc_terms);
static int SaveVsopFile(const VsopModel *model);
static int PositionArcminError(int body, double jd, const double a[3], const double b[3], double *arcmin);
static double VectorLength(const double v[3]);

#define MOON_PERIGEE        0.00238
#define MERCURY_APHELION    0.466697
#define VENUS_APHELION      0.728213
#define EARTH_PERIHELION    0.98327
#define EARTH_APHELION      1.017
#define MARS_PERIHELION     1.382
#define JUPITER_PERIHELION  4.9501
#define SATURN_PERIHELION   9.0412
#define URANUS_PERIHELION   18.33
#define NEPTUNE_PERIHELION  29.81
#define PLUTO_PERIHELION    29.658

/* 
    ErrorRadius[body] is the worst-case (smallest) distance 
    in AU over which a position error affects Earth-based observations. 
    For example, an error in the Earth/Moon Barycenter (EMB) affects
    observations of Venus the most, because that planet comes closest
    to Earth at 0.277 AU.
*/
const double ErrorRadius[] =    
{
    EARTH_PERIHELION - MERCURY_APHELION,    /*  0 = Mercury */
    EARTH_PERIHELION - VENUS_APHELION,      /*  1 = Venus */
    EARTH_PERIHELION - VENUS_APHELION,      /*  2 = Earth/Moon Barycenter (EMB) : Venus has closest approach */
    MARS_PERIHELION    - EARTH_APHELION,    /*  3 = Mars */
    JUPITER_PERIHELION - EARTH_APHELION,    /*  4 = Jupiter */
    SATURN_PERIHELION  - EARTH_APHELION,    /*  5 = Saturn */
    URANUS_PERIHELION  - EARTH_APHELION,    /*  6 = Uranus */
    NEPTUNE_PERIHELION - EARTH_APHELION,    /*  7 = Neptune */
    PLUTO_PERIHELION   - EARTH_APHELION,    /*  8 = Pluto */
    MOON_PERIGEE,                           /*  9 = geocentric Moon */
    EARTH_PERIHELION,                       /* 10 = Sun */
};

int main(int argc, const char *argv[])
{
    if (argc == 2 && !strcmp(argv[1], "all"))
        return GenerateAllSource();

    return PrintUsage();
}

static int PrintUsage(void)
{
    fprintf(stderr, 
        "\n"
        "USAGE:\n"
        "\n"
        "generate all\n"
        "    Generate astronomy source code for all target languages.\n"
        "\n"
    );

    return 1;
}

static int OpenEphem()
{
    const char *filename;
    int error;

    filename = getenv("EPHEM");
    if (filename == NULL)
        filename = "lnxp1600p2200.405";

    error = ephem_open((char *)filename, &jd_begin, &jd_end, &de_number);
    if (error)
        fprintf(stderr, "OpenEphem: Error %d trying to open ephemeris file: %s\n", error, filename);

    return error;
}

static int NovasBodyPos(double jd, int body, double pos[3])
{
    int error, k;
    double jed[2] = { jd, 0.0 };
    double sun_pos[3], moon_pos[3], vel[3];

    if (body == BODY_EARTH || body == BODY_MOON)
    {
        static const double EarthMoonMassRatio = 81.30056;
        double factor;

        /* 
            The caller is asking for the Earth's position or the Moon's position. 
            NOVAS does not directly represent either body.
            Instead, we have to calculate the Earth or Moon using
            the Earth/Moon Barycenter (EMB) and the Geocentric Moon (GM).
        */
        error = state(jed, BODY_EMB, pos, vel);
        if (error)
        {
            fprintf(stderr, "NovasBodyPos: state(%lf, EMB) returned %d\n", jd, error);
            return error;
        }

        error = state(jed, BODY_GM, moon_pos, vel);
        if (error)
        {
            fprintf(stderr, "NovasBodyPos: state(%lf, GM) returned %d\n", jd, error);
            return error;
        }

        if (body == BODY_EARTH)
            /* Calculate the Earth's position away from the EMB, opposite the direction of the Moon. */
            factor = -1.0 / (1.0 + EarthMoonMassRatio);
        else
            /* Calculate the Moon's position away from the EMB, along the vector from the Earth to the Moon. */
            factor = EarthMoonMassRatio / (1.0 + EarthMoonMassRatio);

        for (k=0; k<3; ++k)
            pos[k] += factor * moon_pos[k];
    }
    else
    {
        /* This is a body that NOVAS directly models in its ephemerides. */
        error = state(jed, body, pos, vel);
        if (error)
        {
            fprintf(stderr, "NovasBodyPos: state(%lf, %d) returned %d\n", jd, body, error);
            return error;
        }

        /* Special case: geocentric moon should not be converted to heliocentric coordinates. */
        if (body == BODY_GM) return 0;
    }

    error = state(jed, BODY_SUN, sun_pos, vel);
    if (error)
    {
        fprintf(stderr, "NovasBodyPos: state(%lf, SUN) returned %d\n", jd, error);
        return error;
    }

    pos[0] -= sun_pos[0];
    pos[1] -= sun_pos[1];
    pos[2] -= sun_pos[2];
    return 0;
}

static int LoadVsopFile(VsopModel *model, int body)
{
    static const char * const BodyName[] =
    {
        "mer", "ven", "ear", "mar", "jup", "sat", "ura", "nep"
    //    0      1      2      3      4      5      6      7
    };
    char filename[100];

    if (body < 0 || body > 7)
    {
        fprintf(stderr, "LoadVsopFile: Invalid body=%d\n", body);
        return 1;
    }

    snprintf(filename, sizeof(filename)-1, "vsop/VSOP87B.%s", BodyName[body]);
    return VsopLoadModel(model, filename);
}

static int SearchVsop(int body)
{
    const double arcmin_limit = 0.4;
    int error;
    int trunc_terms, winner_terms = -1;
    double thresh_lo = 1.0e-6;
    double thresh_hi = 1.0e-2;
    double max_arcmin, threshold;
    double winner_threshold = 0.0;
    double winner_arcmin = -1.0;
    VsopModel model;

    error = LoadVsopFile(&model, body);
    if (error) goto fail;

    while (thresh_hi / thresh_lo > 1.000001)
    {
        threshold = sqrt(thresh_lo * thresh_hi);
        error = TestVsopModel(&model, body, threshold, &max_arcmin, &trunc_terms);
        if (error) goto fail;

        if (max_arcmin >= arcmin_limit)
        {
            /* not accurate enough, so decrease threshold next time */
            thresh_hi = threshold;
        }
        else
        {
            /* we could tolerate larger inaccuracy to get smaller data */
            thresh_lo = threshold;

            /* track best-so-far (smallest data size with tolerable error) */
            if (winner_terms < 0 || trunc_terms < winner_terms)
            {
                winner_terms = trunc_terms;
                winner_arcmin = max_arcmin;
                winner_threshold = threshold;
                
                error = SaveVsopFile(&model);
                if (error) goto fail;
            }
        }
    }

    if (winner_terms < 0)
    {
        fprintf(stderr, "SearchVsop: could not find a solution for body %d\n", body);
        error = 1;
        goto fail;
    }

    printf("SearchVsop(WINNER): body=%d, terms=%d, arcmin=%lg, threshold=%lg\n", body, winner_terms, winner_arcmin, winner_threshold);
    fflush(stdout);

fail:
    VsopFreeModel(&model);
    return error;
}

static int SaveVsopFile(const VsopModel *model)
{
    char filename[100];
    
    snprintf(filename, sizeof(filename)-1, "output/vsop_%d.txt", (int)model->body);
    return VsopWriteTrunc(model, filename);
}

static int TestVsopModel(VsopModel *model, int body, double threshold, double *max_arcmin, int *trunc_terms)
{
    int error;
    double jdStart, jdStop, jd, jdDelta;
    double vpos[VSOP_MAX_COORDS];
    double npos[3];
    double arcmin;

    *max_arcmin = -1.0;
    *trunc_terms = -1;

    if (body < 0 || body > 7)
    {
        fprintf(stderr, "TestVsopFile: Invalid body %d\n", body);
        error = 1;
        goto fail;
    }

    jdStart = julian_date(1900, 1, 1, 0.0);
    jdStop = julian_date(2101, 1, 1, 0.0);
    jdDelta = 1.0;

    error = VsopTruncate(model, jdStart, jdStop, threshold);
    if (error) goto fail;
    *trunc_terms = VsopTermCount(model);

    for (jd = jdStart; jd <= jdStop; jd += jdDelta)
    {
        VsopCalc(model, jd, vpos);
        if (body == 2)
            error = NovasEarth(jd, npos);
        else
            error = NovasBodyPos(jd, body, npos);
        if (error) goto fail;
        error = PositionArcminError(body, jd, npos, vpos, &arcmin);
        if (error) goto fail;
        if (arcmin > *max_arcmin) *max_arcmin = arcmin;
    }

    error = 0;

fail:
    return error;
}

static int BuildVsopData(void)
{
    int body;
    int error;

    for (body=0; body < 8; ++body)
        if (0 != (error = SearchVsop(body)))
            break;

    return error;
}

static int GenerateAllSource(void)
{
    int error;
    if (0 != (error = OpenEphem())) goto fail;
    if (0 != (error = BuildVsopData())) goto fail;

fail:
    ephem_close();
    return error;
}

static int PositionArcminError(int body, double jd, const double a[3], const double b[3], double *arcmin)
{
    int error, k;
    double diff, sum, au;
    double epos[3], evel[3], jed[2];
    VectorType adiff, bdiff;

    *arcmin = 99999.0;

    if (body == BODY_SUN)
    {
        /* Trivial case: The Sun as always at heliocentric coordinates (0, 0, 0). */
        /* Let's just check that both positions are very close to the origin. */
        diff =  a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
        diff += b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
        diff = sqrt(diff);
        if (diff > 1.0e-6)
        {
            fprintf(stderr, "PositionArcminError: Excessive error (%lg) in heliocentric Sun.\n", diff);
            return 1;
        }
        *arcmin = 0.0;
        return 0;
    }

    if (body == BODY_MOON) 
    {
        /*
            Measure the position of the moon as seen from the *SURFACE* of the Earth,
            at the point along a line from the geocenter to the object.
            This is closer than the geocenter or the barycenter, so the anglular error is larger.
            For other bodies, there is no significant difference between 
            topocentric error and geocentric error.
        */
        const double scale = (ERAD / AU) / VectorLength(a);
        error = NovasBodyPos(jd, BODY_EARTH, epos);
        if (error) return error;
        for (k=0; k<3; ++k) epos[k] += scale * a[k];
        goto calc;
    }

    if (body == BODY_GM)
    {
        /*
            The given position vectors are geocentric, not heliocentric.
            Adjust them them to be topocentric in the direction of the Moon,
            to maximize the error experienced by an observer on the surface of the Earth.
        */
        const double dilate = (1.0 + EarthMoonMassRatio) / EarthMoonMassRatio;

        adiff.t = jd;   /* doesn't matter, but I don't like leaving undefined memory */
        adiff.v[0] = a[0];
        adiff.v[1] = a[1];
        adiff.v[2] = a[2];
        
        bdiff.t = jd;   /* doesn't matter, but I don't like leaving undefined memory */
        bdiff.v[0] = b[0];
        bdiff.v[1] = b[1];
        bdiff.v[2] = b[2];

        *arcmin = dilate * (RAD2DEG * 60.0) * AngleBetween(adiff, bdiff);
        return 0;
    }

    if (body == BODY_EMB || body == BODY_EARTH)
    {
        /* 
            For the Earth/Moon Barycenter (EMB=2), or Earth=11, we use a special estimate
            of angular error. An error in the Earth's position affects the observed
            position of all other bodies in the solar system.
            It is appropriate to use a much stricter error measure than for other objects.
            The worst case error is for Venus when it is closest to the Earth.
            We approximate that worst case distance by EARTH_PERIHELION - VENUS_APHELION.
        */
        diff = a[0] - b[0];
        sum = diff*diff;

        diff = a[1] - b[1];
        sum += diff*diff;

        diff = a[2] - b[2];
        sum += diff*diff;

        au = sqrt(sum);
        *arcmin = (RAD2DEG * 60.0) * (au / (EARTH_PERIHELION - VENUS_APHELION));
        return 0;
    }

    /* Exclude Sun and Geocentric Moon (GM) from error estimates. */
    /* FIXFIXFIX - Rework to use NovasBodyPos(), so we can support Moon also. */
    /* Can also use Earth position rather than EMB position. */
    if (body < BODY_MERCURY || body > BODY_PLUTO)
    {
        fprintf(stderr, "PositionArcminError(FATAL): Invalid body = %d\n", body);
        return 1;
    }

    /*
        Any planet other than the EMB is assumed to be viewed from the Earth.
        We use the EMB as a good enough estimate of the Earth's position.
        Take two vectors: (a - epos) and (b - epos).
        These are vectors from the EMB to the two calculated positions of the other planet.
        Find the angle between the two vectors using inverse cosine of dot product.
    */
    jed[0] = jd;
    jed[1] = 0.0;
    error = state(jed, BODY_EMB, epos, evel);
    if (error)
    {
        fprintf(stderr, "PositionArcminError: state(%lf, EMB) returned %d\n", jd, error);
        return error;
    }

calc:
    adiff.t = jd;   /* doesn't matter, but I don't like leaving undefined memory */
    adiff.v[0] = a[0] - epos[0];
    adiff.v[1] = a[1] - epos[1];
    adiff.v[2] = a[2] - epos[2];
    
    bdiff.t = jd;   /* doesn't matter, but I don't like leaving undefined memory */
    bdiff.v[0] = b[0] - epos[0];
    bdiff.v[1] = b[1] - epos[1];
    bdiff.v[2] = b[2] - epos[2];

    *arcmin = (RAD2DEG * 60.0) * AngleBetween(adiff, bdiff);
    return 0;
}

static double VectorLength(const double v[3])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
