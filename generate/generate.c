/*
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
#ifdef _WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif

#include "eph_manager.h"
#include "novas.h"

#include "astro_vector.h"
#include "chebyshev.h"
#include "codegen.h"
#include "earth.h"
#include "ephfile.h"
#include "novas_body.h"
#include "vsop.h"
#include "top2013.h"

const double PLUTO_TOLERANCE_ARCMIN = 1.4;

int Verbose;
#define DEBUG(...)  do{if(Verbose)printf(__VA_ARGS__);}while(0)

#define VECTOR_DIM 3
#define MIN_YEAR 1700
#define MAX_YEAR 2200

typedef struct
{
    double rmsPositionError;
    double maxPositionError;
    double rmsArcminError;
    double maxArcminError;
    long dataCount;     /* a measure of how many floating point numbers need to be encoded for ephemeris translation */
}
error_stats_t;

typedef struct
{
    int body;       /* which body we are calculating for */
}
sample_context_t;

static double jd_begin;
static double jd_end;
static short int de_number;

static int ParseDate(const char *text, double *tt);
static int OpenEphem(void);
static int PrintUsage(void);
static int GenerateVsopPlanets(void);
static int TopCalc(const char *name, const char *date);
static int GenerateApsisTestData(void);
static int GenerateSource(void);
static int TestVsopModel(vsop_model_t *model, int body, double threshold, double *max_arcmin, int *trunc_terms);
static int SaveVsopFile(const vsop_model_t *model);
static int PositionArcminError(int body, double jd, const double a[3], const double b[3], double *arcmin);
static double VectorLength(const double v[3]);
static double VectorError(double a[3], double b[3]);
static int CheckTestOutput(const char *filename);
static vsop_body_t LookupBody(const char *name);
static int CheckSkyPos(observer *location, const char *filename, int lnum, const char *line, double *arcmin_equ, double *arcmin_hor, vsop_body_t *body);
static int UnitTestChebyshev(void);
static int DistancePlot(const char *name, double tt1, double tt2);
static int ImproveVsopApsides(vsop_model_t *model);
static int DeltaTimePlot(const char *outFileName);
static int TopFileInfo(const char *filename, const char *name);
static int ValidateTop2013(void);
static int Diff(const char *filename1, const char *filename2, int *nlines);

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
    if (argc > 1 && !strcmp(argv[1], "-v"))
    {
        Verbose = 1;
        --argc;
        ++argv;
    }

    if (argc == 2 && !strcmp(argv[1], "planets"))
        return GenerateVsopPlanets();

    if (argc == 2 && !strcmp(argv[1], "apsis"))
        return GenerateApsisTestData();

    if (argc == 2 && !strcmp(argv[1], "source"))
        return GenerateSource();

    if (argc == 3 && !strcmp(argv[1], "check"))
        return CheckTestOutput(argv[2]);

    if (argc == 2 && !strcmp(argv[1], "chebyshev"))
        return UnitTestChebyshev();

    if (argc == 5 && !strcmp(argv[1], "distplot"))
        return DistancePlot(argv[2], atof(argv[3]), atof(argv[4]));

    if (argc == 3 && !strcmp(argv[1], "dtplot"))
        return DeltaTimePlot(argv[2]);

    if (argc == 4 && !strcmp(argv[1], "topcalc"))
        return TopCalc(argv[2], argv[3]);

    if (argc == 2 && !strcmp(argv[1], "validate_top2013"))
        return ValidateTop2013();

    if (argc == 4 && !strcmp(argv[1], "topinfo"))
        return TopFileInfo(argv[2], argv[3]);

    return PrintUsage();
}

static int PrintUsage(void)
{
    fprintf(stderr,
        "\n"
        "USAGE:\n"
        "\n"
        "generate planets\n"
        "    Generate predictive models for all planets.\n"
        "\n"
        "generate check testfile\n"
        "    Verify the calculations in the testfile generated by a unit test.\n"
        "\n"
        "generate source\n"
        "    Generate source code after running 'generate planets'.\n"
        "\n"
        "generate chebyshev\n"
        "    Run unit test on Chebyshev interpolator.\n"
        "\n"
        "generate distplot planet tt1 tt2\n"
        "    Generate a CSV plot of the given planet's heliocentric\n"
        "    distance as a function of time.\n"
        "\n"
        "generate dtplot outfile.csv\n"
        "   Generate a CSV plot of the Delta-T extrapolator.\n"
        "\n"
        "generate validate_top2013\n"
        "   Validates code for calculating outer planet positions using TOP2013.\n"
        "\n"
        "generate topcalc planet date\n"
        "   Calculate an exact position and velocity of the given planet\n"
        "   (Jupiter, Saturn, Uranus, Neptune, or Pluto)\n"
        "   at the specified date and time.\n"
        "\n"
        "generate topinfo filename planet\n"
        "   Prints summary info about the TOP2013 file.\n"
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


#define NovasBodyPos(jd,body,pos)   NovasBodyPosVel(jd,body,pos,NULL)

/*
 *  Returns heliocentric position and velocity vectors based on NOVAS
 *  interpolation of the ephemeris file. The jd must be a Julian Date
 *  in the range of years 1600..2200.
 */
static int NovasBodyPosVel(double jd, int body, double pos[3], double vel[3])
{
    int error, k;
    double jed[2];
    double sun_pos[3], moon_pos[3];
    double sun_vel[3], moon_vel[3];
    double ignore_vel[3];

    if (vel == NULL)
        vel = ignore_vel;

    jed[0] = jd;
    jed[1] = 0.0;

    if (body == BODY_SSB)
    {
        /*
            The NOVAS state() function does everything with respect to
            the Solar System Barycenter (SSB) already.
            Below we will get the Sun's position relative to the SSB.
            Negating that will result in the SSB position with respect to the Sun.
        */
        pos[0] = pos[1] = pos[2] = 0.0;
        vel[0] = vel[1] = vel[2] = 0.0;
    }
    else if (body == BODY_EARTH || body == BODY_MOON)
    {
        double factor;

        /*
            The caller is asking for the Earth's position or the Moon's position.
            NOVAS does not directly represent either body.
            Instead, we have to calculate the Earth or Moon using
            the Earth/Moon Barycenter (EMB) and the Geocentric Moon (GM).
        */
        error = state(jed, BODY_EMB, pos, vel);
        if (error)
            FAIL("NovasBodyPos: state(%lf, EMB) returned %d\n", jd, error);

        error = state(jed, BODY_GM, moon_pos, moon_vel);
        if (error)
            FAIL("NovasBodyPos: state(%lf, GM) returned %d\n", jd, error);

        if (body == BODY_EARTH)
            /* Calculate the Earth's position away from the EMB, opposite the direction of the Moon. */
            factor = -1.0 / (1.0 + EarthMoonMassRatio);
        else
            /* Calculate the Moon's position away from the EMB, along the vector from the Earth to the Moon. */
            factor = EarthMoonMassRatio / (1.0 + EarthMoonMassRatio);

        for (k=0; k<3; ++k)
        {
            pos[k] += factor * moon_pos[k];
            vel[k] += factor * moon_vel[k];
        }
    }
    else
    {
        /* This is a body that NOVAS directly models in its ephemerides. */
        error = (int)state(jed, (short)body, pos, vel);
        if (error)
            FAIL("NovasBodyPos: state(%lf, %d) returned %d\n", jd, body, error);

        /* Special case: geocentric moon should not be converted to heliocentric coordinates. */
        if (body == BODY_GM) return 0;
    }

    error = (int)state(jed, BODY_SUN, sun_pos, sun_vel);
    if (error)
        FAIL("NovasBodyPos: state(%lf, SUN) returned %d\n", jd, error);

    pos[0] -= sun_pos[0];
    pos[1] -= sun_pos[1];
    pos[2] -= sun_pos[2];
    vel[0] -= sun_vel[0];
    vel[1] -= sun_vel[1];
    vel[2] -= sun_vel[2];
    error = 0;

fail:
    return error;
}

static int LoadVsopFile(vsop_model_t *model, int body)
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

    snprintf(filename, sizeof(filename), "vsop/VSOP87B.%s", BodyName[body]);
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
    vsop_model_t model;

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
                error = ImproveVsopApsides(&model);
                if (error) goto fail;
                VsopTrim(&model);       /* shave off all trailing empty series right before saving */
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

    printf("SearchVsop(WINNER): body=%d, terms=%d, arcmin=%0.6lf, threshold=%0.4le\n", body, winner_terms, winner_arcmin, winner_threshold);
    fflush(stdout);

fail:
    VsopFreeModel(&model);
    return error;
}

static int SaveVsopFile(const vsop_model_t *model)
{
    char filename[100];

    snprintf(filename, sizeof(filename), "output/vsop_%d.txt", (int)model->body);
    return VsopWriteTrunc(model, filename);
}

static int TestVsopModel(vsop_model_t *model, int body, double threshold, double *max_arcmin, int *trunc_terms)
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

    jdStart = julian_date(MIN_YEAR, 1, 1, 0.0);
    jdStop = julian_date(MAX_YEAR, 1, 1, 0.0);
    jdDelta = 1.0;

    error = VsopTruncate(model, jdStart - T0, jdStop - T0, threshold);
    if (error) goto fail;
    *trunc_terms = VsopTermCount(model);

    for (jd = jdStart; jd <= jdStop; jd += jdDelta)
    {
        VsopCalcPos(model, jd - T0, vpos);
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

static int TermContribution(const vsop_term_t *term, int sindex, double *contrib)
{
    switch (sindex)
    {
    case 0:
        /* t^0 term has derivative proprotional to amplitude*frequency */
        *contrib = fabs(term->amplitude * term->frequency);
        return 0;

    default:
        fprintf(stderr, "TermContribution: t^%d derivative rule not yet implemented.\n", sindex);
        *contrib = NAN;
        return 1;
    }
}

static int ExpandSeries(vsop_formula_t *formula, int sindex, int n_top_terms)
{
    int error;
    vsop_series_t *series;
    int t, t_best, n;
    double contrib, max_contrib;

    /*
        Examine the specified truncated series in the given formula.
        Calculate an estimate of the contribution of every term's
        first derivative. Pick n_top_terms of the most important ones.
    */

    if (sindex < 0 || sindex >= formula->nseries_calc)
    {
        fprintf(stderr, "ExpandSeries: Invalid series index %d\n", sindex);
        return 1;
    }

    series = &formula->series[sindex];
    if (series->nterms_calc + n_top_terms >= series->nterms_total)
    {
        fprintf(stderr, "ExpandSeries: series calc=%d, total=%d; cannot expand by %d\n",
            series->nterms_calc, series->nterms_total, n_top_terms);
        return 1;
    }

    for (n=0; n < n_top_terms; ++n)
    {
        /* Iterate through the truncated terms at the end of this series. */
        /* Pick the one whose contribution is greatest and append it next. */
        t_best = -1;
        max_contrib = 0.0;
        for (t = series->nterms_calc; t < series->nterms_total; ++t)
        {
            CHECK(TermContribution(&series->term[t], sindex, &contrib));
            if (contrib > max_contrib)
            {
                t_best = t;
                max_contrib = contrib;
            }
        }
        if (t_best < 0)
        {
            fprintf(stderr, "ExpandSeries: Found only %d of %d top terms\n", n, n_top_terms);
            return 1;
        }
        /* Like bubble sort: swap best with whatever is already just beyond the truncation. */
        if (t_best != series->nterms_calc)
        {
            vsop_term_t swap = series->term[series->nterms_calc];
            series->term[series->nterms_calc] = series->term[t_best];
            series->term[t_best] = swap;
        }
        ++(series->nterms_calc);
    }

    error = 0;
fail:
    return error;
}

static int ImproveVsopApsides(vsop_model_t *model)
{
    int error;
    int n_top_terms;
    vsop_formula_t *radform;        /* formula for calculating radial distance */

    if (model->version != VSOP_HELIO_SPHER_DATE && model->version != VSOP_HELIO_SPHER_J2000)
        FAIL("ImproveVsopApsides: Cannot optimize VSOP version %d\n", model->version);

    if (model->ncoords != 3)
        FAIL("ImproveVsopApsides: INTERNAL ERROR: model has incorrect number of coordinates: %d\n", model->ncoords);

    radform = &model->formula[2];

    switch (model->body)
    {
    case BODY_VENUS:
        n_top_terms = 5;
        break;

    case BODY_EARTH:
        n_top_terms = 5;
        break;

    case BODY_JUPITER:
        n_top_terms = 3;
        break;

    case BODY_SATURN:
        n_top_terms = 7;
        break;

    case BODY_URANUS:
        n_top_terms = 3;
        break;

    case BODY_NEPTUNE:
        n_top_terms = 5;
        break;

    default:
        n_top_terms = 0;
        break;
    }

    CHECK(ExpandSeries(radform, 0, n_top_terms));
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

static int GenerateVsopPlanets(void)
{
    int error;

    CHECK(OpenEphem());
    CHECK(BuildVsopData());

fail:
    ephem_close();
    return error;
}

static double PlanetOrbitalPeriod(int body)
{
    switch (body)
    {
    case BODY_MERCURY:  return     87.969;
    case BODY_VENUS:    return    224.701;
    case BODY_EARTH:    return    365.256;
    case BODY_MARS:     return    686.980;
    case BODY_JUPITER:  return   4332.589;
    case BODY_SATURN:   return  10759.22;
    case BODY_URANUS:   return  30685.4;
    case BODY_NEPTUNE:  return  60189.0;
    case BODY_PLUTO:    return  90560.0;
    default:            return  0.0;        /* invalid body */
    }
}

static int SearchApsis(int body, double direction, double jd1, double jd2, double *jdx, double *xdist)
{
    int error;
    const int npoints = 10;
    int i, best_i;
    double jd;
    double pos[3];
    double interval, dist, max_dist;

    *xdist = *jdx = NAN;

    /* Keep iterating until interval is within 1 second of uncertainty. */
    while (jd2-jd1 > 1.0 / 86400.0)
    {
        interval = (jd2 - jd1) / (npoints - 1.0);
        max_dist = NAN;
        best_i = -1;
        for (i=0; i < npoints; ++i)
        {
            jd = jd1 + (i*interval);
            CHECK(NovasBodyPos(jd, body, pos));
            dist = direction * VectorLength(pos);
            if (i==0 || dist > max_dist)
            {
                best_i = i;
                max_dist = dist;
            }
        }

        /* Narrow in on +/- 1 interval around the extreme point. */
        jd1 += (best_i - 1)*interval;
        jd2 = jd1 + (2 * interval);
    }

    *jdx = (jd1 + jd2) / 2.0;
    CHECK(NovasBodyPos(*jdx, body, pos));
    *xdist = VectorLength(pos);
    error = 0;
fail:
    return error;
}

static void PrintApsis(FILE *outfile, int kind, double jdx, double dist)
{
    short year, month, day;
    int hour, minute, second;
    double frac;

    cal_date(jdx, &year, &month, &day, &frac);
    hour = (int)floor(frac);
    frac = 60.0 * (frac - hour);
    minute = (int)floor(frac);
    frac = 60.0 * (frac - minute);
    second = (int)floor(frac);
    fprintf(outfile, "%d %04d-%02d-%02dT%02d:%02d:%02dZ %12.7lf\n", kind, (int)year, (int)month, (int)day, hour, minute, second, dist);
}

static int GenerateApsisFile(int body, const char *outFileName)
{
    int error = 1;
    FILE *outfile = NULL;
    double pos[3];
    double dist[3];     /* dist[0]=current distance, dist[1]=previous distance, dist[2]=distance before that */
    double jdStart, jdStop, jd;
    double period, interval;
    const double samples_per_orbit = 1000.0;
    struct { double dist; double jd; int kind; } buffer[3];
    int buflen = 0;

    period = PlanetOrbitalPeriod(body);
    if (period <= 0.0)
    {
        fprintf(stderr, "GenerateApsisFile: Invalid body %d\n", body);
        error = 1;
        goto fail;
    }
    interval = period / samples_per_orbit;

    outfile = fopen(outFileName, "wt");
    if (outfile == NULL)
    {
        fprintf(stderr, "GenerateApsisFile: Cannot open output file '%s'\n", outFileName);
        error = 1;
        goto fail;
    }

    jdStart = julian_date(MIN_YEAR, 1, 1, 0.0);
    jdStop = julian_date(MAX_YEAR, 1, 1, 0.0);
    dist[0] = dist[1] = dist[2] = -1.0;
    for (jd=jdStart; jd <= jdStop; jd += interval)
    {
        CHECK(NovasBodyPos(jd, body, pos));
        dist[0] = VectorLength(pos);
        if (dist[2] > 0.0)
        {
            if (dist[1] >= dist[2] && dist[1] >= dist[0])
            {
                buffer[buflen].kind = 1;    /* aphelion */
                CHECK(SearchApsis(body, +1.0, jd-2.0*interval, jd, &buffer[buflen].jd, &buffer[buflen].dist));
                ++buflen;
            }
            else if (dist[1] <= dist[2] && dist[1] <= dist[0])
            {
                buffer[buflen].kind = 0;    /* perihelion */
                CHECK(SearchApsis(body, -1.0, jd-2.0*interval, jd, &buffer[buflen].jd, &buffer[buflen].dist));
                ++buflen;
            }

            /*
                We have to handle weird special cases with Neptune.
                It can have 3 local minima/maxima near perihelion and aphelion,
                due to Sun wobbling around Solar System Barycenter.
            */
            if (buflen > 0)
            {
                if (body == BODY_NEPTUNE)
                {
                    const double mean_dist = 30.1;  /* mean orbital distance of Neptune in AU */

                    if (buflen == 3)
                    {
                        /* If all 3 entries are on the same side of the mean distance, pick the extreme. */
                        if (buffer[0].dist > mean_dist && buffer[1].dist > mean_dist && buffer[2].dist > mean_dist)
                        {
                            int imax = 0;
                            if (buffer[1].dist > buffer[imax].dist) imax = 1;
                            if (buffer[2].dist > buffer[imax].dist) imax = 2;
                            if (buffer[imax].kind != 1)
                            {
                                fprintf(stderr, "FATAL(GenerateApsisFile): expected Neptune aphelion\n");
                                error = 1;
                                goto fail;
                            }
                            PrintApsis(outfile, buffer[imax].kind, buffer[imax].jd, buffer[imax].dist);
                            buflen = 0;     /* empty the buffer */
                        }
                        else if (buffer[0].dist < mean_dist && buffer[1].dist < mean_dist && buffer[2].dist < mean_dist)
                        {
                            int imin = 0;
                            if (buffer[1].dist < buffer[imin].dist) imin = 1;
                            if (buffer[2].dist < buffer[imin].dist) imin = 2;
                            if (buffer[imin].kind != 0)
                            {
                                fprintf(stderr, "FATAL(GenerateApsisFile): expected Neptune perihelion\n");
                                error = 1;
                                goto fail;
                            }
                            PrintApsis(outfile, buffer[imin].kind, buffer[imin].jd, buffer[imin].dist);
                            buflen = 0;     /* empty the buffer */
                        }
                        else
                        {
                            /* We have to fail because the buffer is full and we don't know how to empty it. */
                            fprintf(stderr, "FATAL(GenerateApsisFile): Unhandled special case with Neptune apsis\n");
                            error = 1;
                            goto fail;
                        }
                    }
                }
                else
                {
                    PrintApsis(outfile, buffer[0].kind, buffer[0].jd, buffer[0].dist);
                    buflen = 0;
                }
            }
        }
        dist[2] = dist[1];
        dist[1] = dist[0];
    }

    if (buflen > 0)
    {
        fprintf(stderr, "FATAL(GenerateApsisFile): buffer holds %d unhandled events.\n", buflen);
        error = 1;
        goto fail;
    }

fail:
    if (outfile)
    {
        fclose(outfile);
        if (error) remove(outFileName);
    }
    return error;
}

static int GenerateApsisTestData(void)
{
    int error;

    CHECK(OpenEphem());

    /*
        Tricky: BODY_EARTH=11, but we write to output/apsis_2.txt.
        The BODY_EARTH is my extension here for working around NOVAS,
        but Astronomy Engine uses 2 to represent the Earth.
    */
    CHECK(GenerateApsisFile(BODY_MERCURY, "apsides/apsis_0.txt"));
    CHECK(GenerateApsisFile(BODY_VENUS,   "apsides/apsis_1.txt"));
    CHECK(GenerateApsisFile(BODY_EARTH,   "apsides/apsis_2.txt"));
    CHECK(GenerateApsisFile(BODY_MARS,    "apsides/apsis_3.txt"));
    CHECK(GenerateApsisFile(BODY_JUPITER, "apsides/apsis_4.txt"));
    CHECK(GenerateApsisFile(BODY_SATURN,  "apsides/apsis_5.txt"));
    CHECK(GenerateApsisFile(BODY_URANUS,  "apsides/apsis_6.txt"));
    CHECK(GenerateApsisFile(BODY_NEPTUNE, "apsides/apsis_7.txt"));
    CHECK(GenerateApsisFile(BODY_PLUTO,   "apsides/apsis_8.txt"));

fail:
    ephem_close();
    return error;
}

static int GenerateSource(void)
{
    int error;
    CHECK(GenerateCode(CODEGEN_LANGUAGE_C, "../source/c/astronomy.c", "template/astronomy.c",  "output"));
    CHECK(GenerateCode(CODEGEN_LANGUAGE_CSHARP, "../source/csharp/astronomy.cs", "template/astronomy.cs",  "output"));
    CHECK(GenerateCode(CODEGEN_LANGUAGE_JS, "../source/js/astronomy.js", "template/astronomy.js", "output"));
    CHECK(GenerateCode(CODEGEN_LANGUAGE_PYTHON, "../source/python/astronomy.py", "template/astronomy.py", "output"));
fail:
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

    if ((body != BODY_SSB) && (body < BODY_MERCURY || body > BODY_PLUTO))
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

static double VectorError(double a[3], double b[3])
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

static int ParseObserver(const char *filename, int lnum, const char *line, observer *location)
{
    int error;
    on_surface surface;

    if (3 != sscanf(line, "o %lf %lf %lf", &surface.latitude, &surface.longitude, &surface.height))
    {
        fprintf(stderr, "ParseObserver: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }

    surface.temperature = 20;
    surface.pressure = 1000;

    error = make_observer(1, &surface, NULL, location);
    if (error)
    {
        fprintf(stderr, "ParseObserver: error %d returned by make_observer() for line %d of file %s\n", error, lnum, filename);
    }
    return error;
}

static int EclipticVector(double jd, vsop_body_t body, double ecl_pos[3])
{
    int error;
    double eq_pos[3];

    /* Find equatorial coordinates for the body. */
    CHECK(NovasBodyPos(jd, body, eq_pos));

    /* Convert equatorial cartesian coordinates to ecliptic longitude. */
    error = equ2ecl_vec(jd, 2, 0, eq_pos, ecl_pos);
    if (error)
    {
        fprintf(stderr, "EclipticLongitude: equ2ecl_vec() returned %d\n", error);
        goto fail;
    }

fail:
    return error;
}

static double EclipticLongitude(const double ecl_pos[3])
{
    double lon = RAD2DEG * atan2(ecl_pos[1], ecl_pos[0]);
    if (lon < 0.0)
        lon += 360.0;

    return lon;
}

static double LongitudeOffset(double diff)
{
    double offset = diff;
    while (offset <= -180) offset += 360;
    while (offset > 180) offset -= 360;
    return offset;
}

static int CheckEcliptic(const char *filename, int lnum, const char *line, double *arcmin, vsop_body_t *outbody)
{
    int error;
    vsop_body_t body;
    char name[10];
    char event[10];
    double body_ecl[3];
    double earth_ecl[3];
    double geo[3];
    double tt, jd, lon, dist, blon, elon, diff;
    double tolerance;

    *arcmin = 99999.0;
    *outbody = VSOP_INVALID_BODY;

    /* Example line: "e Jupiter opp <tt> <au>" */
    if (4 != sscanf(line, "e %9[A-Za-z] %9[a-z] %lf %lf", name, event, &tt, &dist))
    {
        fprintf(stderr, "CheckEcliptic: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }

    *outbody = body = LookupBody(name);
    if (body < 0)
    {
        fprintf(stderr, "CheckEcliptic: Unknown body '%s' on line %d of file %s\n", name, lnum, filename);
        return 1;
    }

    jd = tt + T0;
    CHECK(EclipticVector(jd, body, body_ecl));
    CHECK(EclipticVector(jd, VSOP_EARTH, earth_ecl));

    geo[0] = body_ecl[0] - earth_ecl[0];
    geo[1] = body_ecl[1] - earth_ecl[1];
    geo[2] = body_ecl[2] - earth_ecl[2];
    dist = VectorLength(geo);

    if (!strcmp(event, "opp"))
    {
        /* Opposition (superior planets only) */
        lon = 0.0;

        /* Distance could be anything... so ignore it. */
    }
    else if (!strcmp(event, "inf"))
    {
        /* Inferior conjunction (inferior planets only) */
        lon = 0.0;

        /* Distance from Earth must be less than 1 AU for an inferior conjunction. */
        if (dist >= 1.0)
        {
            fprintf(stderr, "CheckEcliptic: Invalid distance %0.3lf AU for inferior conjunction; line %d file %s\n", dist, lnum, filename);
            return 1;
        }
    }
    else if (!strcmp(event, "sup"))
    {
        /* Superior conjunction (inferior or superior planets) */
        lon = 180.0;

        /* Distance from Earth must be greater than 1 AU for a superior conjunction. */
        if (dist <= 1.0)
        {
            fprintf(stderr, "CheckEcliptic: Invalid distance %0.3lf AU for superior conjunction; line %d file %s\n", dist, lnum, filename);
            return 1;
        }
    }
    else
    {
        fprintf(stderr, "CheckEcliptic: Invalid event code '%s' on line %d of file %s\n", event, lnum, filename);
        return 1;
    }

    blon = EclipticLongitude(body_ecl);
    elon = EclipticLongitude(earth_ecl);
    diff = LongitudeOffset((elon - blon) - lon);

    tolerance = (body == BODY_PLUTO) ? 1.3 : 1.0;
    *arcmin = fabs(diff * 60.0);
    if (*arcmin > tolerance)
    {
        fprintf(stderr, "CheckEcliptic: Excessive arcminute error = %0.3lf on line %d of file %s\n", *arcmin, lnum, filename);
        return 1;
    }

fail:
    return error;
}

static int CheckTestVector(const char *filename, int lnum, const char *line, double *arcmin, vsop_body_t *outbody)
{
    int error;
    vsop_body_t body;
    char name[10];
    double tt, jd, xpos[3], npos[3];
    double tolerance;

    *outbody = VSOP_INVALID_BODY;
    *arcmin = 99999.0;

    if (5 != sscanf(line, "v %9[A-Za-z] %lf %lf %lf %lf", name, &tt, &xpos[0], &xpos[1], &xpos[2]))
    {
        fprintf(stderr, "CheckTestVector: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }
    jd = tt + T0;

    *outbody = body = LookupBody(name);
    if (body < 0)
    {
        fprintf(stderr, "CheckTestVector: Unknown body '%s' on line %d of file %s\n", name, lnum, filename);
        return 1;
    }

    error = NovasBodyPos(jd, body, npos);
    if (error)
    {
        fprintf(stderr, "CheckTestVector: NovasBodyPos returned %d on line %d of file %s\n", error, lnum, filename);
        return error;
    }

    error = PositionArcminError(body, jd, npos, xpos, arcmin);
    if (error)
    {
        fprintf(stderr, "CheckTestVector: PositionArcminError returned %d on line %d of file %s\n", error, lnum, filename);
        return error;
    }

    tolerance = (body == BODY_PLUTO) ? PLUTO_TOLERANCE_ARCMIN : 0.4;
    if (*arcmin > tolerance)
    {
        fprintf(stderr, "CheckTestVector: Excessive angular error (%lf arcmin) on line %d of file %s\n", *arcmin, lnum, filename);
        fprintf(stderr, "check   = (%22.16lf, %22.16lf, %22.16lf)\n", xpos[0], xpos[1], xpos[2]);
        fprintf(stderr, "correct = (%22.16lf, %22.16lf, %22.16lf)\n", npos[0], npos[1], npos[2]);
        return 1;
    }

    return 0;   /* success */
}

static int CheckSkyPos(observer *location, const char *filename, int lnum, const char *line, double *arcmin_equ, double *arcmin_hor, vsop_body_t *outbody)
{
    int body, error, bodyIndex;
    char name[10];
    double tt, ut, jd_tt, jd_utc, ra, dec, dist, azimuth, altitude;
    double check_zenith, check_azimuth, check_altitude, check_ra, check_dec;
    object obj;
    sky_pos sky_ast;    /* astrometric (RA,DEC) */
    sky_pos sky_dat;    /* (RA,DEC) using true equator and equinox of date*/
    double delta_ra, delta_dec, delta_az, delta_alt;
    double delta_t_seconds;
    double tolerance;

    *outbody = VSOP_INVALID_BODY;
    *arcmin_equ = *arcmin_hor = 99999.0;

    if (8 != sscanf(line, "s %9[A-Za-z] %lf %lf %lf %lf %lf %lf %lf", name, &tt, &ut, &ra, &dec, &dist, &azimuth, &altitude))
    {
        fprintf(stderr, "CheckSkyPos: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }

    delta_t_seconds = (24.0 * 3600.0) * (tt - ut);
    jd_tt = tt + T0;
    jd_utc = ut + T0;

    *outbody = body = LookupBody(name);
    if (body < 0)
    {
        fprintf(stderr, "CheckSkyPos: Unknown body '%s' on line %d of file %s\n", name, lnum, filename);
        return 1;
    }

    if (body == BODY_EARTH || body == BODY_EMB)
    {
        fprintf(stderr, "CheckSkyPos: Cannot calculate sky position of body '%s'\n", name);
        return 1;
    }

    if (body == BODY_GM)
        bodyIndex = 11;
    else if (body == BODY_SUN)
        bodyIndex = 10;
    else
        bodyIndex = 1 + body;

    error = make_object(0, (short)bodyIndex, name, NULL, &obj);
    if (error)
    {
        fprintf(stderr, "CheckSkyPos: make_object(%d) returned %d on line %d of file %s\n", bodyIndex, error, lnum, filename);
        return error;
    }

    /* First call to place: ask for astrometric coordinates (J2000) : coord_sys=3 */
    error = place(jd_tt, &obj, location, delta_t_seconds, 3, 1, &sky_ast);
    if (error)
    {
        fprintf(stderr, "CheckSkyPos: place(3) returned %d on line %d of file %s\n", error, lnum, filename);
        return error;
    }

    /* Calculate the overall angular error between the two (RA,DEC) pairs. */
    /* Convert both errors to arcminutes. */
    delta_dec = (sky_ast.dec - dec) * 60.0;

    delta_ra = fabs(sky_ast.ra - ra);
    if (delta_ra > 12.0)
    {
        /*
            Sometimes the two RA values can straddle the 24-hour mark.
            For example, one of them is 0.001 and the other 23.999.
            The actual error is then 0.002 hours, not 23.998.
            In general, it is never "fair" to call the error greater than
            12 hours or 180 degrees.
        */
        delta_ra = 24.0 - delta_ra;
    }
    delta_ra *= (15.0 * 60.0);

    /* Right Ascension errors are less significant as the declination approaches the poles. */
    /* For example, a 12-hour RA error for Polaris would not matter very much to its observed position. */
    /* So diminish the error measurement as appropriate for this declination. */
    delta_ra *= cos(sky_ast.dec * DEG2RAD);

    /* Calculate pythagorean error as if both were planar coordinates. */
    *arcmin_equ = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

    tolerance = (body == BODY_PLUTO) ? PLUTO_TOLERANCE_ARCMIN : 0.9;
    if (*arcmin_equ > tolerance)
    {
        fprintf(stderr, "CheckSkyPos: excessive (RA,DEC) error = %lf arcmin at line %d of file %s\n", *arcmin_equ, lnum, filename);
        return 1;
    }

    /* We have tested RA,DEC. Now measure the error in horizontal coordinates. */
    /* We have to call place() again, this time asking for equator and equinox of date (coord_sys=1). */
    error = place(jd_tt, &obj, location, delta_t_seconds, 1, 1, &sky_dat);
    if (error)
    {
        fprintf(stderr, "CheckSkyPos: place(1) returned %d on line %d of file %s\n", error, lnum, filename);
        return error;
    }

    equ2hor(jd_utc, jd_tt-jd_utc, 1, 0.0, 0.0,
        &location->on_surf, sky_dat.ra, sky_dat.dec, 0,
        &check_zenith, &check_azimuth, &check_ra, &check_dec);

    check_altitude = 90.0 - check_zenith;

    delta_az = fabs(azimuth - check_azimuth);
    if (delta_az > 180.0)
        delta_az = 360.0 - delta_az;
    delta_az *= 60.0;
    /* Just like RA, diminish the azimuth error as altitude approaches zenith/nadir... */
    delta_az *= cos(check_altitude * DEG2RAD);

    delta_alt = 60.0 * (altitude - check_altitude);

    *arcmin_hor = sqrt(delta_az*delta_az + delta_alt*delta_alt);

    tolerance = (body == BODY_PLUTO) ? PLUTO_TOLERANCE_ARCMIN : 0.9;
    if (*arcmin_hor > tolerance)
    {
        fprintf(stderr, "CheckSkyPos: excessive (az,alt) error = %lf arcmin for body %d at line %d of file %s\n", *arcmin_hor, body, lnum, filename);
        return 1;
    }

    return 0;
}

static vsop_body_t LookupBody(const char *name)
{
    if (!strcmp(name, "Mercury"))   return BODY_MERCURY;
    if (!strcmp(name, "Venus"))     return BODY_VENUS;
    if (!strcmp(name, "Earth"))     return BODY_EARTH;
    if (!strcmp(name, "Moon"))      return BODY_MOON;
    if (!strcmp(name, "Mars"))      return BODY_MARS;
    if (!strcmp(name, "Jupiter"))   return BODY_JUPITER;
    if (!strcmp(name, "Saturn"))    return BODY_SATURN;
    if (!strcmp(name, "Uranus"))    return BODY_URANUS;
    if (!strcmp(name, "Neptune"))   return BODY_NEPTUNE;
    if (!strcmp(name, "Pluto"))     return BODY_PLUTO;
    if (!strcmp(name, "Sun"))       return BODY_SUN;
    if (!strcmp(name, "EMB"))       return BODY_EMB;
    if (!strcmp(name, "GM"))        return BODY_GM;
    if (!strcmp(name, "SSB"))       return BODY_SSB;
    return BODY_INVALID;
}

typedef struct
{
    int    count;
    double sum;
    double max;
}
error_stat_t;

typedef struct
{
    error_stat_t    helio;
    error_stat_t    equ;
    error_stat_t    hor;
    error_stat_t    eclip;
}
error_bundle_t;

static void UpdateErrorStats(error_stat_t *stats, double arcmin)
{
    ++(stats->count);
    stats->sum += arcmin * arcmin;
    if (arcmin > stats->max)
        stats->max = arcmin;
}

static int PrintErrorStats(vsop_body_t body, const char *tag, const error_stat_t *stats, error_stat_t *tally)
{
    if (stats->count > 0)
    {
        const char *name;
        switch ((int)body)
        {
        case VSOP_MERCURY:  name = "Mercury";   break;
        case VSOP_VENUS:    name = "Venus";     break;
        case VSOP_EARTH:    name = "Earth";     break;
        case VSOP_EMB:      name = "EMB";       break;
        case VSOP_MARS:     name = "Mars";      break;
        case VSOP_JUPITER:  name = "Jupiter";   break;
        case VSOP_SATURN:   name = "Saturn";    break;
        case VSOP_URANUS:   name = "Uranus";    break;
        case VSOP_NEPTUNE:  name = "Neptune";   break;
        case VSOP_PLUTO:    name = "Pluto";     break;
        case VSOP_GM:       name = "GM";        break;
        case VSOP_SUN:      name = "Sun";       break;
        case VSOP_MOON:     name = "Moon";      break;
        case VSOP_SSB:      name = "SSB";       break;
        default:            name = "";          break;
        }
        printf("STATS(%-8s %s):  count= %6d  max= %10.8lf  rms= %10.8lf\n", name, tag, stats->count, stats->max, sqrt(stats->sum / stats->count));

        if (tally != NULL)
        {
            tally->count += stats->count;
            tally->sum += stats->sum;
            if (stats->max > tally->max)
                tally->max = stats->max;
        }

        return 1;   /* printed a line */
    }
    return 0;       /* did not print a line */
}

static int CheckTestOutput(const char *filename)
{
    int error, lnum;
    vsop_body_t body;
    FILE *infile = NULL;
    char line[200];
    double arcmin_helio, arcmin_eclip, arcmin_equ, arcmin_hor;
    error_bundle_t bundle[VSOP_BODY_LIMIT];
    error_stat_t tally;
    observer location;

    memset(&location, 0, sizeof(observer));
    memset(bundle, 0, sizeof(bundle));
    memset(&tally, 0, sizeof(tally));

    CHECK(OpenEphem());

    /* Check input file that contains lines of test output like this: */
    /* v Jupiter 2458455.2291666665 -2.0544644667646907 -4.271606485974493 -1.899398554516329 */
    infile = fopen(filename, "rt");
    if (infile == NULL)
    {
        fprintf(stderr, "CheckTestOutput: Cannot open file: %s\n", filename);
        error = 1;
        goto fail;
    }

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        switch (line[0])
        {
        case '#':
            break;  /* ignore debug output */

        case 'o':
            /* The observer used for all future sky position calculations */
            CHECK(ParseObserver(filename, lnum, line, &location));
            break;

        case 'v':   /* heliocentric cartesian vector represented in J2000 equatorial plane */
            CHECK(CheckTestVector(filename, lnum, line, &arcmin_helio, &body));
            UpdateErrorStats(&bundle[body].helio, arcmin_helio);
            break;

        case 's':   /* sky coordinates: RA, DEC, distance */
            CHECK(CheckSkyPos(&location, filename, lnum, line, &arcmin_equ, &arcmin_hor, &body));
            UpdateErrorStats(&bundle[body].equ, arcmin_equ);
            UpdateErrorStats(&bundle[body].hor, arcmin_hor);
            break;

        case 'e':
            CHECK(CheckEcliptic(filename, lnum, line, &arcmin_eclip, &body));
            UpdateErrorStats(&bundle[body].eclip, arcmin_eclip);
            break;

        default:
            fprintf(stderr, "CheckTestOutput: Invalid first character on line %d of file %s\n", lnum, filename);
            error = 1;
            goto fail;
        }
    }

    printf("CheckTestOutput: Verified %d lines of file %s\n", lnum, filename);

    if (Verbose)
    {
        int nprinted, nbodies=0;
        for (body=0; body < VSOP_BODY_LIMIT; ++body)
        {
            nprinted  = PrintErrorStats(body, "hel", &bundle[body].helio, &tally);
            nprinted += PrintErrorStats(body, "equ", &bundle[body].equ,   &tally);
            nprinted += PrintErrorStats(body, "hor", &bundle[body].hor,   &tally);
            nprinted += PrintErrorStats(body, "ecl", &bundle[body].eclip, &tally);
            if (nprinted > 0)
                ++nbodies;
        }

        if (nbodies > 1)
        {
            printf("---------------------------------------------------------------------\n");
            PrintErrorStats(-1, "ALL", &tally, NULL);
        }

        printf("\n");
    }

    error = 0;

fail:
    ephem_close();
    if (infile != NULL) fclose(infile);
    return error;
}

static int TestChebFunc(const void *context, double t, double f[CHEB_MAX_DIM])
{
    (void)context;
    f[0] = cos(t) - t;
    return 0;
}

static void DumpAlpha(const double alpha[CHEB_MAX_POLYS][CHEB_MAX_POLYS], int npoly)
{
    int j, k;

    for (j=0; j < npoly; ++j)
    {
        printf("alpha[%d][...] = ", j);
        for (k=0; k < npoly; ++k)
        {
            printf(" %15.12lf", alpha[j][k]);
        }
        printf("\n");
    }
}

static int UnitTestChebyshev(void)
{
    const int npoly = 5;
    const int ndimen = 1;
    const double t1 = 0.6;
    const double t2 = 0.8;
    int error;
    int i;
    double x, t, f_exact, f_approx;
    ChebEncoder encoder;
    double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS];

    error = ChebInit(&encoder, npoly, ndimen);
    if (error)
    {
        fprintf(stderr, "UnitTestChebyshev: ChebInit returned error %d\n", error);
        goto fail;
    }

    error = ChebGenerate(&encoder, TestChebFunc, NULL, t1, t2, coeff);
    if (error)
    {
        fprintf(stderr, "UnitTestChebyshev: ChebGenerate returned error %d\n", error);
        goto fail;
    }

    DumpAlpha(encoder.Alpha, encoder.NumPoly);

    /* Verify that the values are exact at all the sample points. */
    for (i=0; i < npoly; ++i)
    {
        x = (2.0 * i)/(npoly-1) - 1.0;
        t = t1 + (t2-t1)*((double)i / (npoly-1));
        TestChebFunc(NULL, t, &f_exact);
        ChebApprox(npoly, ndimen, coeff, x, &f_approx);
        printf("UnitTestChebyshev: i=%d, x=%5.2lf, t=%0.4lf, f_exact=%12.9lf, f_approx=%12.9lf, error=%lg\n", i, x, t, f_exact, f_approx, f_approx - f_exact);
    }

fail:
    return error;
}


static int DistancePlot(const char *name, double tt1, double tt2)
{
    int error;
    vsop_body_t body;
    double pos[3];
    double tt, dist;
    int i;
    const int npoints = 100000;

    CHECK(OpenEphem());

    body = LookupBody(name);
    if (body == BODY_INVALID)
        FAIL("DistancePlot: Invalid body name '%s'\n", name);

    printf("\"tt\",\"distance\",\"x\",\"y\",\"z\"\n");
    for (i=0; i < npoints; ++i)
    {
        tt = tt1 + (((double)i)/((double)(npoints-1)) * (tt2 - tt1));
        CHECK(NovasBodyPos(tt + T0, body, pos));
        dist = VectorLength(pos);
        printf("%0.16lf,%0.16lg,%0.16lg,%0.16lg,%0.16lg\n", tt, dist, pos[0], pos[1], pos[2]);
    }

    error = 0;
fail:
    ephem_close();
    return error;
}


static int DeltaTimePlot(const char *outFileName)
{
    int error = 1;
    FILE *outfile = NULL;
    const int minYear = 1500;
    const int maxYear = 2500;
    int year;
    double dt;

    outfile = fopen(outFileName, "wt");
    if (outfile == NULL)
        FAIL("DeltaTimePlot: Cannot open output file: %s\n", outFileName);

    fprintf(outfile, "\"year\",\"delta_t\"\n");
    for (year=minYear; year <= maxYear; ++year)
    {
        dt = ExtrapolatedDeltaT(year);
        fprintf(outfile, "%d,%lf\n", year, dt);
    }

    printf("DeltaTimePlot: Wrote file %s\n", outFileName);
    error = 0;
fail:
    if (outfile != NULL) fclose(outfile);
    return error;
}

/*-----------------------------------------------------------------------------------------*/

static const char *TopDataFileName = "TOP2013.dat";
static const double TopMillenniaAroundJ2000 = 0.5;   /* optimize for calculations within 500 years of J2000. */
static const double TopThresholdArcmin = 0.5;

static int CalcTop2013(FILE *outfile, const top_model_t *model)
{
    int error = 1;
    double jd, tt;
    int i;
    top_elliptical_t    ellip;
    top_rectangular_t   ecl, equ;
    double pos[3];
    double range, diff, dx, dy, dz;
    double max_diff = 0.0;

    for (i=0; i <= 10; ++i)
    {
        tt = 4000.0 * (i-10);
        jd = tt + 2451545.0;

        fprintf(outfile, "%d %0.1lf\n", model->planet, jd);

        CHECK(TopCalcElliptical(model, tt, &ellip));
        fprintf(outfile, " %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf\n", ellip.a, ellip.lambda, ellip.k, ellip.h, ellip.q, ellip.p);

        CHECK(TopEcliptic(model->planet, &ellip, &ecl));
        fprintf(outfile, " %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf\n", ecl.x, ecl.y, ecl.z, ecl.vx, ecl.vy, ecl.vz);

        CHECK(TopEquatorial(&ecl, &equ));
        fprintf(outfile, " %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf\n", equ.x, equ.y, equ.z, equ.vx, equ.vy, equ.vz);

        /* Compare the equatorial vector against NOVAS calculations. */
        CHECK(NovasBodyPos(jd, model->planet - 1, pos));
        range = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        dx = pos[0] - equ.x;
        dy = pos[1] - equ.y;
        dz = pos[2] - equ.z;
        diff = 60.0 * RAD2DEG * (sqrt(dx*dx + dy*dy + dz*dz) / range);
        if (diff > max_diff)
            max_diff = diff;

        DEBUG("planet %d  jd %10.1lf  arcmin %10.6lf\n", model->planet, jd, diff);
    }

    printf("CalcTop2013: planet=%d, max arcmin error=%lg\n", model->planet, max_diff);

    error = 0;
fail:
    return error;
}

static int ValidateTop2013(void)
{
    int error = 1;
    top_model_t model;
    FILE *mirror = NULL;
    FILE *outfile = NULL;
    int nlines, planet;
    const char *mirrorFileName = "TOP2013.mirror";
    const char *outFileName = "top2013/calc.txt";

    TopInitModel(&model);
    CHECK(OpenEphem());

    outfile = fopen(outFileName, "wt");
    if (outfile == NULL)
        FAIL("ValidateTop2013: Cannot open output file: %s\n", outFileName);

    mirror = fopen(mirrorFileName, "wt");
    if (mirror == NULL)
        FAIL("ValidateTop2013: cannot open output file: %s\n", mirrorFileName);

    for (planet=5; planet <= 9; ++planet)
    {
        CHECK(TopLoadModel(&model, TopDataFileName, planet));
        CHECK(TopWriteModel(&model, mirror));
        CHECK(CalcTop2013(outfile, &model));
        TopFreeModel(&model);
    }

    fclose(mirror);
    mirror = NULL;

    fclose(outfile);
    outfile = NULL;

    /* Verify that the saved model exactly matches what we loaded. */
    CHECK(Diff(TopDataFileName, mirrorFileName, &nlines));
    if (nlines != 336806)
        FAIL("ValidateTop2013(%s): incorrect number of matching lines = %d\n", mirrorFileName, nlines);

    /* Verify that the calculations exactly match those produced by the original FORTRAN code. */
    CHECK(Diff(outFileName, "top2013/correct.txt", NULL));

    /* Clean up after success. */
    remove(mirrorFileName);
    remove(outFileName);
    printf("ValidateTop2013: PASS\n");
    error = 0;

fail:
    ephem_close();
    if (mirror) fclose(mirror);
    if (outfile) fclose(outfile);
    TopFreeModel(&model);
    return error;
}


static int Diff(const char *filename1, const char *filename2, int *nlines)
{
    int lnum = 0;
    int error = 1;
    FILE *infile1 = NULL;
    FILE *infile2 = NULL;
    char line1[200];
    char line2[200];

    infile1 = fopen(filename1, "rt");
    if (infile1 == NULL)
        FAIL("Diff: cannot open input file: %s\n", filename1);

    infile2 = fopen(filename2, "rt");
    if (infile2 == NULL)
        FAIL("Diff: cannot open input file: %s\n", filename2);

    for(;;)
    {
        char *r1 = fgets(line1, sizeof(line1), infile1);
        char *r2 = fgets(line2, sizeof(line2), infile2);
        if (r1 == NULL && r2 == NULL)
            break;
        if (r1 == NULL)
            FAIL("Diff: %s is shorter than %s\n", filename1, filename2);
        if (r2 == NULL)
            FAIL("Diff: %s is shorter than %s\n", filename2, filename1);
        ++lnum;
        if (strcmp(line1, line2))
        {
            fprintf(stderr, "[%s]\n", line1);
            fprintf(stderr, "[%s]\n", line2);
            FAIL("Diff(%d): lines are different.\n", lnum);
        }
    }

    printf("Diff: %d identical lines in files %s and %s\n", lnum, filename1, filename2);
    error = 0;
fail:
    if (infile1 != NULL) fclose(infile1);
    if (infile2 != NULL) fclose(infile2);
    if (nlines != NULL) *nlines = lnum;
    return error;
}


static int TopCalc(const char *name, const char *date)
{
    int error = 1;
    top_model_t model;
    vsop_body_t body;
    top_rectangular_t equ;
    int planet;
    double tt;
    TopInitModel(&model);

    CHECK(ParseDate(date, &tt));

    body = LookupBody(name);
    if (body < 0)
        FAIL("TopCalc: planet name '%s' is not valid.\n", name);

    if (body < BODY_JUPITER || body > BODY_PLUTO)
        FAIL("TopCalc: TOP2013 supports Jupiter through Pluto only.\n");

    planet = body + 1;      /* convert our body ID into TOP2013 planet ID */

    CHECK(TopLoadModel(&model, TopDataFileName, planet));
    CHECK(TopPosition(&model, tt, &equ));
    printf("tt  = %0.8lf\n", tt);
    printf("pos = [%22.16lf, %22.16lf, %22.16lf]\n", equ.x, equ.y, equ.z);
    printf("vel = [%22.16lf, %22.16lf, %22.16lf]\n", equ.vx, equ.vy, equ.vz);

    error = 0;

fail:
    TopFreeModel(&model);
    return error;
}


static int TopFileInfo(const char *filename, const char *name)
{
    int error = 1;
    top_model_t model;
    vsop_body_t body;
    int planet, f;

    TopInitModel(&model);

    body = LookupBody(name);
    if (body < 0)
        FAIL("TopFileInfo: planet name '%s' is not valid.\n", name);

    planet = body + 1;      /* convert our body ID into TOP2013 planet ID */

    CHECK(TopLoadModel(&model, filename, planet));
    printf("%6d [", TopTermCount(&model));
    for (f=0; f < TOP_NCOORDS; ++f)
        printf("%6d", TopTermCountF(&model, f));
    printf("]\n");

    error = 0;
fail:
    TopFreeModel(&model);
    return error;
}

/*------------------------------------------------------------------------------------------------*/

static int ParseDate(const char *text, double *tt)
{
    int error = 1;
    int nscanned;
    int year, month, day, hour, minute;
    double second, float_hours;

    /* Allow passing in an explicit tt value directly. */
    if (text[0] == '@' && 1 == sscanf(text+1, "%lf", tt))
        return 0;

    nscanned = sscanf(text, "%d-%d-%dT%d:%d:%lfZ", &year, &month, &day, &hour, &minute, &second);
    if (nscanned < 5)
        FAIL("ParseDate: text '%s' is not valid. Must be formatted as yyyy-mm-dd[Thh:mm[:ss.sss]Z].\n", text);
    if (nscanned < 6)
        second = 0.0;

    float_hours = hour + minute/60.0 + second/3600.0;
    *tt = julian_date((short)year, (short)month, (short)day, float_hours) - T0;
    error = 0;
fail:
    return error;
}


static const char *QuotedBodyName(int body)
{
    switch (body)
    {
    case BODY_SUN:      return "\"Sun\"    ";
    case BODY_MERCURY:  return "\"Mercury\"";
    case BODY_VENUS:    return "\"Venus\"  ";
    case BODY_EMB:      return "\"Earth\"  ";
    case BODY_MARS:     return "\"Mars\"   ";
    case BODY_JUPITER:  return "\"Jupiter\"";
    case BODY_SATURN:   return "\"Saturn\" ";
    case BODY_URANUS:   return "\"Uranus\" ";
    case BODY_NEPTUNE:  return "\"Neptune\"";
    case BODY_PLUTO:    return "\"Pluto\"  ";
    default:            return "\"\"       ";
    }
}

/*------------------------------------------------------------------------------------------------*/
