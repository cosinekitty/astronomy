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

#define CHECK(x)   do{if(0 != (error = (x))) goto fail;}while(0)

static double jd_begin;
static double jd_end;
static short int de_number;

static int OpenEphem(void);
static int PrintUsage(void);
static int GenerateAllSource(void);
static int FastGenerate(void);
static int TestVsopModel(vsop_model_t *model, int body, double threshold, double *max_arcmin, int *trunc_terms);
static int SaveVsopFile(const vsop_model_t *model);
static int PositionArcminError(int body, double jd, const double a[3], const double b[3], double *arcmin);
static double VectorLength(const double v[3]);
static int MeasureError(const char *inFileName, int nsamples, error_stats_t *stats);
static int Resample(int body, const char *outFileName, int npoly, int startYear, int stopYear, int nsections);
static int SampleFunc(const void *context, double jd, double pos[CHEB_MAX_DIM]);
static double VectorError(double a[3], double b[3]);
static int ManualResample(int body, int npoly, int nsections, int startYear, int stopYear);
static int CheckTestOutput(const char *filename);
static vsop_body_t LookupBody(const char *name);
static int CheckSkyPos(observer *location, const char *filename, int lnum, const char *line, double *arcmin);

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

    if (argc == 2 && !strcmp(argv[1], "fast"))
        return FastGenerate();

    if (argc == 3 && !strcmp(argv[1], "check"))
        return CheckTestOutput(argv[2]);

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
        "generate check testfile\n"
        "    Verify the calculations in the testfile generated by a unit test.\n"
        "\n"
        "generate fast\n"
        "    Generate source code, re-using existing search data from 'generate all'.\n"
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
    double jed[2];
    double sun_pos[3], moon_pos[3], vel[3];

    jed[0] = jd;
    jed[1] = 0.0;

    if (body == BODY_EARTH || body == BODY_MOON)
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
        error = (int)state(jed, (short)body, pos, vel);
        if (error)
        {
            fprintf(stderr, "NovasBodyPos: state(%lf, %d) returned %d\n", jd, body, error);
            return error;
        }

        /* Special case: geocentric moon should not be converted to heliocentric coordinates. */
        if (body == BODY_GM) return 0;
    }

    error = (int)state(jed, BODY_SUN, sun_pos, vel);
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

    CHECK(OpenEphem());
    CHECK(BuildVsopData());
    CHECK(ManualResample(8, 19, 7, MIN_YEAR, MAX_YEAR));
    CHECK(FastGenerate());

fail:
    ephem_close();
    return error;
}

static int FastGenerate(void)
{
    return GenerateCode("../source/js/astronomy.js", "template/astronomy.js", "output");
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

static int ManualResample(int body, int npoly, int nsections, int startYear, int stopYear)
{
    int error, nsamples;
    double daysPerSection;
    error_stats_t stats;
    char filename[80];

    snprintf(filename, sizeof(filename), "output/%02d.eph", body);
    CHECK(Resample(body, filename, npoly, startYear, stopYear, nsections));

    daysPerSection = (365.244 * (stopYear - startYear + 1)) / nsections;
    nsamples = (int)(24.0 * (60.0/10.0) * daysPerSection);    /* approximately one sample per 10 minutes */
    if (nsamples < 100*npoly)
        nsamples = 100*npoly;     /* but always oversample every segment beyond the polynomial count */

    CHECK(MeasureError(filename, nsamples, &stats));
    printf("Chebyshev body=%d, RMS=%0.3lf, max=%0.3lf, data=%ld\n", body, stats.rmsArcminError, stats.maxArcminError, stats.dataCount);

fail:
    return error;
}

static int MeasureError(const char *inFileName, int nsamples, error_stats_t *stats)
{
    int i, nvectors;
    int error = 0;
    double jd, x, sumPosDiff, posDiff, maxPosDiff;
    double arcmin, arcminSum;
    double apos[3];     /* our approximate position vector */
    double npos[3];     /* NOVAS exact position vector */
    eph_file_reader_t reader;
    eph_record_t record;

    memset(stats, 0, sizeof(error_stats_t));
    error = EphFileOpen(&reader, inFileName);
    if (error) 
    {
        fprintf(stderr, "MeasureError: Error %d trying to open file: %s\n", error, inFileName);
        goto fail;
    }

    /* sanity check on the body value, to prevent invalid memory access */
    if (reader.body < 0 || reader.body > 10)
    {
        fprintf(stderr, "MeasureError: Invalid body = %d\n", reader.body);
        error = 1;
        goto fail;
    }

    nvectors = 0;
    sumPosDiff = 0.0;
    maxPosDiff = 0.0;
    arcminSum = 0.0;
    while (EphReadRecord(&reader, &record))
    {
        stats->dataCount += (1 + 3*record.numpoly);     /* tally the beginning Julian Date plus all the Chebyshev coefficients */

        /*
            We sample in such a way as to exclude the endpoints, because
            we know error will be minimal there (and at all the other Chebyshev nodes).
            If nsamples==1, we pick the midpoint.
            If nsamples==2, we pick the 1/3 mark and 2/3 marks.
            Thus we always divide the interval into (1+nsamples) sections.
        */
        for (i=1; i <= nsamples; ++i)
        {
            jd = record.jdStart + (i * record.jdDelta)/(1 + nsamples);

            /* Calculate "correct" position at the time 'jd' (according to NOVAS). */
            error = NovasBodyPos(jd, reader.body, npos);
            if (error) goto fail;
            
            /* Calculate our approximation of the position at the time 'jd'. */
            x = ChebScale(record.jdStart, record.jdStart + record.jdDelta, jd);
            ChebApprox(record.numpoly, 3, record.coeff, x, apos);

            /* Tally the root-mean-square error between the two position vectors. */
            posDiff = VectorError(npos, apos);
            if (posDiff > maxPosDiff)
                maxPosDiff = posDiff;

            /* Estimate angular error */
            error = PositionArcminError(reader.body, jd, npos, apos, &arcmin);
            if (error) goto fail;
            arcminSum += arcmin * arcmin;
            if (arcmin > stats->maxArcminError)
                stats->maxArcminError = arcmin;

            sumPosDiff += posDiff;
            ++nvectors;
        }
    }
    error = record.error;
    if (error)
    {
        fprintf(stderr, "MeasureError: Error %d reading record on line %d of file %s\n", error, reader.lnum, inFileName);
        goto fail;
    }

    if (nvectors < 1)
    {
        fprintf(stderr, "MeasureError(ERROR): nvectors = %d\n", nvectors);
        error = 1;
        goto fail;
    }

    stats->rmsPositionError = sqrt(sumPosDiff / nvectors);
    stats->maxPositionError = sqrt(maxPosDiff);
    stats->rmsArcminError = sqrt(arcminSum / nvectors);

fail:
    EphFileClose(&reader);
    return error;
}

static int Resample(int body, const char *outFileName, int npoly, int startYear, int stopYear, int nsections)
{
    FILE *outfile = NULL;
    ChebEncoder encoder;
    double coeff[CHEB_MAX_DIM][CHEB_MAX_POLYS];    /* independent Chebyshev coefficients for each position coordinate: x, y, z */
    int error, k, i, section;
    double jdStart, jdStop, jdDelta, jd;
    sample_context_t context;

    if (startYear < MIN_YEAR || stopYear < startYear || MAX_YEAR < stopYear)
    {
        fprintf(stderr, "ERROR: Invalid year range %d..%d\n", startYear, stopYear);
        error = 1;
        goto fail;
    }

    if (nsections < 1)
    {
        fprintf(stderr, "ERROR: nsections must be a positive integer.\n");
        error = 2;
        goto fail;
    }

    error = ChebInit(&encoder, npoly, VECTOR_DIM);
    if (error) 
    {
        fprintf(stderr, "ERROR %d returned by ChebInit()\n", error);
        goto fail;
    }

    outfile = fopen(outFileName, "wt");
    if (outfile == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open output file '%s'\n", outFileName);
        error = 3;
        goto fail;
    }

    /* Write header to output file. */
    fprintf(outfile, "body=%d\n", body);
    fprintf(outfile, "\n");     /* blank line terminates the header */

    jdStart = julian_date((short)startYear, 1, 1, 0.0) - 1.0;   /* subtract 1 day for light travel time allowance */
    jdStop = julian_date((short)(stopYear+1), 1, 1, 0.0);
    jdDelta = (jdStop - jdStart) / nsections;
    context.body = body;
    for (section = 0; section < nsections; ++section)
    {
        jd = jdStart + (section * jdDelta);

        /* encode coefficients for vector-valued position function */
        error = ChebGenerate(&encoder, SampleFunc, &context, jd, jd + jdDelta, coeff);
        if (error)
        {
            fprintf(stderr, "ERROR: ChebGenerate() returned %d\n", error);
            goto fail;
        }

        /* Prefix each block of polynomial coefficients with the time range and the number of polynomials. */
        fprintf(outfile, "%0.18lf %0.18lf %d\n", jd, jdDelta, npoly);

        /* Write the Chebyshev coordinates we obtained for the position function. */
        for (k=0; k < npoly; k++)
        {
            for (i=0; i < 3; ++i)
                fprintf(outfile, "%22.18lf", coeff[i][k]);
            fprintf(outfile, "\n");
        }
    }

fail:
    if (outfile != NULL)
        fclose(outfile);

    if (error)
        unlink(outFileName);

    return error;
}

static int SampleFunc(const void *context, double jd, double pos[CHEB_MAX_DIM])
{
    const sample_context_t *c = context;
    int error;
    double jed[2];
    double sun_pos[3];
    double vel[3];      /* we don't care about velocities... ignored */
    
    jed[0] = jd;
    jed[1] = 0.0;
    error = state(jed, (short)c->body, pos, vel);
    if (error) return error;

    error = state(jed, BODY_SUN, sun_pos, vel);
    if (error) return error;

    /* Calculate heliocentric coordinates from barycenric coordinates. */
    pos[0] -= sun_pos[0];
    pos[1] -= sun_pos[1];
    pos[2] -= sun_pos[2];

    return 0;
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

static int CheckTestVector(const char *filename, int lnum, const char *line, double *arcmin)
{
    int body, error;
    char name[10];
    double tt, jd, xpos[3], npos[3];

    *arcmin = 99999.0;

    if (5 != sscanf(line, "v %9[A-Za-z] %lf %lf %lf %lf", name, &tt, &xpos[0], &xpos[1], &xpos[2]))
    {
        fprintf(stderr, "CheckTestVector: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }
    jd = tt + T0;

    body = LookupBody(name);
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

    if (*arcmin > 0.4)
    {
        fprintf(stderr, "CheckTestVector: Excessive angular error (%lf arcmin) on line %d of file %s\n", *arcmin, lnum, filename);
        fprintf(stderr, "check   = (%22.16lf, %22.16lf, %22.16lf)\n", xpos[0], xpos[1], xpos[2]);
        fprintf(stderr, "correct = (%22.16lf, %22.16lf, %22.16lf)\n", npos[0], npos[1], npos[2]);
        return 1;
    }

    return 0;   /* success */
}

static int CheckSkyPos(observer *location, const char *filename, int lnum, const char *line, double *arcmin)
{
    int body, error, bodyIndex;
    char name[10];
    double tt, ut, jd_tt, jd_utc, ra, dec, dist, azimuth, altitude;
    double check_zenith, check_azimuth, check_altitude, check_ra, check_dec;
    object obj;
    sky_pos sky_ast;    /* astrometric (RA,DEC) */
    sky_pos sky_dat;    /* (RA,DEC) using true equator and equinox of date*/
    double delta_ra, delta_dec, delta_az, delta_alt;

    *arcmin = 99999.0;

    if (8 != sscanf(line, "s %9[A-Za-z] %lf %lf %lf %lf %lf %lf %lf", name, &tt, &ut, &ra, &dec, &dist, &azimuth, &altitude))
    {
        fprintf(stderr, "CheckSkyPos: Invalid format on line %d of file %s\n", lnum, filename);
        return 1;
    }
    jd_tt = tt + T0;
    jd_utc = ut + T0;

    body = LookupBody(name);
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
    error = place(jd_tt, &obj, location, jd_tt-jd_utc, 3, 1, &sky_ast);
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
    *arcmin = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

    if (*arcmin > 0.9)
    {
        fprintf(stderr, "CheckSkyPos: excessive (RA,DEC) error = %lf arcmin at line %d of file %s\n", *arcmin, lnum, filename);
        return 1;
    }

    /* We have tested RA,DEC. Now measure the error in horizontal coordinates. */
    /* We have to call place() again, this time asking for equator and equinox of date (coord_sys=1). */
    error = place(jd_tt, &obj, location, jd_tt-jd_utc, 1, 1, &sky_dat);
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

    /* Replace the (RA,DEC) error with the (az,alt) error. */
    *arcmin = sqrt(delta_az*delta_az + delta_alt*delta_alt);   

    if (*arcmin > 0.9)
    {
        fprintf(stderr, "CheckSkyPos: excessive (az,alt) error = %lf arcmin for body %d at line %d of file %s\n", *arcmin, body, lnum, filename);
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
    return BODY_INVALID;
}

static int CheckTestOutput(const char *filename)
{
    int error, lnum;
    FILE *infile = NULL;
    char line[200];
    double arcmin, max_arcmin = 0.0;
    observer location;

    memset(&location, 0, sizeof(observer));

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

        case 'v':   /* cartesian vector represented in J2000 equatorial plane */
            CHECK(CheckTestVector(filename, lnum, line, &arcmin));
            if (arcmin > max_arcmin)
                max_arcmin = arcmin;
            break;

        case 's':   /* sky coordinates: RA, DEC, distance */
            CHECK(CheckSkyPos(&location, filename, lnum, line, &arcmin));
            if (arcmin > max_arcmin)
                max_arcmin = arcmin;
            break;
        
        default:
            fprintf(stderr, "CheckTestOutput: Invalid first character on line %d of file %s\n", lnum, filename);
            error = 1;
            goto fail;
        }
    }

    printf("CheckTestOutput: Verified %d lines of file %s : max error = %lf arcmin\n", lnum, filename, max_arcmin);
    error = 0;

fail:
    ephem_close();
    if (infile != NULL) fclose(infile);
    return error;
}
