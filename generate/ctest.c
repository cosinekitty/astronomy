/*
    ctest.c  -  Don Cross <cosinekitty.com>

    C langauge unit test for Astronomy Engine project.
    https://cosinekitty.github.io/astronomy
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "astronomy.h"

#define PI      3.14159265358979323846

#define CHECK(x)            do{if(0 != (error = (x))) goto fail;}while(0)

static int CheckVector(int lnum, astro_vector_t v)
{
    if (v.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "FAILURE at ctest.c[%d]: vector status = %d\n", lnum, v.status);
        return 1;
    }
    return 0;
}

static int CheckEquator(int lnum, astro_equatorial_t equ)
{
    if (equ.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "FAILURE at ctest.c[%d]: equatorial status = %d\n", lnum, equ.status);
        return 1;
    }
    return 0;
}

#define CHECK_VECTOR(var,expr)   CHECK(CheckVector(__LINE__, ((var) = (expr))))
#define CHECK_EQU(var,expr)      CHECK(CheckEquator(__LINE__, ((var) = (expr))))

static int Test_AstroTime(void);
static int AstroCheck(void);
static int Diff(const char *c_filename, const char *js_filename);
static int DiffLine(int lnum, const char *cline, const char *jline, double *maxdiff, int *worst_lnum);
static int SeasonsTest(const char *filename);
static int MoonPhase(const char *filename);

int main(int argc, const char *argv[])
{
    int error = 1;

    if (argc == 1)
    {
        CHECK(Test_AstroTime());
        CHECK(AstroCheck());
        goto success;
    }

    if (argc == 3)
    {
        const char *verb = argv[1];
        const char *filename = argv[2];

        if (!strcmp(verb, "seasons"))
        {
            CHECK(SeasonsTest(filename));
            goto success;
        }

        if (!strcmp(verb, "moonphase"))
        {
            CHECK(MoonPhase(filename));
            goto success;
        }
    }

    if (argc == 4)
    {
        if (!strcmp(argv[1], "diff"))
        {
            const char *c_filename = argv[2];
            const char *js_filename = argv[3];
            CHECK(Diff(c_filename, js_filename));
            goto success;
        }
    }

    fprintf(stderr, "Invalid command line arguments.\n");
    error = 1;
    goto fail;

success:
    error = 0;

fail:
    fprintf(stderr, "ctest exiting with %d\n", error);
    return error;
}

static int Test_AstroTime(void)
{
    astro_time_t time;
    const double expected_ut = 6910.270978506945;
    const double expected_tt = 6910.271779431480;
    double diff;

    time = Astronomy_MakeTime(2018, 12, 2, 18, 30, 12.543);
    printf("Test_AstroTime: ut=%0.6lf, tt=%0.6lf\n", time.ut, time.tt);

    diff = time.ut - expected_ut;
    if (fabs(diff) > 1.0e-12)
    {
        fprintf(stderr, "Test_AstroTime: excessive UT error %lg\n", diff);
        return 1;
    }

    diff = time.tt - expected_tt;
    if (fabs(diff) > 1.0e-12)
    {
        fprintf(stderr, "Test_AstroTime: excessive TT error %lg\n", diff);
        return 1;
    }

    return 0;
}

static int AstroCheck(void)
{
    int error = 1;
    FILE *outfile = NULL;
    const char *filename = "temp/c_check.txt";
    astro_time_t time;
    astro_time_t stop;
    astro_body_t body;
    astro_vector_t pos;
    astro_equatorial_t j2000;
    astro_equatorial_t ofdate;
    astro_horizon_t hor;
    astro_observer_t observer = Astronomy_MakeObserver(29.0, -81.0, 10.0);
    int b;
    static const astro_body_t bodylist[] =  /* match the order in the JavaScript unit test */
    {
        BODY_SUN, BODY_MERCURY, BODY_VENUS, BODY_EARTH, BODY_MARS, 
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO
    };
    static int nbodies = sizeof(bodylist) / sizeof(bodylist[0]);

    outfile = fopen(filename, "wt");
    if (outfile == NULL)
    {
        fprintf(stderr, "AstroCheck: Cannot open output file: %s\n", filename);
        error = 1;
        goto fail;
    }

    fprintf(outfile, "o %lf %lf %lf\n", observer.latitude, observer.longitude, observer.height);

    time = Astronomy_MakeTime(1700, 1, 1, 0, 0, 0.0);
    stop = Astronomy_MakeTime(2200, 1, 1, 0, 0, 0.0);
    while (time.tt < stop.tt)
    {
        for (b=0; b < nbodies; ++b)
        {
            body = bodylist[b];
            CHECK_VECTOR(pos, Astronomy_HelioVector(body, time));
            fprintf(outfile, "v %s %0.16lf %0.16lf %0.16lf %0.16lf\n", Astronomy_BodyName(body), pos.t.tt, pos.x, pos.y, pos.z);

            if (body != BODY_EARTH)
            {
                CHECK_EQU(j2000, Astronomy_Equator(body, time, observer, 0, 0));
                CHECK_EQU(ofdate, Astronomy_Equator(body, time, observer, 1, 1));
                hor = Astronomy_Horizon(time, observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
                fprintf(outfile, "s %s %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf\n", 
                    Astronomy_BodyName(body), time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude);
            }
        }

        CHECK_VECTOR(pos, Astronomy_GeoVector(BODY_MOON, time, 0));
        fprintf(outfile, "v GM %0.16lf %0.16lf %0.16lf %0.16lf\n", pos.t.tt, pos.x, pos.y, pos.z);

        CHECK_EQU(j2000, Astronomy_Equator(BODY_MOON, time, observer, 0, 0));
        CHECK_EQU(ofdate, Astronomy_Equator(BODY_MOON, time, observer, 1, 1));
        hor = Astronomy_Horizon(time, observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
        fprintf(outfile, "s GM %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf\n", 
            time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude);

        time = Astronomy_AddDays(time, 10.0 + PI/100.0);
    }

fail:
    if (outfile != NULL)
        fclose(outfile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int Diff(const char *c_filename, const char *js_filename)
{
    int error = 1;
    int lnum;
    FILE *cfile = NULL;
    FILE *jfile = NULL;
    char cline[200];
    char jline[200];
    char *cread;
    char *jread;
    double maxdiff = 0.0;
    int worst_lnum = 0;

    cfile = fopen(c_filename, "rt");
    if (cfile == NULL)
    {
        fprintf(stderr, "ctest(Diff): Cannot open input file: %s\n", c_filename);
        error = 1;
        goto fail;
    }

    jfile = fopen(js_filename, "rt");
    if (jfile == NULL)
    {
        fprintf(stderr, "ctest(Diff): Cannot open input file: %s\n", js_filename);
        error = 1;
        goto fail;
    }

    lnum = 0;
    for(;;)
    {
        cread = fgets(cline, sizeof(cline), cfile);
        jread = fgets(jline, sizeof(jline), jfile);
        if (cread==NULL && jread==NULL)
            break;      /* normal end of both files */
        
        if (cread==NULL || jread==NULL)
        {
            fprintf(stderr, "ctest(Diff): Files do not have same number of lines: %s and %s\n", c_filename, js_filename);
            error = 1;
            goto fail;
        }

        ++lnum;
        CHECK(DiffLine(lnum, cline, jline, &maxdiff, &worst_lnum));
    }

    printf("ctest(Diff): Maximum numeric difference = %lg, worst line number = %d\n", maxdiff, worst_lnum);
    if (maxdiff > 1.8e-12)
    {
        fprintf(stderr, "ERROR: Excessive error comparing files %s and %s\n", c_filename, js_filename);
        error = 1;
        goto fail;
    }

    error = 0;

fail:
    if (cfile != NULL) fclose(cfile);
    if (jfile != NULL) fclose(jfile);
    return error;
}

static int DiffLine(int lnum, const char *cline, const char *jline, double *maxdiff, int *worst_lnum)
{
    int error = 1;
    char cbody[10];
    char jbody[10];
    double cdata[7];
    double jdata[7];
    double diff;
    int i, nc, nj, nrequired = -1;

    /* be paranoid: make sure we can't possibly have a fake match. */
    memset(cdata, 0xdc, sizeof(cdata));
    memset(jdata, 0xce, sizeof(jdata));

    /* Make sure the two data records are the same type. */
    if (cline[0] != jline[0])
    {
        fprintf(stderr, "ctest(DiffLine): Line %d mismatch record type: '%c' vs '%c'.\n", lnum, cline[0], jline[0]);
        error = 1;
        goto fail;
    }

    switch (cline[0])
    {
    case 'o':       /* observer */
        nc = sscanf(cline, "o %lf %lf %lf", &cdata[0], &cdata[1], &cdata[2]);
        nj = sscanf(jline, "o %lf %lf %lf", &jdata[0], &jdata[1], &jdata[2]);
        cbody[0] = jbody[0] = '\0';
        nrequired = 3;
        break;

    case 'v':       /* heliocentric vector */
        nc = sscanf(cline, "v %9[A-Za-z] %lf %lf %lf", cbody, &cdata[0], &cdata[1], &cdata[2]);
        nj = sscanf(jline, "v %9[A-Za-z] %lf %lf %lf", jbody, &jdata[0], &jdata[1], &jdata[2]);
        nrequired = 4;
        break;

    case 's':       /* sky coords: ecliptic and horizontal */
        nc = sscanf(cline, "s %9[A-Za-z] %lf %lf %lf %lf %lf %lf %lf", cbody, &cdata[0], &cdata[1], &cdata[2], &cdata[3], &cdata[4], &cdata[5], &cdata[6]);
        nj = sscanf(jline, "s %9[A-Za-z] %lf %lf %lf %lf %lf %lf %lf", jbody, &jdata[0], &jdata[1], &jdata[2], &jdata[3], &jdata[4], &jdata[5], &jdata[6]);
        nrequired = 8;
        break;

    default:
        fprintf(stderr, "ctest(DiffLine): Line %d type '%c' is not a valid record type.\n", lnum, cline[0]);
        error = 1;
        goto fail;
    }

    if (nc != nj)
    {
        fprintf(stderr, "ctest(DiffLine): Line %d mismatch data counts: %d vs %d\n", lnum, nc, nj);
        error = 1;
        goto fail;
    }

    if (nc != nrequired)
    {
        fprintf(stderr, "ctest(DiffLine): Line %d incorrect number of scanned arguments: %d\n", lnum, nc);
        error = 1;
        goto fail;
    }

    if (strcmp(cbody, jbody))
    {
        fprintf(stderr, "ctest(DiffLine): Line %d body mismatch: '%s' vs '%s'\n.", lnum, cbody, jbody);
        error = 1;
        goto fail;
    }

    if (cbody[0])
    {
        /* This is one of the record types that contains a body name. */
        /* Therefore, we need to correct the number of numeric data. */
        --nrequired;
    }

    /* Verify all the numeric data are very close. */
    for (i=0; i < nrequired; ++i)
    {
        diff = fabs(cdata[i] - jdata[i]);
        if (diff > *maxdiff)
        {
            *maxdiff = diff;
            *worst_lnum = lnum;
        }
    }
    error = 0;

fail:
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int SeasonsTest(const char *filename)
{
    int error = 1;
    int lnum;
    FILE *infile = NULL;
    char line[200];
    int nscanned, year, month, day, hour, minute;
    int current_year = 0;
    char name[20];
    astro_time_t correct_time;
    astro_time_t calc_time;
    astro_seasons_t seasons;
    double diff_minutes, max_minutes = 0.0;
    int mar_count=0, jun_count=0, sep_count=0, dec_count=0;

    memset(&seasons, 0, sizeof(seasons));    

    infile = fopen(filename, "rt");
    if (infile == NULL)
    {
        fprintf(stderr, "SeasonsTest: Cannot open input file: %s\n", filename);
        error = 1;
        goto fail;
    }

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        /*
            2019-01-03T05:20Z Perihelion
            2019-03-20T21:58Z Equinox
            2019-06-21T15:54Z Solstice
            2019-07-04T22:11Z Aphelion
            2019-09-23T07:50Z Equinox
            2019-12-22T04:19Z Solstice
        */
        nscanned = sscanf(line, "%d-%d-%dT%d:%dZ %10[A-Za-z]", &year, &month, &day, &hour, &minute, name);
        if (nscanned != 6)
        {
            fprintf(stderr, "SeasonsTest: %s line %d : scanned %d, expected 6\n", filename, lnum, nscanned);
            error = 1;
            goto fail;
        }

        if (year != current_year)
        {
            current_year = year;
            seasons = Astronomy_Seasons(year);
            if (seasons.status != ASTRO_SUCCESS)
            {
                fprintf(stderr, "SeasonsTest: Astronomy_Seasons(%d) returned %d\n", year, seasons.status);
                error = 1;
                goto fail;
            }
        }

        memset(&calc_time, 0xcd, sizeof(calc_time));
        correct_time = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);
        if (!strcmp(name, "Equinox"))
        {
            switch (month)
            {
            case 3:
                calc_time = seasons.mar_equinox;
                ++mar_count;
                break;
            case 9:
                calc_time = seasons.sep_equinox;
                ++sep_count;
                break;
            default:
                fprintf(stderr, "SeasonsTest: Invalid equinox date in test data: %s line %d\n", filename, lnum);
                error = 1;
                goto fail;
            }
        }
        else if (!strcmp(name, "Solstice"))
        {
            switch (month)
            {
            case 6:
                calc_time = seasons.jun_solstice;
                ++jun_count;
                break;
            case 12:
                calc_time = seasons.dec_solstice;
                ++dec_count;
                break;
            default:
                fprintf(stderr, "SeasonsTest: Invalid solstice date in test data: %s line %d\n", filename, lnum);
                error = 1;
                goto fail;
            }
        }
        else if (!strcmp(name, "Aphelion"))
        {
            /* not yet calculated */
            continue;
        }
        else if (!strcmp(name, "Perihelion"))
        {
            /* not yet calculated */
            continue;
        }
        else
        {
            fprintf(stderr, "SeasonsTest: %s line %d: unknown event type '%s'\n", filename, lnum, name);
            error = 1;
            goto fail;
        }

        /* Verify that the calculated time matches the correct time for this event. */
        diff_minutes = (24.0 * 60.0) * fabs(calc_time.tt - correct_time.tt);
        if (diff_minutes > max_minutes)
            max_minutes = diff_minutes;

        if (diff_minutes > 1.7)
        {
            fprintf(stderr, "SeasonsTest: %s line %d: excessive error (%s): %lf minutes.\n", filename, lnum, name, diff_minutes);
            error = 1;
            goto fail;
        }
    }

    printf("SeasonsTest: verified %d lines from file %s : max error minutes = %0.3lf\n", lnum, filename, max_minutes);
    printf("SeasonsTest: Event counts: mar=%d, jun=%d, sep=%d, dec=%d\n", mar_count, jun_count, sep_count, dec_count);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int MoonPhase(const char *filename)
{
    int error = 1;
    FILE *infile = NULL;
    int lnum, nscanned;
    int quarter, year, month, day, hour, minute;
    double second;
    char line[200];    

    infile = fopen(filename, "rt");
    if (infile == NULL)
    {
        fprintf(stderr, "MoonPhase: Cannot open input file '%s'\n", filename);
        error = 1;
        goto fail;
    }

    /*
        0 1800-01-25T03:21:00.000Z
        1 1800-02-01T20:40:00.000Z
        2 1800-02-09T17:26:00.000Z
        3 1800-02-16T15:49:00.000Z
    */
    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        nscanned = sscanf(line, "%d %d-%d-%dT%d:%d:%lfZ", &quarter, &year, &month, &day, &hour, &minute, &second);
        if (nscanned != 7)
        {
            fprintf(stderr, "MoonPhase(%s line %d): Invalid data format\n", filename, lnum);
            error = 1;
            goto fail;
        }

        if (quarter < 0 || quarter > 3)
        {
            fprintf(stderr, "MoonPhase(%s line %d): Invalid quarter %d\n", filename, lnum, quarter);
            error = 1;
            goto fail;
        }
    }

    printf("MoonPhase: passed %d lines for file %s\n", lnum, filename);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/
