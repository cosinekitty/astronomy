/*
    ctest.c  -  Don Cross <cosinekitty.com>

    C langauge unit test for Astronomy Engine project.
    https://github.com/cosinekitty/astronomy
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "astronomy.h"

#define PI      3.14159265358979323846
static const double DEG2RAD = 0.017453292519943296;

#define CHECK(x)                do{if(0 != (error = (x))) goto fail;}while(0)
#define FAIL(format, ...)       do{fprintf(stderr, format, __VA_ARGS__); error = 1; goto fail;}while(0)
#define FAILSTR(text)           do{fprintf(stderr, text); error = 1; goto fail;}while(0)
#define FAILRET(format, ...)    do{fprintf(stderr, format, __VA_ARGS__); return 1;}while(0)
#define FAILRETSTR(text)        do{fprintf(stderr, text); return 1;}while(0)

static int CheckStatus(int lnum, const char *varname, astro_status_t status)
{
    if (status != ASTRO_SUCCESS)
        FAILRET("FAILURE at ctest.c[%d]: %s.status = %d\n", lnum, varname, status);
    return 0;
}

static int CheckVector(int lnum, astro_vector_t v)
{
    if (v.status != ASTRO_SUCCESS)
        FAILRET("FAILURE at ctest.c[%d]: vector status = %d\n", lnum, v.status);
    return 0;
}

static int CheckEquator(int lnum, astro_equatorial_t equ)
{
    if (equ.status != ASTRO_SUCCESS)
        FAILRET("FAILURE at ctest.c[%d]: equatorial status = %d\n", lnum, equ.status);
    return 0;
}

#define CHECK_VECTOR(var,expr)   CHECK(CheckVector(__LINE__, ((var) = (expr))))
#define CHECK_EQU(var,expr)      CHECK(CheckEquator(__LINE__, ((var) = (expr))))
#define CHECK_STATUS(expr)       CHECK(CheckStatus(__LINE__, #expr, (expr).status))

static int Issue46(void);
static int Issue48(void);
static int Test_AstroTime(void);
static int AstroCheck(void);
static int Diff(const char *c_filename, const char *js_filename);
static int DiffLine(int lnum, const char *cline, const char *jline, double *maxdiff, int *worst_lnum);
static int SeasonsTest(const char *filename);
static int MoonPhase(const char *filename);
static int RiseSet(const char *filename);
static int LunarApsis(const char *filename);
static int EarthApsis(const char *filename);
static int PlanetApsis(void);
static int ElongationTest(void);
static int MagnitudeTest(void);
static int MoonTest(void);
static int RotationTest(void);
static int TestMaxMag(astro_body_t body, const char *filename);
static const char *ParseJplHorizonsDateTime(const char *text, astro_time_t *time);
static int VectorDiff(astro_vector_t a, astro_vector_t b, double *diff);
static int RefractionTest(void);
static int ConstellationTest(void);
static int LunarEclipseTest(void);

int main(int argc, const char *argv[])
{
    int error = 1;

    if (argc == 2)
    {
        const char *verb = argv[1];
        if (!strcmp(verb, "check"))
        {
            CHECK(Test_AstroTime());
            CHECK(AstroCheck());
            goto success;
        }

        if (!strcmp(verb, "elongation"))
        {
            CHECK(ElongationTest());
            goto success;
        }

        if (!strcmp(verb, "magnitude"))
        {
            CHECK(MagnitudeTest());
            goto success;
        }

        if (!strcmp(verb, "moon"))
        {
            CHECK(MoonTest());
            goto success;
        }

        if (!strcmp(verb, "rotation"))
        {
            CHECK(RotationTest());
            goto success;
        }

        if (!strcmp(verb, "refraction"))
        {
            CHECK(RefractionTest());
            goto success;
        }

        if (!strcmp(verb, "planet_apsis"))
        {
            CHECK(PlanetApsis());
            goto success;
        }

        if (!strcmp(verb, "issue46"))
        {
            CHECK(Issue46());
            goto success;
        }

        if (!strcmp(verb, "issue48"))
        {
            CHECK(Issue48());
            goto success;
        }

        if (!strcmp(verb, "constellation"))
        {
            CHECK(ConstellationTest());
            goto success;
        }

        if (!strcmp(verb, "lunar_eclipse"))
        {
            CHECK(LunarEclipseTest());
            goto success;
        }
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

        if (!strcmp(verb, "riseset"))
        {
            CHECK(RiseSet(filename));
            goto success;
        }

        if (!strcmp(verb, "moon_apsis"))
        {
            CHECK(LunarApsis(filename));
            goto success;
        }

        if (!strcmp(verb, "earth_apsis"))
        {
            CHECK(EarthApsis(filename));
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

    FAILSTR("ctest: Invalid command line arguments.\n");

success:
    error = 0;

fail:
    return error;
}

static int Test_AstroTime(void)
{
    astro_time_t time;
    astro_utc_t utc;
    const double expected_ut = 6910.270978506945;
    const double expected_tt = 6910.271779431480;
    double diff;

    const int year = 2018;
    const int month = 12;
    const int day = 2;
    const int hour = 18;
    const int minute = 30;
    const double second = 12.543;

    time = Astronomy_MakeTime(year, month, day, hour, minute, second);
    printf("C Test_AstroTime: ut=%0.6lf, tt=%0.6lf\n", time.ut, time.tt);

    diff = time.ut - expected_ut;
    if (fabs(diff) > 1.0e-12)
        FAILRET("C Test_AstroTime: excessive UT error %lg\n", diff);

    diff = time.tt - expected_tt;
    if (fabs(diff) > 1.0e-12)
        FAILRET("C Test_AstroTime: excessive TT error %lg\n", diff);

    utc = Astronomy_UtcFromTime(time);
    if (utc.year != year || utc.month != month || utc.day != day || utc.hour != hour || utc.minute != minute)
    {
        FAILRET("C Test_AstroTime: UtcFromTime FAILURE - Expected %04d-%02d-%02dT%02d:%02dZ, found %04d-%02d-%02dT%02d:%02dZ\n",
            year, month, day, hour, minute,
            utc.year, utc.month, utc.day, utc.hour, utc.minute);
    }

    diff = utc.second - second;
    if (fabs(diff) > 2.0e-5)
        FAILRET("C Test_AstroTime: excessive UTC second error %lg\n", diff);

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
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO,
        BODY_SSB, BODY_EMB
    };
    static int nbodies = sizeof(bodylist) / sizeof(bodylist[0]);

    outfile = fopen(filename, "wt");
    if (outfile == NULL)
        FAIL("C AstroCheck: Cannot open output file: %s\n", filename);

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

            if (body != BODY_EARTH && body != BODY_EMB && body != BODY_SSB)
            {
                CHECK_EQU(j2000, Astronomy_Equator(body, &time, observer, EQUATOR_J2000, NO_ABERRATION));
                CHECK_EQU(ofdate, Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION));
                hor = Astronomy_Horizon(&time, observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
                fprintf(outfile, "s %s %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf %0.16lf\n",
                    Astronomy_BodyName(body), time.tt, time.ut, j2000.ra, j2000.dec, j2000.dist, hor.azimuth, hor.altitude);
            }
        }

        CHECK_VECTOR(pos, Astronomy_GeoVector(BODY_MOON, time, NO_ABERRATION));
        fprintf(outfile, "v GM %0.16lf %0.16lf %0.16lf %0.16lf\n", pos.t.tt, pos.x, pos.y, pos.z);

        CHECK_EQU(j2000, Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_J2000, NO_ABERRATION));
        CHECK_EQU(ofdate, Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_OF_DATE, ABERRATION));
        hor = Astronomy_Horizon(&time, observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
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
        FAIL("ctest(Diff): Cannot open input file: %s\n", c_filename);

    jfile = fopen(js_filename, "rt");
    if (jfile == NULL)
        FAIL("ctest(Diff): Cannot open input file: %s\n", js_filename);

    lnum = 0;
    for(;;)
    {
        cread = fgets(cline, sizeof(cline), cfile);
        jread = fgets(jline, sizeof(jline), jfile);
        if (cread==NULL && jread==NULL)
            break;      /* normal end of both files */

        if (cread==NULL || jread==NULL)
            FAIL("ctest(Diff): Files do not have same number of lines: %s and %s\n", c_filename, js_filename);

        ++lnum;
        CHECK(DiffLine(lnum, cline, jline, &maxdiff, &worst_lnum));
    }

    printf("ctest(Diff): Maximum numeric difference = %lg, worst line number = %d\n", maxdiff, worst_lnum);
    if (maxdiff > 1.819e-12)
        FAIL("C ERROR: Excessive error comparing files %s and %s\n", c_filename, js_filename);

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
        FAIL("ctest(DiffLine): Line %d mismatch record type: '%c' vs '%c'.\n", lnum, cline[0], jline[0]);

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
        FAIL("ctest(DiffLine): Line %d type '%c' is not a valid record type.\n", lnum, cline[0]);
    }

    if (nc != nj)
        FAIL("ctest(DiffLine): Line %d mismatch data counts: %d vs %d\n", lnum, nc, nj);

    if (nc != nrequired)
        FAIL("ctest(DiffLine): Line %d incorrect number of scanned arguments: %d\n", lnum, nc);

    if (strcmp(cbody, jbody))
        FAIL("ctest(DiffLine): Line %d body mismatch: '%s' vs '%s'\n.", lnum, cbody, jbody);

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
        FAIL("C SeasonsTest: Cannot open input file: %s\n", filename);

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
            FAIL("C SeasonsTest: %s line %d : scanned %d, expected 6\n", filename, lnum, nscanned);

        if (year != current_year)
        {
            current_year = year;
            seasons = Astronomy_Seasons(year);
            if (seasons.status != ASTRO_SUCCESS)
                FAIL("C SeasonsTest: Astronomy_Seasons(%d) returned %d\n", year, seasons.status);
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
                FAIL("C SeasonsTest: Invalid equinox date in test data: %s line %d\n", filename, lnum);
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
                FAIL("C SeasonsTest: Invalid solstice date in test data: %s line %d\n", filename, lnum);
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
            FAIL("C SeasonsTest: %s line %d: unknown event type '%s'\n", filename, lnum, name);

        /* Verify that the calculated time matches the correct time for this event. */
        diff_minutes = (24.0 * 60.0) * fabs(calc_time.tt - correct_time.tt);
        if (diff_minutes > max_minutes)
            max_minutes = diff_minutes;

        if (diff_minutes > 1.7)
            FAIL("C SeasonsTest: %s line %d: excessive error (%s): %lf minutes.\n", filename, lnum, name, diff_minutes);
    }

    printf("C SeasonsTest: verified %d lines from file %s : max error minutes = %0.3lf\n", lnum, filename, max_minutes);
    printf("C SeasonsTest: Event counts: mar=%d, jun=%d, sep=%d, dec=%d\n", mar_count, jun_count, sep_count, dec_count);
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
    int expected_quarter, quarter_count = 0;
    int prev_year = 0;
    double second, expected_elong;
    astro_time_t expected_time, start_time;
    astro_angle_result_t result;
    double degree_error, arcmin, max_arcmin = 0.0;
    double diff_seconds, maxdiff = 0.0;
    const double threshold_seconds = 120.0; /* max tolerable prediction error in seconds */
    astro_moon_quarter_t mq;
    char line[200];

    memset(&mq, 0xcd, sizeof(mq));

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C MoonPhase: Cannot open input file '%s'\n", filename);

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
            FAIL("C MoonPhase(%s line %d): Invalid data format\n", filename, lnum);

        if (quarter < 0 || quarter > 3)
            FAIL("C MoonPhase(%s line %d): Invalid quarter %d\n", filename, lnum, quarter);

        expected_elong = 90.0 * quarter;
        expected_time = Astronomy_MakeTime(year, month, day, hour, minute, second);
        result = Astronomy_MoonPhase(expected_time);
        degree_error = fabs(result.angle - expected_elong);
        if (degree_error > 180.0)
            degree_error = 360 - degree_error;
        arcmin = 60.0 * degree_error;

        if (arcmin > 1.0)
            FAIL("C MoonPhase(%s line %d): EXCESSIVE ANGULAR ERROR: %lg arcmin\n", filename, lnum, arcmin);

        if (arcmin > max_arcmin)
            max_arcmin = arcmin;

        if (year != prev_year)
        {
            prev_year = year;
            /* The test data contains a single year's worth of data for every 10 years. */
            /* Every time we see the year value change, it breaks continuity of the phases. */
            /* Start the search over again. */
            start_time = Astronomy_MakeTime(year, 1, 1, 0, 0, 0.0);
            mq = Astronomy_SearchMoonQuarter(start_time);
            expected_quarter = -1;  /* we have no idea what the quarter should be */
        }
        else
        {
            /* Yet another lunar quarter in the same year. */
            expected_quarter = (1 + mq.quarter) % 4;
            mq = Astronomy_NextMoonQuarter(mq);

            /* Make sure we find the next expected quarter. */
            if (expected_quarter != mq.quarter)
                FAIL("C MoonPhase(%s line %d): Astronomy_SearchMoonQuarter returned quarter %d, but expected %d\n", filename, lnum, mq.quarter, expected_quarter);
        }

        if (mq.status != ASTRO_SUCCESS)
            FAIL("C MoonPhase(%s line %d): Astronomy_SearchMoonQuarter returned %d\n", filename, lnum, mq.status);

        ++quarter_count;

        /* Make sure the time matches what we expect. */
        diff_seconds = fabs(mq.time.tt - expected_time.tt) * (24.0 * 3600.0);
        if (diff_seconds > threshold_seconds)
            FAIL("C MoonPhase(%s line %d): excessive time error %0.3lf seconds\n", filename, lnum, diff_seconds);

        if (diff_seconds > maxdiff)
            maxdiff = diff_seconds;
    }

    printf("C MoonPhase: passed %d lines for file %s : max_arcmin = %0.6lf, maxdiff = %0.3lf seconds, %d quarters\n", lnum, filename, max_arcmin, maxdiff, quarter_count);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int TestElongFile(const char *filename, double targetRelLon)
{
    int error = 1;
    FILE *infile = NULL;
    int lnum;
    char line[100];
    char name[20];
    int year, month, day, hour, minute;
    int nscanned;
    astro_time_t search_date, expected_time;
    astro_body_t body;
    astro_search_result_t search_result;
    double diff_minutes;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C TestElongFile: Cannot open input file: %s\n", filename);

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        /* 2018-05-09T00:28Z Jupiter */
        nscanned = sscanf(line, "%d-%d-%dT%d:%dZ %9[A-Za-z]", &year, &month, &day, &hour, &minute, name);
        if (nscanned != 6)
            FAIL("C TestElongFile(%s line %d): Invalid data format.\n", filename, lnum);

        body = Astronomy_BodyCode(name);
        if (body == BODY_INVALID)
            FAIL("C TestElongFile(%s line %d): Invalid body name '%s'\n", filename, lnum, name);

        search_date = Astronomy_MakeTime(year, 1, 1, 0, 0, 0.0);
        expected_time = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);
        search_result = Astronomy_SearchRelativeLongitude(body, targetRelLon, search_date);
        if (search_result.status != ASTRO_SUCCESS)
            FAIL("C TestElongFile(%s line %d): SearchRelativeLongitude returned %d\n", filename, lnum, search_result.status);

        diff_minutes = (24.0 * 60.0) * (search_result.time.tt - expected_time.tt);
        printf("C TestElongFile: %-7s error = %6.3lf minutes\n", name, diff_minutes);
        if (fabs(diff_minutes) > 15.0)
            FAIL("C TestElongFile(%s line %d): EXCESSIVE ERROR\n", filename, lnum);
    }

    printf("C TestElongFile: passed %d rows of data\n", lnum);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

typedef struct
{
    astro_body_t        body;
    const char         *searchDate;
    const char         *eventDate;
    double              angle;
    astro_visibility_t  visibility;
}
elong_test_t;

static const elong_test_t ElongTestData[] =
{
    /* Max elongation data obtained from: */
    /* http://www.skycaramba.com/greatest_elongations.shtml */
    { BODY_MERCURY, "2010-01-17T05:22Z", "2010-01-27T05:22Z", 24.80, VISIBLE_MORNING },
    { BODY_MERCURY, "2010-05-16T02:15Z", "2010-05-26T02:15Z", 25.10, VISIBLE_MORNING },
    { BODY_MERCURY, "2010-09-09T17:24Z", "2010-09-19T17:24Z", 17.90, VISIBLE_MORNING },
    { BODY_MERCURY, "2010-12-30T14:33Z", "2011-01-09T14:33Z", 23.30, VISIBLE_MORNING },
    { BODY_MERCURY, "2011-04-27T19:03Z", "2011-05-07T19:03Z", 26.60, VISIBLE_MORNING },
    { BODY_MERCURY, "2011-08-24T05:52Z", "2011-09-03T05:52Z", 18.10, VISIBLE_MORNING },
    { BODY_MERCURY, "2011-12-13T02:56Z", "2011-12-23T02:56Z", 21.80, VISIBLE_MORNING },
    { BODY_MERCURY, "2012-04-08T17:22Z", "2012-04-18T17:22Z", 27.50, VISIBLE_MORNING },
    { BODY_MERCURY, "2012-08-06T12:04Z", "2012-08-16T12:04Z", 18.70, VISIBLE_MORNING },
    { BODY_MERCURY, "2012-11-24T22:55Z", "2012-12-04T22:55Z", 20.60, VISIBLE_MORNING },
    { BODY_MERCURY, "2013-03-21T22:02Z", "2013-03-31T22:02Z", 27.80, VISIBLE_MORNING },
    { BODY_MERCURY, "2013-07-20T08:51Z", "2013-07-30T08:51Z", 19.60, VISIBLE_MORNING },
    { BODY_MERCURY, "2013-11-08T02:28Z", "2013-11-18T02:28Z", 19.50, VISIBLE_MORNING },
    { BODY_MERCURY, "2014-03-04T06:38Z", "2014-03-14T06:38Z", 27.60, VISIBLE_MORNING },
    { BODY_MERCURY, "2014-07-02T18:22Z", "2014-07-12T18:22Z", 20.90, VISIBLE_MORNING },
    { BODY_MERCURY, "2014-10-22T12:36Z", "2014-11-01T12:36Z", 18.70, VISIBLE_MORNING },
    { BODY_MERCURY, "2015-02-14T16:20Z", "2015-02-24T16:20Z", 26.70, VISIBLE_MORNING },
    { BODY_MERCURY, "2015-06-14T17:10Z", "2015-06-24T17:10Z", 22.50, VISIBLE_MORNING },
    { BODY_MERCURY, "2015-10-06T03:20Z", "2015-10-16T03:20Z", 18.10, VISIBLE_MORNING },
    { BODY_MERCURY, "2016-01-28T01:22Z", "2016-02-07T01:22Z", 25.60, VISIBLE_MORNING },
    { BODY_MERCURY, "2016-05-26T08:45Z", "2016-06-05T08:45Z", 24.20, VISIBLE_MORNING },
    { BODY_MERCURY, "2016-09-18T19:27Z", "2016-09-28T19:27Z", 17.90, VISIBLE_MORNING },
    { BODY_MERCURY, "2017-01-09T09:42Z", "2017-01-19T09:42Z", 24.10, VISIBLE_MORNING },
    { BODY_MERCURY, "2017-05-07T23:19Z", "2017-05-17T23:19Z", 25.80, VISIBLE_MORNING },
    { BODY_MERCURY, "2017-09-02T10:14Z", "2017-09-12T10:14Z", 17.90, VISIBLE_MORNING },
    { BODY_MERCURY, "2017-12-22T19:48Z", "2018-01-01T19:48Z", 22.70, VISIBLE_MORNING },
    { BODY_MERCURY, "2018-04-19T18:17Z", "2018-04-29T18:17Z", 27.00, VISIBLE_MORNING },
    { BODY_MERCURY, "2018-08-16T20:35Z", "2018-08-26T20:35Z", 18.30, VISIBLE_MORNING },
    { BODY_MERCURY, "2018-12-05T11:34Z", "2018-12-15T11:34Z", 21.30, VISIBLE_MORNING },
    { BODY_MERCURY, "2019-04-01T19:40Z", "2019-04-11T19:40Z", 27.70, VISIBLE_MORNING },
    { BODY_MERCURY, "2019-07-30T23:08Z", "2019-08-09T23:08Z", 19.00, VISIBLE_MORNING },
    { BODY_MERCURY, "2019-11-18T10:31Z", "2019-11-28T10:31Z", 20.10, VISIBLE_MORNING },
    { BODY_MERCURY, "2010-03-29T23:32Z", "2010-04-08T23:32Z", 19.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2010-07-28T01:03Z", "2010-08-07T01:03Z", 27.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2010-11-21T15:42Z", "2010-12-01T15:42Z", 21.50, VISIBLE_EVENING },
    { BODY_MERCURY, "2011-03-13T01:07Z", "2011-03-23T01:07Z", 18.60, VISIBLE_EVENING },
    { BODY_MERCURY, "2011-07-10T04:56Z", "2011-07-20T04:56Z", 26.80, VISIBLE_EVENING },
    { BODY_MERCURY, "2011-11-04T08:40Z", "2011-11-14T08:40Z", 22.70, VISIBLE_EVENING },
    { BODY_MERCURY, "2012-02-24T09:39Z", "2012-03-05T09:39Z", 18.20, VISIBLE_EVENING },
    { BODY_MERCURY, "2012-06-21T02:00Z", "2012-07-01T02:00Z", 25.70, VISIBLE_EVENING },
    { BODY_MERCURY, "2012-10-16T21:59Z", "2012-10-26T21:59Z", 24.10, VISIBLE_EVENING },
    { BODY_MERCURY, "2013-02-06T21:24Z", "2013-02-16T21:24Z", 18.10, VISIBLE_EVENING },
    { BODY_MERCURY, "2013-06-02T16:45Z", "2013-06-12T16:45Z", 24.30, VISIBLE_EVENING },
    { BODY_MERCURY, "2013-09-29T09:59Z", "2013-10-09T09:59Z", 25.30, VISIBLE_EVENING },
    { BODY_MERCURY, "2014-01-21T10:00Z", "2014-01-31T10:00Z", 18.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2014-05-15T07:06Z", "2014-05-25T07:06Z", 22.70, VISIBLE_EVENING },
    { BODY_MERCURY, "2014-09-11T22:20Z", "2014-09-21T22:20Z", 26.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2015-01-04T20:26Z", "2015-01-14T20:26Z", 18.90, VISIBLE_EVENING },
    { BODY_MERCURY, "2015-04-27T04:46Z", "2015-05-07T04:46Z", 21.20, VISIBLE_EVENING },
    { BODY_MERCURY, "2015-08-25T10:20Z", "2015-09-04T10:20Z", 27.10, VISIBLE_EVENING },
    { BODY_MERCURY, "2015-12-19T03:11Z", "2015-12-29T03:11Z", 19.70, VISIBLE_EVENING },
    { BODY_MERCURY, "2016-04-08T14:00Z", "2016-04-18T14:00Z", 19.90, VISIBLE_EVENING },
    { BODY_MERCURY, "2016-08-06T21:24Z", "2016-08-16T21:24Z", 27.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2016-12-01T04:36Z", "2016-12-11T04:36Z", 20.80, VISIBLE_EVENING },
    { BODY_MERCURY, "2017-03-22T10:24Z", "2017-04-01T10:24Z", 19.00, VISIBLE_EVENING },
    { BODY_MERCURY, "2017-07-20T04:34Z", "2017-07-30T04:34Z", 27.20, VISIBLE_EVENING },
    { BODY_MERCURY, "2017-11-14T00:32Z", "2017-11-24T00:32Z", 22.00, VISIBLE_EVENING },
    { BODY_MERCURY, "2018-03-05T15:07Z", "2018-03-15T15:07Z", 18.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2018-07-02T05:24Z", "2018-07-12T05:24Z", 26.40, VISIBLE_EVENING },
    { BODY_MERCURY, "2018-10-27T15:25Z", "2018-11-06T15:25Z", 23.30, VISIBLE_EVENING },
    { BODY_MERCURY, "2019-02-17T01:23Z", "2019-02-27T01:23Z", 18.10, VISIBLE_EVENING },
    { BODY_MERCURY, "2019-06-13T23:14Z", "2019-06-23T23:14Z", 25.20, VISIBLE_EVENING },
    { BODY_MERCURY, "2019-10-10T04:00Z", "2019-10-20T04:00Z", 24.60, VISIBLE_EVENING },
    { BODY_VENUS,   "2010-12-29T15:57Z", "2011-01-08T15:57Z", 47.00, VISIBLE_MORNING },
    { BODY_VENUS,   "2012-08-05T08:59Z", "2012-08-15T08:59Z", 45.80, VISIBLE_MORNING },
    { BODY_VENUS,   "2014-03-12T19:25Z", "2014-03-22T19:25Z", 46.60, VISIBLE_MORNING },
    { BODY_VENUS,   "2015-10-16T06:57Z", "2015-10-26T06:57Z", 46.40, VISIBLE_MORNING },
    { BODY_VENUS,   "2017-05-24T13:09Z", "2017-06-03T13:09Z", 45.90, VISIBLE_MORNING },
    { BODY_VENUS,   "2018-12-27T04:24Z", "2019-01-06T04:24Z", 47.00, VISIBLE_MORNING },
    { BODY_VENUS,   "2010-08-10T03:19Z", "2010-08-20T03:19Z", 46.00, VISIBLE_EVENING },
    { BODY_VENUS,   "2012-03-17T08:03Z", "2012-03-27T08:03Z", 46.00, VISIBLE_EVENING },
    { BODY_VENUS,   "2013-10-22T08:00Z", "2013-11-01T08:00Z", 47.10, VISIBLE_EVENING },
    { BODY_VENUS,   "2015-05-27T18:46Z", "2015-06-06T18:46Z", 45.40, VISIBLE_EVENING },
    { BODY_VENUS,   "2017-01-02T13:19Z", "2017-01-12T13:19Z", 47.10, VISIBLE_EVENING },
    { BODY_VENUS,   "2018-08-07T17:02Z", "2018-08-17T17:02Z", 45.90, VISIBLE_EVENING }
};

static const int ElongTestCount = sizeof(ElongTestData) / sizeof(ElongTestData[0]);

static int ParseDate(const char *text, astro_time_t *time)
{
    int year, month, day, hour, minute, nscanned;
    double second = 0.0;

    nscanned = sscanf(text, "%d-%d-%dT%d:%dZ", &year, &month, &day, &hour, &minute);
    if (nscanned != 5)
    {
        nscanned = sscanf(text, "%d-%d-%dT%d:%d:%lfZ", &year, &month, &day, &hour, &minute, &second);
        if (nscanned != 6)
        {
            fprintf(stderr, "C ParseDate: Invalid date text '%s'\n", text);
            time->ut = time->tt = NAN;
            return 1;
        }
    }

    *time = Astronomy_MakeTime(year, month, day, hour, minute, second);
    return 0;
}

static int TestMaxElong(const elong_test_t *test)
{
    int error;
    astro_time_t searchTime, eventTime;
    astro_elongation_t evt;
    double hour_diff, arcmin_diff;
    const char *name = NULL;
    const char *vis = NULL;

    switch (test->body)
    {
    case BODY_MERCURY:  name = "Mercury";   break;
    case BODY_VENUS:    name = "Venus";     break;
    default:
        FAIL("C TestMaxElong: invalid body %d in test data.\n", test->body);
    }

    switch (test->visibility)
    {
    case VISIBLE_MORNING:   vis = "morning";    break;
    case VISIBLE_EVENING:   vis = "evening";    break;
    default:
        FAIL("C TestMaxElong: invalid visibility %d in test data.\n", test->visibility);
    }

    CHECK(ParseDate(test->searchDate, &searchTime));
    CHECK(ParseDate(test->eventDate,  &eventTime));

    evt = Astronomy_SearchMaxElongation(test->body, searchTime);
    if (evt.status != ASTRO_SUCCESS)
        FAIL("C TestMaxElong(%s %s): SearchMaxElongation returned %d\n", name, test->searchDate, evt.status);

    hour_diff = 24.0 * fabs(evt.time.tt - eventTime.tt);
    arcmin_diff = 60.0 * fabs(evt.elongation - test->angle);

    printf("C TestMaxElong: %-7s %-7s elong=%5.2lf (%4.2lf arcmin, %5.3lf hours)\n", name, vis, evt.elongation, arcmin_diff, hour_diff);

    if (hour_diff > 0.603)
        FAIL("C TestMaxElong(%s %s): excessive hour error.\n", name, test->searchDate);

    if (arcmin_diff > 3.4)
        FAIL("C TestMaxElong(%s %s): excessive arcmin error.\n", name, test->searchDate);

fail:
    return error;
}

static int SearchElongTest()
{
    int error = 1;
    int i;

    for (i=0; i < ElongTestCount; ++i)
        CHECK(TestMaxElong(&ElongTestData[i]));

    printf("C SearchElongTest: Passed %d rows\n", ElongTestCount);
    error = 0;

fail:
    return error;
}

static int TestPlanetLongitudes(
    astro_body_t body,
    const char *outFileName,
    const char *zeroLonEventName)
{
    int error = 1;
    const int startYear = 1700;
    const int stopYear  = 2200;
    astro_time_t time, stopTime;
    double rlon = 0.0;
    const char *event;
    astro_search_result_t search_result;
    int count = 0;
    double day_diff, min_diff = 1.0e+99, max_diff = 1.0e+99, sum_diff = 0.0;
    astro_vector_t geo;
    double dist;
    FILE *outfile = NULL;
    const char *name;
    double ratio, thresh;

    name = Astronomy_BodyName(body);
    if (!name[0])
        FAIL("C TestPlanetLongitudes: Invalid body code %d\n", body);

    outfile = fopen(outFileName, "wt");
    if (outfile == NULL)
        FAIL("C TestPlanetLongitudes: Cannot open output file: %s\n", outFileName);

    time = Astronomy_MakeTime(startYear, 1, 1, 0, 0, 0.0);
    stopTime = Astronomy_MakeTime(stopYear, 1, 1, 0, 0, 0.0);
    while (time.tt < stopTime.tt)
    {
        ++count;
        event = (rlon == 0.0) ? zeroLonEventName : "sup";
        search_result = Astronomy_SearchRelativeLongitude(body, rlon, time);
        if (search_result.status != ASTRO_SUCCESS)
            FAIL("C TestPlanetLongitudes(%s): SearchRelativeLongitude returned %d\n", name, search_result.status);

        if (count >= 2)
        {
            /* Check for consistent intervals. */
            /* Mainly I don't want to skip over an event! */
            day_diff = search_result.time.tt - time.tt;
            sum_diff += day_diff;
            if (count == 2)
            {
                min_diff = max_diff = day_diff;
            }
            else
            {
                if (day_diff < min_diff)
                    min_diff = day_diff;

                if (day_diff > max_diff)
                    max_diff = day_diff;
            }
        }

        geo = Astronomy_GeoVector(body, search_result.time, ABERRATION);
        if (geo.status != ASTRO_SUCCESS)
            FAIL("C TestPlanetLongitudes(%s): GeoVector returned %d\n", name, geo.status);

        dist = Astronomy_VectorLength(geo);
        fprintf(outfile, "e %s %s %0.16lf %0.16lf\n", name, event, search_result.time.tt, dist);

        /* Search for the opposite longitude event next time. */
        time = search_result.time;
        rlon = 180.0 - rlon;
    }

    switch (body)
    {
    case BODY_MERCURY:  thresh = 1.65;  break;
    case BODY_MARS:     thresh = 1.30;  break;
    default:            thresh = 1.07;  break;
    }

    ratio = max_diff / min_diff;
    printf("C TestPlanetLongitudes(%-7s): %5d events, ratio=%5.3lf, file: %s\n", name, count, ratio, outFileName);

    if (ratio > thresh)
        FAIL("C TestPlanetLongitudes(%s): excessive event interval ratio.\n", name);

    error = 0;
fail:
    if (outfile != NULL) fclose(outfile);
    return error;
}

static int ElongationTest(void)
{
    int error;

    CHECK(TestElongFile("longitude/opposition_2018.txt", 0.0));

    CHECK(TestPlanetLongitudes(BODY_MERCURY, "temp/c_longitude_Mercury.txt", "inf"));
    CHECK(TestPlanetLongitudes(BODY_VENUS,   "temp/c_longitude_Venus.txt",   "inf"));
    CHECK(TestPlanetLongitudes(BODY_MARS,    "temp/c_longitude_Mars.txt",    "opp"));
    CHECK(TestPlanetLongitudes(BODY_JUPITER, "temp/c_longitude_Jupiter.txt", "opp"));
    CHECK(TestPlanetLongitudes(BODY_SATURN,  "temp/c_longitude_Saturn.txt",  "opp"));
    CHECK(TestPlanetLongitudes(BODY_URANUS,  "temp/c_longitude_Uranus.txt",  "opp"));
    CHECK(TestPlanetLongitudes(BODY_NEPTUNE, "temp/c_longitude_Neptune.txt", "opp"));
    CHECK(TestPlanetLongitudes(BODY_PLUTO,   "temp/c_longitude_Pluto.txt",   "opp"));

    CHECK(SearchElongTest());

fail:
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int RiseSet(const char *filename)
{
    int error = 1;
    FILE *infile = NULL;
    char line[100];
    char name[20];
    double longitude, latitude;
    int year, month, day, hour, minute;
    char kind[2];       /* "r" or "s" for rise or set */
    int direction;      /* +1 for rise, -1 for set */
    int lnum, nscanned;
    astro_time_t correct_date;
    astro_body_t body, current_body;
    astro_observer_t observer;
    astro_time_t r_search_date;
    astro_search_result_t r_evt, s_evt;     /* rise event, set event: search results */
    astro_search_result_t a_evt, b_evt;     /* chronologically first and second events */
    astro_time_t s_search_date;
    int a_dir = 0, b_dir = 0;
    double error_minutes, rms_minutes;
    double sum_minutes = 0.0;
    double max_minutes = 0.0;
    const double nudge_days = 0.01;

    observer.latitude = observer.longitude = observer.height = NAN;
    current_body = BODY_INVALID;
    a_evt.status = b_evt.status = r_evt.status = s_evt.status = ASTRO_NOT_INITIALIZED;
    b_evt.time.tt = b_evt.time.ut = a_evt.time.tt = a_evt.time.ut = NAN;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C RiseSet: cannot open input file: %s\n", filename);

    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        /* Moon  103 -61 1944-01-02T17:08Z s */
        /* Moon  103 -61 1944-01-03T05:47Z r */
        nscanned = sscanf(line, "%9[A-Za-z] %lf %lf %d-%d-%dT%d:%dZ %1[rs]",
            name, &longitude, &latitude, &year, &month, &day, &hour, &minute, kind);

        if (nscanned != 9)
            FAIL("C RiseSet(%s line %d): invalid format\n", filename, lnum);

        correct_date = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);

        if (!strcmp(kind, "r"))
            direction = +1;
        else if (!strcmp(kind, "s"))
            direction = -1;
        else
            FAIL("C RiseSet(%s line %d): invalid kind '%s'\n", filename, lnum, kind);

        body = Astronomy_BodyCode(name);
        if (body == BODY_INVALID)
            FAIL("C RiseSet(%s line %d): invalid body name '%s'", filename, lnum, name);

        /* Every time we see a new geographic location, start a new iteration */
        /* of finding all rise/set times for that UTC calendar year. */
        if (observer.latitude != latitude || observer.longitude != longitude || current_body != body)
        {
            current_body = body;
            observer = Astronomy_MakeObserver(latitude, longitude, 0.0);
            r_search_date = s_search_date = Astronomy_MakeTime(year, 1, 1, 0, 0, 0.0);
            b_evt.time.tt = b_evt.time.ut = NAN;
            b_evt.status = ASTRO_NOT_INITIALIZED;
            printf("C RiseSet: %-7s lat=%0.1lf lon=%0.1lf\n", name, latitude, longitude);
        }

        if (b_evt.status == ASTRO_SUCCESS)      /* has b_evt been initialized? (does it contain a valid event?) */
        {
            /* Recycle the second event from the previous iteration as the first event. */
            a_evt = b_evt;
            a_dir = b_dir;
            b_evt.status = ASTRO_NOT_INITIALIZED;   /* invalidate b_evt for the next time around the loop */
        }
        else
        {
            r_evt = Astronomy_SearchRiseSet(body, observer, DIRECTION_RISE, r_search_date, 366.0);
            if (r_evt.status != ASTRO_SUCCESS)
                FAIL("C RiseSet(%s line %d): did not find %s rise event.\n", filename, lnum, name);

            s_evt = Astronomy_SearchRiseSet(body, observer, DIRECTION_SET, s_search_date, 366.0);
            if (s_evt.status != ASTRO_SUCCESS)
                FAIL("C RiseSet(%s line %d): did not find %s set event.\n", filename, lnum, name);

            /* Expect the current event to match the earlier of the found dates. */
            if (r_evt.time.tt < s_evt.time.tt)
            {
                a_evt = r_evt;
                b_evt = s_evt;
                a_dir = +1;
                b_dir = -1;
            }
            else
            {
                a_evt = s_evt;
                b_evt = r_evt;
                a_dir = -1;
                b_dir = +1;
            }

            /* Nudge the event times forward a tiny amount. */
            r_search_date = Astronomy_AddDays(r_evt.time, nudge_days);
            s_search_date = Astronomy_AddDays(s_evt.time, nudge_days);
        }

        if (a_dir != direction)
            FAIL("C RiseSet(%s line %d): expected dir=%d but found %d\n", filename, lnum, a_dir, direction);

        error_minutes = (24.0 * 60.0) * fabs(a_evt.time.tt - correct_date.tt);
        sum_minutes += error_minutes * error_minutes;
        if (error_minutes > max_minutes)
            max_minutes = error_minutes;

        if (error_minutes > 0.57)
            FAIL("C RiseSet(%s line %d): excessive prediction time error = %lg minutes.\n", filename, lnum, error_minutes);
    }

    rms_minutes = sqrt(sum_minutes / lnum);
    printf("C RiseSet: passed %d lines: time errors in minutes: rms=%0.4lf, max=%0.4lf\n", lnum, rms_minutes, max_minutes);
    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int CheckMagnitudeData(astro_body_t body, const char *filename)
{
    int error = 1;
    int lnum, count;
    FILE *infile = NULL;
    char line[200];
    const char *rest;
    astro_time_t time;
    astro_illum_t illum;
    int nscanned;
    double mag, sbrt, dist, rdot, delta, deldot, phase_angle;
    double diff, diff_lo = NAN, diff_hi = NAN, sum_squared_diff = 0.0, rms;
    const double limit = 0.012;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C CheckMagnitudeData: cannot open input file: %s\n", filename);

    count = lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        rest = ParseJplHorizonsDateTime(line, &time);
        /* Ignore non-data rows and data rows that contain "n.a." */
        if (rest && !strstr(rest, "n.a."))
        {
            nscanned = sscanf(rest, "%lf %lf %lf %lf %lf %lf %lf",
                &mag, &sbrt, &dist, &rdot, &delta, &deldot, &phase_angle);

            if (nscanned != 7)
                FAIL("C CheckMagnitudeData(%s line %d): invalid data format\n", filename, lnum);

            illum = Astronomy_Illumination(body, time);
            if (illum.status != ASTRO_SUCCESS)
                FAIL("C CheckMagnitudeData(%s line %d): Astronomy_Illumination returned %d\n", filename, lnum, illum.status);

            diff = illum.mag - mag;
            if (fabs(diff) > limit)
                FAIL("C CheckMagnitudeData(%s line %d): EXCESSIVE ERROR: correct mag=%lf, calc mag=%lf, diff=%lf\n", filename, lnum, mag, illum.mag, diff);

            sum_squared_diff += diff * diff;
            if (count == 0)
            {
                diff_lo = diff_hi = diff;
            }
            else
            {
                if (diff < diff_lo)
                    diff_lo = diff;

                if (diff > diff_hi)
                    diff_hi = diff;
            }

            ++count;
        }
    }

    if (count == 0)
        FAIL("C CheckMagnitudeData: Did not find any data in file: %s\n", filename);

    rms = sqrt(sum_squared_diff / count);
    printf("C CheckMagnitudeData: %-21s %5d rows diff_lo=%0.4lf diff_hi=%0.4lf rms=%0.4lf\n", filename, count, diff_lo, diff_hi, rms);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

static int CheckSaturn()
{
    /* JPL Horizons does not include Saturn's rings in its magnitude models. */
    /* I still don't have authoritative test data for Saturn's magnitude. */
    /* For now, I just test for consistency with Paul Schlyter's formulas at: */
    /* http://www.stjarnhimlen.se/comp/ppcomp.html#15 */

    static const struct test_case
    {
        const char *date;
        double mag;
        double tilt;
    }
    data[] =
    {
        { "1972-01-01T00:00Z", -0.31904865,  +24.50061220 },
        { "1980-01-01T00:00Z", +0.85213663,   -1.85761461 },
        { "2009-09-04T00:00Z", +1.01626809,   +0.08380716 },
        { "2017-06-15T00:00Z", -0.12318790,  -26.60871409 },
        { "2019-05-01T00:00Z", +0.32954097,  -23.53880802 },
        { "2025-09-25T00:00Z", +0.51286575,   +1.52327932 },
        { "2032-05-15T00:00Z", -0.04652109,  +26.95717765 }
    };

    static const int ncases = sizeof(data) / sizeof(data[0]);

    int i;
    astro_illum_t illum;
    astro_time_t time;
    double mag_diff, tilt_diff;

    for (i=0; i < ncases; ++i)
    {
        int error = ParseDate(data[i].date, &time);
        if (error)
            return 1;

        illum = Astronomy_Illumination(BODY_SATURN, time);
        if (illum.status != ASTRO_SUCCESS)
            FAILRET("C CheckSaturn(%d): Illumination returned %d\n", i, illum.status);

        printf("Saturn: date=%s  calc mag=%12.8lf  ring_tilt=%12.8lf\n", data[i].date, illum.mag, illum.ring_tilt);

        mag_diff = fabs(illum.mag - data[i].mag);
        if (mag_diff > 1.0e-4)
            FAILRET("C ERROR: Excessive magnitude error %lg\n", mag_diff);

        tilt_diff = fabs(illum.ring_tilt - data[i].tilt);
        if (tilt_diff > 3.0e-5)
            FAILRET("C ERROR: Excessive ring tilt error %lg\n", tilt_diff);
    }

    return 0;
}

static int MoonTest(void)
{
    astro_time_t time = Astronomy_MakeTime(2019, 6, 24, 15, 45, 37.0);
    astro_vector_t vec = Astronomy_GeoMoon(time);
    double dx, dy, dz, diff;

    if (vec.status != ASTRO_SUCCESS)
        FAILRET("C MoonTest: ERROR: vec.status = %d\n", vec.status);

    printf("C MoonTest: %0.16lg %0.16lg %0.16lg\n", vec.x, vec.y, vec.z);


    dx = vec.x - (+0.002674036155459549);
    dy = vec.y - (-0.0001531716308218381);
    dz = vec.z - (-0.0003150201604895409);
    diff = sqrt(dx*dx + dy*dy + dz*dz);
    printf("C MoonTest: diff = %lg\n", diff);
    if (diff > 4.34e-19)
    {
        fprintf(stderr, "C MoonTest: EXCESSIVE ERROR\n");
        return 1;
    }

    return 0;
}

static int MagnitudeTest(void)
{
    int nfailed = 0;

    nfailed += CheckMagnitudeData(BODY_SUN,     "magnitude/Sun.txt");
    nfailed += CheckMagnitudeData(BODY_MOON,    "magnitude/Moon.txt");
    nfailed += CheckMagnitudeData(BODY_MERCURY, "magnitude/Mercury.txt");
    nfailed += CheckMagnitudeData(BODY_VENUS,   "magnitude/Venus.txt");
    nfailed += CheckMagnitudeData(BODY_MARS,    "magnitude/Mars.txt");
    nfailed += CheckMagnitudeData(BODY_JUPITER, "magnitude/Jupiter.txt");
    nfailed += CheckSaturn();
    nfailed += CheckMagnitudeData(BODY_URANUS,  "magnitude/Uranus.txt");
    nfailed += CheckMagnitudeData(BODY_NEPTUNE, "magnitude/Neptune.txt");
    nfailed += CheckMagnitudeData(BODY_PLUTO,   "magnitude/Pluto.txt");

    nfailed += TestMaxMag(BODY_VENUS, "magnitude/maxmag_Venus.txt");

    if (nfailed > 0)
        fprintf(stderr, "C MagnitudeTest: FAILED %d test(s).\n", nfailed);

    return nfailed;
}

static int TestMaxMag(astro_body_t body, const char *filename)
{
    int error = 1;
    FILE *infile = NULL;
    int lnum, nscanned;
    int year1, month1, day1, hour1, minute1;
    int year2, month2, day2, hour2, minute2;
    double correct_mag, correct_angle1, correct_angle2;
    char line[100];
    astro_time_t search_time, time1, time2, center_time;
    astro_illum_t illum;
    double mag_diff, hours_diff;

    /*
        Example of input data:

        2001-02-21T08:00Z 2001-02-27T08:00Z 23.17 19.53 -4.84

        JPL Horizons test data has limited floating point precision in the magnitude values.
        There is a pair of dates for the beginning and end of the max magnitude period,
        given the limited precision. We pick the point halfway between as the supposed max magnitude time.
    */

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C TestMaxMag: Cannot open input file: %s\n", filename);

    lnum = 0;
    search_time = Astronomy_MakeTime(2001, 1, 1, 0, 0, 0.0);
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        nscanned = sscanf(line, "%d-%d-%dT%d:%dZ %d-%d-%dT%d:%dZ %lf %lf %lf\n",
            &year1, &month1, &day1, &hour1, &minute1,
            &year2, &month2, &day2, &hour2, &minute2,
            &correct_angle1, &correct_angle2, &correct_mag);

        if (nscanned != 13)
            FAIL("C TestMaxMag(%s line %d): invalid data format.\n", filename, lnum);

        time1 = Astronomy_MakeTime(year1, month1, day1, hour1, minute1, 0.0);
        time2 = Astronomy_MakeTime(year2, month2, day2, hour2, minute2, 0.0);
        center_time = Astronomy_AddDays(time1, 0.5*(time2.ut - time1.ut));

        illum = Astronomy_SearchPeakMagnitude(body, search_time);
        if (illum.status != ASTRO_SUCCESS)
            FAIL("C TestMaxMag(%s line %d): SearchPeakMagnitude returned %d\n", filename, lnum, illum.status);

        mag_diff = fabs(illum.mag - correct_mag);
        hours_diff = 24.0 * fabs(illum.time.ut - center_time.ut);
        printf("C TestMaxMag: mag_diff=%0.3lf, hours_diff=%0.3lf\n", mag_diff, hours_diff);
        if (hours_diff > 7.1)
            FAIL("C TestMaxMag(%s line %d): EXCESSIVE TIME DIFFERENCE.\n", filename, lnum);

        if (mag_diff > 0.005)
            FAIL("C TestMaxMag(%s line %d): EXCESSIVE MAGNITUDE DIFFERENCE.\n", filename, lnum);

        search_time = time2;
    }

    printf("C TestMaxMag: processed %d lines from file %s\n", lnum, filename);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static const char *ParseJplHorizonsDateTime(const char *text, astro_time_t *time)
{
    int year, month, day, hour, minute;
    int nscanned;
    char mtext[4];
    char verify[30];
    size_t length;

    time->ut = time->tt = NAN;

    while (*text == ' ' && *text != '\0')
        ++text;

    nscanned = sscanf(text, "%d-%3[A-Za-z]-%d %d:%d", &year, mtext, &day, &hour, &minute);
    if (nscanned != 5)
        return NULL;

    if (year < 1000 || year > 9999)
        return NULL;       /* if not a 4-digit year, we are skipping wrong number of chars */

    if (day < 1 || day > 31)
        return NULL;

    if (hour < 0 || hour > 23)
        return NULL;

    if (minute < 0 || minute > 59)
        return NULL;

    if (!strcmp(mtext, "Jan"))
        month = 1;
    else if (!strcmp(mtext, "Feb"))
        month = 2;
    else if (!strcmp(mtext, "Mar"))
        month = 3;
    else if (!strcmp(mtext, "Apr"))
        month = 4;
    else if (!strcmp(mtext, "May"))
        month = 5;
    else if (!strcmp(mtext, "Jun"))
        month = 6;
    else if (!strcmp(mtext, "Jul"))
        month = 7;
    else if (!strcmp(mtext, "Aug"))
        month = 8;
    else if (!strcmp(mtext, "Sep"))
        month = 9;
    else if (!strcmp(mtext, "Oct"))
        month = 10;
    else if (!strcmp(mtext, "Nov"))
        month = 11;
    else if (!strcmp(mtext, "Dec"))
        month = 12;
    else
        return NULL;

    /* Make absolutely sure we know how many characters the date text is. */
    /* If anything is fishy, fail! */
    snprintf(verify, sizeof(verify), "%04d-%s-%02d %02d:%02d", year, mtext, day, hour, minute);
    length = strlen(verify);
    if (length != 17)
        return NULL;

    if (memcmp(verify, text, length))
        return NULL;

    *time = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);
    return text + length;   /* return remaining portion of string */
}

/*-----------------------------------------------------------------------------------------------------------*/

static int LunarApsis(const char *filename)
{
    int error = 1;
    FILE *infile = NULL;
    int lnum, nscanned;
    char line[100];
    int kind, year, month, day, hour, minute;
    double dist_km;
    astro_time_t start_time, correct_time;
    astro_apsis_t apsis;
    double diff_minutes, diff_km;
    double max_minutes = 0.0, max_km = 0.0;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C LunarApsis: Cannot open input file: %s\n", filename);

    /*
        0 2001-01-10T08:59Z 357132
        1 2001-01-24T19:02Z 406565
    */

    start_time = Astronomy_MakeTime(2001, 1, 1, 0, 0, 0.0);
    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        if (lnum == 1)
            apsis = Astronomy_SearchLunarApsis(start_time);
        else
            apsis = Astronomy_NextLunarApsis(apsis);

        if (apsis.status != ASTRO_SUCCESS)
            FAIL("C LunarApsis(%s line %d): Failed to find apsis.\n", filename, lnum);

        nscanned = sscanf(line, "%d %d-%d-%dT%d:%dZ %lf", &kind, &year, &month, &day, &hour, &minute, &dist_km);
        if (nscanned != 7)
            FAIL("C LunarApsis(%s line %d): invalid data format\n", filename, lnum);

        if (kind != apsis.kind)
            FAIL("C LunarApsis(%s line %d): expected apsis kind %d but found %d\n", filename, lnum, kind, apsis.kind);

        correct_time = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);

        diff_minutes = (24.0 * 60.0) * fabs(apsis.time.ut - correct_time.ut);
        diff_km = fabs(apsis.dist_km - dist_km);

        if (diff_minutes > 35.0)
            FAIL("C LunarApsis(%s line %d): Excessive time error: %lf minutes.\n", filename, lnum, diff_minutes);

        if (diff_km > 25.0)
            FAIL("C LunarApsis(%s line %d): Excessive distance error: %lf km.\n", filename, lnum, diff_km);

        if (diff_minutes > max_minutes)
            max_minutes = diff_minutes;

        if (diff_km > max_km)
            max_km = diff_km;
    }

    printf("C LunarApsis: found %d events, max time error = %0.3lf minutes, max distance error = %0.3lf km.\n", lnum, max_minutes, max_km);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int EarthApsis(const char *filename)
{
    int error = 1;
    FILE *infile = NULL;
    int lnum, nscanned;
    char line[100];
    int kind, year, month, day, hour, minute;
    double dist_au;
    astro_time_t start_time, correct_time;
    astro_apsis_t apsis;
    double diff_minutes, diff_au;
    double max_minutes = 0.0, max_au = 0.0;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C EarthApsis: Cannot open input file: %s\n", filename);

    /*
        0 2001-01-04T08:52Z 0.9832860
        1 2001-07-04T13:37Z 1.0166426
    */

    start_time = Astronomy_MakeTime(2001, 1, 1, 0, 0, 0.0);
    lnum = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        if (lnum == 1)
            apsis = Astronomy_SearchPlanetApsis(BODY_EARTH, start_time);
        else
            apsis = Astronomy_NextPlanetApsis(BODY_EARTH, apsis);

        if (apsis.status != ASTRO_SUCCESS)
            FAIL("C EarthApsis(%s line %d): Failed to find apsis.\n", filename, lnum);

        nscanned = sscanf(line, "%d %d-%d-%dT%d:%dZ %lf", &kind, &year, &month, &day, &hour, &minute, &dist_au);
        if (nscanned != 7)
            FAIL("C EarthApsis(%s line %d): invalid data format\n", filename, lnum);

        if (kind != apsis.kind)
            FAIL("C EarthApsis(%s line %d): expected apsis kind %d but found %d\n", filename, lnum, kind, apsis.kind);

        correct_time = Astronomy_MakeTime(year, month, day, hour, minute, 0.0);

        diff_minutes = (24.0 * 60.0) * fabs(apsis.time.ut - correct_time.ut);
        diff_au = fabs(apsis.dist_au - dist_au);

        if (diff_minutes > 120.5)
            FAIL("C EarthApsis(%s line %d): Excessive time error: %lf minutes.\n", filename, lnum, diff_minutes);

        if (diff_au > 1.2e-5)
            FAIL("C EarthApsis(%s line %d): Excessive distance error: %lg AU.\n", filename, lnum, diff_au);

        if (diff_minutes > max_minutes)
            max_minutes = diff_minutes;

        if (diff_au > max_au)
            max_au = diff_au;
    }

    printf("C EarthApsis: found %d events, max time error = %0.3lf minutes, max distance error = %0.3lg AU.\n", lnum, max_minutes, max_au);
    error = 0;

fail:
    if (infile != NULL) fclose(infile);
    return error;
}

static double PlanetOrbitalPeriod(astro_body_t body)
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

static int PlanetApsis(void)
{
    int error;
    const double degree_threshold = 0.1;
    astro_body_t body;
    astro_time_t start_time, prev_time;
    astro_apsis_t apsis;
    astro_utc_t utc;
    int count;
    int bad_planets_found = 0;
    double interval, min_interval, max_interval;
    FILE *infile = NULL;
    char filename[100];
    char line[100];
    char expected_time_text[100];
    astro_time_t expected_time;
    int expected_kind;
    double expected_distance;
    double period;
    double diff_days, diff_degrees, diff_dist_ratio;
    double max_diff_days, max_dist_ratio;

    start_time = Astronomy_MakeTime(MIN_YEAR, 1, 1, 0, 0, 0.0);

    for (body = BODY_MERCURY; body <= BODY_PLUTO; ++body)
    {
        period = PlanetOrbitalPeriod(body);
        max_dist_ratio = 0.0;
        max_diff_days = 0.0;
        snprintf(filename, sizeof(filename)-1, "apsides/apsis_%d.txt", (int)body);
        if (infile) fclose(infile);
        infile = fopen(filename, "rt");
        if (infile == NULL)
            FAIL("C PlanetApsis: ERROR - cannot open input file: %s\n", filename);

        min_interval = max_interval = -1.0;
        apsis = Astronomy_SearchPlanetApsis(body, start_time);
        if (apsis.status != ASTRO_SUCCESS)
            FAIL("C PlanetApsis: ERROR %d finding first apsis for %s\n", apsis.status, Astronomy_BodyName(body));

        count = 1;
        for(;;)
        {
            if (NULL == fgets(line, sizeof(line), infile))
                break;  /* normal end of test data */

            /* Parse the line of test data. */
            if (   (3 != sscanf(line, "%d %s %lf", &expected_kind, expected_time_text, &expected_distance))
                || (expected_kind & ~1)     /* must be either 0=perihelion or 1=aphelion */
                || (expected_distance <= 0.0)
                || ParseDate(expected_time_text, &expected_time))
            {
                FAIL("C PlanetApsis: INPUT SYNTAX ERROR (%s line %d): '%s'\n", filename, count, line);
            }

            /* Compare computed values against expected values. */
            if (apsis.kind != expected_kind)
                FAIL("C PlanetApsis: WRONG APSIS KIND (%s line %d)\n", filename, count);

            diff_days = fabs(expected_time.tt - apsis.time.tt);
            if (diff_days > max_diff_days) max_diff_days = diff_days;
            diff_degrees = (diff_days / period) * 360.0;
            if (diff_degrees > degree_threshold)
                bad_planets_found = 1;

            diff_dist_ratio = fabs(expected_distance - apsis.dist_au) / expected_distance;
            if (diff_dist_ratio > max_dist_ratio) max_dist_ratio = diff_dist_ratio;
            if (diff_dist_ratio > 1.0e-4)
            {
                FAIL("C PlanetApsis: EXCESSIVE DISTANCE ERROR for %s (%s line %d): expected=%0.16lf, calculated=%0.16lf, error ratio=%lg\n",
                    Astronomy_BodyName(body), filename, count, expected_distance, apsis.dist_au, diff_dist_ratio);
            }

            /* Calculate the next apsis. */
            prev_time = apsis.time;
            utc = Astronomy_UtcFromTime(apsis.time);
            apsis = Astronomy_NextPlanetApsis(body, apsis);
            if (apsis.status == ASTRO_BAD_TIME && body == BODY_PLUTO)
                break;      /* Pluto is limited by MAX_YEAR; OK for it to fail with this error. */

            if (apsis.status != ASTRO_SUCCESS)
            {
                FAIL("C PlanetApsis: ERROR %d finding apsis for %s after %04d-%02d-%02d\n",
                    apsis.status, Astronomy_BodyName(body), utc.year, utc.month, utc.day);
            }

            /* Update statistics. */
            ++count;
            interval = apsis.time.tt - prev_time.tt;
            if (min_interval < 0.0)
            {
                min_interval = max_interval = interval;
            }
            else
            {
                if (interval < min_interval)
                    min_interval = interval;
                if (interval > max_interval)
                    max_interval = interval;
            }
        }

        if (count < 2)
            FAIL("C PlanetApsis: FAILED to find apsides for %s\n", Astronomy_BodyName(body));

        printf("C PlanetApsis: %5d apsides for %-9s -- intervals: min=%9.2lf, max=%9.2lf, ratio=%8.6lf; max day=%lg, degrees=%0.3lf, dist ratio=%lg\n",
            count, Astronomy_BodyName(body),
            min_interval, max_interval, max_interval / min_interval,
            max_diff_days,
            (max_diff_days / period) * 360.0,
            max_dist_ratio);
    }

    if (bad_planets_found)
        FAIL("C PlanetApsis: FAIL - planet(s) exceeded angular threshold (%lg degrees)\n", degree_threshold);

    printf("C PlanetApsis: PASS\n");
    error = 0;
fail:
    if (infile) fclose(infile);
    return error;
}


/*-----------------------------------------------------------------------------------------------------------*/

static int MathCheck(astro_body_t body, double ut)
{
    astro_observer_t observer = Astronomy_MakeObserver(29, -81, 10);
    astro_time_t time = Astronomy_TimeFromDays(ut);
    astro_equatorial_t j2000;
    astro_equatorial_t ofdate;
    astro_horizon_t hor;
    int error = 1;

    printf("C time.ut     = %0.16lf\n", time.ut);
    printf("C time.tt     = %0.16lf\n", time.tt);

    CHECK_EQU(j2000, Astronomy_Equator(body, &time, observer, EQUATOR_J2000, NO_ABERRATION));
    printf("C j2000  ra   = %0.16lf\n", j2000.ra);
    printf("C j2000  dec  = %0.16lf\n", j2000.dec);

    CHECK_EQU(ofdate, Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION));
    printf("C ofdate ra   = %0.16lf\n", ofdate.ra);
    printf("C ofdate dec  = %0.16lf\n", ofdate.dec);

    hor = Astronomy_Horizon(&time, observer, ofdate.ra, ofdate.dec, REFRACTION_NONE);
    printf("C azimuth     = %0.16lf\n", hor.azimuth);
    printf("C altitude    = %0.16lf\n", hor.altitude);

    error = 0;
fail:
    return error;
}

static int Issue46(void)
{
    /* https://github.com/cosinekitty/astronomy/issues/46 */
    return MathCheck(BODY_SUN, -93692.7685882873047376);
}

static int Issue48(void)
{
    /* https://github.com/cosinekitty/astronomy/issues/48 */
    return MathCheck(BODY_VENUS, -39864.1907264927140204);
}

/*-----------------------------------------------------------------------------------------------------------*/

static int CheckUnitVector(int lnum, const char *name, astro_rotation_t r, int i0, int j0, int di, int dj)
{
    int k;
    double x;
    double sum = 0.0;

    for (k=0; k<3; ++k)
    {
        x = r.rot[i0+k*di][j0+k*dj];
        sum += x*x;
    }
    x = fabs(sum - 1.0);
    if (x > 1.0e-15)
        FAILRET("C CheckUnitVector ERROR(%s line %d): unit error = %lg for i0=%d, j0=%d, di=%d, dj=%d\n", name, lnum, x, i0, j0, di, dj);

    return 0;
}

static int CheckRotationMatrix(int lnum, const char *name, astro_rotation_t r)
{
    if (r.status != ASTRO_SUCCESS)
        FAILRET("C CheckRotationMatrix ERROR(%s line %d): status = %d\n", name, lnum, r.status);

    /* Verify that every row and every column is a unit vector. */
    if (CheckUnitVector(lnum, name, r, 0, 0, 1, 0)) return 1;
    if (CheckUnitVector(lnum, name, r, 0, 1, 1, 0)) return 1;
    if (CheckUnitVector(lnum, name, r, 0, 2, 1, 0)) return 1;
    if (CheckUnitVector(lnum, name, r, 0, 0, 0, 1)) return 1;
    if (CheckUnitVector(lnum, name, r, 1, 0, 0, 1)) return 1;
    if (CheckUnitVector(lnum, name, r, 2, 0, 0, 1)) return 1;
    return 0;
}

#define CHECK_ROTMAT(r)   CHECK(CheckRotationMatrix(__LINE__, #r, (r)))

static int CompareMatrices(const char *caller, astro_rotation_t a, astro_rotation_t b, double tolerance)
{
    int i, j;

    if (a.status != ASTRO_SUCCESS)
        FAILRET("C CompareMatrices ERROR(%s): a.status = %d\n", caller, a.status);

    if (b.status != ASTRO_SUCCESS)
        FAILRET("C CompareMatrices ERROR(%s): b.status = %d\n", caller, b.status);

    for (i=0; i < 3; ++i)
    {
        for (j=0; j < 3; ++j)
        {
            double diff = fabs(a.rot[i][j] - b.rot[i][j]);
            if (diff > tolerance)
                FAILRET("C CompareMatrices ERROR(%s): matrix[%d][%d]=%lg, expected %lg, diff %lg\n", caller, i, j, a.rot[i][j], b.rot[i][j], diff);
        }
    }

    return 0;
}

static int Rotation_MatrixInverse(void)
{
    astro_rotation_t a, b, v;
    int error;

    a.status = ASTRO_SUCCESS;
    a.rot[0][0] = 1.0; a.rot[1][0] = 2.0; a.rot[2][0] = 3.0;
    a.rot[0][1] = 4.0; a.rot[1][1] = 5.0; a.rot[2][1] = 6.0;
    a.rot[0][2] = 7.0; a.rot[1][2] = 8.0; a.rot[2][2] = 9.0;

    v.status = ASTRO_SUCCESS;
    v.rot[0][0] = 1.0; v.rot[1][0] = 4.0; v.rot[2][0] = 7.0;
    v.rot[0][1] = 2.0; v.rot[1][1] = 5.0; v.rot[2][1] = 8.0;
    v.rot[0][2] = 3.0; v.rot[1][2] = 6.0; v.rot[2][2] = 9.0;

    b = Astronomy_InverseRotation(a);
    CHECK(CompareMatrices("Rotation_MatrixInverse", b, v, 0.0));

    printf("Rotation_MatrixInverse: PASS\n");
    error = 0;

fail:
    return error;
}

static int Rotation_MatrixMultiply(void)
{
    astro_rotation_t  a, b, c, v;
    int error;

    a.status = ASTRO_SUCCESS;
    a.rot[0][0] = 1.0; a.rot[1][0] = 2.0; a.rot[2][0] = 3.0;
    a.rot[0][1] = 4.0; a.rot[1][1] = 5.0; a.rot[2][1] = 6.0;
    a.rot[0][2] = 7.0; a.rot[1][2] = 8.0; a.rot[2][2] = 9.0;

    b.status = ASTRO_SUCCESS;
    b.rot[0][0] = 10.0; b.rot[1][0] = 11.0; b.rot[2][0] = 12.0;
    b.rot[0][1] = 13.0; b.rot[1][1] = 14.0; b.rot[2][1] = 15.0;
    b.rot[0][2] = 16.0; b.rot[1][2] = 17.0; b.rot[2][2] = 18.0;

    v.status = ASTRO_SUCCESS;
    v.rot[0][0] =  84.0; v.rot[1][0] =  90.0; v.rot[2][0] =  96.0;
    v.rot[0][1] = 201.0; v.rot[1][1] = 216.0; v.rot[2][1] = 231.0;
    v.rot[0][2] = 318.0; v.rot[1][2] = 342.0; v.rot[2][2] = 366.0;

    /* Let c = a*b. The order looks backwards in the call to make rotation semantics more intuitive. */
    c = Astronomy_CombineRotation(b, a);

    /* Verify that c = v. */
    CHECK(CompareMatrices("Rotation_MatrixMultiply", c, v, 0.0));

    printf("C Rotation_MatrixMultiply: PASS\n");
    error = 0;

fail:
    return error;
}

static int TestVectorFromAngles(double lat, double lon, double x, double y, double z)
{
    astro_spherical_t sphere;
    astro_vector_t vector;
    astro_time_t time;
    double diff, dx, dy, dz;

    /* Confirm the expected vector really is a unit vector. */
    diff = fabs((x*x + y*y + z*z) - 1.0);
    if (diff > 1.0e-16)
        FAILRET("C TestVectorFromAngles: EXCESSIVE unit error = %lg\n", diff);

    sphere.status = ASTRO_SUCCESS;
    sphere.lat = lat;
    sphere.lon = lon;
    sphere.dist = 1.0;

    time = Astronomy_MakeTime(2015, 3, 5, 12, 0, 0.0);
    vector = Astronomy_VectorFromSphere(sphere, time);

    if (vector.status != ASTRO_SUCCESS)
        FAILRET("C ERROR(TestVectorFromAngles): vector.status = %d\n", vector.status);

    dx = x - vector.x;
    dy = y - vector.y;
    dz = z - vector.z;
    diff = sqrt(dx*dx + dy*dy + dz*dz);

    printf("TestVectorFromAngles(%lf, %lf): diff = %lg\n", lat, lon, diff);
    if (diff > 2.0e-16)
        FAILRETSTR("C TestVectorFromAngles: EXCESSIVE ERROR.\n");

    return 0;
}

static int TestAnglesFromVector(double lat, double lon, double x, double y, double z)
{
    astro_vector_t vector;
    astro_spherical_t sphere;
    double diff, latdiff, londiff;

    /* Confirm the expected vector really is a unit vector. */
    diff = fabs((x*x + y*y + z*z) - 1.0);
    if (diff > 1.0e-16)
        FAILRET("C TestAnglesFromVector(lat=%lf, lon=%lf, x=%lf, y=%lf, z=%lf): EXCESSIVE unit error = %lg\n", lat, lon, x, y, z, diff);

    vector.status = ASTRO_SUCCESS;
    vector.t = Astronomy_MakeTime(2015, 3, 5, 12, 0, 0.0);
    vector.x = x;
    vector.y = y;
    vector.z = z;

    sphere = Astronomy_SphereFromVector(vector);
    if (sphere.status != ASTRO_SUCCESS)
        FAILRET("C ERROR TestAnglesFromVector(lat=%lf, lon=%lf, x=%lf, y=%lf, z=%lf): sphere.status = %d\n", lat, lon, x, y, z, sphere.status);

    latdiff = fabs(sphere.lat - lat);
    londiff = fabs(sphere.lon - lon);
    printf("TestAnglesFromVector(x=%lf, y=%lf, z=%lf): latdiff=%lg, londiff=%lg\n", x, y, z, latdiff, londiff);
    if (latdiff > 8.0e-15 || londiff > 8.0e-15)
        FAILRETSTR("C TestAnglesFromVector: EXCESSIVE ERROR\n");

    return 0;
}

static astro_rotation_t RotateX(double rot)
{
    astro_rotation_t m;
    double ca, sa;

    /* https://en.wikipedia.org/wiki/Rotation_matrix */

    ca = cos(rot * DEG2RAD);
    sa = sin(rot * DEG2RAD);
    m.status = ASTRO_SUCCESS;
    m.rot[0][0] = 1.0;  m.rot[1][0] = 0.0;  m.rot[2][0] = 0.0;
    m.rot[0][1] = 0.0;  m.rot[1][1] = +ca;  m.rot[2][1] = -sa;
    m.rot[0][2] = 0.0;  m.rot[1][2] = +sa;  m.rot[2][2] = +ca;

    return m;
}

static astro_rotation_t RotateY(double rot)
{
    astro_rotation_t m;
    double ca, sa;

    /* https://en.wikipedia.org/wiki/Rotation_matrix */

    ca = cos(rot * DEG2RAD);
    sa = sin(rot * DEG2RAD);
    m.status = ASTRO_SUCCESS;
    m.rot[0][0] = +ca;  m.rot[1][0] = 0.0;  m.rot[2][0] = +sa;
    m.rot[0][1] = 0.0;  m.rot[1][1] = 1.0;  m.rot[2][1] = 0.0;
    m.rot[0][2] = -sa;  m.rot[1][2] = 0.0;  m.rot[2][2] = +ca;

    return m;
}

static astro_rotation_t RotateZ(double rot)
{
    astro_rotation_t m;
    double ca, sa;

    /* https://en.wikipedia.org/wiki/Rotation_matrix */

    ca = cos(rot * DEG2RAD);
    sa = sin(rot * DEG2RAD);
    m.status = ASTRO_SUCCESS;
    m.rot[0][0] = +ca;  m.rot[1][0] = -sa;  m.rot[2][0] = 0.0;
    m.rot[0][1] = +sa;  m.rot[1][1] = +ca;  m.rot[2][1] = 0.0;
    m.rot[0][2] = 0.0;  m.rot[1][2] = 0.0;  m.rot[2][2] = 1.0;

    return m;
}

static int TestSpin(
    double xrot, double yrot, double zrot,
    double sx, double sy, double sz,
    double tx, double ty, double tz)
{
    astro_rotation_t mx, my, mz, m;
    astro_vector_t sv, tv;
    double diff, dx, dy, dz;
    int error;

    /* https://en.wikipedia.org/wiki/Rotation_matrix */

    mx = RotateX(xrot);
    CHECK_ROTMAT(mx);
    my = RotateY(yrot);
    CHECK_ROTMAT(my);
    mz = RotateZ(zrot);
    CHECK_ROTMAT(mz);
    m = Astronomy_CombineRotation(mx, my);
    CHECK_ROTMAT(m);
    m = Astronomy_CombineRotation(m, mz);
    CHECK_ROTMAT(m);

    sv.status = ASTRO_SUCCESS;
    sv.x = sx;
    sv.y = sy;
    sv.z = sz;
    sv.t = Astronomy_MakeTime(2019, 5, 5, 12, 30, 0.0);

    tv = Astronomy_RotateVector(m, sv);
    if (tv.status != ASTRO_SUCCESS)
        FAIL("C ERROR(TestSpin): tv.status = %d\n", tv.status);

    if (tv.t.ut != sv.t.ut)
        FAILSTR("C ERROR(TestSpin): tv time != sv time\n");

    dx = tx - tv.x;
    dy = ty - tv.y;
    dz = tz - tv.z;
    diff = sqrt(dx*dx + dy*dy + dz*dz);
    printf("C TestSpin(xrot=%0.0lf, yrot=%0.0lf, zrot=%0.0lf, sx=%0.0lf, sy=%0.0lf, sz=%0.0lf): diff = %lg\n", xrot, yrot, zrot, sx, sy, sz, diff);
    if (diff > 1.0e-15)
        FAILSTR("C TestSpin: EXCESSIVE ERROR\n");

    error = 0;
fail:
    return error;
}

static int Test_EQJ_ECL(void)
{
    astro_rotation_t r;
    astro_time_t time;
    astro_vector_t ev;      /* Earth vector in equatorial J2000 */
    astro_ecliptic_t ecl;
    astro_vector_t ee;      /* Earth vector in ecliptic J2000 */
    astro_vector_t et;      /* Test of reconstructing original Earth vector in equatorial J2000 */
    double diff, dx, dy, dz;
    int error;

    r = Astronomy_Rotation_EQJ_ECL();
    CHECK_ROTMAT(r);

    printf("C Test_EQJ_ECL:\n[%0.18lf  %0.18lf  %0.18lf]\n[%0.18lf  %0.18lf  %0.18lf]\n[%0.18lf  %0.18lf  %0.18lf]\n",
        r.rot[0][0], r.rot[1][0], r.rot[2][0],
        r.rot[0][1], r.rot[1][1], r.rot[2][1],
        r.rot[0][2], r.rot[1][2], r.rot[2][2]);

    /* Calculate heliocentric Earth position at a test time (when I wrote this code). */
    time = Astronomy_MakeTime(2019, 12, 8, 19, 39, 15.0);
    ev = Astronomy_HelioVector(BODY_EARTH, time);
    if (ev.status != ASTRO_SUCCESS)
        FAIL("C Test_EQJ_ECL: Astronomy_HelioVector returned error %d\n", ev.status);

    /* Use the existing Astronomy_Ecliptic() to calculate ecliptic vector and angles. */
    ecl = Astronomy_Ecliptic(ev);
    if (ecl.status != ASTRO_SUCCESS)
        FAIL("C Test_EQJ_ECL: Astronomy_Ecliptic returned error %d\n", ecl.status);

    printf("C Test_EQJ_ECL ecl = (%0.18lf, %0.18lf,%0.18lf)\n", ecl.ex, ecl.ey, ecl.ez);

    /* Now compute the same vector via rotation matrix. */
    ee = Astronomy_RotateVector(r, ev);
    if (ee.status != ASTRO_SUCCESS)
        FAIL("C Test_EQJ_ECL: Astronomy_RotateVector returned error %d\n", ee.status);

    dx = ee.x - ecl.ex;
    dy = ee.y - ecl.ey;
    dz = ee.z - ecl.ez;
    diff = sqrt(dx*dx + dy*dy + dz*dz);
    printf("C Test_EQJ_ECL  ee = (%0.18lf, %0.18lf,%0.18lf);  diff=%lg\n", ee.x, ee.y, ee.z, diff);
    if (diff > 1.0e-16)
        FAILSTR("C Test_EQJ_ECL: EXCESSIVE VECTOR ERROR\n");

    /* Reverse the test: go from ecliptic back to equatorial. */
    r = Astronomy_Rotation_ECL_EQJ();
    CHECK_ROTMAT(r);
    et = Astronomy_RotateVector(r, ee);
    CHECK(VectorDiff(et, ev, &diff));
    printf("C Test_EQJ_ECL  ev diff=%lg\n", diff);
    if (diff > 2.0e-16)
        FAILSTR("C Test_EQJ_ECL: EXCESSIVE REVERSE ROTATION ERROR\n");

    error = 0;
fail:
    return error;
}

static int Test_EQJ_EQD(astro_body_t body)
{
    astro_time_t time;
    astro_observer_t observer;
    astro_equatorial_t eq2000, eqdate, eqcheck;
    astro_vector_t v2000, vdate, t2000;
    astro_rotation_t r;
    double ra_diff, dec_diff, dist_diff, diff;
    int error;

    /* Verify conversion of equatorial J2000 to equatorial of-date, and back. */

    /* Use established functions to calculate spherical coordinates for the body, in both EQJ and EQD. */
    time = Astronomy_MakeTime(2019, 12, 8, 20, 50, 0.0);
    observer = Astronomy_MakeObserver(35, -85, 0);
    eq2000 = Astronomy_Equator(body, &time, observer, EQUATOR_J2000, ABERRATION);
    CHECK_STATUS(eq2000);
    eqdate = Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    CHECK_STATUS(eqdate);

    /* Convert EQJ angular coordinates to vector. */
    v2000 = Astronomy_VectorFromEquator(eq2000, time);
    CHECK_STATUS(v2000);

    /* Find rotation matrix. */
    r = Astronomy_Rotation_EQJ_EQD(time);
    CHECK_ROTMAT(r);

    /* Rotate EQJ vector to EQD vector. */
    vdate = Astronomy_RotateVector(r, v2000);
    CHECK_STATUS(vdate);

    /* Convert vector back to angular coordinates. */
    eqcheck = Astronomy_EquatorFromVector(vdate);
    CHECK_STATUS(eqcheck);

    /* Compare the result with the eqdate. */
    ra_diff = fabs(eqcheck.ra - eqdate.ra);
    dec_diff = fabs(eqcheck.dec - eqdate.dec);
    dist_diff = fabs(eqcheck.dist - eqdate.dist);
    printf("C Test_EQJ_EQD: %s ra=%0.3lf, dec=%0.3lf, dist=%0.3lf, ra_diff=%lg, dec_diff=%lg, dist_diff=%lg\n", Astronomy_BodyName(body), eqdate.ra, eqdate.dec, eqdate.dist, ra_diff, dec_diff, dist_diff);
    if (ra_diff > 1.0e-14 || dec_diff > 1.0e-14 || dist_diff > 4.0e-15)
        FAILSTR("C Test_EQJ_EQD: EXCESSIVE ERROR\n");

    /* Perform the inverse conversion back to equatorial J2000 coordinates. */
    r = Astronomy_Rotation_EQD_EQJ(time);
    CHECK_ROTMAT(r);
    t2000 = Astronomy_RotateVector(r, vdate);
    CHECK_STATUS(t2000);
    CHECK(VectorDiff(t2000, v2000, &diff));
    printf("C Test_EQJ_EQD: %s inverse diff = %lg\n", Astronomy_BodyName(body), diff);
    if (diff > 3.0e-15)
        FAILSTR("C Test_EQJ_EQD: EXCESSIVE INVERSE ERROR\n");

    error = 0;
fail:
    return error;
}

static int Test_EQD_HOR(astro_body_t body)
{
    astro_time_t time;
    astro_observer_t observer;
    astro_equatorial_t eqd, eqj;
    astro_horizon_t hor;
    astro_rotation_t rot;
    astro_spherical_t sphere;
    astro_vector_t vec_eqd, vec_eqj, vec_hor, check_eqd, check_eqj, check_hor;
    double diff_az, diff_alt, diff;
    int error;

    /* Use existing functions to calculate horizontal coordinates of the body for the time+observer. */
    time = Astronomy_MakeTime(1970, 12, 13, 5, 15, 0.0);
    observer = Astronomy_MakeObserver(-37.0, +45.0, 0.0);
    CHECK_EQU(eqd, Astronomy_Equator(body, &time, observer, EQUATOR_OF_DATE, ABERRATION));
    printf("C Test_EQD_HOR %s: OFDATE ra=%0.3lf, dec=%0.3lf\n", Astronomy_BodyName(body), eqd.ra, eqd.dec);
    hor = Astronomy_Horizon(&time, observer, eqd.ra, eqd.dec, REFRACTION_NORMAL);

    /* Calculate the position of the body as an equatorial vector of date. */
    CHECK_VECTOR(vec_eqd, Astronomy_VectorFromEquator(eqd, time));

    /* Calculate rotation matrix to convert equatorial J2000 vector to horizontal vector. */
    rot = Astronomy_Rotation_EQD_HOR(time, observer);
    CHECK_ROTMAT(rot);

    /* Rotate the equator of date vector to a horizontal vector. */
    CHECK_VECTOR(vec_hor, Astronomy_RotateVector(rot, vec_eqd));

    /* Convert the horizontal vector to horizontal angular coordinates. */
    sphere = Astronomy_HorizonFromVector(vec_hor, REFRACTION_NORMAL);
    CHECK_STATUS(sphere);

    diff_alt = fabs(sphere.lat - hor.altitude);
    diff_az = fabs(sphere.lon - hor.azimuth);

    printf("C Test_EQD_HOR %s: trusted alt=%0.3lf, az=%0.3lf; test alt=%0.3lf, az=%0.3lf; diff_alt=%lg, diff_az=%lg\n",
        Astronomy_BodyName(body), hor.altitude, hor.azimuth, sphere.lat, sphere.lon, diff_alt, diff_az);

    if (diff_alt > 2.0e-14 || diff_az > 4e-14)
        FAILSTR("C Test_EQD_HOR: EXCESSIVE HORIZONTAL ERROR.\n");

    /* Confirm that we can convert back to horizontal vector. */
    CHECK_VECTOR(check_hor, Astronomy_VectorFromHorizon(sphere, time, REFRACTION_NORMAL));
    CHECK(VectorDiff(check_hor, vec_hor, &diff));
    printf("C Test_EQD_HOR %s: horizontal recovery: diff = %lg\n", Astronomy_BodyName(body), diff);
    if (diff > 2.0e-15)
        FAILSTR("C Test_EQD_HOR: EXCESSIVE ERROR IN HORIZONTAL RECOVERY.\n");

    /* Verify the inverse translation from horizontal vector to equatorial of-date vector. */
    rot = Astronomy_Rotation_HOR_EQD(time, observer);
    CHECK_ROTMAT(rot);
    CHECK_VECTOR(check_eqd, Astronomy_RotateVector(rot, vec_hor));
    CHECK(VectorDiff(check_eqd, vec_eqd, &diff));
    printf("C Test_EQD_HOR %s: OFDATE inverse rotation diff = %lg\n", Astronomy_BodyName(body), diff);
    if (diff > 2.0e-15)
        FAILSTR("C Test_EQD_HOR: EXCESSIVE OFDATE INVERSE HORIZONTAL ERROR.\n");

    /* Exercise HOR to EQJ translation. */
    CHECK_EQU(eqj, Astronomy_Equator(body, &time, observer, EQUATOR_J2000, ABERRATION));
    CHECK_VECTOR(vec_eqj, Astronomy_VectorFromEquator(eqj, time));

    rot = Astronomy_Rotation_HOR_EQJ(time, observer);
    CHECK_ROTMAT(rot);
    CHECK_VECTOR(check_eqj, Astronomy_RotateVector(rot, vec_hor));
    CHECK(VectorDiff(check_eqj, vec_eqj, &diff));
    printf("C Test_EQD_HOR %s: J2000 inverse rotation diff = %lg\n", Astronomy_BodyName(body), diff);
    if (diff > 4.0e-15)
        FAILSTR("C Test_EQD_HOR: EXCESSIVE J2000 INVERSE HORIZONTAL ERROR.\n");

    /* Verify the inverse translation: EQJ to HOR. */
    rot = Astronomy_Rotation_EQJ_HOR(time, observer);
    CHECK_ROTMAT(rot);
    CHECK_VECTOR(check_hor, Astronomy_RotateVector(rot, vec_eqj));
    CHECK(VectorDiff(check_hor, vec_hor, &diff));
    printf("C Test_EQD_HOR %s: EQJ inverse rotation diff = %lg\n", Astronomy_BodyName(body), diff);
    if (diff > 3.0e-15)
        FAILSTR("C Test_EQD_HOR: EXCESSIVE EQJ INVERSE HORIZONTAL ERROR.\n");

    error = 0;
fail:
    return error;
}

static int CheckInverse(const char *aname, const char *bname, astro_rotation_t arot, astro_rotation_t brot)
{
    int error;
    astro_rotation_t crot;
    astro_rotation_t identity;

    CHECK(CheckRotationMatrix(__LINE__, aname, arot));
    CHECK(CheckRotationMatrix(__LINE__, bname, brot));

    crot = Astronomy_CombineRotation(arot, brot);
    CHECK_ROTMAT(crot);

    identity.status = ASTRO_SUCCESS;
    identity.rot[0][0] = 1; identity.rot[1][0] = 0; identity.rot[2][0] = 0;
    identity.rot[0][1] = 0; identity.rot[1][1] = 1; identity.rot[2][1] = 0;
    identity.rot[0][2] = 0; identity.rot[1][2] = 0; identity.rot[2][2] = 1;

    CHECK(CompareMatrices("CheckInverse", crot, identity, 2.0e-15));

    error = 0;
fail:
    return error;
}

#define CHECK_INVERSE(a,b)   CHECK(CheckInverse(#a, #b, a, b))

static int CheckCycle(
    const char *cyclename,
    astro_rotation_t arot, astro_rotation_t brot, astro_rotation_t crot)
{
    int error;
    astro_rotation_t xrot;

    xrot = Astronomy_CombineRotation(arot, brot);
    CHECK_ROTMAT(xrot);
    xrot = Astronomy_InverseRotation(xrot);
    CHECK(CompareMatrices(cyclename, crot, xrot, 2.0e-15));
    error = 0;
fail:
    return error;
}

#define CHECK_CYCLE(a,b,c)  CHECK(CheckCycle(("CheckCycle: " #a ", " #b ", " #c), a, b, c))

static int Test_RotRoundTrip(void)
{
    int error;
    astro_time_t time;
    astro_observer_t observer;
    astro_rotation_t eqj_ecl, ecl_eqj;
    astro_rotation_t eqj_hor, hor_eqj;
    astro_rotation_t eqj_eqd, eqd_eqj;
    astro_rotation_t hor_eqd, eqd_hor;
    astro_rotation_t eqd_ecl, ecl_eqd;
    astro_rotation_t hor_ecl, ecl_hor;

    time = Astronomy_MakeTime(2067, 5, 30, 14, 45, 0.0);
    observer = Astronomy_MakeObserver(+28.0, -82.0, 0.0);

    /*
        In each round trip, calculate a forward rotation and a backward rotation.
        Verify the two are inverse matrices.
    */

    /* Round trip #1: EQJ <==> EQD. */
    eqj_eqd = Astronomy_Rotation_EQJ_EQD(time);
    eqd_eqj = Astronomy_Rotation_EQD_EQJ(time);
    CHECK_INVERSE(eqj_eqd, eqd_eqj);

    /* Round trip #2: EQJ <==> ECL. */
    eqj_ecl = Astronomy_Rotation_EQJ_ECL();
    ecl_eqj = Astronomy_Rotation_ECL_EQJ();
    CHECK_INVERSE(eqj_ecl, ecl_eqj);

    /* Round trip #3: EQJ <==> HOR. */
    eqj_hor = Astronomy_Rotation_EQJ_HOR(time, observer);
    hor_eqj = Astronomy_Rotation_HOR_EQJ(time, observer);
    CHECK_INVERSE(eqj_hor, hor_eqj);

    /* Round trip #4: EQD <==> HOR. */
    eqd_hor = Astronomy_Rotation_EQD_HOR(time, observer);
    hor_eqd = Astronomy_Rotation_HOR_EQD(time, observer);
    CHECK_INVERSE(eqd_hor, hor_eqd);

    /* Round trip #5: EQD <==> ECL. */
    eqd_ecl = Astronomy_Rotation_EQD_ECL(time);
    ecl_eqd = Astronomy_Rotation_ECL_EQD(time);
    CHECK_INVERSE(eqd_ecl, ecl_eqd);

    /* Round trip #6: HOR <==> ECL. */
    hor_ecl = Astronomy_Rotation_HOR_ECL(time, observer);
    ecl_hor = Astronomy_Rotation_ECL_HOR(time, observer);
    CHECK_INVERSE(hor_ecl, ecl_hor);

    /*
        Verify that combining different sequences of rotations result
        in the expected combination.
        For example, (EQJ ==> HOR ==> ECL) must be the same matrix as (EQJ ==> ECL).
        Each of these is a "triangle" of relationships between 3 orientations.
        There are 4 possible ways to pick 3 orientations from the 4 to form a triangle.
        Because we have just proved that each transformation is reversible,
        we only need to verify the triangle in one cyclic direction.
    */
    CHECK_CYCLE(eqj_ecl, ecl_eqd, eqd_eqj);     /* excluded corner = HOR */
    CHECK_CYCLE(eqj_hor, hor_ecl, ecl_eqj);     /* excluded corner = EQD */
    CHECK_CYCLE(eqj_hor, hor_eqd, eqd_eqj);     /* excluded corner = ECL */
    CHECK_CYCLE(ecl_eqd, eqd_hor, hor_ecl);     /* excluded corner = EQJ */

    printf("C Test_RotRoundTrip: PASS\n");
    error = 0;
fail:
    return error;
}

static int RotationTest(void)
{
    int error;
    CHECK(Rotation_MatrixInverse());
    CHECK(Rotation_MatrixMultiply());

    /* Verify conversion of spherical coordinates to vector. */
    CHECK(TestVectorFromAngles(0.0, 0.0, 1.0, 0.0, 0.0));
    CHECK(TestVectorFromAngles(0.0, 90.0, 0.0, 1.0, 0.0));
    CHECK(TestVectorFromAngles(0.0, 180.0, -1.0, 0.0, 0.0));
    CHECK(TestVectorFromAngles(0.0, 270.0, 0.0, -1.0, 0.0));
    CHECK(TestVectorFromAngles(+90.0, 0.0, 0.0, 0.0, 1.0));
    CHECK(TestVectorFromAngles(-90.0, 0.0, 0.0, 0.0, -1.0));
    CHECK(TestVectorFromAngles(-30.0, +60.0, 0.43301270189221946, 0.75, -0.5));

    /* Verify conversion of vector to spherical coordinates. */
    CHECK(TestAnglesFromVector(0.0, 0.0, 1.0, 0.0, 0.0));
    CHECK(TestAnglesFromVector(0.0, 90.0, 0.0, 1.0, 0.0));
    CHECK(TestAnglesFromVector(0.0, 180.0, -1.0, 0.0, 0.0));
    CHECK(TestAnglesFromVector(0.0, 270.0, 0.0, -1.0, 0.0));
    CHECK(TestAnglesFromVector(+90.0, 0.0, 0.0, 0.0, 1.0));
    CHECK(TestAnglesFromVector(-90.0, 0.0, 0.0, 0.0, -1.0));
    CHECK(TestAnglesFromVector(-30.0, +60.0, 0.43301270189221946, 0.75, -0.5));

    /* Verify rotation of vectors */
    CHECK(TestSpin(0.0, 0.0, 90.0, +1, +2, +3, -2, +1, +3));
    CHECK(TestSpin(0.0, 0.0, 0.0, 1.0, 2.0, -3.0, 1.0, 2.0, -3.0));
    CHECK(TestSpin(90.0, 0.0, 0.0, +1, +2, +3, +1, -3, +2));
    CHECK(TestSpin(0.0, 90.0, 0.0, +1, +2, +3, +3, +2, -1));
    CHECK(TestSpin(180.0, 180.0, 180.0, +1, +2, +3, +1, +2, +3));
    CHECK(TestSpin(0.0, 0.0, -45.0, +1, 0, 0, +0.7071067811865476, -0.7071067811865476, 0));

    CHECK(Test_EQJ_ECL());

    CHECK(Test_EQJ_EQD(BODY_MERCURY));
    CHECK(Test_EQJ_EQD(BODY_VENUS));
    CHECK(Test_EQJ_EQD(BODY_MARS));
    CHECK(Test_EQJ_EQD(BODY_JUPITER));
    CHECK(Test_EQJ_EQD(BODY_SATURN));

    CHECK(Test_EQD_HOR(BODY_MERCURY));
    CHECK(Test_EQD_HOR(BODY_VENUS));
    CHECK(Test_EQD_HOR(BODY_MARS));
    CHECK(Test_EQD_HOR(BODY_JUPITER));
    CHECK(Test_EQD_HOR(BODY_SATURN));

    CHECK(Test_RotRoundTrip());

    error = 0;
fail:
    return error;
}

static int VectorDiff(astro_vector_t a, astro_vector_t b, double *diff)
{
    double dx, dy, dz;

    *diff = 1.0e+99;

    if (a.status != ASTRO_SUCCESS)
        FAILRET("C VectorDiff: ERROR - first vector has status %d\n", a.status);

    if (b.status != ASTRO_SUCCESS)
        FAILRET("C VectorDiff: ERROR - second vector has status %d\n", b.status);

    dx = a.x - b.x;
    dy = a.y - b.y;
    dz = a.z - b.z;
    *diff = sqrt(dx*dx + dy*dy + dz*dz);
    return 0;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int RefractionTest(void)
{
    double alt, corrected, refr, inv_refr, check_alt, diff;

    for (alt = -90.1; alt <= +90.1; alt += 0.001)
    {
        refr = Astronomy_Refraction(REFRACTION_NORMAL, alt);
        corrected = alt + refr;
        inv_refr = Astronomy_InverseRefraction(REFRACTION_NORMAL, corrected);
        check_alt = corrected + inv_refr;
        diff = fabs(check_alt - alt);
        if (diff > 2.0e-14)
        {
            printf("C ERROR(RefractionTest): alt=%8.3lf, refr=%10.6lf, diff=%lg\n", alt, refr, diff);
            return 1;
        }
    }

    printf("C RefractionTest: PASS\n");
    return 0;
}

/*-----------------------------------------------------------------------------------------------------------*/

static int ConstellationTest(void)
{
    int error = 1;
    FILE *infile = NULL;
    const char *inFileName = "constellation/test_input.txt";
    char line[100];
    double ra, dec;
    char symbol[4];
    int lnum, failcount, id;
    astro_constellation_t constel;

    infile = fopen(inFileName, "rt");
    if (infile == NULL)
        FAIL("C ConstellationTest: Cannot open input file: %s\n", inFileName);

    lnum = 0;
    failcount = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (4 != sscanf(line, "%d %lf %lf %3s", &id, &ra, &dec, symbol) || 3 != strlen(symbol))
            FAIL("C ConstellationTest: Invalid test data in %s line %d\n", inFileName, lnum);

        constel = Astronomy_Constellation(ra, dec);
        if (constel.status != ASTRO_SUCCESS)
            FAIL("C ConstellationTest: FAILED star %d with status %d: %s line %d\n", id, constel.status, inFileName, lnum);

        if (constel.symbol == NULL || constel.name == NULL)
            FAIL("C ConstellationTest: UNEXPECTED NULL name/symbol: star %d, %s line %d\n", id, inFileName, lnum);

        if (strcmp(constel.symbol, symbol))
        {
            fprintf(stderr, "Star %6d: expected %s, found %s at B1875 RA=%10.6lf, DEC=%10.6lf\n", id, symbol, constel.symbol, constel.ra_1875, constel.dec_1875);
            ++failcount;
        }
    }

    if (failcount > 0)
        FAIL("C ConstellationTest: %d failures\n", failcount);

    printf("C ConstellationTest: PASS (verified %d)\n", lnum);
    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/

static void PrintTime(astro_time_t time)
{
    astro_utc_t utc = Astronomy_UtcFromTime(time);
    printf("%04d-%02d-%02dT%02d:%02d:%06.3lfZ", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
}

static int LunarEclipseTest(void)
{
    const char *filename = "eclipse/lunar_eclipse.txt";
    FILE *infile = NULL;
    int lnum = 0;
    int error = 1;
    char line[100];
    astro_time_t peak_time;
    astro_lunar_eclipse_t eclipse;
    double partial_minutes, total_minutes;
    double diff_minutes, max_diff_minutes = 0.0, sum_diff_minutes = 0.0;
    int diff_count = 0;
    double diff_days;
    int valid = 0;
    int skip_count = 0;
    const double diff_limit = 6.853;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("C LunarEclipseTest: Cannot open input file: %s\n", filename);

    eclipse = Astronomy_SearchLunarEclipse(Astronomy_MakeTime(1701, 1, 1, 0, 0, 0.0));
    CHECK_STATUS(eclipse);

    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;

        /* scan test data */

        /* 2021-05-26T11:19Z  94   9 */
        if (strlen(line) < 17)
            FAIL("C LunarEclipseTest(%s line %d): line is too short.\n", filename, lnum);

        CHECK(ParseDate(line, &peak_time));
        if (2 != sscanf(line+17, "%lf %lf", &partial_minutes, &total_minutes))
            FAIL("C LunarEclipseTest(%s line %d): invalid data format.\n", filename, lnum);

        /* verify that the calculated eclipse semi-durations are consistent with the kind */

        switch (eclipse.kind)
        {
        case ECLIPSE_PENUMBRAL:
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial == 0.0) && (eclipse.sd_total == 0.0);
            break;

        case ECLIPSE_PARTIAL:
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total == 0.0);
            break;

        case ECLIPSE_TOTAL:
            valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total > 0.0);
            break;

        default:
            FAIL("C LunarEclipseTest(%s line %d): invalid eclipse kind %d.\n", filename, lnum, eclipse.kind);
        }

        if (!valid)
        {
            FAIL("C LunarEclipseTest(%s line %d): invalid semiduration(s) for kind %d: penum=%lf, partial=%lf, total=%lf.\n",
                filename, lnum, eclipse.kind, eclipse.sd_penum, eclipse.sd_partial, eclipse.sd_total);
        }

        /* check eclipse center */

        diff_days = eclipse.center.ut - peak_time.ut;
        /* tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse. */
        if (partial_minutes == 0.0 && diff_days > 20.0)
        {
            ++skip_count;
            continue;
        }

        diff_minutes = (24.0 * 60.0) * fabs(diff_days);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit)
        {
            printf("C LunarEclipseTest expected center: ");
            PrintTime(peak_time);
            printf("\n");
            printf("C LunarEclipseTest found    center: ");
            PrintTime(eclipse.center);
            printf("\n");
            FAIL("C LunarEclipseTest(%s line %d): EXCESSIVE center time error = %lf minutes (%lf days).\n", filename, lnum, diff_minutes, diff_days);
        }

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* check partial eclipse duration */

        diff_minutes = fabs(partial_minutes - eclipse.sd_partial);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit)
            FAIL("C LunarEclipseTest(%s line %d): EXCESSIVE partial eclipse semiduration error: %lf minutes\n", filename, lnum, diff_minutes);

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* check total eclipse duration */

        diff_minutes = fabs(total_minutes - eclipse.sd_total);
        sum_diff_minutes += diff_minutes;
        ++diff_count;

        if (diff_minutes > diff_limit)
            FAIL("C LunarEclipseTest(%s line %d): EXCESSIVE total eclipse semiduration error: %lf minutes\n", filename, lnum, diff_minutes);

        if (diff_minutes > max_diff_minutes)
            max_diff_minutes = diff_minutes;

        /* calculate for next iteration */

        eclipse = Astronomy_NextLunarEclipse(eclipse.center);
        if (eclipse.status != ASTRO_SUCCESS)
            FAIL("C LunarEclipseTest(%s line %d): Astronomy_NextLunarEclipse returned status %d\n", filename, lnum, eclipse.status);
    }

    printf("C LunarEclipseTest: PASS (verified %d, skipped %d, max_diff_minutes = %0.3lf, avg_diff_minutes = %0.3lf)\n", lnum, skip_count, max_diff_minutes, (sum_diff_minutes / diff_count));
    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    return error;
}

/*-----------------------------------------------------------------------------------------------------------*/
