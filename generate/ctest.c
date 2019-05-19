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

static int AdHoc(void);
static int Test_AstroTime(void);
static int AstroCheck(void);

int main(int argc, const char *argv[])
{
    int error = 0;

    if (argc == 2 && !strcmp(argv[1], "adhoc"))
    {
        CHECK(AdHoc());      /* ad hoc test for debugging */
    }
    else
    {
        CHECK(Test_AstroTime());
        CHECK(AstroCheck());
    }

fail:
    printf("ctest exiting with %d\n", error);
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

static int AdHoc(void)
{
    /* Put ad-hoc stuff here to test. */
    return 0;
}
