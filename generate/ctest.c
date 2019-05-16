/*
    ctest.c  -  Don Cross <cosinekitty.com>

    C langauge unit test for Astronomy Engine project.
    https://cosinekitty.github.io/astronomy
*/

#include <stdio.h>
#include <math.h>
#include "astronomy.h"

#define CHECK(x)   do{if(0 != (error = (x))) goto fail;}while(0)

static int Test_AstroTime(void);

int main()
{
    int error;

    CHECK(Test_AstroTime());

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
