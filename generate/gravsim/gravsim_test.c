/*
    gravsim_test.c  -  Don Cross <cosinekitty@gmail.com>  -  2021-11-09.

    Code to measure accuracy of different numeric integrator approaches
    for calculating the orbit of Pluto.
*/

#include <stdio.h>
#include "astronomy.h"

int TestYear(int year)
{
    astro_time_t time;
    astro_vector_t pos;

    time = Astronomy_MakeTime(year, 1, 1, 0, 0, 0.0);
    pos = Astronomy_HelioVector(BODY_PLUTO, time);
    if (pos.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "TestYear(%d): HelioVector returned %d\n", year, (int)pos.status);
        return 1;
    }
    printf("year=%d, tt=%0.6lf, pos=(%0.16lf, %0.16lf, %0.16lf)\n", year, pos.t.tt, pos.x, pos.y, pos.z);
    return 0;
}

int main(void)
{
    if (TestYear(1950)) return 1;
    if (TestYear(2050)) return 1;
    if (TestYear(2150)) return 1;
    return 0;
}
