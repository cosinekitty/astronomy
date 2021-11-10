/*
    gravsim_test.c  -  Don Cross <cosinekitty@gmail.com>  -  2021-11-09.

    Code to measure accuracy of different numeric integrator approaches
    for calculating the orbit of Pluto.
*/

#include <stdio.h>
#include "astronomy.h"

int main(void)
{
    astro_time_t time;

    time = Astronomy_MakeTime(2020, 1, 1, 0, 0, 0.0);
    Astronomy_HelioVector(BODY_PLUTO, time);

    return 0;
}
