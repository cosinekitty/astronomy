/*
    gravsim_test.c  -  Don Cross <cosinekitty@gmail.com>  -  2021-11-09.

    Code to measure accuracy of different numeric integrator approaches
    for calculating the orbit of Pluto.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "astronomy.h"
#include "top2013.h"
#include "codegen.h"        // for CHECK macro
#include "pluto_gravsim.h"

static int DebugMode;

int TestPluto(double tt, const top_model_t *model, double *arcmin)
{
    astro_time_t time;
    astro_vector_t pos;
    top_rectangular_t equ;
    double dx, dy, dz;
    int error;

    *arcmin = 999.0;    /* bogus value in case of failure */

    if (DebugMode) printf("tt = %0.1lf\n", tt);

    time = Astronomy_TerrestrialTime(tt);
    pos = Astronomy_HelioVector(BODY_PLUTO, time);
    if (pos.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "TestPluto(%0.0lf): HelioVector returned %d\n", tt, (int)pos.status);
        return 1;
    }

    if (DebugMode) printf("HelioVector : pos=(%23.16lf, %23.16lf, %23.16lf)\n", pos.x, pos.y, pos.z);

    /* Compare with untruncated TOP2013 model. */
    CHECK(TopPosition(model, tt, &equ));
    if (DebugMode) printf("TOP2013     : pos=(%23.16lf, %23.16lf, %23.16lf)\n", equ.x, equ.y, equ.z);

    /* Calculate relative error and scale in arcminute units (as seen from the Sun). */
    dx = equ.x - pos.x;
    dy = equ.y - pos.y;
    dz = equ.z - pos.z;
    *arcmin = (RAD2DEG * 60.0) * sqrt((dx*dx + dy*dy + dz*dz) / (equ.x*equ.x + equ.y*equ.y + equ.z*equ.z));
    if (DebugMode) printf("arcmin_error = %0.6lf\n", *arcmin);
    if (DebugMode) printf("\n");

    error = 0;
fail:
    return error;
}

int main(int argc, const char *argv[])
{
    top_model_t model;
    double arcmin, tt, dev;
    int error, n, count;

    DebugMode = (argc > 1) && !strcmp(argv[1], "-v");

    printf("PLUTO_TIME_STEP = %d\n", PLUTO_TIME_STEP);
    printf("PLUTO_NSTEPS    = %d\n", PLUTO_NSTEPS);

    TopInitModel(&model);
    CHECK(TopLoadModel(&model, "../TOP2013.dat", 9));

    /* Verify that we are in sync at exact checkpoints. */
    if (DebugMode) printf("(Expect exact matches...)\n");
    for (n = -5; n < 5; ++n)
    {
        tt = n * PLUTO_TIME_STEP;
        CHECK(TestPluto(tt, &model, &arcmin));
        if (arcmin > 1.0e-6)
            FAIL("EXCESSIVE ERROR at node %d: %lf arcmin\n", n, arcmin);
    }

    /* Now check worst-case at near-midpoints. */
    if (DebugMode) printf("(Expect simulation differences...)\n");
    count = 0;
    dev = 0.0;
    for (n = -5; n < 5; ++n)
    {
        /* Calculate a worst case value: the midpoint of a step, near the midpoint of a segment. */
        /* Each segment is 36500 days, and each step is 250 days. */
        tt = ((n + 0.5) * PLUTO_TIME_STEP) + (PLUTO_DT / 2.0);
        CHECK(TestPluto(tt, &model, &arcmin));
        dev += arcmin * arcmin;
        ++count;
    }

    dev = sqrt(dev / count);
    printf("gravsim_test.c: Pluto score = %0.6lf arcmin, over %d data points.\n", dev, count);
    if (dev > 0.058327)
        FAIL("gravsim_test.c: EXCESSIVE ERROR\n");

    error = 0;
fail:
    TopFreeModel(&model);
    return error;
}
