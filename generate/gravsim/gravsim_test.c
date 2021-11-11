/*
    gravsim_test.c  -  Don Cross <cosinekitty@gmail.com>  -  2021-11-09.

    Code to measure accuracy of different numeric integrator approaches
    for calculating the orbit of Pluto.
*/

#include <stdio.h>
#include <math.h>
#include "astronomy.h"
#include "top2013.h"

int TestPluto(double tt, const top_model_t *model)
{
    astro_time_t time;
    astro_vector_t pos;
    top_elliptical_t ellip;
    top_rectangular_t ecl, equ;
    double dx, dy, dz, arcmin;

    printf("tt = %0.1lf\n", tt);

    time = Astronomy_TerrestrialTime(tt);
    pos = Astronomy_HelioVector(BODY_PLUTO, time);
    if (pos.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "TestPluto(%0.0lf): HelioVector returned %d\n", tt, (int)pos.status);
        return 1;
    }

    printf("HelioVector : pos=(%23.16lf, %23.16lf, %23.16lf)\n", pos.x, pos.y, pos.z);

    /* Compare with untruncated TOP2013 model. */
    if (TopCalcElliptical(model, tt, &ellip)) return 1;
    if (TopEcliptic(model->planet, &ellip, &ecl)) return 1;
    if (TopEquatorial(&ecl, &equ)) return 1;
    printf("TOP2013     : pos=(%23.16lf, %23.16lf, %23.16lf)\n", equ.x, equ.y, equ.z);

    /* Calculate relative error and scale in arcminute units (as seen from the Sun). */
    dx = equ.x - pos.x;
    dy = equ.y - pos.y;
    dz = equ.z - pos.z;
    arcmin = (RAD2DEG * 60.0) * sqrt((dx*dx + dy*dy + dz*dz) / (equ.x*equ.x + equ.y*equ.y + equ.z*equ.z));
    printf("arcmin_error = %0.6lf\n", arcmin);
    printf("\n");

    return 0;
}

int main(void)
{
    top_model_t model;
    const double delta_tt = 18250.0;
    int error;

    TopInitModel(&model);
    error = TopLoadModel(&model, "../TOP2013.dat", 9);
    if (error) goto fail;

    /*
        The PlutoStateTable has exact values at tt = { ... , -36500, 0, +36500, ... }.
        I want to exercise the calculator at the midpoints tt = { ..., -18250, +18250, ... }.
        We also exercise at the exact points where we coincide with TOP2013, to verify we get
        exact match.
    */
    for (int n = -5; n <= 5; ++n)
    {
        double tt = n * delta_tt;
        if (TestPluto(tt, &model)) goto fail;
    }

fail:
    TopFreeModel(&model);
    return error;
}
