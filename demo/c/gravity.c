/*
    gravity.c  -  Don Cross  -  2021-07-19

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

#include <stdio.h>
#include <math.h>
#include "astronomy.h"

const char * const UsageText =
"\n"
"    USAGE:\n"
"\n"
"    gravity latitude height\n"
"\n"
"    Calculates the gravitational acceleration experienced\n"
"    by an observer on the surface of the Earth at the specified\n"
"    latitude (degrees north of the equator) and height\n"
"    (meters above sea level).\n"
"    The output is the gravitational acceleration in m/s^2.\n"
"\n";

int main(int argc, const char *argv[])
{
    const double MAX_HEIGHT_METERS = 100000.0;
    double latitude, height, gravity;

    if (argc != 3)
    {
        fprintf(stderr, "%s", UsageText);
        return 1;
    }

    if (1 != sscanf(argv[1], "%lf", &latitude) || !isfinite(latitude) || latitude < -90.0 || latitude > +90.0)
    {
        fprintf(stderr, "ERROR: Invalid latitude '%s'. Must be a number between -90 and +90.\n", argv[1]);
        return 1;
    }

    if (1 != sscanf(argv[2], "%lf", &height) || !isfinite(height) || height < 0.0 || height > MAX_HEIGHT_METERS)
    {
        fprintf(stderr, "ERROR: Invalid height '%s'. Must be a number of meters between 0 and %0.0lf.\n", argv[2], MAX_HEIGHT_METERS);
        return 1;
    }

    gravity = Astronomy_ObserverGravity(latitude, height);
    printf("latitude = %8.4lf,  height = %6.0lf,  gravity = %8.6lf\n", latitude, height, gravity);
    return 0;
}
