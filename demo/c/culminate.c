/*
    culminate.c  -  by Don Cross - 2019-06-18

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This example program shows how to calculate the time
    the Sun, Moon, and planets will next reach their highest point in the sky
    as seen by an observer at a given location on the Earth.
    This is called culmination, and is found by finding when
    each body's "hour angle" is 0.

    Having an hour angle of 0 is another way of saying that the body is
    crossing the meridian, the imaginary semicircle in the sky that passes
    from due north on the horizon, through the zenith (straight up),
    toward due south on the horizon. At this moment the body appears to
    have an azimuth of either 180 degrees (due south) or 0 (due north).
*/

#include <stdio.h>
#include "astro_demo_common.h"

int PrintEvent(const char *name, astro_hour_angle_t evt)
{
    printf("%-8s : ", name);
    if (evt.status == ASTRO_SUCCESS)
    {
        PrintTime(evt.time);
        printf("  altitude=%6.2lf  azimuth=%6.2lf\n", evt.hor.altitude, evt.hor.azimuth);
        return 0;
    }
    else
    {
        printf("ERROR %d\n", evt.status);
        return 1;
    }
}

int main(int argc, const char *argv[])
{
    static const astro_body_t bodies[] =
    {
        BODY_SUN, BODY_MOON, BODY_MERCURY, BODY_VENUS, BODY_MARS,
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO
    };
    static const int nbodies = sizeof(bodies) / sizeof(bodies[0]);

    astro_observer_t observer;
    astro_time_t time;
    astro_hour_angle_t evt;
    int i;

    if (ParseArgs(argc, argv, &observer, &time))
        return 1;

    printf("search   : ");
    PrintTime(time);
    printf("\n");

    for (i=0; i < nbodies; ++i)
    {
        evt = Astronomy_SearchHourAngleEx(bodies[i], observer, 0.0, time, +1);
        PrintEvent(Astronomy_BodyName(bodies[i]), evt);
    }

    return 0;
}
