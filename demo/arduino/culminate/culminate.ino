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
char printBuffer[80];

int PrintEvent(const char *name, astro_hour_angle_t evt)
{
    snprintf(printBuffer, sizeof printBuffer,"%-8s : ", name);
    Serial.println(printBuffer);
    if (evt.status == ASTRO_SUCCESS)
    {
        PrintTime(evt.time);
        snprintf(printBuffer, sizeof printBuffer,"  altitude=%6.2lf  azimuth=%6.2lf\n", evt.hor.altitude, evt.hor.azimuth);
        Serial.println(printBuffer);
        return 0;
    }
    else
    {
        snprintf(printBuffer, sizeof printBuffer,"ERROR %d\n", evt.status);
        Serial.println(printBuffer);
        return 1;
    }
}


void setup() 
{
  Serial.begin(9600);
}

void loop() 
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

    // if (ParseArgs(argc, argv, &observer, &time))
    //     return 1;
        // const char *argvs[] = {"22.2","11.8","1988-02-02T12:55:33Z"};
    const char *argvs[] = {"22.2","33.8","2023-10-07T00:38:57+03:30"};
    int argcs = 3;
        int error;
    error = ParseArgs(argcs, argvs, &observer, &time);
    if (error)
        Serial.println("error.......");

    snprintf(printBuffer, sizeof printBuffer,"search   : ");
    Serial.println(printBuffer);
    PrintTime(time);
    snprintf(printBuffer, sizeof printBuffer,"\n");
    Serial.println(printBuffer);

    for (i=0; i < nbodies; ++i)
    {
        evt = Astronomy_SearchHourAngleEx(bodies[i], observer, 0.0, time, +1);
        PrintEvent(Astronomy_BodyName(bodies[i]), evt);
    }
    delay(1000);
    // return 0;
}
