/*
    solar_time.c  -  by Don Cross - 2019-06-11

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Given an observer's geographic latitude and longitude,
    and an optional date and time, this program displays
    true solar time for that observer and time.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"

int main(int argc, const char *argv[])
{
    int error;
    astro_observer_t observer;
    astro_time_t time;
    astro_func_result_t ha;
    double solarTimeHours;
    int hour, minute, second, milli;

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    ha = Astronomy_HourAngle(BODY_SUN, &time, observer);
    if (ha.status != ASTRO_SUCCESS)
    {
        printf("ERROR %d in Astronomy_HourAngle().\n", ha.status);
        return 1;
    }

    solarTimeHours = fmod(ha.value + 12.0, 24.0);

    milli = (int)round(solarTimeHours * 3600000.0);
    second = milli / 1000;
    milli %= 1000;
    minute = second / 60;
    second %= 60;
    hour = minute / 60;
    minute %= 60;
    hour %= 24;

    printf("True solar time = %0.4lf hours (%02d:%02d:%02d.%03d)\n", solarTimeHours, hour, minute, second, milli);
    return 0;
}
