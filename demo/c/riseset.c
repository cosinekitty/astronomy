/*
    riseset.c  -  by Don Cross - 2019-06-14

    Example C program for Astronomy Engine:
    https://cosinekitty.github.io/astronomy/

    This program calculates sunrise, sunset, moonrise, and moonset
    times for an observer at a given latitude and longitude.
*/

#include <stdio.h>
#include "astro_demo_common.h"

void PrintTime(astro_time_t time)
{
    astro_utc_t utc;

    utc = Astronomy_UtcFromTime(time);
    printf("%04d-%02d-%02dT%02d:%02d:%06.3lfZ", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
}

int PrintEvent(const char *name, astro_search_result_t evt)
{
    printf("%-8s : ", name);
    if (evt.status == ASTRO_SUCCESS)
    {
        PrintTime(evt.time);
        printf("\n");
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
    astro_observer_t observer;
    astro_time_t time;
    astro_search_result_t sunrise, sunset, moonrise, moonset;

    if (ParseArgs(argc, argv, &observer, &time))
        return 1;

    printf("search   : ");
    PrintTime(time);
    printf("\n");

    sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 300.0);
    sunset   = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET,  time, 300.0);
    moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 300.0);
    moonset  = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET,  time, 300.0);

    if (PrintEvent("sunrise", sunrise))
        return 1;

    if (PrintEvent("sunset", sunset))
        return 1;

    if (PrintEvent("moonrise", moonrise))
        return 1;

    if (PrintEvent("moonset", moonset))
        return 1;

    return 0;
}
