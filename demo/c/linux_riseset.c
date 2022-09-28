/*
    linux_riseset.c  -  by Don Cross - 2020-09-08

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates sunrise, sunset, moonrise, and moonset
    times for an observer at a given latitude and longitude.

    It uses timezone functions available in Linux to print
    the events the system's configured local time, instead of UTC.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "astro_demo_common.h"

void PrintLocalTime(astro_time_t time)
{
    astro_time_t linux_epoch;
    time_t seconds;
    struct tm local;
    const char *dow;

    /* Calculate the Linux Epoch as an astronomy time. */
    linux_epoch = Astronomy_MakeTime(1970, 1, 1, 0, 0, 0);

    /* Subtract the two UT values to find the number of seconds elapsed. */
    seconds = (time_t) round((time.ut - linux_epoch.ut) * 86400.0);

    /* Obtain the corresponding local time for this system's timezone. */
    localtime_r(&seconds, &local);

    switch (local.tm_wday)
    {
    case 0:  dow = "Sun"; break;
    case 1:  dow = "Mon"; break;
    case 2:  dow = "Tue"; break;
    case 3:  dow = "Wed"; break;
    case 4:  dow = "Thu"; break;
    case 5:  dow = "Fri"; break;
    case 6:  dow = "Sat"; break;
    default: dow = "---"; break;
    }

    printf("%04d-%02d-%02d %s %02d:%02d:%02d %s",
        local.tm_year+1900, local.tm_mon+1, local.tm_mday,
        dow,
        local.tm_hour, local.tm_min, local.tm_sec,
        (local.tm_isdst ? "DT" : "ST"));
}


typedef struct
{
    const char *name;
    astro_time_t time;
    double altitude;    /* angle of culmination above horizon; otherwise NAN */
}
event_t;


void AppendCulm(
    event_t evtlist[],
    int *evtcount,
    const char *name,
    astro_body_t body,
    astro_observer_t observer,
    astro_time_t search_time)
{
    astro_hour_angle_t culm;

    culm = Astronomy_SearchHourAngleEx(body, observer, 0.0, search_time, +1);
    if (culm.status == ASTRO_SUCCESS)
    {
        evtlist[*evtcount].name = name;
        evtlist[*evtcount].time = culm.time;
        evtlist[*evtcount].altitude = culm.hor.altitude;
        ++(*evtcount);
    }
    else
    {
        printf("ERROR %d finding culmination of %s\n", culm.status, name);
    }
}


void AppendEvent(const char *name, event_t evtlist[], int *evtcount, astro_search_result_t result)
{
    if (result.status == ASTRO_SUCCESS)
    {
        evtlist[*evtcount].name = name;
        evtlist[*evtcount].time = result.time;
        evtlist[*evtcount].altitude = NAN;
        ++(*evtcount);
    }
    else
    {
        printf("ERROR %d finding %s\n", result.status, name);
    }
}


int EventCompare(const void *aptr, const void *bptr)
{
    const event_t *a = aptr;
    const event_t *b = bptr;
    if (a->time.ut < b->time.ut)
        return -1;
    if (a->time.ut > b->time.ut)
        return +1;
    return 0;
}


int main(int argc, const char *argv[])
{
    astro_observer_t observer;
    astro_time_t time;
    astro_search_result_t result;
    event_t evtlist[6];
    int evtcount = 0;
    int i;

    if (ParseArgs(argc, argv, &observer, &time))
        return 1;

    result = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 300.0);
    AppendEvent("sunrise", evtlist, &evtcount, result);

    result = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET,  time, 300.0);
    AppendEvent("sunset", evtlist, &evtcount, result);

    result = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 300.0);
    AppendEvent("moonrise", evtlist, &evtcount, result);

    result = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET,  time, 300.0);
    AppendEvent("moonset", evtlist, &evtcount, result);

    AppendCulm(evtlist, &evtcount, "moonculm", BODY_MOON, observer, time);
    AppendCulm(evtlist, &evtcount, "sunculm",  BODY_SUN, observer, time);

    /* Sort the events in chronological order. */
    qsort(evtlist, (size_t)evtcount, sizeof(event_t), EventCompare);

    /* Print the sorted events. */
    for (i=0; i < evtcount; ++i)
    {
        printf("%-8s : ", evtlist[i].name);
        PrintLocalTime(evtlist[i].time);
        if (isfinite(evtlist[i].altitude))
            printf("   alt = %0.2lf", evtlist[i].altitude);
        printf("\n");
    }

    return 0;
}
