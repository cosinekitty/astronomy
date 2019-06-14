/*
    moonphase.c  -  by Don Cross - 2019-05-25

    Example C program for Astronomy Engine:
    https://cosinekitty.github.io/astronomy/

    This program calculates the Moon's phase for a given date and time,
    or the computer's current date and time if none is given.
    It also finds the dates and times of the subsequent 10 quarter phase changes.
*/

#include <stdio.h>
#include "astro_demo_common.h"

static const char *QuarterName(int quarter)
{
    switch (quarter)
    {
    case 0:   return "New Moon";
    case 1:   return "First Quarter";
    case 2:   return "Full Moon";
    case 3:   return "Third Quarter";
    default:  return "INVALID QUARTER";
    }
}

static void PrintTime(astro_time_t time)
{
    astro_utc_t utc = Astronomy_UtcFromTime(time);
    printf("%04d-%02d-%02d %02d:%02d:%02.0lf UTC", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
}

int main(int argc, const char *argv[])
{
    astro_time_t time;
    astro_angle_result_t phase;
    astro_moon_quarter_t mq;
    int i;

    /* 
        Calculate the Moon's current phase angle,
        which ranges from 0 to 360 degrees.

          0 = new moon,
         90 = first quarter,
        180 = full moon,
        270 = third quarter.
    */

    switch (argc)
    {
    case 1:
        time = Astronomy_CurrentTime();
        break;

    case 2:
        if (ParseTime(argv[1], &time))
            return 1;
        break;

    default:
        fprintf(stderr, "USAGE: moonphase [date]\n");
        return 1;
    }

    phase = Astronomy_MoonPhase(time);
    if (phase.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_MoonPhase error %d\n", phase.status);
        return 1;
    }

    PrintTime(time);
    printf(" : Moon's phase angle = %0.6lf degrees.\n", phase.angle);

    /* Find the next 10 lunar quarter phases. */
    printf("\nThe next 10 lunar quarters are:\n");
    for (i=0; i < 10; ++i)
    {
        if (i == 0)
            mq = Astronomy_SearchMoonQuarter(time);
        else
            mq = Astronomy_NextMoonQuarter(mq);

        if (mq.status != ASTRO_SUCCESS)
        {
            printf("Error %d trying to find moon quarter.\n", mq.status);
            return 1;
        }

        PrintTime(mq.time);
        printf(" : %s\n", QuarterName(mq.quarter));
    }

    return 0;
}
